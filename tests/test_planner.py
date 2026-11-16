"""Experiment planner: power maths, variance-component recovery, interval coverage, API speed, outputs."""
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kerncount import planner
from kerncount import settings as st_mod


def _settings(**over):
    s = st_mod.load(Path("/nonexistent"))
    for sec, kv in over.items():
        s[sec].update(kv)
    return s


def simulate(rng, n=6, treatments=6, wells=2, fields=6, s_field=0.3, s_pos=0.2, s_well=0.1, s_bio=0.2,
             s_shift=0.1, level=80.0, effect=0.6, drop=()):
    """Nested log-normal data: image -> (well, position) consecutive, as the default [design].
    log value = log level + log effect (treatments) + shift_n (treatments) + bio_{n,t} + well + pos_{n,p} + field.
    The per-n log ratio treatment / Control then has variance s_shift^2 + 2 s_bio^2 + measurement."""
    rows = []
    names = ["Control", "LyoVec"] + [f"miR-{i}" for i in range(treatments)]
    for i in range(n):
        pos = rng.normal(0, s_pos, fields)
        shift = rng.normal(0, s_shift)
        for t in names:
            mu = math.log(level) + rng.normal(0, s_bio)
            if t.startswith("miR"):
                mu += math.log(effect) + shift
            for w in range(wells):
                wl = rng.normal(0, s_well)
                for p in range(fields):
                    img = w * fields + p + 1
                    if (f"n{i + 1}", t, img) in drop:
                        continue
                    rows.append(dict(n=f"n{i + 1}", treatment=t, image=img,
                                     value=math.exp(mu + wl + pos[p] + rng.normal(0, s_field))))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- power maths

def test_power_matches_known_value_and_simulation():
    # one-sample / paired t, d = 1, n = 10, alpha 0.05: power 0.803 (standard tables)
    assert planner.power(1.0, 1.0, 10) == pytest.approx(0.8031, abs=0.002)
    assert planner.power(0.0, 1.0, 10) == pytest.approx(0.05, abs=1e-6)
    rng = np.random.default_rng(0)
    x = rng.normal(0.4, 0.5, (20000, 6))
    sim = (stats.ttest_1samp(x, 0, axis=1).pvalue < 0.05).mean()
    assert planner.power(0.4, 0.5, 6) == pytest.approx(sim, abs=0.015)
    # nc for 80% power reproduces 80% power, also for df = 1 where scipy's nct is fragile
    for d in (1, 2, 3, 10):
        nc = planner.nc_for_power(d, 0.8)
        assert planner.power(nc / math.sqrt(d + 1), 1.0, d + 1) == pytest.approx(0.8, abs=1e-4)
    assert np.isnan(planner.power(0.5, 0.5, 1))


def test_n_needed_monotone():
    a = planner.n_needed(math.log(1 / 0.7), [0.05, 0.1, 0.2, 0.4])
    assert list(a) == sorted(a) and a[0] >= 2


# ---------------------------------------------------------------- variance model

def test_component_recovery_on_simulated_nested_data():
    rng = np.random.default_rng(1)
    df = simulate(rng, n=12, treatments=8, wells=3, fields=8, s_field=0.3, s_pos=0.2, s_well=0.12, s_bio=0.15,
                  s_shift=0.1)
    plan = planner.build(df, _settings(design={"wells": 3}))
    c = plan["components_point"]
    cv2 = lambda s: math.exp(s ** 2) - 1  # relative-scale variance of a log-normal
    assert c["field"] == pytest.approx(cv2(0.3), rel=0.2)
    assert c["position"] == pytest.approx(cv2(0.2), rel=0.6)
    assert c["well"] == pytest.approx(cv2(0.12), rel=0.6)
    b = plan["bio"]["Control|treatments"]
    true_bio = 0.1 ** 2 + 2 * 0.15 ** 2
    assert b["var"] == pytest.approx(true_bio, rel=0.6)
    assert b["ci95"][0] <= true_bio <= b["ci95"][1]
    assert plan["position_share_of_field_variance"] > 0.15


def test_component_intervals_cover_truth():
    """95% GPQ intervals of field, well and the per-n ratio variance bracket the truth in most simulations."""
    rng = np.random.default_rng(2)
    s_field, s_well, s_bio, s_shift = 0.3, 0.15, 0.2, 0.1
    hit = dict(field=0, well=0, V=0)
    reps = 40
    for _ in range(reps):
        df = simulate(rng, n=4, treatments=5, wells=2, fields=6, s_field=s_field, s_pos=0.0, s_well=s_well,
                      s_bio=s_bio, s_shift=s_shift)
        df = df.assign(well=(df["image"] - 1) // 6 + 1, position=(df["image"] - 1) % 6 + 1)
        wc = planner.within_components(df, positions=False)
        d = wc.pop("data")
        wd = planner.within_draws(wc, rng, 1000)
        hit["field"] += np.quantile(wd["field"], 0.025) <= math.exp(s_field ** 2) - 1 <= np.quantile(wd["field"], 0.975)
        tw = math.exp(s_well ** 2) - 1
        hit["well"] += np.quantile(wd["well"], 0.025) <= tw <= np.quantile(wd["well"], 0.975)
        lr = planner.log_ratios(d, "Control")
        b = planner.bio_fit(lr[[c for c in lr if c.startswith("miR")]], planner._meas_coef(d), "Control")
        Vd, _ = planner.bio_draws(b, wd, rng, 1000)
        true_V = s_shift ** 2 + 2 * s_bio ** 2 + 2 * (tw / 2 + (math.exp(s_field ** 2) - 1) / 12)
        hit["V"] += np.quantile(Vd, 0.025) <= true_V <= np.quantile(Vd, 0.975)
    for k, v in hit.items():
        assert v / reps >= 0.85, (k, v)


def test_classes_from_settings_and_default():
    s = _settings()
    cls, src, notes = planner.treatment_classes(s, ["Control", "LyoVec", "miR-1", "CL075"], ["Control", "LyoVec"])
    assert cls == {"Control": "controls", "LyoVec": "controls", "miR-1": "treatments", "CL075": "treatments"}
    s["analysis"]["classes"] = {"mirna": ["miR-*"], "positive": "CL075"}
    cls, src, notes = planner.treatment_classes(s, ["Control", "miR-1", "miR-2", "CL075"], ["Control"])
    assert cls["miR-2"] == "mirna" and cls["CL075"] == "positive" and cls["Control"] == "controls"
    assert src == "analysis.classes"


# ---------------------------------------------------------------- API

@pytest.fixture(scope="module")
def plan():
    return planner.build(simulate(np.random.default_rng(3), n=4), _settings())


def test_evaluate_bands_ordered_and_fast(plan):
    r = planner.evaluate(plan, 4, 2, 6, 0.7)
    for k in ("power", "ci_halfwidth_pct", "mde_ratio", "sd_log_ratio"):
        assert r[k]["low"] <= r[k]["mid"] <= r[k]["high"], k
    assert 0 < r["mde_ratio"]["mid"] < 1 and r["cost"] > 0 and r["notes"]
    # more n's: more power, narrower CI
    r8 = planner.evaluate(plan, 8, 2, 6, 0.7)
    assert r8["power"]["mid"] > r["power"]["mid"] and r8["ci_halfwidth_pct"]["mid"] < r["ci_halfwidth_pct"]["mid"]
    # random positions are noisier than standardised ones
    assert planner.evaluate(plan, 4, 2, 6, 0.7, positions="random")["sd_log_ratio"]["mid"] > r["sd_log_ratio"]["mid"]
    assert planner.evaluate(plan, 4, 2, 6, 0.7, reference_wells="auto")["reference_wells"] > 2
    obs = planner.evaluate(plan, 4, 2, 6, "miR-0")
    assert obs["effect_label"] == "miR-0" and obs["effect_ratio"] < 1
    t = time.perf_counter()
    for _ in range(50):
        planner.evaluate(plan, 5, 3, 7, 0.6)
    assert (time.perf_counter() - t) / 50 < 0.05
    with pytest.raises(ValueError):
        planner.evaluate(plan, 1, 2, 6, 0.7)


def test_curve(plan):
    pts = planner.curve(plan, "n", {"wells": 2, "fields": 6}, 0.7)
    assert [p["x"] for p in pts] == list(range(2, 13))
    assert all(a["power"]["mid"] <= b["power"]["mid"] for a, b in zip(pts, pts[1:]))
    pts = planner.curve(plan, "fields", {"n": 4}, 0.7, values=[1, 6, 18])
    assert pts[0]["ci_halfwidth_pct"]["mid"] > pts[-1]["ci_halfwidth_pct"]["mid"]


def test_more_n_sane():
    rng = np.random.default_rng(4)
    p = planner.build(simulate(rng, n=3, effect=0.5, s_bio=0.1, s_shift=0.05), _settings())
    tab, summ = planner.more_n(p, sims=2000)
    g = tab[(tab["reference"] == "Control") & (tab["treatment"] == "miR-0")].sort_values("k_added")
    assert g["prob_p_below_alpha"].iloc[-1] > 0.8  # a strong effect is almost surely confirmed
    assert g["prob_p_below_alpha"].iloc[-1] >= g["prob_p_below_alpha"].iloc[0] - 0.05
    assert tab["prob_p_below_alpha"].between(0, 1).all()
    assert summ["false_positive_optional_stopping"] > 0.08  # optional stopping inflates the 5%
    lv = tab[(tab["reference"] == "Control") & (tab["treatment"] == "LyoVec")]
    # no true effect: clearly less promising than the strong effect (it is centred on the observed mean of 3 n's,
    # so it need not be small; see docs/planner.md)
    assert lv["prob_p_below_alpha"].mean() < g["prob_p_below_alpha"].mean() - 0.2


# ---------------------------------------------------------------- outputs and odd designs

def test_run_writes_everything_and_load(tmp_path):
    df = simulate(np.random.default_rng(5), n=4, drop={("n2", "LyoVec", i) for i in range(1, 13)})
    out = tmp_path / "prj" / "sampling" / "g1"
    (out).mkdir(parents=True)
    (out / "REPORT.txt").write_text("How many images are needed?\nold report\n")
    planner.run(df, _settings(), out, sims=500)
    for f in ("planner.json", "design_forward.csv", "mde.csv", "cost_frontier.csv", "allocation.csv", "more_n.csv",
              "manifest.json", "variance_components.png", "power_vs_n.png", "ci_width_vs_fields.png",
              "cost_frontier.png", "more_n.png"):
        assert (out / f).exists(), f
    rep = (out / "REPORT.txt").read_text()
    assert rep.startswith("Recommendations for the next experiment") and "old report" in rep
    assert "n2 lacks LyoVec" in rep and "NOT a stopping rule" in rep
    planner.run(df, _settings(), out, sims=500)  # re-run replaces, does not stack, the recommendations
    rep2 = (out / "REPORT.txt").read_text()
    assert rep2.count("Recommendations for the next experiment") == 1 and "old report" in rep2
    man = json.loads((out / "manifest.json").read_text())["files"]
    assert {"path", "kind", "title", "description"} <= set(man[0])
    assert any(m["path"] == "more_n.png" for m in man)
    p = planner.load(tmp_path / "prj")
    r = planner.evaluate(p, 4, 2, 6, 0.7)
    assert 0 <= r["power"]["mid"] <= 1
    fwd = pd.read_csv(out / "design_forward.csv")
    assert {"power_low", "power_mid", "power_high", "ci_halfwidth_pct_mid", "cost"} <= set(fwd.columns)
    assert fwd["fields_per_well"].max() == 18 and fwd["n"].max() == 12 and fwd["wells"].max() == 6
    al = pd.read_csv(out / "allocation.csv")
    assert al["recommended"].sum() == 1 and (al.loc[al["recommended"], "wells_used"] <= 24).all()


def test_unbalanced_single_well_and_random_positions_do_not_crash(tmp_path):
    rng = np.random.default_rng(6)
    drop = {("n1", "miR-0", 12), ("n2", "Control", 3), ("n3", "miR-2", 7), ("n3", "miR-2", 8)}
    planner.run(simulate(rng, n=3, drop=drop), _settings(), tmp_path / "unb", sims=200)
    one = simulate(rng, n=3, wells=1, fields=6)
    plan = planner.run(one, _settings(design={"wells": 1}), tmp_path / "single", sims=200)
    assert plan["components"]["well"]["var"] is None
    assert any("single well" in n for n in plan["notes"])
    planner.evaluate(plan, 4, 2, 6, 0.7)
    plan = planner.run(simulate(rng, n=3), _settings(design={"field_positions": "random"}), tmp_path / "rnd", sims=200)
    assert plan["components"]["position"]["var"] is None
    one_field = simulate(rng, n=3, wells=2, fields=1)
    planner.run(one_field, _settings(design={"fields_per_well": 1}), tmp_path / "onefield", sims=200)


def test_planner_settings_and_form():
    from kerncount import settings_form

    assert list(st_mod.DEFAULTS)[-1] == "planner"
    assert set(settings_form.FIELDS["planner"]) == set(st_mod.DEFAULTS["planner"])
    import tomllib

    parsed = tomllib.loads(st_mod.TEMPLATE.format(raw_dir="x"))
    assert list(parsed)[-1] == "planner" and parsed["planner"] == st_mod.DEFAULTS["planner"]
    o = planner.options(_settings(planner={"effects": ["0.7", "0,5"], "plate_format": "96", "reference_wells": "4"}))
    assert o["effects"] == [0.7, 0.5] and o["plate_format"] == 96 and o["reference_wells"] == 4
    with pytest.raises(ValueError):
        planner.options(_settings(planner={"plate_format": 30}))
    with pytest.raises(ValueError):
        planner.options(_settings(planner={"effects": [1.0]}))
