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
    plan = planner.build(df, _settings(design={"wells": 3}, analysis={"classes": {}}))  # pooled class (controls vs the rest)
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
    # default "auto" = the statistics' rule: one class per treatment, each control its own
    assert cls == {"Control": "Control", "LyoVec": "LyoVec", "miR-1": "miR-1", "CL075": "CL075"}
    assert src.startswith("auto")
    s["analysis"]["classes"] = {}  # "two": controls vs the rest
    cls, src, notes = planner.treatment_classes(s, ["Control", "LyoVec", "miR-1", "CL075"], ["Control", "LyoVec"])
    assert cls == {"Control": "controls", "LyoVec": "controls", "miR-1": "treatments", "CL075": "treatments"}
    s["analysis"]["classes"] = '{ mirna = ["miR-*"] }'  # TOML text form
    cls, src, notes = planner.treatment_classes(s, ["Control", "miR-1", "CL075"], ["Control"])
    assert cls == {"Control": "controls", "miR-1": "mirna", "CL075": "other"}
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
    assert set(settings_form.FORM_ITEMS["planner"]) == set(st_mod.DEFAULTS["planner"])
    import tomllib

    parsed = tomllib.loads(st_mod.TEMPLATE.format(raw_dir="x"))
    assert list(parsed)[-1] == "planner" and parsed["planner"] == st_mod.DEFAULTS["planner"]
    o = planner.options(_settings(planner={"effects": ["0.7", "0,5"], "plate_format": "96", "reference_wells": "4"}))
    assert o["effects"] == [0.7, 0.5] and o["plate_format"] == 96 and o["reference_wells"] == 4
    with pytest.raises(ValueError):
        planner.options(_settings(planner={"plate_format": 30}))
    with pytest.raises(ValueError):
        planner.options(_settings(planner={"effects": [1.0]}))


# ---------------------------------------------------------------- I1: power for a star (family-adjusted)

def test_family_critical_matches_the_statistics():
    from kerncount import models

    for k, rho, df in ((7, 0.5, 3.0), (7, 0.3, 6.0), (4, 0.0, 10.0), (7, 0.9, 4.5)):
        R = np.full((k, k), rho)
        np.fill_diagonal(R, 1)
        _, crit = models._mvt_adjust(np.full(k, 2.0), np.full(k, df), R)
        assert planner.family_critical(k, rho, df) == pytest.approx(crit[0], abs=0.01), (k, rho, df)
    assert planner.family_critical(1, 0.5, 3) == pytest.approx(stats.t.ppf(0.975, 3))


def test_evaluate_uses_the_family_adjusted_level_by_default(plan):
    adj = planner.evaluate(plan, 4, 2, 6, 0.7)
    raw = planner.evaluate(plan, 4, 2, 6, 0.7, adjust_for_family=False)
    k = plan["design"]["conditions"] - 1
    assert adj["adjust_for_family"] and adj["family_size"] == k and 0 <= adj["family_correlation"] < 1
    assert adj["t_critical"] == pytest.approx(planner.family_critical(k, adj["family_correlation"], adj["df"]),
                                              abs=0.01)
    assert 3 <= adj["df"] <= 6  # Satterthwaite between n - 1 and 2 (n - 1): one class per treatment
    assert not raw["adjust_for_family"] and raw["df"] == 3 and raw["t_critical"] == pytest.approx(stats.t.ppf(0.975, 3))
    assert adj["power"]["mid"] < raw["power"]["mid"]
    assert adj["mde_ratio"]["mid"] < raw["mde_ratio"]["mid"] < 1  # a star needs a bigger effect
    assert "star" in adj["test"] and "unadjusted" in raw["test"]
    assert any("share one adjustment" in x for x in adj["notes"])
    off = dict(plan, options=dict(plan["options"], adjust_for_family=False))
    assert planner.evaluate(off, 4, 2, 6, 0.7)["t_critical"] == raw["t_critical"]  # [planner] adjust_for_family
    pts = planner.curve(plan, "n", {"wells": 2, "fields": 6}, 0.7, adjust_for_family=False)
    assert pts[2]["t_critical"] == pytest.approx(stats.t.ppf(0.975, pts[2]["x"] - 1))


def test_tables_and_report_use_the_same_test(tmp_path):
    df = simulate(np.random.default_rng(7), n=4)
    out = tmp_path / "g1"
    plan = planner.run(df, _settings(), out, sims=300)
    mde = pd.read_csv(out / "mde.csv")
    key = (plan["design"]["main_reference"], plan["default_class"])
    row = mde[(mde["reference"] == key[0]) & (mde["class"] == key[1]) & (mde["n"] == 5) & (mde["wells"] == 2)
              & (mde["fields_per_well"] == 6)].iloc[0]
    ev = planner.evaluate(plan, 5, 2, 6, 0.7)
    assert row["mde_ratio_mid"] == pytest.approx(ev["mde_ratio"]["mid"], abs=0.01)  # coarse grid, same test
    rep = (out / "REPORT.txt").read_text()
    assert "chance of a STAR" in rep and "single-step multivariate t" in rep
    plan2 = planner.run(df, _settings(planner={"adjust_for_family": False}), tmp_path / "g2", sims=300)
    assert plan2["options"]["adjust_for_family"] is False
    rep2 = (tmp_path / "g2" / "REPORT.txt").read_text()
    assert "unadjusted ratio paired t-test" in rep2
    m2 = pd.read_csv(tmp_path / "g2" / "mde.csv")
    row2 = m2[(m2["reference"] == key[0]) & (m2["class"] == key[1]) & (m2["n"] == 5) & (m2["wells"] == 2)
              & (m2["fields_per_well"] == 6)].iloc[0]
    assert row2["mde_ratio_mid"] > row["mde_ratio_mid"]  # unadjusted: smaller effects look detectable


def _star_rate(seed):
    from kerncount import metrics, models

    d = simulate(np.random.default_rng(seed), n=4, **STAR_KW)
    d = d.assign(well=(d["image"] - 1) // 6 + 1, position=(d["image"] - 1) % 6 + 1, timepoint="")
    s = _settings()
    r = models.analyse(d, metrics.Metric("value", "value", "v", "a.u.", "positive"), s, positions_ok=False,
                       leave_one_out=False, density=False)
    c = r.contrasts
    c = c[(c["family"] == f"vs {s['analysis']['reference']}") & c["group"].str.startswith("miR")]
    return float((c["p_adj"] < 0.05).mean())


STAR_KW = dict(s_field=0.3, s_pos=0.0, s_well=0.1, s_bio=0.15, s_shift=0.05, effect=0.75)


@pytest.mark.slow
def test_adjusted_power_predicts_the_share_of_stars():
    """The planner's power (plan from 40 simulated n's) vs the share of stars the mixed model actually gives at
    4 n's: the family-adjusted prediction is close; the unadjusted paired-t power (the old default) is not."""
    plan = planner.build(simulate(np.random.default_rng(99), n=40, **STAR_KW), _settings())
    adj = planner.evaluate(plan, 4, 2, 6, STAR_KW["effect"], positions="random")["power"]
    raw = planner.evaluate(plan, 4, 2, 6, STAR_KW["effect"], positions="random", adjust_for_family=False)["power"]
    rate = float(np.mean([_star_rate(s) for s in range(80)]))
    assert adj["low"] - 0.05 <= rate <= adj["high"] + 0.05, (rate, adj)
    assert abs(rate - adj["mid"]) < abs(rate - raw["mid"])


def test_bio_fit_keeps_the_pooled_variance_when_the_split_is_clipped():
    """M9: an interaction mean square bigger than the pooled variance must not inflate V."""
    ns = [f"n{i}" for i in range(4)]
    # unbalanced (two missing cells): the additive fit's interaction mean square exceeds the pooled variance
    lr = pd.DataFrame(np.random.default_rng(24).normal(0, 0.3, (4, 3)), index=ns, columns=list("abc"))
    lr.iloc[0, 1] = lr.iloc[1, 2] = np.nan
    coef = pd.DataFrame([dict(n=n, t=t, a_well=0.1, a_field=0.01) for n in ns for t in ("a", "b", "c", "ref")]
                        ).set_index(["n", "t"])
    b = planner.bio_fit(lr, coef, "ref")
    pooled = float(((lr - lr.mean()) ** 2).sum().sum() / (lr.notna().sum() - 1).sum())
    assert b["kind"] == "two-way" and b["split_clipped"] and b["V"] == pytest.approx(pooled)
    assert b["ms_e_raw"] * (1 - 1 / b["tbar"]) > pooled  # what used to inflate V
    assert b["ms_n"] / b["tbar"] + (1 - 1 / b["tbar"]) * b["ms_e"] == pytest.approx(pooled, rel=1e-6)


# ---------------------------------------------------------------- I2: one plan per cell line, enforced acceptance

def test_sampling_blocks_per_cell_line_without_enforced_failures(tmp_path):
    from kerncount import pipeline

    used = pd.concat([simulate(np.random.default_rng(1), n=3).assign(cell_line="U251"),
                      simulate(np.random.default_rng(2), n=3).assign(cell_line="LN229",
                                                                     n=lambda d: d["n"] + "b")])
    g = tmp_path / "g1"
    g.mkdir()
    pd.DataFrame([dict(n="n1", criterion="overall", value="", limit="", result="fail", detail=""),
                  dict(n="n2", criterion="overall", value="", limit="", result="pass", detail="")]
                 ).to_csv(g / "acceptance.csv", index=False)
    blocks, dropped = pipeline.sampling_blocks(used, g, _settings())
    assert dropped == [] and [b[0] for b in blocks] == ["LN229", "U251"]
    blocks, dropped = pipeline.sampling_blocks(used, g, _settings(acceptance={"enforce": True}))
    assert dropped == ["n1"] and set(dict(blocks)["U251"]["n"]) == {"n2", "n3"}
    assert set(dict(blocks)["LN229"]["n"]) == {"n1b", "n2b", "n3b"}
    one = used[used["cell_line"] == "U251"].copy()
    one.loc[one["n"] == "n3", "cell_line"] = ""
    blocks, _ = pipeline.sampling_blocks(one, g, _settings())
    assert [b[0] for b in blocks] == ["U251 + unlabelled"]


def test_plan_per_cell_line_load_and_page(tmp_path):
    from kerncount import gui

    prj = tmp_path / "ws" / "analysis" / "E"
    root = prj / "sampling" / "g1"
    for line, seed in (("U251", 1), ("LN229", 2)):
        planner.run(simulate(np.random.default_rng(seed), n=3), _settings(), root / line, sims=200, cell_line=line,
                    excluded_ns=["n9"])
    (root / "cell_lines.json").write_text(json.dumps({"cell_lines": [
        {"cell_line": "U251", "path": "U251", "n_count": 3}, {"cell_line": "LN229", "path": "LN229", "n_count": 3}],
        "excluded_ns": ["n9"]}))
    assert planner.load(prj, "g1", "LN229")["cell_line"] == "LN229"
    assert planner.load(prj)["cell_line"] == "U251"  # default: the first line
    rep = (root / "LN229" / "REPORT.txt").read_text()
    assert "cell line LN229 only" in rep and "n9 left out" in rep
    info = gui.planner_info(prj, "g1", "LN229")
    assert info["cell_line"] == "LN229" and [c["cell_line"] for c in info["cell_lines"]] == ["U251", "LN229"]
    assert info["base"] == "sampling/g1/LN229/" and info["recommendations"].startswith("Recommendations")
    r = gui.planner_evaluate(prj, {"gen": "g1", "cell_line": "LN229", "n": 4, "adjust_for_family": "0"})
    assert r["result"]["adjust_for_family"] is False and r["curve"]
    assert gui.planner_evaluate(prj, {"gen": "g1", "cell_line": "LN229", "n": 4})["result"]["adjust_for_family"]
