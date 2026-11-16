"""Hierarchical model, comparison families, robustness, acceptance and the generation output contracts."""
import json

import numpy as np
import pandas as pd
import pytest

from kerncount import acceptance, generations, manifest, metrics, models, settings_form
from kerncount import settings as st_mod

TR = ["Control", "LyoVec", "CL075", "A", "B", "C"]


def _settings(**analysis):
    a = dict(reference="LyoVec", normalise_to=["Control", "LyoVec"], treatment_order=TR) | analysis
    return st_mod._merge(st_mod.DEFAULTS, {"analysis": a})


def simulate(rng, n_n=4, effects=None, var_n=0.05, var_ctrl=0.0, var_trt=0.1, var_well=0.03, var_pos=0.1,
             var_res=0.2, wells=2, fields=6, base=4.0, treatments=TR, drop=()):
    """Log-normal field values with every level of the hierarchy. effects: {treatment: log ratio vs baseline}."""
    effects = effects or {}
    rows = []
    for i in range(n_n):
        n = f"n{i + 1}"
        bn = rng.normal(0, np.sqrt(var_n))
        posfx = rng.normal(0, np.sqrt(var_pos), fields)
        for t in treatments:
            if (n, t) in drop:
                continue
            v = var_ctrl if t in ("Control", "LyoVec", "CL075") else var_trt
            bnt = rng.normal(0, np.sqrt(v)) if v > 0 else 0.0
            for w in range(1, wells + 1):
                bw = rng.normal(0, np.sqrt(var_well))
                for k in range(fields):
                    y = base + effects.get(t, 0) + bn + bnt + bw + posfx[k] + rng.normal(0, np.sqrt(var_res))
                    rows.append(dict(n=n, treatment=t, well=w, position=k + 1, image=(w - 1) * fields + k + 1,
                                     timepoint="", value=float(np.exp(y))))
    return pd.DataFrame(rows)


POS = metrics.Metric("x", "x", "x", "a.u.", "positive")


def _row(res, family, comparison):
    c = res.contrasts
    return c[(c["family"] == family) & (c["comparison"] == comparison)].iloc[0]


# ---------------------------------------------------------------- model

def test_model_recovers_known_effects_and_absolute_means():
    rng = np.random.default_rng(3)
    d = simulate(rng, n_n=6, effects={"A": np.log(0.5), "B": np.log(0.3)}, base=np.log(60))
    res = models.analyse(d, POS, _settings(), positions_ok=True, leave_one_out=False)
    assert res.ok and res.status == "ok"
    a, b, c = (_row(res, "vs LyoVec", f"{t} vs LyoVec") for t in ("A", "B", "C"))
    assert a.ci_low < 0.5 < a.ci_high and abs(np.log(a.estimate / 0.5)) < 0.35
    assert b.ci_low < 0.3 < b.ci_high and b.p_adj < 0.05
    assert c.ci_low < 1 < c.ci_high
    assert a.estimate_type == "ratio" and np.isfinite(a.df) and a.df >= 1
    assert a.estimand == "ratio of mean a.u." and res.unit == "well"
    e = res.emm.set_index("treatment")
    assert 30 < e.loc["LyoVec", "estimate"] < 120  # absolute units (mean per field), not log
    assert e.loc["LyoVec", "ci_low"] < e.loc["LyoVec", "estimate"] < e.loc["LyoVec", "ci_high"]
    assert set(res.varcomp["component"]) >= {"n (experiment: day, plate, passage)",
                                             "well (residual: well-to-well within n × treatment)"}
    # one n x treatment variance per non-control treatment (classes = "auto"), controls pooled
    assert {"n × treatment [controls]", "n × treatment [A]", "n × treatment [CL075]"} <= set(res.varcomp["component"])
    res = models.analyse(d, POS, _settings(model_unit="field"), positions_ok=True, leave_one_out=False)
    assert res.ok and res.unit == "field"
    assert set(res.varcomp["component"]) >= {"n (experiment: day, plate, passage)", "well (within n × treatment)",
                                             "field position (within n)", "field (residual)"}


def test_class_specific_variances_are_recovered():
    rng = np.random.default_rng(5)
    d = simulate(rng, n_n=10, var_ctrl=0.0, var_trt=0.25, var_n=0.02)
    s = _settings(classes={"mirna": ["A", "B", "C"]})
    res = models.analyse(d, POS, s, positions_ok=True, leave_one_out=False)
    assert "class-specific" in res.model
    vc = res.varcomp.set_index("component")["variance"]
    assert vc["n × treatment [mirna]"] > 0.08
    assert vc["n × treatment [controls]"] < 0.05
    # comparisons that rely on the variable class get fewer df / wider CIs than control-vs-control
    assert _row(res, "vs LyoVec", "A vs LyoVec").se > 2 * _row(res, "vs LyoVec", "Control vs LyoVec").se


def test_null_type_one_error_is_roughly_nominal_at_four_n():
    """Quick null simulation (the docs report a larger one): per-comparison rejection rate near 5%."""
    rng = np.random.default_rng(11)
    ps = []
    for _ in range(25):
        res = models.analyse(simulate(rng, n_n=4), POS, _settings(vehicle_check=""), positions_ok=True,
                             leave_one_out=False)
        assert res.ok
        ps += res.contrasts.loc[res.contrasts["family"] == "vs LyoVec", "p"].tolist()
    rate = float(np.mean(np.asarray(ps) < 0.05))
    assert rate <= 0.12, rate  # 125 correlated tests: nominal 0.05; a Wald z would give ~0.15-0.25 here


def test_fallback_ladder_and_rm_anova_fallback(monkeypatch):
    rng = np.random.default_rng(2)
    d = simulate(rng)
    real = models._fit_spec

    def fail_full(p, cfg, spec, *a, **k):
        if spec == "full":
            return None, "forced failure", []
        return real(p, cfg, spec, *a, **k)

    monkeypatch.setattr(models, "_fit_spec", fail_full)
    res = models.analyse(d, POS, _settings(), positions_ok=True, leave_one_out=False)
    assert res.ok and res.spec == "common" and "forced failure" in res.attempts[0]

    monkeypatch.setattr(models, "_fit_spec", lambda *a, **k: (None, "forced failure", []))
    res = models.analyse(d, POS, _settings(), positions_ok=True)
    assert not res.ok and res.status == "fallback_rm_anova" and res.primary_test == "rm_anova"
    assert len(res.loo) and set(res.loo["estimate_type"].dropna()) == {"difference"}  # LOO with the RM ANOVA


def test_equivalence_tost():
    rng = np.random.default_rng(8)
    same = simulate(rng, n_n=8, var_n=0.02, var_res=0.05, var_pos=0.02, var_well=0.005)
    res = models.analyse(same, POS, _settings(reference="Control", vehicle_check="LyoVec vs Control"),
                         positions_ok=True, leave_one_out=False)
    r = _row(res, "vehicle check (equivalence)", "LyoVec vs Control")
    assert r.conclusion == "equivalent" and 1 / 1.2 < r.ci_low and r.ci_high < 1.2 and r.p < 0.05
    far = simulate(rng, n_n=8, effects={"LyoVec": np.log(2.0)}, var_n=0.02, var_res=0.05, var_pos=0.02)
    res = models.analyse(far, POS, _settings(reference="Control"), positions_ok=True, leave_one_out=False)
    assert _row(res, "vehicle check (equivalence)", "LyoVec vs Control").conclusion == "not equivalent"
    few = simulate(rng, n_n=3, var_res=0.4)
    res = models.analyse(few, POS, _settings(reference="Control"), positions_ok=True, leave_one_out=False)
    assert _row(res, "vehicle check (equivalence)", "LyoVec vs Control").conclusion == "inconclusive"


def test_families_positive_control_and_among_classes():
    rng = np.random.default_rng(4)
    d = simulate(rng, effects={"A": -1.0})
    s = _settings(positive_control="CL075", classes={"mirna": ["A", "B", "C"]}, among_classes=["mirna"],
                  primary_family="Control")
    res = models.analyse(d, POS, s, positions_ok=True, leave_one_out=False)
    fams = set(res.contrasts["family"])
    assert {"vs Control", "vs LyoVec", "vs positive control CL075", "within class mirna (all pairs)",
            "vehicle check (equivalence)"} <= fams
    assert set(res.contrasts.loc[res.contrasts["family_kind"] == "positive_control", "group"]) == {"A", "B", "C"}
    assert len(res.contrasts[res.contrasts["family_kind"] == "among"]) == 3
    adj = res.contrasts[res.contrasts["family"] == "vs Control"]
    assert (adj["p_adj"] >= adj["p"] - 1e-9).all()


def test_leave_one_out_and_per_n_ratios():
    rng = np.random.default_rng(6)
    d = simulate(rng, effects={"A": np.log(0.4)})
    res = models.analyse(d, POS, _settings(), positions_ok=True)
    assert set(res.loo["dropped_n"]) == {"n1", "n2", "n3", "n4"}
    s = res.loo_summary.set_index("comparison")
    assert set(s.columns) >= {"loo_min", "loo_max", "stable", "same_direction", "same_significance"}
    assert s.loc["A vs LyoVec", "same_direction"]
    assert set(s["stable"]) <= {"yes", "no"}
    pn = res.per_n[(res.per_n["control"] == "LyoVec") & (res.per_n["treatment"] == "A")]
    assert len(pn) == 4 and (pn["ratio"] < 1).all()


def test_unbalanced_and_missing_reference():
    rng = np.random.default_rng(9)
    d = simulate(rng, drop={("n2", "LyoVec"), ("n3", "C"), ("n4", "C")}, effects={"A": -0.7})
    d = d[~((d["n"] == "n1") & (d["treatment"] == "B") & (d["image"] > 9))]  # unbalanced wells
    res = models.analyse(d, POS, _settings(), positions_ok=True)
    assert res.ok
    assert _row(res, "vs LyoVec", "A vs LyoVec").N == 3
    assert _row(res, "vs Control", "C vs Control").N == 2
    # a comparison present in only one n is reported but not tested
    d1 = d[~((d["treatment"] == "C") & (d["n"] != "n1"))]
    res = models.analyse(d1, POS, _settings(), positions_ok=True, leave_one_out=False)
    r = _row(res, "vs LyoVec", "C vs LyoVec")
    assert np.isnan(r.p) and "not estimable" in r.note


def test_one_n_and_one_well_do_not_crash():
    rng = np.random.default_rng(1)
    one_n = simulate(rng, n_n=1)
    res = models.analyse(one_n, POS, _settings(), positions_ok=True)
    assert res.status == "not_estimable" and "one n" in res.model and res.descriptive and res.has_estimates
    one_well = simulate(rng, wells=1)
    res = models.analyse(one_well, POS, _settings(), positions_ok=True, leave_one_out=False)
    assert res.ok
    w = res.varcomp.set_index("component").loc["well (residual: well-to-well within n × treatment)"]
    assert np.isnan(w["variance"]) and "not identifiable" in w["note"]
    res = models.analyse(one_well, POS, _settings(model_unit="field"), positions_ok=True, leave_one_out=False)
    w = res.varcomp.set_index("component").loc["well (within n × treatment)"]
    assert np.isnan(w["variance"]) and "not estimable" in w["note"]


def test_fraction_and_count_scales():
    y, sc = models.transform(pd.Series([0.0, 0.5, 1.0]), "fraction", pd.Series([10, 10, 10]))
    assert np.all(np.isfinite(y)) and sc.estimate_type == "odds_ratio"
    assert sc.back(0.0) == pytest.approx(0.5)
    y, sc = models.transform(pd.Series([0, 9.5]), "count")
    assert y.iloc[1] == pytest.approx(np.log(10)) and sc.back(np.log(10)) == pytest.approx(9.5)
    y, sc = models.transform(pd.Series([0.0, 2.0, 4.0]), "positive")
    assert np.all(np.isfinite(y)) and "zero" in sc.note


# ---------------------------------------------------------------- metrics, settings, acceptance

def test_metric_registry():
    reg = metrics.registry(["n", "count", "intensity_DAPI_nucleus", "nuclei_per_mm2", "readout_frac_dead",
                            "readout_mystery", "qc_level"], {})
    assert list(reg) == ["count", "intensity_DAPI_nucleus", "nuclei_per_mm2", "readout_frac_dead", "readout_mystery"]
    assert reg["count"].kind == "count" and reg["intensity_DAPI_nucleus"].kind == "positive"
    assert reg["readout_frac_dead"].kind == "fraction" and reg["readout_mystery"].kind == "positive"
    legacy = metrics.normalise({"intensity_C1_cell": "C1 cell mean intensity"})
    assert legacy["intensity_C1_cell"].label == "C1 cell mean intensity" and "count" in legacy


def test_settings_form_new_fields():
    for sec in ("analysis", "acceptance"):
        assert set(settings_form.FIELDS[sec]) == set(st_mod.DEFAULTS[sec])
    assert list(st_mod.DEFAULTS).index("acceptance") == list(st_mod.DEFAULTS).index("analysis") + 1
    assert list(settings_form.FIELDS).index("acceptance") == list(settings_form.FIELDS).index("analysis") + 1
    typed = settings_form.coerce({"analysis": dict(classes='{ mirna = ["miR-*"] }', primary_test="rm_anova",
                                                   among_classes="mirna"),
                                  "acceptance": dict(vehicle_within="", positive_control_max_ratio="0.8", enforce=True)})
    assert typed["analysis"]["classes"] == {"mirna": ["miR-*"]}
    assert typed["acceptance"] == dict(vehicle_within="", positive_control_max_ratio=0.8, enforce=True)
    text = settings_form.update_text(st_mod.TEMPLATE.format(raw_dir="x"), typed)
    assert 'classes = { mirna = ["miR-*"] }' in text
    with pytest.raises(ValueError):
        settings_form.coerce({"analysis": dict(classes="{ mirna = 3 }")})


def test_acceptance_criteria_and_enforce():
    rows = []
    for n, lyo in (("n1", 60), ("n2", 100)):
        for t, v in (("Control", 60), ("LyoVec", lyo), ("CL075", 40)):
            for i in range(1, 7):
                rows.append(dict(n=n, treatment=t, image=i, count=v, excluded=False, qc_level="high", source="auto"))
    used = pd.DataFrame(rows)
    s = _settings(reference="Control", positive_control="CL075")
    s["acceptance"].update(vehicle_within=0.3, positive_control_max_ratio=0.8, enforce=True)
    acc = acceptance.evaluate(used, s)
    assert acceptance.failing(acc) == ["n2"]
    v = acc[(acc["n"] == "n2") & (acc["criterion"] == "vehicle_within")].iloc[0]
    assert v["result"] == "fail" and v["value"] == pytest.approx(100 / 60)
    assert acceptance.summary(acc, s)["excluded"] == ["n2"]


# ---------------------------------------------------------------- generation contracts

def _gen_table(rng, ns=("n1", "n2", "n3"), drop=(), wells_fields=12):
    rows = []
    for k, n in enumerate(ns):
        for t, base in (("Control", 60), ("LyoVec", 62), ("CL075", 50), ("miR-1", 25), ("miR-2", 30)):
            if (n, t) in drop:
                continue
            for i in range(1, wells_fields + 1):
                rows.append(dict(n=n, treatment=t, image=i, count=int(max(0, rng.normal(base * (1 + 0.1 * k), 8))),
                                 source="auto", excluded=False))
    return pd.DataFrame(rows)


def _build(tmp_path, tab, **analysis):
    s = st_mod._merge(st_mod.DEFAULTS, {"analysis": dict(reference="LyoVec", normalise_to=["Control", "LyoVec"]) | analysis})
    img = tab[["n", "treatment", "image"]].assign(qc_level="high", intensity_C1_nucleus=lambda d: 1000 + d["image"],
                                                  readout_frac_x=0.2)
    gen = dict(name="g1", table=tab, source="auto", changes=0, fingerprint="x")
    reg = metrics.registry(img.columns, s)
    return generations.build(tmp_path, gen, img, {k: m for k, m in reg.items() if k != "count"}, s), s


def test_manifest_and_stats_long_contract(tmp_path):
    rng = np.random.default_rng(0)
    out, _ = _build(tmp_path, _gen_table(rng, drop={("n2", "LyoVec")}), positive_control="CL075",
                    classes={"mirna": ["miR-*"]})
    man = json.loads((out / "manifest.json").read_text())
    listed = {f["path"] for f in man["files"]}
    on_disk = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.suffix in (".png", ".csv")}
    assert on_disk <= listed  # every plot and table is registered
    for f in man["files"]:
        assert (out / f["path"]).exists()
        assert f["level"] in manifest.LEVELS and f["kind"] in manifest.KINDS
        assert f["title"] and f["description"]
        assert set(f) >= {"path", "level", "n", "treatment", "metric", "kind", "title", "description"}
        if f["level"] in ("n", "treatment"):
            assert f["n"]
        if f.get("pdf"):
            assert (out / f["pdf"]).exists()
    assert [m["id"] for m in man["metrics"]] == ["count", "intensity_C1_nucleus", "readout_frac_x"]
    assert man["primary"]["test"] == "mixed" and man["primary"]["reference"] == "LyoVec"
    assert set(man["model_status"]) == {"count", "intensity_C1_nucleus", "readout_frac_x"}
    assert man["acceptance"]["failing"] == [] or isinstance(man["acceptance"]["failing"], list)
    kinds = {f["kind"] for f in man["files"]}
    assert {"absolute", "relative", "model", "contrast", "robustness", "diagnostic", "design", "acceptance", "data"} <= kinds
    sl = pd.read_csv(out / "stats_long.csv")
    assert list(sl.columns) == manifest.STATS_LONG_COLUMNS
    assert set(sl["level"]) == {"n", "all_n"}
    assert set(sl["estimate_type"].dropna()) <= set(manifest.ESTIMATE_TYPES)
    assert sl.loc[sl["level"] == "all_n", "estimand"].notna().any()
    frac = sl[(sl["metric"] == "readout_frac_x") & (sl["model"] != "classical") & (sl["level"] == "all_n")]
    assert (frac["estimate_type"] == "pp_difference").any()  # the headline for fractions
    assert (sl["family"] == "vehicle check (equivalence)").any()
    assert (sl.loc[(sl["metric"] == "readout_frac_x") & (sl["model"] != "classical"), "estimate_type"] == "odds_ratio").any()
    xl = pd.ExcelFile(out / "results.xlsx").sheet_names
    for sheet in ("About", "Counts used", "Count - all n", "Model - means (absolute)", "Model - ratios vs controls",
                  "Comparison families", "Variance components", "Robustness - leave-one-out", "Acceptance"):
        assert sheet in xl
    assert (out / "all_n" / "count_forest_vs_LyoVec.png").exists()
    assert (out / "all_n" / "count_prism_pct_of_Control.csv").exists()  # existing outputs kept


def test_generation_survives_one_n_and_rm_anova_primary(tmp_path):
    rng = np.random.default_rng(1)
    out, _ = _build(tmp_path, _gen_table(rng, ns=("n1",)))
    man = json.loads((out / "manifest.json").read_text())
    assert man["model_status"]["count"]["status"] == "not_estimable"
    out, _ = _build(tmp_path, _gen_table(rng), primary_test="rm_anova")
    man = json.loads((out / "manifest.json").read_text())
    assert man["primary"]["test"] == "rm_anova" and man["model_status"]["count"]["stars_from"] == "rm_anova"


def test_generation_enforced_acceptance_excludes_n(tmp_path):
    rng = np.random.default_rng(2)
    tab = _gen_table(rng)
    tab.loc[(tab["n"] == "n3") & (tab["treatment"] == "LyoVec"), "count"] = 200
    s_over = dict(reference="Control")
    out, _ = _build(tmp_path, tab, **s_over)
    acc = pd.read_csv(out / "acceptance.csv")
    assert acc.loc[(acc["n"] == "n3") & (acc["criterion"] == "overall"), "result"].iloc[0] == "fail"
    import tomllib  # noqa: F401  (settings come from DEFAULTS here; enforce via a direct build)
    s = st_mod._merge(st_mod.DEFAULTS, {"analysis": dict(reference="Control"), "acceptance": dict(enforce=True)})
    img = tab[["n", "treatment", "image"]].assign(qc_level="high")
    out = generations.build(tmp_path, dict(name="g1", table=tab, source="auto", changes=0, fingerprint="x"), img, {}, s)
    man = json.loads((out / "manifest.json").read_text())
    assert man["acceptance"]["excluded"] == ["n3"]
    assert "n3" not in set(pd.read_csv(out / "all_n" / "count_per_n_ratios.csv")["n"])
