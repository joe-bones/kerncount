"""Well-unit estimands, conservative default classes, fraction effects in percentage points, one n, cell lines,
negative control, density adjustment and duplicate raw images."""
import json

import numpy as np
import pandas as pd
import pytest

from kerncount import acceptance, duplicates, generations, manifest, metrics, models
from kerncount import settings as st_mod

TR = ["Control", "LyoVec", "CL075", "miR-1", "miR-2", "miR-3"]
CTRL = ("Control", "LyoVec")
COUNT = metrics.COUNT
FRAC = metrics.Metric("readout_frac_x", "readout_frac_x", "x fraction", "fraction", "fraction")


def _settings(**analysis):
    a = dict(reference="LyoVec", normalise_to=["Control", "LyoVec"], treatment_order=TR, vehicle_check="") | analysis
    return st_mod._merge(st_mod.DEFAULTS, {"analysis": a})


def sim_counts(rng, n_n=4, effects=None, var_nt=None, field_sd=None, var_n=0.05, var_well=0.02, pos=None,
               base=np.log(70), treatments=TR):
    """Poisson-lognormal field counts. effects: {t: ratio of ARITHMETIC means}; field_sd: {t: field log-SD} (the
    lognormal mean correction keeps the arithmetic mean fixed whatever the field dispersion)."""
    effects, var_nt, field_sd = effects or {}, var_nt or {}, field_sd or {}
    pos = np.zeros(6) if pos is None else np.asarray(pos)
    rows = []
    for i in range(n_n):
        n = f"n{i + 1}"
        bn = rng.normal(0, np.sqrt(var_n))
        for t in treatments:
            bnt = rng.normal(0, np.sqrt(var_nt[t])) if var_nt.get(t, 0) > 0 else 0.0
            sd = field_sd.get(t, 0.4)
            for w in (1, 2):
                bw = rng.normal(0, np.sqrt(var_well))
                for k in range(6):
                    mu = base + np.log(effects.get(t, 1.0)) + bn + bnt + bw + pos[k] + rng.normal(0, sd) - sd ** 2 / 2
                    rows.append(dict(n=n, treatment=t, well=w, position=k + 1, image=(w - 1) * 6 + k + 1,
                                     timepoint="", value=float(rng.poisson(np.exp(mu)))))
    d = pd.DataFrame(rows)
    return d.assign(count=d["value"])


def _row(c, family, comparison):
    return c[(c["family"] == family) & (c["comparison"] == comparison)].iloc[0]


# ---------------------------------------------------------------- 1. estimand: ratio of mean nuclei per field

def test_well_unit_recovers_ratio_of_arithmetic_means_with_uneven_fields():
    """miRNA fields far more uneven than control fields: the field-level model compares geometric means (biased
    low), the well unit recovers the ratio of mean nuclei per field."""
    rng = np.random.default_rng(1)
    d = sim_counts(rng, n_n=6, effects={"miR-1": 0.5}, field_sd={"miR-1": 1.1}, var_n=0.02, var_well=0.005)
    s = _settings()
    well = models.analyse(d, COUNT, s, positions_ok=True, leave_one_out=False, density=False)
    field = models.analyse(d, COUNT, _settings(model_unit="field"), positions_ok=True, leave_one_out=False,
                           density=False)
    rw = _row(well.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    rf = _row(field.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    m = d.groupby(["n", "treatment"])["value"].mean().unstack()
    empirical = float(np.exp(np.log(m["miR-1"] / m["LyoVec"]).mean()))
    # well unit = ratio of mean nuclei per field (log of a mean of 6 very uneven fields is a little low: ~4%)
    assert abs(np.log(rw.estimate / empirical)) < 0.1
    assert rw.ci_low < 0.5 < rw.ci_high
    assert rf.estimate < 0.8 * rw.estimate  # geometric means of uneven fields are dragged down (~exp(-(1.1²-0.4²)/2))
    assert rw.estimand == "ratio of mean nuclei per field"
    e = well.emm.set_index("treatment")
    assert e.loc["LyoVec", "mean_type"].startswith("mean per field")


def test_position_adjusted_well_mean_with_missing_fields():
    """Strong position pattern; the far positions of the treated wells are missing. A plain mean of the fields
    present is biased; the position-adjusted mean (ratio estimator) is not."""
    rng = np.random.default_rng(2)
    pos = np.log([0.5, 0.7, 0.9, 1.1, 1.4, 1.8])
    d = sim_counts(rng, n_n=4, effects={"miR-1": 0.6}, pos=pos, var_n=0.01, var_well=0.0005,
                   field_sd={t: 0.1 for t in TR})
    d = d[~((d["treatment"] == "miR-1") & (d["position"] >= 5))]
    w, notes = models.well_means(d, "value", positions_ok=True)
    assert w.loc[w["treatment"] == "miR-1", "adjusted"].all() and not w.loc[w["treatment"] == "LyoVec", "adjusted"].any()
    assert any("position-adjusted" in x for x in notes)
    res = models.analyse(d, COUNT, _settings(), positions_ok=True, leave_one_out=False, density=False)
    r = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    true_ratio = 0.6
    plain = d.groupby(["n", "treatment"])["value"].mean().unstack()
    plain_ratio = float(np.exp(np.log(plain["miR-1"] / plain["LyoVec"]).mean()))
    assert abs(np.log(plain_ratio / true_ratio)) > 0.2  # the naive mean is clearly biased (low positions only)
    assert abs(np.log(r.estimate / true_ratio)) < 0.08
    res_r = models.analyse(d, COUNT, _settings(), positions_ok=False, leave_one_out=False, density=False)
    assert abs(np.log(_row(res_r.contrasts, "vs LyoVec", "miR-1 vs LyoVec").estimate / plain_ratio)) < 0.05


# ---------------------------------------------------------------- 2. conservative classes + type I error

def test_auto_classes_one_per_treatment_each_control_included():
    cfg = models.config(_settings(), TR)
    assert cfg.classes_mode == "auto"
    assert cfg.classes["Control"] == "Control" and cfg.classes["LyoVec"] == "LyoVec"
    assert {cfg.classes[t] for t in ("CL075", "miR-1", "miR-2")} == {"CL075", "miR-1", "miR-2"}
    assert models.config(_settings(classes={}), TR).classes["miR-1"] == "treatments"
    assert models.config(_settings(classes='{ mirna = ["miR-*"] }'), TR).classes["miR-2"] == "mirna"


def test_vehicle_comparisons_rest_on_the_vehicles_own_variance():
    """C1: the vehicle (LyoVec) varies between n's, Control does not. With the controls pooled in one class the
    vehicle's n-to-n variance was diluted by Control's zero and "vs LyoVec" got too many df / too small SEs."""
    rng = np.random.default_rng(31)
    d = sim_counts(rng, n_n=4, var_nt={"LyoVec": 0.3})
    res = models.analyse(d, COUNT, _settings(), positions_ok=True, leave_one_out=False, density=False)
    assert res.ok and "class-specific" in res.model
    vc = res.varcomp.set_index("component")["variance"]
    assert "n × treatment [LyoVec]" in vc.index and "n × treatment [Control]" in vc.index
    assert "n × treatment [controls]" not in vc.index
    r = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    assert r.df <= 6  # (4 - 1) for miR-1 + (4 - 1) for LyoVec


@pytest.mark.parametrize("scenario", ["misspecified", "uneven_fields", "vehicle_varies"])
@pytest.mark.slow
def test_type_one_error_default_including_misspecification(scenario):
    """Quick null simulation under the default (well unit, classes = auto). misspecified: one miRNA varies much
    more between n's than the others (pooled classes gave 0.16 for it). uneven_fields: miRNA fields are far more
    dispersed, equal arithmetic means (the field-level model would see a 'difference'). vehicle_varies: LyoVec varies
    between n's, Control does not (pooled controls gave 0.13 for "vs LyoVec"). Full tables:
    docs/statistics.md."""
    rng = np.random.default_rng({"misspecified": 21, "uneven_fields": 22, "vehicle_varies": 23}[scenario])
    ps, p1 = [], []
    for _ in range(30):
        if scenario == "misspecified":
            d = sim_counts(rng, n_n=4, var_nt={"miR-1": 0.3, "miR-2": 0.02, "miR-3": 0.02, "CL075": 0.02})
        elif scenario == "vehicle_varies":  # LyoVec varies between n's, Control and the treatments do not (C1)
            d = sim_counts(rng, n_n=4, var_nt={"LyoVec": 0.3})
        else:
            d = sim_counts(rng, n_n=4, var_nt={t: 0.1 for t in TR if t not in CTRL},
                           field_sd={t: (1.0 if t.startswith("miR") else 0.4) for t in TR})
        res = models.analyse(d, COUNT, _settings(), positions_ok=True, leave_one_out=False, density=False)
        assert res.ok
        c = res.contrasts[res.contrasts["family"] == "vs LyoVec"]
        ps += c["p"].tolist()
        p1.append(float(c.loc[c["group"] == "miR-1", "p"].iloc[0]))
    rate = float(np.mean(np.asarray(ps) < 0.05))
    assert rate <= 0.10, rate  # 150 correlated tests; nominal 0.05 (simulation in the docs: <= 0.06)
    assert float(np.mean(np.asarray(p1) < 0.05)) <= 0.17  # 30 tests of the odd one out: nominal 0.05


# ---------------------------------------------------------------- 3. fractions: percentage points

def sim_fraction(rng, n_n=4, p=None, nuclei=80, var_n=0.03, var_nt=0.02, var_well=0.02, cov=None):
    """Binomial positives per field with logit-normal n / n x treatment / well effects. p: {t: proportion}.
    cov: {t: coverage level} adds a readout_coverage column (field coverage around that level)."""
    p = p or {}
    rows = []
    for i in range(n_n):
        n = f"n{i + 1}"
        bn = rng.normal(0, np.sqrt(var_n))
        for t in TR:
            bnt = rng.normal(0, np.sqrt(var_nt))
            for w in (1, 2):
                bw = rng.normal(0, np.sqrt(var_well))
                lp = np.log(p.get(t, 0.1) / (1 - p.get(t, 0.1))) + bn + bnt + bw
                for k in range(6):
                    N = int(rng.poisson(nuclei))
                    K = rng.binomial(N, 1 / (1 + np.exp(-lp)))
                    r = dict(n=n, treatment=t, well=w, position=k + 1, image=(w - 1) * 6 + k + 1, timepoint="",
                             count=N, value=K / N if N else np.nan)
                    if cov is not None:
                        r["readout_coverage"] = float(cov.get(t, 0.1) * np.exp(rng.normal(0, 0.15)))
                    rows.append(r)
    return pd.DataFrame(rows)


def test_fraction_effect_is_percentage_points_with_ci():
    rng = np.random.default_rng(3)
    d = sim_fraction(rng, n_n=6, p={"miR-1": 0.25, "LyoVec": 0.10, "Control": 0.10})
    res = models.analyse(d, FRAC, _settings(), positions_ok=True, leave_one_out=False, density=False)
    r = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    assert r.estimate_type == "odds_ratio" and r.estimand.startswith("difference in percentage points")
    assert 10 < r.pp_difference < 20 and r.pp_ci_low < 15 < r.pp_ci_high  # true: +15 pp
    assert r.pp_ci_low < r.pp_difference < r.pp_ci_high and r.p_adj < 0.05
    assert r.pct_group == pytest.approx(r.pct_comparator + r.pp_difference)
    assert r.pct_of_comparator == pytest.approx(100 * r.pct_group / r.pct_comparator)
    assert r.estimate == pytest.approx(np.exp(r.log_estimate))  # odds ratio kept as secondary
    assert np.isfinite(r.pp_ci_low_adj) and r.pp_ci_low_adj <= r.pp_ci_low
    e = res.emm.set_index("treatment")
    assert 7 < e.loc["LyoVec", "percent"] < 14 and e.loc["LyoVec", "percent_ci_low"] < e.loc["LyoVec", "percent"]
    # near-zero reference: the pp CI stays in the possible range and is asymmetric (no clipping at a boundary)
    z = sim_fraction(rng, n_n=4, p={"LyoVec": 0.004, "Control": 0.004, "miR-1": 0.05})
    res = models.analyse(z, FRAC, _settings(), positions_ok=True, leave_one_out=False, density=False)
    r = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    assert -100 <= r.pp_ci_low < r.pp_difference < r.pp_ci_high <= 100
    assert r.pp_ci_high - r.pp_difference > r.pp_difference - r.pp_ci_low  # skewed away from 0 %


# ---------------------------------------------------------------- 3b. background-subtracted intensities (I3)

INT = metrics.Metric("intensity_C1_nucleus", "intensity_C1_nucleus", "C1 nucleus", "intensity (a.u.)", "positive")


def sim_bgsub(rng, n_n=4, mean=1.0, sd=3.0, diff=None):
    """Background-subtracted intensity per field: true mean `mean` (+ diff[t]), field SD `sd`, so many fields and
    some well means are <= 0 (log impossible)."""
    diff = diff or {}
    rows = []
    for i in range(n_n):
        bn = rng.normal(0, 0.3)
        for t in TR:
            for w in (1, 2):
                for k in range(6):
                    rows.append(dict(n=f"n{i + 1}", treatment=t, well=w, position=k + 1, image=(w - 1) * 6 + k + 1,
                                     timepoint="", value=mean + diff.get(t, 0.0) + bn + rng.normal(0, sd)))
    return pd.DataFrame(rows)


def test_positive_metric_with_well_means_at_or_below_zero_is_a_raw_difference():
    """I3: true mean 1, SD 3, 4 n. Dropping wells with mean <= 0 (log impossible) biased the ratio; the metric is
    analysed on the raw scale as a difference (identity link), flagged everywhere."""
    rng = np.random.default_rng(8)
    d = sim_bgsub(rng, diff={"miR-1": 2.0})
    assert models.needs_raw_scale(d, INT) > 0
    res = models.analyse(d, INT, _settings(), positions_ok=True, leave_one_out=False, density=False)
    assert res.ok and res.scale.kind == "difference" and res.scale.estimate_type == "difference"
    st = res.status_dict()
    assert st["raw_scale"] and st["effect_type"] == "difference" and "RAW SCALE" in st["raw_scale_warning"]
    assert any("RAW SCALE" in x for x in res.notes)
    r = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    assert r.estimate_type == "difference" and r.estimand.startswith("DIFFERENCE")
    assert r.ci_low < 2.0 < r.ci_high and 0.5 < r.estimate < 3.5  # true difference +2 intensity units
    n0 = _row(res.contrasts, "vs LyoVec", "miR-2 vs LyoVec")
    assert n0.ci_low < 0 < n0.ci_high
    assert res.emm["estimate"].between(-3, 6).all()  # means on the raw scale, not exp() of them
    vc = res.varcomp
    assert vc["cv_pct"].isna().all()  # no CV on the raw scale
    pn = res.per_n
    assert pn["difference"].notna().all() and pn["ratio"].isna().all()
    # all positive: the usual log scale and ratio
    pos = models.analyse(d.assign(value=d["value"] + 30), INT, _settings(), positions_ok=True, leave_one_out=False,
                         density=False)
    assert pos.scale.kind == "positive" and not pos.status_dict()["raw_scale"]


@pytest.mark.slow
def test_raw_scale_metric_flagged_in_generation_outputs(tmp_path):
    rng = np.random.default_rng(9)
    tab = _gen_table(rng, ns=("n1", "n2", "n3", "n4"))
    out, _ = _build(tmp_path, tab, img_extra={"intensity_C1_nucleus": rng.normal(1.0, 3.0, len(tab))})
    man = json.loads((out / "manifest.json").read_text())
    ms = man["model_status"]["intensity_C1_nucleus"]
    assert ms["raw_scale"] and ms["effect_type"] == "difference"
    assert man["estimands"]["intensity_C1_nucleus"]["headline"] == "difference"
    assert man["model_status"]["count"]["raw_scale"] is False
    sl = pd.read_csv(out / "stats_long.csv")
    m = sl[(sl["level"] == "all_n") & (sl["metric"] == "intensity_C1_nucleus") & sl["family"].eq("vs LyoVec")
           & sl["model"].str.contains("n×treatment", regex=False, na=False)]
    assert len(m) and (m["estimate_type"] == "difference").all() and m["note"].str.contains("RAW SCALE").all()
    info = pd.read_csv(out / "all_n" / "intensity_C1_nucleus_model_info.csv")
    assert info.astype(str).apply(lambda c: c.str.contains("RAW SCALE")).any().any()


# ---------------------------------------------------------------- 4. one n

def _gen_table(rng, ns=("n1", "n2", "n3"), extra=()):
    rows = []
    for k, n in enumerate(ns):
        for t, base in (("Control", 60), ("LyoVec", 62), ("CL075", 50), ("miR-1", 25), ("miR-2", 30)) + tuple(extra):
            for i in range(1, 13):
                rows.append(dict(n=n, treatment=t, image=i, count=int(max(0, rng.normal(base * (1 + 0.1 * k), 8))),
                                 source="auto", excluded=False))
    return pd.DataFrame(rows)


def _build(tmp_path, tab, img_extra=None, settings_over=None, **analysis):
    s = st_mod._merge(st_mod.DEFAULTS, {"analysis": dict(reference="LyoVec", normalise_to=["Control", "LyoVec"]) | analysis})
    if settings_over:
        s = st_mod._merge(s, settings_over)
    rng = np.random.default_rng(5)
    img = tab[["n", "treatment", "image"]].assign(qc_level="high", intensity_C1_nucleus=lambda d: 1000 + d["image"],
                                                  readout_frac_x=rng.uniform(0.05, 0.3, len(tab)),
                                                  readout_coverage=rng.uniform(0.05, 0.12, len(tab)))
    for k, v in (img_extra or {}).items():
        img[k] = v
    gen = dict(name="g1", table=tab, source="auto", changes=0, fingerprint="x")
    reg = metrics.registry(img.columns, s)
    return generations.build(tmp_path, gen, img, {k: m for k, m in reg.items() if k != "count"}, s), s


def _shape(out):
    """Every all-n file with its CSV columns (the 'layout' of the all-n level); the one-column-per-n tables
    (summary, Prism) naturally have one column per n and are compared without those columns."""
    files = {}
    ns = set(pd.read_csv(out / "counts_used.csv")["n"])
    for p in sorted((out / "all_n").rglob("*")):
        if p.suffix == ".csv":
            files[p.relative_to(out).as_posix()] = [c for c in pd.read_csv(p).columns if c not in ns]
        elif p.suffix == ".png":
            files[p.relative_to(out).as_posix()] = None
    return files


@pytest.mark.slow
def test_one_n_keeps_all_n_files_columns_and_has_no_stars(tmp_path):
    rng = np.random.default_rng(1)
    many, _ = _build(tmp_path / "a", _gen_table(rng))
    one, _ = _build(tmp_path / "b", _gen_table(rng, ns=("n1",)))
    a, b = _shape(many), _shape(one)
    assert set(a) == set(b)  # same files
    assert all(a[k] == b[k] for k in a)  # same columns
    sl = pd.read_csv(one / "stats_long.csv")
    alln = sl[sl["level"] == "all_n"]
    assert alln["significance"].fillna("").isin(["", "inconclusive"]).all()  # no stars anywhere
    assert alln["p"].isna().all() and alln["ci_low"].isna().all()
    model = alln[alln["model"] != "classical"]
    assert model["note"].str.contains("needs ≥ 2 n").any()
    means = alln[(alln["family"] == "estimated means") & (alln["metric"] == "count")]
    assert means["estimate"].notna().all() and (means["N"] == 1).all()  # estimates filled in from the one n
    c = pd.read_csv(one / "all_n" / "count_contrasts.csv")
    assert c["estimate"].notna().all() and c["p"].isna().all() and (c["note"] == "needs ≥ 2 n").all()
    st = pd.read_csv(one / "all_n" / "count_stats.csv")
    assert st["p"].isna().all() and st["note"].str.contains("needs ≥ 2 n").all()
    man = json.loads((one / "manifest.json").read_text())
    assert man["model_status"]["count"]["across_n_inference"] is False
    assert any("only one n" in x for x in man["notes"])
    # per-n level: technical intervals, labelled descriptive
    summ = pd.read_csv(one / "n1" / "count_summary.csv")
    assert summ["ci_low"].notna().all() and summ["interval"].str.startswith("descriptive: technical").all()
    pn = pd.read_csv(one / "all_n" / "count_per_n_ratios.csv")
    assert pn["tech_ci_low"].notna().all() and (pn["interval"] == models.TECH).all()


@pytest.mark.slow
def test_comparison_present_in_one_n_only_gets_no_stars(tmp_path):
    rng = np.random.default_rng(2)
    tab = _gen_table(rng)
    tab = tab[~((tab["treatment"] == "miR-2") & (tab["n"] != "n1"))]
    out, _ = _build(tmp_path, tab)
    c = pd.read_csv(out / "all_n" / "count_contrasts.csv")
    r = c[c["comparison"] == "miR-2 vs LyoVec"].iloc[0]
    assert np.isfinite(r["estimate"]) and np.isnan(r["p"]) and not isinstance(r["significance"], str)
    assert "needs ≥ 2 n" in r["note"]


# ---------------------------------------------------------------- 5. cell lines, 6. negative control

@pytest.mark.slow
def test_cell_lines_are_never_pooled(tmp_path):
    rng = np.random.default_rng(3)
    tab = pd.concat([_gen_table(rng, ns=("n1", "n2")), _gen_table(rng, ns=("n3", "n4"))])
    line = np.where(tab["n"].isin(["n1", "n2"]), "U251", "LN229")
    out, _ = _build(tmp_path, tab, img_extra={"cell_line": line})
    assert (out / "all_n" / "U251" / "count_contrasts.csv").exists()
    assert (out / "all_n" / "LN229" / "count_contrasts.csv").exists()
    assert not (out / "all_n" / "count_contrasts.csv").exists()
    pn = pd.read_csv(out / "all_n" / "U251" / "count_per_n_ratios.csv")
    assert set(pn["n"]) == {"n1", "n2"}
    man = json.loads((out / "manifest.json").read_text())
    assert man["cell_lines"] == ["LN229", "U251"]
    assert {f["path"] for f in man["all_n_folders"]} == {"all_n/LN229", "all_n/U251"}
    assert set(man["model_status"]) >= {"U251/count", "LN229/count"}
    e = next(f for f in man["files"] if f["path"] == "all_n/U251/count_contrasts.csv")
    assert e["cell_line"] == "U251" and e["level"] == "all_n" and "(U251)" in e["title"]
    assert next(f for f in man["files"] if f["path"] == "n3/count_summary.csv")["cell_line"] == "LN229"
    sl = pd.read_csv(out / "stats_long.csv")
    assert set(sl.loc[sl["level"] == "all_n", "cell_line"]) == {"U251", "LN229"}
    used = pd.read_csv(out / "counts_used.csv")
    assert "cell_line" in used


@pytest.mark.slow
def test_one_cell_line_plus_unlabelled_ns_is_labelled_so(tmp_path):
    """M6: U251 n's plus n's without a cell line form one block; it must not be called plain "U251"."""
    rng = np.random.default_rng(4)
    tab = _gen_table(rng)
    out, _ = _build(tmp_path, tab, img_extra={"cell_line": np.where(tab["n"] == "n3", "", "U251")})
    man = json.loads((out / "manifest.json").read_text())
    assert man["all_n_folders"][0]["cell_line"] == "U251 + unlabelled"
    assert any("n3 have no cell line set" in x for x in man["notes"])
    sl = pd.read_csv(out / "stats_long.csv")
    assert set(sl.loc[sl["level"] == "all_n", "cell_line"]) == {"U251 + unlabelled"}


@pytest.mark.slow
def test_negative_control_is_excluded_everywhere(tmp_path):
    rng = np.random.default_rng(4)
    tab = _gen_table(rng, extra=(("2nd only", 5),))
    out, _ = _build(tmp_path, tab, settings_over={"experiment": {"negative_control": "2nd only"}})
    sl = pd.read_csv(out / "stats_long.csv")
    assert not sl["comparison"].astype(str).str.contains("2nd only").any()
    assert not (out / "n1" / "2nd_only").exists()
    assert "2nd only" not in set(pd.read_csv(out / "all_n" / "count_summary.csv")["treatment"])
    assert "2nd only" not in pd.read_csv(out / "all_n" / "count_prism.csv").columns
    used = pd.read_csv(out / "counts_used.csv")
    neg = used[used["treatment"] == "2nd only"]
    assert len(neg) and neg["excluded"].all() and neg["source"].str.contains("negative control").all()
    man = json.loads((out / "manifest.json").read_text())
    assert man["negative_control"] == {"treatment": "2nd only", "present": True}
    assert any("negative control" in x for x in man["notes"])
    about = pd.read_excel(out / "results.xlsx", sheet_name="About", header=None)[0].astype(str)
    assert about.str.contains("negative control '2nd only'").any()


# ---------------------------------------------------------------- 7. density adjustment

def test_density_adjustment_and_overlap_flag():
    rng = np.random.default_rng(6)
    # out-of-focus fraction driven only by crowding; miR-1 wells are much sparser (no density overlap)
    rows = []
    for i in range(5):
        for t in TR:
            for w in (1, 2):
                for k in range(6):
                    cov = (0.03 if t == "miR-1" else 0.10) * np.exp(rng.normal(0, 0.1))
                    if t == "miR-2":
                        cov = rng.uniform(0.07, 0.13)
                    p = 1 / (1 + np.exp(-(-1.0 - 25 * cov)))  # sparser -> more out of focus
                    N = 80
                    rows.append(dict(n=f"n{i + 1}", treatment=t, well=w, position=k + 1, image=(w - 1) * 6 + k + 1,
                                     timepoint="", count=N, value=rng.binomial(N, p) / N, readout_coverage=cov))
    d = pd.DataFrame(rows)
    res = models.analyse(d, FRAC, _settings(), positions_ok=True, leave_one_out=False)
    raw = _row(res.contrasts, "vs LyoVec", "miR-1 vs LyoVec")
    assert raw.pp_difference > 5  # crowding difference shows up as an effect in the raw estimate
    dens = res.density
    assert dens is not None and dens.covariate == "readout_coverage"
    adj = dens.contrasts.set_index("comparison")
    assert adj.loc["miR-1 vs LyoVec", "density_flag"] == models.DENSITY_FLAG
    assert adj.loc["miR-1 vs LyoVec", "density_overlap"] < models.OVERLAP_MIN
    assert adj.loc["miR-2 vs LyoVec", "density_flag"] == ""  # overlapping densities: no flag
    assert abs(adj.loc["miR-2 vs LyoVec", "pp_difference"]) < 3
    assert dens.slope["slope"] < 0 and dens.slope["ci_high"] < 0
    pts, curves = models.density_curves(dens)
    assert len(pts) == len(dens.prepared.d) and set(curves["treatment"]) == set(TR)
    off = models.analyse(d, FRAC, _settings(density_covariate=""), positions_ok=True, leave_one_out=False)
    assert off.density is None


@pytest.mark.slow
def test_density_outputs_in_generation(tmp_path):
    rng = np.random.default_rng(7)
    out, _ = _build(tmp_path, _gen_table(rng))
    for f in ("readout_frac_x_density_adjusted_means.csv", "readout_frac_x_density_adjusted_contrasts.csv",
              "readout_frac_x_density.png"):
        assert (out / "all_n" / f).exists(), f
    assert (out / "n1" / "readout_frac_x_density.png").exists()
    ct = pd.read_csv(out / "all_n" / "readout_frac_x_density_adjusted_contrasts.csv")
    assert {"density_overlap", "density_flag", "raw_pp_difference", "slope_slope"} <= set(ct.columns)
    man = json.loads((out / "manifest.json").read_text())
    kinds = {f["path"]: f["kind"] for f in man["files"]}
    assert kinds["all_n/readout_frac_x_density.png"] == "model"
    sl = pd.read_csv(out / "stats_long.csv")
    assert sl["family"].str.endswith("(density-adjusted)").any()


# ---------------------------------------------------------------- 8. duplicate raw images

def test_duplicate_raw_images_are_found_cached_and_fail_acceptance(tmp_path):
    raw = tmp_path / "raw"
    for folder in ("ctrl", "lyo", "mir"):
        (raw / folder).mkdir(parents=True)
    for i in (1, 2):
        (raw / "ctrl" / f"{i}.tif").write_bytes(b"image %d" % i)
        (raw / "lyo" / f"{i}.tif").write_bytes(b"image %d" % i)  # copied folder
        (raw / "mir" / f"{i}.tif").write_bytes(b"other %d" % i)
    table = pd.DataFrame([dict(n="n3", treatment=t, image=i, source_file=str(raw / f / f"{i}.tif"))
                          for t, f in (("Control", "ctrl"), ("LyoVec", "lyo"), ("miR-1", "mir")) for i in (1, 2)])
    logged = []
    dup = duplicates.check(tmp_path, table, log=logged.append)
    assert len(dup) == 4 and dup["cross_condition"].all()
    assert set(dup["treatment"]) == {"Control", "LyoVec"}
    assert (tmp_path / "duplicates.csv").exists() and (tmp_path / "_work" / "file_hashes.json").exists()
    assert any("WARNING" in x for x in logged)
    assert duplicates.summary_lines(dup) == [
        "n3 / Control images 1-2 = n3 / LyoVec images 1-2 (2 identical file(s); different conditions)"]
    # cached by size + mtime: a second run does not read the files again
    cache = json.loads((tmp_path / "_work" / "file_hashes.json").read_text())
    k = str(raw / "ctrl" / "1.tif")
    cache[k]["sha256"] = "pretend"
    (tmp_path / "_work" / "file_hashes.json").write_text(json.dumps(cache))
    assert duplicates.file_hashes([k], tmp_path / "_work" / "file_hashes.json")[k] == "pretend"
    # acceptance: the n with images shared across conditions fails
    used = table.assign(count=50, excluded=False, qc_level="high")
    s = _settings()
    acc = acceptance.evaluate(used, s, dup)
    r = acc[acc["criterion"] == "no_duplicate_images"].iloc[0]
    assert r["result"] == "fail" and r["value"] == 4 and "LyoVec" in r["detail"]
    assert acceptance.failing(acc) == ["n3"]
    clean = acceptance.evaluate(used, s, duplicates.find(table, {}))
    assert (clean.loc[clean["criterion"] == "no_duplicate_images", "result"] == "pass").all()
    nochk = acceptance.evaluate(used, s, None)
    assert (nochk.loc[nochk["criterion"] == "no_duplicate_images", "result"] == "not applicable").all()


def test_duplicates_in_generation_about_and_manifest(tmp_path):
    rng = np.random.default_rng(8)
    tab = _gen_table(rng)
    dup = pd.DataFrame([dict(group=1, n="n3", treatment=t, image=1, source_file=f"/x/{t}/1.tif", hash="h", copies=2,
                             cross_condition=True, shared_with=f"n3 / {o} / 1")
                        for t, o in (("Control", "LyoVec"), ("LyoVec", "Control"))], columns=duplicates.COLUMNS)
    s = st_mod._merge(st_mod.DEFAULTS, {"analysis": dict(reference="LyoVec")})
    gen = dict(name="g1", table=tab, source="auto", changes=0, fingerprint="x")
    out = generations.build(tmp_path, gen, tab[["n", "treatment", "image"]].assign(qc_level="high"), {}, s,
                            duplicates=dup)
    man = json.loads((out / "manifest.json").read_text())
    assert man["duplicates"]["cross_condition"] == 2 and man["duplicates"]["affected"] == {"n3": ["Control", "LyoVec"]}
    assert any("DUPLICATE RAW IMAGES" in x for x in man["notes"])
    assert any(f["path"] == "duplicates.csv" for f in man["files"])
    about = pd.read_excel(out / "results.xlsx", sheet_name="About", header=None)[0].astype(str)
    assert about.str.contains("DUPLICATE RAW IMAGES").any()
    acc = pd.read_csv(out / "acceptance.csv")
    assert acc.loc[(acc["n"] == "n3") & (acc["criterion"] == "no_duplicate_images"), "result"].iloc[0] == "fail"


# ---------------------------------------------------------------- contracts

def test_manifest_lists_every_file_with_new_outputs(tmp_path):
    rng = np.random.default_rng(9)
    out, _ = _build(tmp_path, _gen_table(rng))
    man = json.loads((out / "manifest.json").read_text())
    listed = {f["path"] for f in man["files"]}
    on_disk = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.suffix in (".png", ".csv", ".xlsx")}
    assert on_disk == listed
    for f in man["files"]:
        assert f["level"] in manifest.LEVELS and f["kind"] in manifest.KINDS and f["description"]
        assert "cell_line" in f
    assert man["model_unit"] == "well" and man["primary"]["classes_mode"] == "auto"
    assert man["estimands"]["count"]["effect"] == "ratio of mean nuclei per field"
    assert man["estimands"]["readout_frac_x"]["headline"] == "pp_difference"
    sl = pd.read_csv(out / "stats_long.csv")
    assert list(sl.columns) == manifest.STATS_LONG_COLUMNS
