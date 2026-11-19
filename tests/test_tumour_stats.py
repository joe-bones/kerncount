"""Tumour factorial statistics (ticket 04, ADR 0003 amendment): contrasts from a declarative family spec, one N-way
log-scale model with pooled variance, Holm within family, confounding with imaging session, missing levels.

All tests use synthetic readout tables (no images)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kerncount.tumour.config import Factor, Family

GEN = Factor("genotype", ("WT", "TLR7 KO"), "WT", "row", "session")
LINE = Factor("tumour_line", ("cntrl", "HERV", "N4"), "cntrl", "x")
CLOD = Factor("clodronate", ("-clodro", "+clodro"), "-clodro", "panel")
FACS = [GEN, LINE, CLOD]

F1 = Family("F1", "simple", ("tumour_line",), "confirmatory")
F2 = Family("F2", "simple", ("clodronate",), "confirmatory")
F3 = Family("F3", "interaction", ("tumour_line", "clodronate"), "exploratory")
F4 = Family("F4", "interaction", ("genotype", "tumour_line"), "confirmatory")
F6 = Family("F6", "simple", ("genotype",), "descriptive")
F7 = Family("F7", "interaction", ("genotype", "tumour_line", "clodronate"), "exploratory")


# ---------------------------------------------------------------- seam 1: factors + family spec -> contrasts

def _w(c) -> dict:
    return {k: v for k, v in c.weights.items() if v != 0}


def test_simple_family_compares_each_level_with_the_reference_within_every_cell_of_the_other_factors():
    from kerncount.tumour.families import contrasts

    cs = contrasts(FACS, F1)
    assert len(cs) == 2 * 2 * 2  # (HERV, N4) x genotype x clodronate
    c = next(c for c in cs if c.within == ("TLR7 KO", "+clodro") and c.level == "N4")
    assert _w(c) == {("TLR7 KO", "N4", "+clodro"): 1.0, ("TLR7 KO", "cntrl", "+clodro"): -1.0}
    assert c.comparison == "N4 vs cntrl | TLR7 KO, +clodro"


def test_two_way_interaction_is_a_ratio_of_ratios_of_the_last_factor_within_the_rest():
    from kerncount.tumour.families import contrasts

    cs = contrasts(FACS, F4)
    assert len(cs) == 2 * 2  # (HERV, N4) x clodronate
    c = next(c for c in cs if c.within == ("-clodro",) and c.level == "TLR7 KO, HERV")
    assert _w(c) == {("TLR7 KO", "HERV", "-clodro"): 1.0, ("TLR7 KO", "cntrl", "-clodro"): -1.0,
                     ("WT", "HERV", "-clodro"): -1.0, ("WT", "cntrl", "-clodro"): 1.0}
    assert c.comparison == "(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | -clodro"


def test_three_way_interaction_has_eight_signed_cells_and_no_within():
    from kerncount.tumour.families import contrasts

    cs = contrasts(FACS, F7)
    assert len(cs) == 2
    c = next(c for c in cs if "N4" in c.level)
    w = _w(c)
    assert len(w) == 8 and sum(w.values()) == 0
    assert w[("TLR7 KO", "N4", "+clodro")] == 1.0 and w[("WT", "cntrl", "-clodro")] == -1.0
    assert w[("WT", "N4", "+clodro")] == -1.0
    assert c.within == ()


def test_families_are_not_tied_to_factor_names():
    from kerncount.tumour.families import contrasts

    dose = Factor("dose", ("0", "1", "10"), "0")
    drug = Factor("drug", ("veh", "A"), "veh")
    cs = contrasts([dose, drug], Family("D", "simple", ("dose",), "confirmatory"))
    assert {c.comparison for c in cs} == {"1 vs 0 | veh", "10 vs 0 | veh", "1 vs 0 | A", "10 vs 0 | A"}


# ---------------------------------------------------------------- synthetic readout tables

def _truth_log(cell, line_fx=None, clod_fx=0.0, gxl=None, ko_shift=0.0) -> float:
    """log mean of a cell: base 0 (=1 mm³), line effects, clodronate effect, genotype x line effect, KO session shift."""
    g, line, c = cell
    line_fx = line_fx or {"cntrl": 0.0, "HERV": np.log(2.0), "N4": np.log(0.7)}
    gxl = gxl or {"cntrl": 0.0, "HERV": np.log(1.5), "N4": 0.0}
    mu = line_fx[line] + (clod_fx if c == "+clodro" else 0.0)
    if g == "TLR7 KO":
        mu += gxl[line] + ko_shift
    return mu


def synth_table(rng, n=None, sigma=0.3, genotypes=("WT", "TLR7 KO"), **truth) -> pd.DataFrame:
    """One row per tumour; n(cell) -> group size (default: 12 for -clodro cntrl/HERV, else 6, as in the KO set)."""
    rows = []
    for g in genotypes:
        for line in LINE.levels:
            for c in CLOD.levels:
                k = n((g, line, c)) if callable(n) else (12 if c == "-clodro" and line != "N4" else 6)
                v = np.exp(_truth_log((g, line, c), **truth) + rng.normal(0, sigma, k))
                for i, x in enumerate(v):
                    rows.append({"tumour_id": f"{g}/{line}/{c}/{i}", "imaging_session": f"session {g}",
                                 "genotype": g, "tumour_line": line, "clodronate": c, "volume": x})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- seam 2: table -> model -> estimates and CIs

def test_planted_effects_are_recovered_with_about_95_percent_ci_coverage():
    from kerncount.tumour import stats
    from kerncount.tumour.families import contrasts

    rng = np.random.default_rng(1)
    targets = {  # (family, contrast label) -> true ratio
        ("F1", "HERV vs cntrl | WT, -clodro"): 2.0,
        ("F1", "N4 vs cntrl | TLR7 KO, +clodro"): 0.7,
        ("F4", "(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | +clodro"): 1.5,
        ("F2", "+clodro vs -clodro | TLR7 KO, HERV"): 0.8,
    }
    fams = {"F1": F1, "F2": F2, "F4": F4}
    hits, ests = {k: 0 for k in targets}, {k: [] for k in targets}
    sims = 300
    for _ in range(sims):
        t = synth_table(rng, clod_fx=np.log(0.8))
        m = stats.fit(t, FACS, "volume")
        for (fid, label), true in targets.items():
            c = next(c for c in contrasts(FACS, fams[fid]) if c.comparison == label)
            e = stats.estimate(m, c)
            ests[(fid, label)].append(np.log(e.ratio))
            hits[(fid, label)] += e.ci_low <= true <= e.ci_high
    for k, true in targets.items():
        assert np.mean(ests[k]) == pytest.approx(np.log(true), abs=0.03), k
        assert 0.91 <= hits[k] / sims <= 0.985, (k, hits[k] / sims)


def test_estimate_matches_ratio_of_geometric_means_and_pooled_variance_with_unbalanced_groups():
    from kerncount.tumour import stats
    from kerncount.tumour.families import contrasts

    rng = np.random.default_rng(2)
    t = synth_table(rng, n=lambda cell: {"cntrl": 12, "HERV": 6, "N4": 4}[cell[1]])
    m = stats.fit(t, FACS, "volume")
    assert m.df == len(t) - 12
    c = next(c for c in contrasts(FACS, F1) if c.comparison == "N4 vs cntrl | WT, -clodro")
    e = stats.estimate(m, c)
    lv = np.log(t["volume"])
    sel = lambda line: (t.genotype == "WT") & (t.tumour_line == line) & (t.clodronate == "-clodro")  # noqa: E731
    assert e.ratio == pytest.approx(np.exp(lv[sel("N4")].mean() - lv[sel("cntrl")].mean()))
    resid = lv - t.groupby(["genotype", "tumour_line", "clodronate"])["volume"].transform(lambda v: np.log(v).mean())
    s2 = (resid ** 2).sum() / (len(t) - 12)
    assert e.se == pytest.approx(np.sqrt(s2 * (1 / 4 + 1 / 12)))
    assert e.n == 16
    # p from the pooled t: matches an independent computation
    from scipy import stats as sps
    assert e.p == pytest.approx(2 * sps.t.sf(abs(np.log(e.ratio)) / e.se, len(t) - 12))


def test_model_fits_with_only_the_levels_present_and_marks_contrasts_needing_absent_cells():
    from kerncount.tumour import stats
    from kerncount.tumour.families import contrasts

    t = synth_table(np.random.default_rng(3), genotypes=("TLR7 KO",))
    m = stats.fit(t, FACS, "volume")
    assert m.df == len(t) - 6
    ko = next(c for c in contrasts(FACS, F1) if c.comparison == "HERV vs cntrl | TLR7 KO, -clodro")
    assert np.isfinite(stats.estimate(m, ko).ratio)
    wt = next(c for c in contrasts(FACS, F4))
    e = stats.estimate(m, wt)
    assert not np.isfinite(e.ratio) and e.missing and all(cell[0] == "WT" for cell in e.missing)


# ---------------------------------------------------------------- seam 3: family -> Holm-adjusted, labelled results

from kerncount.tumour.config import Readout  # noqa: E402

VOL = Readout("volume", "volume", "Volume (mm³)", "confirmatory")
SHAPE = Readout("shape", "volume", "Shape", "exploratory")
ALL = [F1, F2, F3, F4, Family("F5", "interaction", ("genotype", "clodronate"), "exploratory"), F6, F7]


def test_holm_is_applied_within_each_family_and_family_membership_is_visible():
    from kerncount.tumour import stats

    t = synth_table(np.random.default_rng(4))
    r = stats.readout_results(t, FACS, ALL, VOL)
    assert set(r["family"]) == {"F1", "F2", "F3", "F4", "F5", "F6", "F7"}
    for fid, g in r[r["family_role"] != "descriptive"].groupby("family"):
        m = len(g)
        assert (g["p_adj"] >= g["p"] - 1e-15).all()
        assert g["p_adj"].min() == pytest.approx(min(1.0, m * g["p"].min())), fid
    assert len(r[r.family == "F1"]) == 8 and len(r[r.family == "F4"]) == 4 and len(r[r.family == "F7"]) == 2


def test_stars_only_on_confirmatory_families_of_confirmatory_readouts():
    from kerncount.tumour import stats

    t = synth_table(np.random.default_rng(5), sigma=0.1)  # strong, clear effects
    conf = stats.readout_results(t, FACS, ALL, VOL)
    f1 = conf[conf.family == "F1"]
    assert (f1["significance"] != "").all() and "****" in set(f1["significance"])
    f3 = conf[conf.family == "F3"]
    assert (f3["significance"] == "").all() and f3["note"].str.contains("exploratory").all()
    assert f3["p_adj"].notna().all()
    expl = stats.readout_results(t, FACS, ALL, SHAPE)
    assert (expl["significance"] == "").all()
    assert expl[expl.family == "F1"]["note"].str.contains("exploratory readout").all()


def test_ratio_contrast_is_unaffected_by_a_session_shift_while_direct_genotype_contrasts_are_flagged():
    from kerncount.tumour import stats

    a = synth_table(np.random.default_rng(6))
    b = synth_table(np.random.default_rng(6), ko_shift=np.log(3.0))  # KO session imaged 3x brighter/bigger
    ra, rb = (stats.readout_results(t, FACS, ALL, VOL).set_index("comparison") for t in (a, b))
    f4 = ra[ra.family == "F4"].index
    np.testing.assert_allclose(ra.loc[f4, "estimate"], rb.loc[f4, "estimate"], rtol=1e-10)
    herv = "(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | -clodro"
    assert ra.loc[herv, "ci_low"] < 1.5 < ra.loc[herv, "ci_high"]
    f6 = rb[rb.family == "F6"]
    assert f6["confounded"].all() and f6["note"].str.contains("confounded with imaging session").all()
    assert (f6["significance"] == "").all() and f6["p"].isna().all()
    assert np.allclose(f6["estimate"] / ra[ra.family == "F6"]["estimate"], 3.0)
    assert not rb[rb.family.isin(["F1", "F2", "F3", "F4", "F5", "F7"])]["confounded"].any()


def test_a_confirmatory_contrast_spanning_sessions_is_never_starred():
    from kerncount.tumour import stats

    t = synth_table(np.random.default_rng(7), sigma=0.05, ko_shift=np.log(3.0))
    direct = Family("G", "simple", ("genotype",), "confirmatory")
    r = stats.readout_results(t, FACS, [direct], VOL)
    assert r["confounded"].all() and (r["significance"] == "").all() and r["p_adj"].isna().all()


def test_with_one_genotype_families_needing_the_other_are_reported_as_needing_its_data():
    from kerncount.tumour import stats

    t = synth_table(np.random.default_rng(8), genotypes=("TLR7 KO",))
    r = stats.readout_results(t, FACS, ALL, VOL)
    f4 = r[r.family == "F4"]
    assert len(f4) == 4 and f4["note"].str.contains("needs WT data").all() and f4["estimate"].isna().all()
    assert r[r.family == "F6"]["note"].str.contains("needs WT data").all()
    f1 = r[r.family == "F1"]
    ko = f1[f1.comparison.str.contains("TLR7 KO")]
    assert ko["estimate"].notna().all() and (ko["significance"] != "").all()
    wt = f1[f1.comparison.str.contains("WT")]
    assert len(wt) == 4 and wt["note"].str.contains("needs WT data").all()
    # Holm counts only the estimable contrasts of the family
    assert ko["p_adj"].min() == pytest.approx(min(1.0, 4 * ko["p"].min()))


# ---------------------------------------------------------------- normalised values

def test_normalised_controls_scatter_around_100_and_normalised_cis_agree_with_the_forest_fold_changes():
    from kerncount.tumour import stats

    t = synth_table(np.random.default_rng(9))
    m = stats.fit(t, FACS, "volume")
    pct = stats.normalised(t, m, LINE)
    ctrl = pct[t.tumour_line == "cntrl"]
    assert ctrl.std() > 5  # not pinned to 100
    for (g, c), v in ctrl.groupby([t.genotype, t.clodronate]):
        assert np.exp(np.log(v).mean()) == pytest.approx(100.0), (g, c)  # model estimate of cntrl = 100 %
    herv = t[(t.genotype == "TLR7 KO") & (t.tumour_line == "HERV") & (t.clodronate == "+clodro")].index
    cntrl = t[(t.genotype == "TLR7 KO") & (t.tumour_line == "cntrl") & (t.clodronate == "+clodro")].index
    assert pct[herv].tolist() == pytest.approx((100 * t.volume[herv] / np.exp(np.log(t.volume[cntrl]).mean())).tolist())
    groups = stats.normalised_groups(m, LINE).set_index(["genotype", "tumour_line", "clodronate"])
    forest = stats.readout_results(t, FACS, [F1], VOL).set_index("comparison")
    g = groups.loc[("TLR7 KO", "HERV", "+clodro")]
    f = forest.loc["HERV vs cntrl | TLR7 KO, +clodro"]
    assert (g["pct"], g["ci_low"], g["ci_high"]) == pytest.approx((100 * f.estimate, 100 * f.ci_low, 100 * f.ci_high))
    assert groups.loc[("TLR7 KO", "cntrl", "+clodro"), "pct"] == pytest.approx(100.0)


# ---------------------------------------------------------------- seam 4: end-to-end analyse on a readout table

def _settings(**tumour) -> dict:
    from kerncount import settings as settings_mod

    return settings_mod._merge(settings_mod.DEFAULTS, {"experiment": {"type": "tumour_volume"}, "tumour": tumour})


def _vol_table(rng, **kw) -> pd.DataFrame:
    t = synth_table(rng, **kw).rename(columns={"volume": "main_mass_volume_mm3"})
    t["main_mass_volume_um3"] = t["main_mass_volume_mm3"] * 1e9
    return t


def test_analyse_writes_stats_long_results_workbook_prism_csvs_plots_and_a_locked_report(tmp_path):
    from kerncount.experiment_types import RunOptions, get
    from kerncount.manifest import STATS_LONG_COLUMNS

    t = _vol_table(np.random.default_rng(10), n=lambda cell: 4 if cell == ("WT", "N4", "+clodro") else 6)
    get("tumour_volume").analyse(tmp_path, t, _settings(), RunOptions())
    st = pd.read_csv(tmp_path / "stats_long.csv")
    assert list(st.columns) == STATS_LONG_COLUMNS
    fams = set(st["family"])
    assert {"F1", "F2", "F3", "F4", "F5", "F6", "F7", "estimated means"} <= fams
    assert (st.loc[st.family == "F1", "estimate_type"] == "ratio").all()
    xl = pd.ExcelFile(tmp_path / "results.xlsx")
    assert {"About", "Tumours", "Statistics", "Estimated means", "Families"} <= set(xl.sheet_names)
    tum = pd.read_excel(xl, "Tumours")
    assert "main_mass_volume_pct" in tum.columns and len(tum) == len(t)
    prism = pd.read_csv(tmp_path / "prism" / "main_mass_volume_prism.csv")
    assert prism.shape[1] == 12 and prism.count().sum() == len(t)
    assert "WT / N4 / +clodro" in prism.columns and prism["WT / N4 / +clodro"].count() == 4
    assert (tmp_path / "prism" / "main_mass_volume_pct_prism.csv").exists()
    for kind in ("absolute", "normalised", "forest", "interaction"):
        assert (tmp_path / "plots" / f"main_mass_volume_{kind}.png").exists(), kind
    report = (tmp_path / "stats_report.md").read_text()
    lock = st["note"].dropna().str.extract(r"lock (\w{12})")[0].dropna().unique()
    assert len(lock) == 1 and lock[0] in report
    assert "main_mass_volume (confirmatory)" in report and "F6" in report and "descriptive" in report
    about = " ".join(pd.read_excel(xl, "About", header=None)[0].astype(str))
    assert lock[0] in about


def test_lock_hash_changes_when_the_family_spec_or_readout_split_changes():
    from kerncount.tumour import results

    base = _settings()
    other = _settings(readouts=[{"name": "main_mass_volume", "column": "main_mass_volume_mm3",
                                 "role": "exploratory"}])
    assert results.lock(base)[0] != results.lock(other)[0]
    assert results.lock(base)[0] == results.lock(_settings())[0]


def test_analyse_with_one_genotype_and_groups_of_one_does_not_fail(tmp_path):
    from kerncount.experiment_types import RunOptions, get

    t = _vol_table(np.random.default_rng(11), genotypes=("TLR7 KO",), n=lambda cell: 1)
    get("tumour_volume").analyse(tmp_path, t, _settings(), RunOptions())
    st = pd.read_csv(tmp_path / "stats_long.csv")
    assert st["note"].str.contains("needs WT data").any()
    assert st["note"].str.contains("no residual df").any()
    assert (tmp_path / "plots" / "main_mass_volume_forest.png").exists()


def test_a_declared_readout_missing_from_the_table_is_reported_not_fatal(tmp_path, capsys):
    from kerncount.experiment_types import RunOptions, get

    t = _vol_table(np.random.default_rng(12))
    s = _settings(readouts=[{"name": "main_mass_volume", "column": "main_mass_volume_mm3", "role": "confirmatory"},
                            {"name": "surface", "column": "main_mass_surface_mm2", "role": "confirmatory"}])
    get("tumour_volume").analyse(tmp_path, t, s, RunOptions())
    assert "main_mass_surface_mm2" in (tmp_path / "stats_report.md").read_text()
