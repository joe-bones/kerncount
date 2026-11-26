"""Tumour sensitivity models (ticket 08, ADR 0003): every contrast re-estimated with per-group variances (Welch-type,
Satterthwaite df) and tested by permutation; a concordance flag marks where conclusions differ from the primary
pooled-variance model. All tests use synthetic readout tables (no images)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from test_tumour_stats import ALL, F1, F2, F4, F6, FACS, VOL, synth_table


def _contrast(fam, label):
    from kerncount.tumour.families import contrasts

    return next(c for c in contrasts(FACS, fam) if c.comparison == label)


def _logs(t, cell):
    g, line, c = cell
    sel = (t.genotype == g) & (t.tumour_line == line) & (t.clodronate == c)
    return np.log(t.loc[sel, "volume"].to_numpy(float))


# ---------------------------------------------------------------- seam 1: Welch-type estimate of one contrast

def test_welch_simple_contrast_equals_the_two_sample_welch_t_test():
    from kerncount.tumour import sensitivity

    t = synth_table(np.random.default_rng(1), n=lambda cell: {"cntrl": 12, "HERV": 5, "N4": 7}[cell[1]])
    groups = sensitivity.groups(t, FACS, "volume")
    c = _contrast(F1, "HERV vs cntrl | WT, -clodro")
    e = sensitivity.welch(groups, c)
    a, b = _logs(t, ("WT", "HERV", "-clodro")), _logs(t, ("WT", "cntrl", "-clodro"))
    ref = sps.ttest_ind(a, b, equal_var=False)
    assert e.ratio == pytest.approx(np.exp(a.mean() - b.mean()))
    assert e.p == pytest.approx(ref.pvalue)
    assert e.df == pytest.approx(ref.df)
    lo, hi = ref.confidence_interval()
    assert (e.ci_low, e.ci_high) == pytest.approx((np.exp(lo), np.exp(hi)))


def hetero_table(rng, sd, n, genotypes=("WT", "TLR7 KO"), shift=None) -> pd.DataFrame:
    """No true effects; per-cell SD sd(cell) and size n(cell); optional log shift(cell)."""
    rows = []
    for g in genotypes:
        for line in ("cntrl", "HERV", "N4"):
            for c in ("-clodro", "+clodro"):
                cell = (g, line, c)
                mu = shift(cell) if shift else 0.0
                for i, x in enumerate(np.exp(mu + rng.normal(0, sd(cell), n(cell)))):
                    rows.append({"tumour_id": f"{g}/{line}/{c}/{i}", "imaging_session": f"session {g}",
                                 "genotype": g, "tumour_line": line, "clodronate": c, "volume": x})
    return pd.DataFrame(rows)


def test_welch_interaction_ci_keeps_coverage_when_one_small_group_is_much_noisier():
    from kerncount.tumour import sensitivity, stats

    noisy = ("TLR7 KO", "HERV", "-clodro")
    sd = lambda cell: 0.9 if cell == noisy else 0.15  # noqa: E731
    n = lambda cell: 4 if cell == noisy else 12  # noqa: E731
    c = _contrast(F4, "(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | -clodro")
    rng = np.random.default_rng(2)
    sims, welch_hit, pooled_hit = 400, 0, 0
    for _ in range(sims):
        t = hetero_table(rng, sd, n)
        e = sensitivity.welch(sensitivity.groups(t, FACS, "volume"), c)
        welch_hit += e.ci_low <= 1.0 <= e.ci_high
        p = stats.estimate(stats.fit(t, FACS, "volume"), c)
        pooled_hit += p.ci_low <= 1.0 <= p.ci_high
    assert 0.92 <= welch_hit / sims <= 0.985
    assert pooled_hit / sims < 0.85  # the primary model's pooled variance is far too small for this contrast


# ---------------------------------------------------------------- seam 2: permutation test of one contrast

def test_simple_contrast_permutation_is_exact_and_matches_scipy_permutation_test():
    from kerncount.tumour import sensitivity

    t = synth_table(np.random.default_rng(3), n=lambda cell: 6, sigma=0.4)
    g = sensitivity.groups(t, FACS, "volume")
    c = _contrast(F1, "HERV vs cntrl | TLR7 KO, +clodro")
    r = sensitivity.permutation(g, c, n_perm=9999, seed=1)
    assert r.exact and r.n_perm == 924  # C(12, 6) splits of 6 + 6 tumours, all enumerated
    a, b = _logs(t, ("TLR7 KO", "HERV", "+clodro")), _logs(t, ("TLR7 KO", "cntrl", "+clodro"))
    ref = sps.permutation_test((a, b), lambda x, y, axis: x.mean(axis=axis) - y.mean(axis=axis),
                               permutation_type="independent", n_resamples=np.inf, vectorized=True)
    assert r.p == pytest.approx(ref.pvalue)
    assert r.p >= 2 / 924  # smallest attainable two-sided p for 6 vs 6


def test_large_groups_use_seeded_monte_carlo_permutations_that_are_reproducible():
    from kerncount.tumour import sensitivity

    t = synth_table(np.random.default_rng(4), n=lambda cell: 12, sigma=0.4)  # C(24, 12) = 2.7 million splits
    g = sensitivity.groups(t, FACS, "volume")
    c = _contrast(F1, "N4 vs cntrl | WT, -clodro")
    r1, r2 = (sensitivity.permutation(g, c, n_perm=1999, seed=7) for _ in range(2))
    r3 = sensitivity.permutation(g, c, n_perm=1999, seed=8)
    assert not r1.exact and r1.n_perm == 1999 and r1.seed == 7
    assert r1.p == r2.p and r1.p != r3.p
    assert r1.p >= 1 / 2000
    a, b = _logs(t, ("WT", "N4", "-clodro")), _logs(t, ("WT", "cntrl", "-clodro"))
    assert r1.p == pytest.approx(sps.ttest_ind(a, b).pvalue, abs=0.03)


def test_interaction_permutation_freedman_lane_holds_its_level_and_detects_a_planted_interaction():
    from kerncount.tumour import sensitivity

    c = _contrast(F4, "(N4 / cntrl in TLR7 KO) ÷ (N4 / cntrl in WT) | +clodro")
    rng = np.random.default_rng(5)
    # main effects present (line and genotype), no interaction: the reduced model keeps them, so p is ~uniform
    main = lambda cell: (0.7 if cell[1] == "N4" else 0.0) + (0.5 if cell[0] == "TLR7 KO" else 0.0)  # noqa: E731
    sims = 300
    ps = []
    for i in range(sims):
        t = hetero_table(rng, lambda cell: 0.3, lambda cell: 6, shift=main)
        ps.append(sensitivity.permutation(sensitivity.groups(t, FACS, "volume"), c, n_perm=199, seed=i).p)
    assert 0.02 <= np.mean(np.array(ps) <= 0.05) <= 0.09
    assert 0.4 <= np.mean(ps) <= 0.6
    big = lambda cell: main(cell) + (np.log(2.5) if cell == ("TLR7 KO", "N4", "+clodro") else 0)  # noqa: E731
    t = hetero_table(np.random.default_rng(6), lambda cell: 0.3, lambda cell: 6, shift=big)
    r = sensitivity.permutation(sensitivity.groups(t, FACS, "volume"), c, n_perm=4999, seed=1)
    assert not r.exact and r.p < 0.01 and "Freedman" in r.scheme


def test_freedman_lane_p_agrees_with_the_t_test_over_the_contrasts_cells_for_normal_equal_variance_data():
    from kerncount.tumour import sensitivity

    shift = lambda cell: np.log(1.4) if cell == ("WT", "HERV", "-clodro") else 0.0  # noqa: E731
    t = hetero_table(np.random.default_rng(7), lambda cell: 0.3, lambda cell: 8, shift=shift)
    c = _contrast(F4, "(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | -clodro")
    r = sensitivity.permutation(sensitivity.groups(t, FACS, "volume"), c, n_perm=9999, seed=3)
    # independent: classical t of the ratio of ratios, variance pooled over its 4 groups of 8 (df 28)
    ys = {cell: _logs(t, cell) for cell in c.weights}
    est = sum(w * ys[cell].mean() for cell, w in c.weights.items())
    s2 = sum(((y - y.mean()) ** 2).sum() for y in ys.values()) / 28
    p = 2 * sps.t.sf(abs(est) / np.sqrt(s2 * 4 / 8), 28)
    assert r.p == pytest.approx(p, abs=0.02)


# ------------------------------------------------------------ seam 3: one readout -> sensitivity table + concordance

def test_every_tested_contrast_gets_welch_and_permutation_holm_within_family_and_untested_stay_untested():
    from kerncount.tumour import sensitivity, stats

    t = synth_table(np.random.default_rng(8), ko_shift=np.log(2.0))
    t = t[~((t.genotype == "WT") & (t.tumour_line == "N4") & (t.clodronate == "+clodro"))]  # one empty group
    prim = stats.readout_results(t, FACS, ALL, VOL)
    s = sensitivity.readout_sensitivity(t, FACS, ALL, VOL, prim, n_perm=999, seed=11)
    assert len(s) == len(prim) and list(s["comparison"]) == list(prim["comparison"])
    tested = prim["p"].notna().to_numpy()
    assert s.loc[tested, ["welch_p", "perm_p", "concordant"]].notna().all().all()
    assert s.loc[~tested, ["welch_p", "welch_p_adj", "perm_p", "perm_p_adj"]].isna().all().all()
    assert s.loc[~tested, "concordant"].isna().all()
    assert s[s.family == "F6"]["note"].str.contains("not tested").all()
    for fid, g in s[tested].groupby("family"):
        for col in ("welch", "perm"):
            assert g[f"{col}_p_adj"].min() == pytest.approx(min(1.0, len(g) * g[f"{col}_p"].min())), (fid, col)
    np.testing.assert_allclose(s["estimate"], prim["estimate"])  # same point estimates (cell means)
    assert (s["perm_seed"] == 11).all()


def test_concordance_flag_marks_a_contrast_whose_significance_rests_on_the_pooled_variance():
    from kerncount.tumour import sensitivity, stats

    t = hetero_table(np.random.default_rng(9), lambda cell: 0.15, lambda cell: 10)
    noisy = (t.genotype == "TLR7 KO") & (t.tumour_line == "HERV") & (t.clodronate == "-clodro")
    t = t[~noisy | (t.groupby(noisy).cumcount() < 4)].copy()  # the noisy group has 4 tumours
    noisy = (t.genotype == "TLR7 KO") & (t.tumour_line == "HERV") & (t.clodronate == "-clodro")
    t.loc[noisy, "volume"] = np.exp(np.log(2.5) + np.array([-1.0, -0.3, 0.3, 1.0]))  # big shift, very noisy
    prim = stats.readout_results(t, FACS, [F1], VOL)
    s = sensitivity.readout_sensitivity(t, FACS, [F1], VOL, prim, n_perm=999, seed=1).set_index("comparison")
    key = "HERV vs cntrl | TLR7 KO, -clodro"
    assert prim.set_index("comparison").loc[key, "p_adj"] < 0.05  # the primary model calls it
    assert s.loc[key, "welch_p_adj"] > 0.05
    assert s.loc[key, "concordant"] == False  # noqa: E712
    assert "Welch" in s.loc[key, "concordance"]
    others = s.drop(index=key)
    assert others["concordant"].all()  # nothing else is significant anywhere (no true effects)


def test_clear_effects_with_equal_variances_are_concordant():
    from kerncount.tumour import sensitivity, stats

    t = synth_table(np.random.default_rng(10), sigma=0.1)
    prim = stats.readout_results(t, FACS, [F1, F2, F4], VOL)
    s = sensitivity.readout_sensitivity(t, FACS, [F1, F2, F4], VOL, prim, n_perm=999, seed=1)
    f1 = s[s.family == "F1"]
    assert f1["concordant"].all() and (f1["concordance"] == "agrees").all()
    assert (f1["perm_p_adj"] < 0.05).all() and (f1["welch_p_adj"] < 0.05).all()


# ---------------------------------------------------------------- seam 4: end-to-end analyse

def test_analyse_adds_sensitivity_sheet_long_rows_report_section_and_residual_plots(tmp_path):
    from kerncount.experiment_types import RunOptions, get
    from test_tumour_stats import _settings, _vol_table

    t = _vol_table(np.random.default_rng(12))
    s = _settings(sensitivity_permutations=499, sensitivity_seed=42)
    get("tumour_volume").analyse(tmp_path, t, s, RunOptions())
    xl = pd.ExcelFile(tmp_path / "results.xlsx")
    assert "Sensitivity" in xl.sheet_names
    sens = pd.read_excel(xl, "Sensitivity")
    assert {"welch_p_adj", "perm_p_adj", "concordant", "concordance"} <= set(sens.columns)
    assert (sens["perm_seed"] == 42).all()
    st = pd.read_csv(tmp_path / "stats_long.csv")
    w = st[st.family == "F1 (Welch-type)"]
    p = st[st.family == "F1 (permutation)"]
    assert len(w) == len(p) == 8 and (w["significance"].fillna("") == "").all()
    assert w["note"].str.contains("concordant with primary").all() and p["note"].str.contains("seed 42").all()
    assert not st["family"].str.startswith("F6 (").any()  # untested contrasts get no sensitivity rows
    assert st.loc[st.family == "F1", "p"].notna().all()  # the primary rows are unchanged
    report = (tmp_path / "stats_report.md").read_text()
    assert "## Sensitivity" in report and "Freedman–Lane" in report and "seed 42" in report
    assert (tmp_path / "plots" / "main_mass_volume_residuals.png").exists()
    # identical reruns (seeded)
    get("tumour_volume").analyse(tmp_path / "again", t, s, RunOptions())
    again = pd.read_excel(tmp_path / "again" / "results.xlsx", "Sensitivity")
    pd.testing.assert_series_equal(sens["perm_p"], again["perm_p"])
