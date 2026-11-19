"""Tumour power planner (ticket 09, ADR 0003): minimum detectable fold-change per family at the current group sizes
and tumours per group for a chosen fold-change, from the pooled residual SD of the current data (log scale).
Bonferroni within family (conservative for Holm). Analytic power validated by simulation. Synthetic tables only."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from test_tumour_stats import ALL, F1, F4, FACS, VOL, synth_table


def _row(plan, label):
    return plan.set_index("comparison").loc[label]


def _simulated_power(rng, sims, n, sigma, contrast, log_effect, m, alpha=0.05):
    """Share of simulated experiments in which the contrast is significant after Bonferroni over m contrasts, data
    drawn with the true effect on the contrast's first positive cell and residual SD sigma."""
    from kerncount.tumour import stats

    target = next(cell for cell, w in contrast.weights.items() if w > 0)
    hits = 0
    for _ in range(sims):
        rows = []
        for cell in [(g, li, c) for g in FACS[0].levels for li in FACS[1].levels for c in FACS[2].levels]:
            k = n(cell)
            if not k:
                continue
            mu = log_effect if cell == target else 0.0
            for x in np.exp(mu + rng.normal(0, sigma, k)):
                rows.append(dict(zip(("genotype", "tumour_line", "clodronate"), cell), volume=x,
                                 imaging_session=f"session {cell[0]}"))
        e = stats.estimate(stats.fit(pd.DataFrame(rows), FACS, "volume"), contrast)
        hits += min(1.0, m * e.p) < alpha
    return hits / sims


def test_minimum_detectable_fold_change_has_80_percent_power_in_simulation():
    from kerncount.tumour import planner, stats
    from kerncount.tumour.families import contrasts

    t = synth_table(np.random.default_rng(1))  # groups 12 / 6, sigma 0.3
    m = stats.fit(t, FACS, "volume")
    plan = planner.plan_readout(m, FACS, [F1], VOL)
    label = "N4 vs cntrl | WT, +clodro"  # 6 vs 6
    r = _row(plan, label)
    assert r["tests_in_family"] == 8 and r["alpha_per_test"] == pytest.approx(0.05 / 8)
    assert 1.2 < r["mdfc"] < 2.5
    c = next(c for c in contrasts(FACS, F1) if c.comparison == label)
    sizes = {cell: n for cell, n in m.n.items()}
    pw = _simulated_power(np.random.default_rng(2), 1500, lambda cell: sizes.get(cell, 0), np.sqrt(m.s2), c,
                          np.log(r["mdfc"]), 8)
    assert pw == pytest.approx(0.80, abs=0.035)


def test_large_df_matches_the_normal_approximation():
    from kerncount.tumour import planner, stats

    t = synth_table(np.random.default_rng(3), n=lambda cell: 60)
    m = stats.fit(t, FACS, "volume")
    plan = planner.plan_readout(m, FACS, [F1], VOL, alpha=0.05)
    r = _row(plan, "HERV vs cntrl | TLR7 KO, -clodro")
    from scipy.stats import norm

    se = np.sqrt(m.s2 * 2 / 60)
    z = norm.ppf(1 - 0.05 / 8 / 2) + norm.ppf(0.8)
    assert np.log(r["mdfc"]) == pytest.approx(z * se, rel=0.01)


def test_tumours_per_group_for_the_chosen_fold_change_reach_80_percent_power_in_simulation():
    from kerncount.tumour import planner, stats
    from kerncount.tumour.families import contrasts

    t = synth_table(np.random.default_rng(4))
    m = stats.fit(t, FACS, "volume")
    plan = planner.plan_readout(m, FACS, [F1], VOL, fold_change=1.5)
    label = "HERV vs cntrl | TLR7 KO, +clodro"
    k = int(_row(plan, label)["n_balanced"])
    assert 5 <= k <= 40
    c = next(c for c in contrasts(FACS, F1) if c.comparison == label)
    s = np.sqrt(m.s2)
    rng = np.random.default_rng(5)
    at_k = _simulated_power(rng, 1000, lambda cell: k, s, c, np.log(1.5), 8)
    below = _simulated_power(rng, 1000, lambda cell: k - 2, s, c, np.log(1.5), 8)
    assert 0.77 <= at_k <= 0.90 and below < 0.80


def test_unbalanced_design_names_the_limiting_group_and_a_top_up_that_reaches_the_target():
    from kerncount.tumour import planner, stats

    small = ("TLR7 KO", "N4", "+clodro")
    t = synth_table(np.random.default_rng(6), n=lambda cell: 3 if cell == small else 12)
    m = stats.fit(t, FACS, "volume")
    plan = planner.plan_readout(m, FACS, [F1], VOL, fold_change=1.5)
    r = _row(plan, "N4 vs cntrl | TLR7 KO, +clodro")
    assert r["limiting_group"] == "TLR7 KO / N4 / +clodro (n 3)"
    other = _row(plan, "HERV vs cntrl | TLR7 KO, +clodro")
    assert r["mdfc"] > other["mdfc"]  # the small group makes this contrast the hardest to detect
    assert r["power_at_target"] < other["power_at_target"]
    # top-up: raising only the 3-tumour group (cntrl stays at 12) — needs at least the balanced k for that group,
    # and with the cntrl side already at 12 never more than ~ the balanced plan in the limiting group
    assert r["n_topup"] >= 3 and r["n_topup"] <= 2 * r["n_balanced"]


def test_only_tested_families_and_contrasts_are_planned_with_bonferroni_within_family():
    from kerncount.tumour import planner, stats

    t = synth_table(np.random.default_rng(7))
    plan = planner.plan_readout(stats.fit(t, FACS, "volume"), FACS, ALL, VOL)
    assert "F6" not in set(plan["family"])  # descriptive (and confounded with imaging session)
    counts = plan.groupby("family")["comparison"].count().to_dict()
    assert counts == {"F1": 8, "F2": 6, "F3": 4, "F4": 4, "F5": 3, "F7": 2}
    assert (plan["alpha_per_test"] * plan["tests_in_family"]).round(12).eq(0.05).all()


def test_genotype_by_line_is_planned_from_ko_variance_when_wt_is_absent():
    from kerncount.tumour import planner, stats

    t = synth_table(np.random.default_rng(8), genotypes=("TLR7 KO",))
    plan = planner.plan_readout(stats.fit(t, FACS, "volume"), FACS, ALL, VOL)
    f4 = plan[plan.family == "F4"]
    assert len(f4) == 4
    assert f4["mdfc"].isna().all()  # nothing detectable with no WT tumours at all
    assert f4["n_balanced"].notna().all() and f4["n_topup"].notna().all()
    assert f4["note"].str.contains("assumes WT variance equals that of TLR7 KO").all()
    # a ratio of ratios has four groups: it needs more tumours per group than a simple effect at the same fold-change
    f1 = plan[plan.family == "F1"].dropna(subset=["mdfc"])
    assert f4["n_balanced"].min() > f1["n_balanced"].max()


def test_family_summary_takes_the_hardest_contrast_as_limiting():
    from kerncount.tumour import planner, stats

    small = ("TLR7 KO", "N4", "+clodro")
    t = synth_table(np.random.default_rng(9), n=lambda cell: 3 if cell == small else 12)
    plan = planner.plan_readout(stats.fit(t, FACS, "volume"), FACS, [F1, F4], VOL)
    s = planner.summarise(plan).set_index("family")
    f1 = plan[plan.family == "F1"]
    assert s.loc["F1", "mdfc_worst"] == pytest.approx(f1["mdfc"].max())
    assert s.loc["F1", "mdfc_best"] == pytest.approx(f1["mdfc"].min())
    assert s.loc["F1", "n_balanced"] == f1["n_balanced"].max()
    assert s.loc["F1", "limiting_contrast"] == "N4 vs cntrl | TLR7 KO, +clodro"
    assert "TLR7 KO / N4 / +clodro (n 3)" in s.loc["F1", "limiting_group"]


# ---------------------------------------------------------------- end-to-end: separated from confirmatory results

def test_analyse_writes_a_separate_planning_aid_for_confirmatory_readouts_only(tmp_path):
    from kerncount.experiment_types import RunOptions, get
    from test_tumour_stats import _settings, _vol_table

    t = _vol_table(np.random.default_rng(10), genotypes=("TLR7 KO",))
    t["shape_index"] = np.exp(np.random.default_rng(1).normal(0, 0.2, len(t)))
    s = _settings(planner_fold_change=2.0, sensitivity_permutations=199,
                  readouts=[{"name": "main_mass_volume", "column": "main_mass_volume_mm3", "role": "confirmatory"},
                            {"name": "shape", "column": "shape_index", "role": "exploratory"}])
    get("tumour_volume").analyse(tmp_path, t, s, RunOptions())
    xl = pd.ExcelFile(tmp_path / "results.xlsx")
    assert {"Planning aid", "Planning aid (contrasts)"} <= set(xl.sheet_names)
    summ = pd.read_excel(xl, "Planning aid")
    rows = pd.read_excel(xl, "Planning aid (contrasts)")
    assert set(rows["readout"]) == {"main_mass_volume"}  # confirmatory readouts only
    assert (rows["target_fold_change"] == 2.0).all()
    assert "F4" in set(summ["family"]) and summ.set_index("family").loc["F4", "note"].count("assumes WT variance")
    st = pd.read_csv(tmp_path / "stats_long.csv")
    assert not st["family"].astype(str).str.contains("plan", case=False).any()
    report = (tmp_path / "stats_report.md").read_text()
    sec = report[report.index("## Planning aid"):]
    assert "not part of the confirmatory results" in sec and "Bonferroni" in sec and "2.0" in sec
