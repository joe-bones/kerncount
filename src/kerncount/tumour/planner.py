"""Tumour power planner (ticket 09, ADR 0003): a **planning aid**, never part of the confirmatory results.

From the pooled residual SD of the current data (log scale, `stats.fit`), per declared family (descriptive ones are
never tested, so not planned) and confirmatory readout:

- **minimum detectable fold-change** (MDFC) of every contrast at the current group sizes: the fold-change (either
  direction, symmetric on the log scale) a two-sided t test detects with the target power (default 80 %);
- **tumours per group** needed for a chosen fold-change (default 1.5×): balanced (every group of the contrast at k)
  and as a top-up of the current design (raise only the groups below k), naming the **limiting group** (the group
  of the contrast contributing most to its variance, i.e. the smallest).

Contrast SE = s·√(Σ w²/n) (as in `stats.estimate`); power = P(|T| > t_crit), T noncentral t (df, δ/SE).
**Multiplicity:** Holm within family is approximated by **Bonferroni within family** (α / m, m = contrasts the
family tests): Holm's first step *is* Bonferroni and later steps are less strict, so the power to reject any one
contrast under Holm is at least the Bonferroni power — the numbers are conservative.

Contrasts needing cells without tumours (e.g. F4 genotype × line with only KO imaged) are planned from the variance
of the present data, with the note "assumes <absent level> variance equals that of <present level>".
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from .config import Factor
from .families import Cell, contrasts
from .stats import CellModel

PLAN_COLUMNS = ["readout", "family", "family_role", "comparison", "groups", "se", "df", "tests_in_family",
                "alpha_per_test", "mdfc", "target_fold_change", "power_at_target", "n_balanced", "limiting_group",
                "n_topup", "note"]
CAP = 500  # largest group size searched


def power(log_effect: float, se: float, df: float, alpha: float) -> float:
    """Two-sided t-test power at level alpha: P(|T| > t_{1-alpha/2, df}), T ~ noncentral t(df, |log_effect| / se)."""
    if not (df > 0 and se > 0 and np.isfinite(se)):
        return float("nan")
    tc = float(sps.t.ppf(1 - alpha / 2, df))
    nc = abs(log_effect) / se
    up = float(sps.nct.sf(tc, df, nc))
    if not np.isfinite(up):
        up = 1 - float(sps.nct.cdf(tc, df, nc))
    lo = float(np.nan_to_num(sps.nct.cdf(-tc, df, nc), nan=0.0))
    return float(np.clip(up + lo, 0, 1))


def detectable(se: float, df: float, alpha: float, target: float = 0.8) -> float:
    """Smallest |log effect| with `target` power."""
    from ..planner import nc_for_power

    if not (df > 0 and se > 0 and np.isfinite(se)):
        return float("nan")
    return nc_for_power(df, target, alpha) * se


def _session_confounded(facs: list[Factor], weights: dict[Cell, float]) -> bool:
    """Design-level version of `stats.confounded`: factors read from the imaging-session folder (source = session)
    are confounded with session; a contrast whose weights do not cancel within each of their level combinations is
    never tested, so it is not planned."""
    idx = [i for i, f in enumerate(facs) if f.source == "session"]
    per: dict[tuple, float] = {}
    for cell, w in weights.items():
        key = tuple(cell[i] for i in idx)
        per[key] = per.get(key, 0.0) + w
    return any(abs(v) > 1e-9 for v in per.values())


def _absent_note(model: CellModel, missing: list[Cell]) -> str:
    absent, present = [], []
    for i, f in enumerate(model.factors):
        have = [lv for lv in f.levels if any(c[i] == lv for c in model.n)]
        gone = [lv for lv in f.levels if lv not in have and any(c[i] == lv for c in missing)]
        if gone:
            absent += gone
            present += have
    if absent:
        return (f"needs {', '.join(absent)} data; planned from the present data only: assumes {', '.join(absent)} "
                f"variance equals that of {', '.join(present) or 'the present groups'}")
    return ("needs data for " + "; ".join(" / ".join(c) for c in missing)
            + "; planned assuming their variance equals the pooled variance")


def plan_readout(model: CellModel, facs: list[Factor], fams: list, readout, fold_change: float = 1.5,
                 target: float = 0.8, alpha: float = 0.05) -> pd.DataFrame:
    """One row per tested contrast of every non-descriptive family (see module docstring). The Bonferroni count m is
    every contrast the family tests in the full declared design (also those still needing data), so the plan holds
    once every level is imaged; today's Holm may count fewer."""
    s = float(np.sqrt(model.s2)) if model.df > 0 else float("nan")
    n_cells_design = int(np.prod([len(f.levels) for f in facs]))
    delta = abs(float(np.log(fold_change)))
    rows = []
    for fam in fams:
        if fam.role == "descriptive":
            continue
        cs = [c for c in contrasts(facs, fam) if not _session_confounded(facs, c.weights)]
        m = len(cs)
        a = alpha / m if m else alpha
        for c in cs:
            n = {cell: model.n.get(cell, 0) for cell in c.weights}
            missing = [cell for cell, k in n.items() if k == 0]
            notes = [_absent_note(model, missing)] if missing else []
            if not np.isfinite(s):
                notes.append("no residual df in the current data (every group n = 1): cannot plan")
            v = {cell: c.weights[cell] ** 2 for cell in c.weights}  # w² per cell
            se = s * float(np.sqrt(sum(v[cell] / n[cell] for cell in n))) if not missing else float("nan")
            mdfc = float(np.exp(detectable(se, model.df, a, target))) if not missing else float("nan")
            pw = power(delta, se, model.df, a) if not missing else float("nan")
            # balanced: every group of the design at k tumours (df = cells × (k − 1))
            k_bal = _search(lambda k: power(delta, s * float(np.sqrt(sum(v.values()) / k)),
                                            n_cells_design * (k - 1), a) >= target)
            # top-up: raise only the contrast's groups below k (other groups unchanged; each tumour adds 1 df)
            def topup(k: int) -> bool:
                nk = {cell: max(n[cell], k) for cell in n}
                df = model.df + sum(nk[cell] - n[cell] for cell in n) - len(missing)
                return power(delta, s * float(np.sqrt(sum(v[cell] / nk[cell] for cell in n))), df, a) >= target
            k_top = _search(topup)
            worst = max(v[cell] / n[cell] if n[cell] else np.inf for cell in n)
            limiting = [cell for cell in n if (v[cell] / n[cell] if n[cell] else np.inf) == worst]
            rows.append(dict(
                readout=readout.name, family=fam.id, family_role=fam.role, comparison=c.comparison,
                groups="; ".join(f"{' / '.join(cell)}: {n[cell]}" for cell in n), se=se, df=model.df,
                tests_in_family=m, alpha_per_test=a, mdfc=mdfc, target_fold_change=fold_change, power_at_target=pw,
                n_balanced=k_bal, limiting_group="; ".join(f"{' / '.join(cell)} (n {n[cell]})" for cell in limiting),
                n_topup=k_top, note="; ".join(notes)))
    return pd.DataFrame(rows, columns=PLAN_COLUMNS)


def _search(ok, lo: int = 2, hi: int = CAP) -> float:
    """Smallest k in [lo, hi] with ok(k) (power is monotone in k); NaN when not reached by hi."""
    if not ok(hi):
        return float("nan")
    while lo < hi:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return float(lo)


SUMMARY_COLUMNS = ["readout", "family", "family_role", "tests_in_family", "alpha_per_test", "mdfc_best", "mdfc_worst",
                   "target_fold_change", "power_at_target_worst", "n_balanced", "limiting_contrast", "limiting_group",
                   "n_topup", "note"]


def summarise(plan: pd.DataFrame) -> pd.DataFrame:
    """One row per readout and family: MDFC range over the contrasts, and the hardest contrast (largest MDFC, or the
    largest balanced n when no contrast is estimable yet) with its limiting group and the n that covers every
    contrast of the family."""
    rows = []
    for (ro, fid), g in plan.groupby(["readout", "family"], sort=False):
        key = g["mdfc"] if g["mdfc"].notna().any() else g["n_balanced"]
        hard = g.loc[key.idxmax()] if key.notna().any() else g.iloc[0]
        counts = g["note"][g["note"] != ""].value_counts(sort=False)
        notes = [x if k == len(g) else f"{k} of {len(g)} contrasts: {x}" for x, k in counts.items()]
        rows.append(dict(readout=ro, family=fid, family_role=g["family_role"].iloc[0],
                         tests_in_family=int(g["tests_in_family"].iloc[0]), alpha_per_test=g["alpha_per_test"].iloc[0],
                         mdfc_best=g["mdfc"].min(), mdfc_worst=g["mdfc"].max(),
                         target_fold_change=g["target_fold_change"].iloc[0],
                         power_at_target_worst=g["power_at_target"].min(), n_balanced=g["n_balanced"].max(),
                         limiting_contrast=hard["comparison"], limiting_group=hard["limiting_group"],
                         n_topup=hard["n_topup"], note=" | ".join(notes)))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def options(settings: dict) -> tuple[float, float, float]:
    """(fold-change, power, alpha) from [tumour] planner_fold_change / planner_power / planner_alpha."""
    from .config import tumour_settings

    ts = tumour_settings(settings)
    return float(ts["planner_fold_change"]), float(ts["planner_power"]), float(ts["planner_alpha"])


def report_lines(summ: pd.DataFrame, s_by_readout: dict[str, float], fold_change: float, target: float,
                 alpha: float) -> list[str]:
    """`## Planning aid` section of stats_report.md (kept apart from the confirmatory results)."""
    def f(x, fmt="{:.2f}"):
        return fmt.format(x) if pd.notna(x) else "—"

    md = ["## Planning aid (power planner; not part of the confirmatory results)", "",
          "- From the pooled residual SD of the current data (log scale) of each confirmatory readout; tumours "
          "treated as independent, as in the analysis.",
          f"- Minimum detectable fold-change (MDFC) = the fold-change (up or down) a two-sided t test detects with "
          f"{target:.0%} power at the current group sizes; α = {alpha} per family.",
          "- Multiplicity: Holm within family is approximated by **Bonferroni within family** (α / number of tested "
          "contrasts). Holm is never stricter than Bonferroni for any one contrast, so these numbers are "
          "conservative.",
          f"- Tumours per group for a fold-change of {fold_change} (`[tumour] planner_fold_change`): *balanced* = "
          "every group at k; *top-up* = raise only the groups of the hardest contrast below k (the limiting group "
          "is its smallest group). Contrasts confounded with imaging session and descriptive families are not "
          "planned (they are never tested).",
          "- Families needing absent levels (e.g. genotype × tumour line with one genotype imaged) are planned from "
          "the present data's variance: see the note (assumes the absent level's variance equals the present one).",
          ""]
    for ro, g in summ.groupby("readout", sort=False):
        md += [f"### {ro} (pooled SD {f(s_by_readout.get(ro, float('nan')))} on log scale)", "",
               f"| family | tests | α per test | MDFC (best – worst) | power for {fold_change}× (worst) | "
               f"n per group for {fold_change}× (balanced) | limiting contrast | limiting group | top-up to | note |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for _, r in g.iterrows():
            from .sensitivity import md_cell

            md.append(f"| {r['family']} [{r['family_role']}] | {r['tests_in_family']} | {r['alpha_per_test']:.4g} | "
                      f"{f(r['mdfc_best'])} – {f(r['mdfc_worst'])} | {f(r['power_at_target_worst'], '{:.0%}')} | "
                      f"{f(r['n_balanced'], '{:.0f}')} | {md_cell(r['limiting_contrast'])} | {r['limiting_group']} | "
                      f"{f(r['n_topup'], '{:.0f}')} | {r['note']} |")
        md.append("")
    return md
