"""Sensitivity models for the tumour statistics (ticket 08, ADR 0003): every contrast of every family re-estimated
with **per-group variances** (Welch-type) and tested by **permutation**; a concordance flag marks contrasts where the
conclusion (Holm-adjusted significance at 0.05, direction) differs from the primary pooled-variance model.

- Welch-type: the same point estimate Σ w·mean (cell means of log readout), SE² = Σ w²·s²_cell / n_cell with each
  cell's own variance, Satterthwaite df = SE⁴ / Σ (w²·s²_cell / n_cell)² / (n_cell − 1). Needs n ≥ 2 in every cell the
  contrast uses.
- Permutation (see `permutation`): simple effects permute the focus-factor labels between the two compared groups
  (they lie in one cell of the other factors); interactions use Freedman–Lane residual permutation under the reduced
  model (the cell-means model with this contrast constrained to 0), within the cells the contrast uses.

Nothing here knows readout or factor names; contrasts come from `families.contrasts` (weights on cells).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as sps

from .config import Factor
from .families import Cell, Contrast
from .stats import cell_of, usable

Groups = dict[Cell, np.ndarray]  # cell -> log(readout) of its tumours


def groups(table: pd.DataFrame, facs: list[Factor], column: str) -> Groups:
    """log(readout) per present cell (same usable rows as `stats.fit`)."""
    ok = usable(table, facs, column)
    t = table[ok]
    y = np.log(pd.to_numeric(t[column]).to_numpy(float))
    keys = cell_of(t, facs)
    return {cell: y[t.index.get_indexer(idx)] for cell, idx in keys.groupby(keys).groups.items()}


@dataclass
class WelchEstimate:
    ratio: float
    ci_low: float
    ci_high: float
    p: float
    se: float  # log scale
    df: float  # Satterthwaite
    note: str = ""


def welch(g: Groups, c: Contrast, level: float = 0.95) -> WelchEstimate:
    nan = float("nan")
    if any(cell not in g for cell in c.weights):
        return WelchEstimate(nan, nan, nan, nan, nan, nan, "needs data")
    est = float(sum(w * g[cell].mean() for cell, w in c.weights.items()))
    if any(len(g[cell]) < 2 for cell in c.weights):
        return WelchEstimate(float(np.exp(est)), nan, nan, nan, nan, nan, "a group of one has no own variance")
    parts = np.array([w * w * g[cell].var(ddof=1) / len(g[cell]) for cell, w in c.weights.items()])
    dfs = np.array([len(g[cell]) - 1 for cell in c.weights], float)
    v = float(parts.sum())
    if v <= 0:
        return WelchEstimate(float(np.exp(est)), nan, nan, nan, 0.0, nan, "no within-group variation")
    df = v * v / float((parts ** 2 / dfs).sum())
    se = float(np.sqrt(v))
    q = float(sps.t.ppf(0.5 + level / 2, df))
    p = float(2 * sps.t.sf(abs(est) / se, df))
    return WelchEstimate(float(np.exp(est)), float(np.exp(est - q * se)), float(np.exp(est + q * se)), p, se, df)


SIMPLE_SCHEME = ("label permutation: the focus-factor labels of the two compared groups are permuted (both lie in "
                 "one cell of the other factors; other levels of the focus factor are not exchangeable under this "
                 "contrast's null)")
FL_SCHEME = ("Freedman–Lane: residuals of the reduced model (cell means with this contrast constrained to 0) "
             "permuted across the tumours of the cells the contrast uses, added back to the reduced fit; statistic "
             "= t of the contrast")


@dataclass
class Permutation:
    p: float
    n_perm: int  # permutations evaluated (all distinct splits when exact)
    exact: bool
    seed: int
    scheme: str
    note: str = ""


def _rng(seed: int, c: Contrast) -> np.random.Generator:
    """Per-contrast stream: reproducible from the run seed and the contrast, independent of evaluation order."""
    import zlib

    return np.random.default_rng([int(seed), zlib.crc32(c.comparison.encode())])


def permutation(g: Groups, c: Contrast, n_perm: int = 9999, seed: int = 0) -> Permutation:
    """Two-sided permutation p of one contrast. Simple effects (two cells): exact enumeration of every split when
    there are at most `n_perm` of them (p = share of splits at least as extreme, observed included), else `n_perm`
    random label permutations (p = (1 + #extreme) / (n_perm + 1)). Interactions: Freedman–Lane with `n_perm` random
    residual permutations. Both use |t| of the contrast (for two groups equivalent to |mean difference|)."""
    nan = float("nan")
    cells = list(c.weights)
    scheme = SIMPLE_SCHEME if len(cells) == 2 else FL_SCHEME
    if any(cell not in g for cell in cells):
        return Permutation(nan, 0, False, seed, scheme, "needs data")
    ys = [g[cell] for cell in cells]
    n = np.array([len(y) for y in ys])
    w = np.array([c.weights[cell] for cell in cells], float)
    if n.sum() - len(cells) <= 0:
        return Permutation(nan, 0, False, seed, scheme, "no residual df")
    y = np.concatenate(ys)
    lab = np.repeat(np.arange(len(cells)), n)
    obs = _t(y[None, :], lab[None, :], w, n)[0]
    if not np.isfinite(obs):
        return Permutation(nan, 0, False, seed, scheme, "no within-group variation")
    tol = abs(obs) * (1 - 1e-9)  # ties (e.g. the mirrored split of equal groups) count as extreme
    rng = _rng(seed, c)
    if len(cells) == 2:
        from math import comb

        n_split = comb(int(n.sum()), int(n[0]))
        if n_split <= n_perm:
            import itertools

            first = np.zeros((n_split, int(n.sum())), bool)
            for i, pick in enumerate(itertools.combinations(range(int(n.sum())), int(n[0]))):
                first[i, list(pick)] = True
            labs = np.where(first, 0, 1)
            tt = _t(np.broadcast_to(y, labs.shape), labs, w, n)
            return Permutation(float(np.mean(np.abs(tt) >= tol)), n_split, True, seed, scheme)
        perms = np.argsort(rng.random((n_perm, len(y))), axis=1)
        tt = _t(y[perms], np.broadcast_to(lab, perms.shape), w, n)
    else:
        means = np.array([yy.mean() for yy in ys])
        lam = float(w @ means) / float((w * w / n).sum())
        fit_red = (means - lam * w / n)[lab]
        resid = y - fit_red
        perms = np.argsort(rng.random((n_perm, len(y))), axis=1)
        tt = _t(fit_red[None, :] + resid[perms], np.broadcast_to(lab, perms.shape), w, n)
    hits = int(np.sum(np.abs(tt) >= tol))
    return Permutation((1 + hits) / (n_perm + 1), n_perm, False, seed, scheme)


def _t(y: np.ndarray, lab: np.ndarray, w: np.ndarray, n: np.ndarray) -> np.ndarray:
    """t of the contrast Σ w·mean per row of y (rows = permutations), variance pooled over the contrast's cells.
    `lab` gives each column's cell; cell sizes `n` are fixed under permutation."""
    k = len(w)
    onehot = lab[..., None] == np.arange(k)  # (B, N, k)
    sums = np.einsum("bn,bnk->bk", y, onehot)
    means = sums / n
    fitted = np.take_along_axis(means, lab, axis=1)
    ss = ((y - fitted) ** 2).sum(axis=1)
    s2 = ss / (n.sum() - k)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (means @ w) / np.sqrt(s2 * float((w * w / n).sum()))


SENSITIVITY_COLUMNS = ["readout", "readout_role", "family", "family_role", "comparison", "estimate", "N",
                       "primary_p_adj", "welch_ci_low", "welch_ci_high", "welch_df", "welch_p", "welch_p_adj",
                       "perm_p", "perm_p_adj", "perm_n", "perm_exact", "perm_seed", "concordant", "concordance",
                       "note"]
ALPHA = 0.05


def _sig(p: float) -> bool:
    return bool(np.isfinite(p) and p < ALPHA)


def readout_sensitivity(table: pd.DataFrame, facs: list[Factor], fams: list, readout, primary: pd.DataFrame,
                        n_perm: int = 9999, seed: int = 0) -> pd.DataFrame:
    """Welch-type and permutation results for every contrast of `primary` (stats.readout_results, same families and
    order). Only contrasts the primary model tests are tested here (confounded, descriptive and not-estimable ones
    stay untested); Holm within each family, per method. Concordance (tested contrasts): the Holm-adjusted
    significance at 0.05 of each sensitivity method equals the primary model's (the point estimate, hence the
    direction, is the same in all three)."""
    from ..analysis import holm
    from .families import contrasts

    g = groups(table, facs, readout.column)
    prim = primary.set_index(["family", "comparison"])
    rows = []
    for fam in fams:
        fam_rows = []
        for c in contrasts(facs, fam):
            pr = prim.loc[(fam.id, c.comparison)]
            tested = bool(np.isfinite(pr["p"]))
            r = dict(readout=readout.name, readout_role=readout.role, family=fam.id, family_role=fam.role,
                     comparison=c.comparison, estimate=pr["estimate"], N=pr["N"], primary_p_adj=pr["p_adj"],
                     welch_ci_low=np.nan, welch_ci_high=np.nan, welch_df=np.nan, welch_p=np.nan,
                     welch_p_adj=np.nan, perm_p=np.nan, perm_p_adj=np.nan, perm_n=0, perm_exact=False,
                     perm_seed=seed, concordant=None, concordance="", note="")
            notes = []
            if not tested:
                notes.append("not tested (as in the primary model" + (f": {pr['note']})" if pr["note"] else ")"))
            else:
                w = welch(g, c)
                r.update(welch_ci_low=w.ci_low, welch_ci_high=w.ci_high, welch_df=w.df, welch_p=w.p)
                if w.note:
                    notes.append(f"Welch: {w.note}")
                pm = permutation(g, c, n_perm, seed)
                r.update(perm_p=pm.p, perm_n=pm.n_perm, perm_exact=pm.exact)
                if pm.note:
                    notes.append(f"permutation: {pm.note}")
            r["note"] = "; ".join(notes)
            fam_rows.append(r)
        for col in ("welch", "perm"):
            idx = [i for i, r in enumerate(fam_rows) if np.isfinite(r[f"{col}_p"])]
            for i, pa in zip(idx, holm([fam_rows[i][f"{col}_p"] for i in idx])):
                fam_rows[i][f"{col}_p_adj"] = pa
        for r in fam_rows:
            if not np.isfinite(r["primary_p_adj"]):
                continue
            ref = _sig(r["primary_p_adj"])
            diffs = [f"{name} {'' if _sig(r[col]) else 'not '}significant (Holm p {r[col]:.3g})"
                     for name, col in (("Welch", "welch_p_adj"), ("permutation", "perm_p_adj"))
                     if np.isfinite(r[col]) and _sig(r[col]) != ref]
            r["concordant"] = not diffs
            r["concordance"] = ("agrees" if not diffs else "differs from primary "
                                f"({'significant' if ref else 'not significant'}): " + "; ".join(diffs))
        rows += fam_rows
    return pd.DataFrame(rows, columns=SENSITIVITY_COLUMNS)


# ---------------------------------------------------------------- outputs

WELCH_MODEL = "Welch-type: cell means of log(readout), per-group variances, Satterthwaite df"
PERM_MODEL = "permutation: label permutation (simple effects) / Freedman–Lane (interactions)"


def options(settings: dict) -> tuple[int, int]:
    """(permutations, seed) from [tumour] sensitivity_permutations / sensitivity_seed."""
    from .config import tumour_settings

    ts = tumour_settings(settings)
    return int(ts["sensitivity_permutations"]), int(ts["sensitivity_seed"])


def long_rows(sens: pd.DataFrame, readout, lock_hash: str) -> list[dict]:
    """stats_long rows (kerncount long format) for the tested contrasts: family `<id> (Welch-type)` and
    `<id> (permutation)`, like `<family> (leave-one-n-out)` in the cell-count statistics. Never starred."""
    rows = []
    for _, r in sens[sens["welch_p"].notna() | sens["perm_p"].notna()].iterrows():
        conc = ("concordant with primary: " + ("yes" if r["concordant"] else "no — " + r["concordance"])
                if pd.notna(r["concordant"]) else "")
        base = dict(level="all_n", n=None, metric=readout.name, comparison=r["comparison"], estimate=r["estimate"],
                    estimate_type="ratio", significance="", N=r["N"], estimand=f"fold-change of {readout.label}",
                    cell_line="")
        tail = [x for x in (f"sensitivity ({r['family_role']} family)", conc, r["note"], f"lock {lock_hash}") if x]
        rows.append(dict(base, family=f"{r['family']} (Welch-type)",
                         test=f"sensitivity: Welch-type t (Satterthwaite df {r['welch_df']:.1f}), Holm within family",
                         ci_low=r["welch_ci_low"], ci_high=r["welch_ci_high"], p=r["welch_p"],
                         p_adj=r["welch_p_adj"], model=WELCH_MODEL, note="; ".join(tail)))
        how = (f"exact, all {int(r['perm_n'])} splits" if r["perm_exact"] else
               f"{int(r['perm_n'])} random permutations, seed {int(r['perm_seed'])}")
        rows.append(dict(base, family=f"{r['family']} (permutation)",
                         test=f"sensitivity: permutation test ({how}), Holm within family", ci_low=np.nan,
                         ci_high=np.nan, p=r["perm_p"], p_adj=r["perm_p_adj"], model=PERM_MODEL,
                         note="; ".join([*tail[:-1], f"seed {int(r['perm_seed'])}", tail[-1]])))
    return rows


def residual_plot(g: Groups, facs: list[Factor], label: str, out) -> object:
    """Residual diagnostics of the primary model: log residuals per group (with each group's SD) and a normal QQ plot
    of the studentised residuals (pooled SD; groups of one have no residual and are left out)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import DOT, EST, INK, MUTED, _save

    order = sorted(g, key=lambda c: tuple(f.levels.index(x) for f, x in zip(facs, c)))
    res = {c: g[c] - g[c].mean() for c in order}
    df = sum(len(v) for v in res.values()) - len(res)
    allr = np.concatenate(list(res.values())) if res else np.array([])
    s = float(np.sqrt((allr ** 2).sum() / df)) if df > 0 else float("nan")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(max(6.0, 2.5 + 0.45 * len(order)) + 3.4, 3.8),
                                 gridspec_kw=dict(width_ratios=[max(2.0, 0.45 * len(order)), 3.0]))
    rng = np.random.default_rng(0)
    for i, c in enumerate(order):
        r = res[c]
        a1.scatter(i + rng.uniform(-0.15, 0.15, len(r)), r, s=12, color=DOT, alpha=0.8, linewidths=0)
        if len(r) > 1:
            a1.text(i, 1.02, f"SD {r.std(ddof=1):.2f}", transform=a1.get_xaxis_transform(), ha="center",
                    fontsize=6, color=MUTED, rotation=90, va="bottom")
    a1.axhline(0, color=MUTED, linewidth=0.8)
    a1.set_xticks(range(len(order)), [" / ".join(c) for c in order], rotation=60, ha="right", fontsize=7)
    a1.set_ylabel("residual, log scale")
    a1.set_title("residuals by group", fontsize=9, color=INK, pad=34)
    if np.isfinite(s) and s > 0 and sum(len(v) > 1 for v in res.values()):
        z = np.sort(np.concatenate([v / (s * np.sqrt(1 - 1 / len(v))) for v in res.values() if len(v) > 1]))
        q = sps.norm.ppf((np.arange(1, len(z) + 1) - 0.375) / (len(z) + 0.25))
        a2.scatter(q, z, s=10, color=DOT, linewidths=0)
        lim = [min(q.min(), z.min()), max(q.max(), z.max())]
        a2.plot(lim, lim, color=EST, linewidth=1)
    else:
        a2.text(0.5, 0.5, "no residual df", ha="center", va="center", transform=a2.transAxes, color=MUTED)
    a2.set_xlabel("normal quantile")
    a2.set_ylabel("studentised residual")
    a2.set_title(f"normal QQ (pooled SD {s:.2f} on log scale)", fontsize=9, color=INK)
    for ax in (a1, a2):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(f"{label}: residual diagnostics of the primary model (log scale)", fontsize=9, color=MUTED)
    fig.tight_layout()
    return _save(fig, out)


def md_cell(text: str) -> str:
    return str(text).replace("|", "\\|")  # comparisons contain "|" (within …), which would split a table cell


def report_lines(sens: list[pd.DataFrame], n_perm: int, seed: int) -> list[str]:
    """`## Sensitivity` section of stats_report.md: method, then per readout and family the contrasts with Welch and
    permutation Holm p and the concordance flag."""
    md = ["## Sensitivity (per-group variance and permutation; ADR 0003)", "",
          f"- {WELCH_MODEL}; same point estimates as the primary model.",
          f"- Permutation, simple effects: {SIMPLE_SCHEME}. Exact enumeration of every split when there are at most "
          f"{n_perm}, else {n_perm} random permutations.",
          f"- Permutation, interactions: {FL_SCHEME}; {n_perm} random permutations.",
          f"- Random permutations are seeded: seed {seed} (`[tumour] sensitivity_seed`) combined with each contrast's "
          "label, so reruns give identical p-values.",
          "- Holm within each family and readout, per method. Contrasts the primary model does not test (confounded "
          "with imaging session, descriptive, needing absent groups) are not tested here either.",
          "- **Concordant** = the Holm-adjusted significance at 0.05 of both sensitivity methods agrees with the "
          "primary model (direction always agrees: all three share the point estimate). Sensitivity results never "
          "carry stars.",
          ""]
    for s in sens:
        tested = s[s["concordant"].notna()]
        if not len(tested):
            continue
        n_diff = int((tested["concordant"] == False).sum())  # noqa: E712
        md += [f"### {s['readout'].iloc[0]} ({s['readout_role'].iloc[0]}): "
               f"{len(tested) - n_diff} of {len(tested)} tested contrasts concordant", "",
               "| family | comparison | p (Holm), primary | Welch | permutation | concordant |",
               "|---|---|---|---|---|---|"]
        for _, r in tested.iterrows():
            f = lambda x: f"{x:.3g}" if np.isfinite(x) else ""  # noqa: E731
            flag = "yes" if r["concordant"] else f"**no** — {r['concordance']}"
            md.append(f"| {r['family']} | {md_cell(r['comparison'])} | {f(r['primary_p_adj'])} | "
                      f"{f(r['welch_p_adj'])} | {f(r['perm_p_adj'])} | {flag} |")
        md.append("")
    return md
