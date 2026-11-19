"""Factorial statistics for tumour readouts (ADR 0003 + amendment): one N-way linear model on log(readout).

The model has one mean per **cell** (combination of all declared factor levels = full interaction) and one residual
variance pooled across all cells; it is fitted with whatever cells are present. A contrast (families.py) is a set of
weights on cell means: estimate = Σ w·mean, SE = s·√(Σ w²/n), t with the pooled residual df; exp() gives a fold-change
(ratio of geometric means) or a ratio of ratios. Nothing here knows the readout's or the factors' names.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats as sps

from .config import Factor
from .families import Cell, Contrast

CONFOUNDED = "confounded with imaging session"
SESSION_COLUMN = "imaging_session"


@dataclass
class CellModel:
    factors: list[Factor]
    column: str
    n: dict[Cell, int]  # tumours per present cell
    mean: dict[Cell, float]  # mean log(readout) per present cell
    s2: float  # pooled residual variance (log scale); nan when df = 0
    df: int  # residual degrees of freedom = tumours − present cells
    sessions: dict[Cell, dict[str, float]] = field(default_factory=dict)  # cell -> {session: share of its tumours}
    dropped: int = 0  # rows left out (missing, non-positive or outside the declared levels)

    @property
    def N(self) -> int:
        return sum(self.n.values())


@dataclass
class Estimate:
    ratio: float  # exp(contrast): fold-change or ratio of ratios
    ci_low: float
    ci_high: float
    p: float
    se: float  # log scale
    df: int
    n: int  # tumours in the cells the contrast uses
    missing: list[Cell] = field(default_factory=list)  # cells the contrast needs that have no tumours
    confounded: bool = False  # a multiplicative imaging-session shift does not cancel


def cell_of(table: pd.DataFrame, facs: list[Factor]) -> pd.Series:
    return pd.Series(list(zip(*(table[f.name].astype(str) for f in facs))), index=table.index, dtype=object)


def usable(table: pd.DataFrame, facs: list[Factor], column: str) -> pd.Series:
    """Rows the model can use: positive finite readout and every factor at a declared level."""
    v = pd.to_numeric(table[column], errors="coerce")
    ok = np.isfinite(v) & (v > 0)
    for f in facs:
        ok &= table[f.name].astype(str).isin(f.levels)
    return ok


def fit(table: pd.DataFrame, facs: list[Factor], column: str, session: str = SESSION_COLUMN) -> CellModel:
    ok = usable(table, facs, column)
    t = table[ok]
    y = np.log(pd.to_numeric(t[column]).to_numpy(float))
    keys = cell_of(t, facs)
    n: dict[Cell, int] = {}
    mean: dict[Cell, float] = {}
    sessions: dict[Cell, dict[str, float]] = {}
    ss = 0.0
    for cell, idx in keys.groupby(keys).groups.items():
        pos = t.index.get_indexer(idx)
        yc = y[pos]
        n[cell], mean[cell] = len(yc), float(yc.mean())
        ss += float(((yc - yc.mean()) ** 2).sum())
        if session in t.columns:
            sessions[cell] = t.loc[idx, session].astype(str).value_counts(normalize=True).to_dict()
    df = len(y) - len(n)
    return CellModel(list(facs), column, n, mean, ss / df if df > 0 else float("nan"), df, sessions,
                     int((~ok).sum()))


def confounded(model: CellModel, weights: dict[Cell, float]) -> bool:
    """True when the contrast's weight does not cancel within every imaging session: a multiplicative shift of one
    session (day effect) would then move the estimate. Within-session ratios (and ratios of them) cancel."""
    per: dict[str, float] = {}
    for cell, w in weights.items():
        for s, share in model.sessions.get(cell, {}).items():
            per[s] = per.get(s, 0.0) + w * share
    return any(abs(v) > 1e-9 for v in per.values())


def estimate(model: CellModel, c: Contrast, level: float = 0.95) -> Estimate:
    missing = [cell for cell in c.weights if cell not in model.n]
    n = sum(model.n.get(cell, 0) for cell in c.weights)
    if missing:
        nan = float("nan")
        return Estimate(nan, nan, nan, nan, nan, model.df, n, missing)
    est = sum(w * model.mean[cell] for cell, w in c.weights.items())
    se = float(np.sqrt(model.s2 * sum(w * w / model.n[cell] for cell, w in c.weights.items())))
    conf = confounded(model, c.weights)
    if model.df > 0 and se > 0:
        q = float(sps.t.ppf(0.5 + level / 2, model.df))
        p = float(2 * sps.t.sf(abs(est) / se, model.df))
        lo, hi = float(np.exp(est - q * se)), float(np.exp(est + q * se))
    else:
        p = lo = hi = float("nan")
    return Estimate(float(np.exp(est)), lo, hi, p, se, model.df, n, [], conf)


RESULT_COLUMNS = ["readout", "readout_role", "family", "family_role", "family_label", "comparison", "estimate",
                  "ci_low", "ci_high", "p", "p_adj", "significance", "N", "df", "confounded", "note"]


def needs_text(model: CellModel, missing: list[Cell]) -> str:
    """'needs WT data' when whole levels are absent, else the empty groups by name."""
    absent: list[str] = []
    for i, f in enumerate(model.factors):
        have = {cell[i] for cell in model.n}
        absent += [lv for lv in f.levels if lv not in have and any(c[i] == lv for c in missing) and lv not in absent]
    if absent:
        return f"needs {', '.join(absent)} data"
    return "needs data for " + "; ".join(" / ".join(c) for c in missing)


def readout_results(table: pd.DataFrame, facs: list[Factor], fams: list, readout, model: CellModel | None = None
                    ) -> pd.DataFrame:
    """Every contrast of every family for one readout, Holm within each family, labelled (ADR 0003):
    stars only for confirmatory families of confirmatory readouts; descriptive families and contrasts confounded with
    imaging session are shown with estimate + CI but never tested; contrasts needing absent cells are kept as rows
    saying which data they need."""
    from ..analysis import holm, stars
    from .families import contrasts

    m = model or fit(table, facs, readout.column)
    rows = []
    for fam in fams:
        fam_rows = []
        for c in contrasts(facs, fam):
            e = estimate(m, c)
            notes = []
            tested = True
            if e.missing:
                notes.append(needs_text(m, e.missing))
                tested = False
            else:
                if e.confounded:
                    notes.append(CONFOUNDED + " (not tested)")
                    tested = False
                if fam.role == "descriptive":
                    notes.append("descriptive only (not tested)")
                    tested = False
                if m.df <= 0:
                    notes.append("no residual df (every group n = 1): not tested")
                    tested = False
            if fam.role == "exploratory":
                notes.append("exploratory family")
            if readout.role != "confirmatory":
                notes.append("exploratory readout")
            fam_rows.append(dict(readout=readout.name, readout_role=readout.role, family=fam.id,
                                 family_role=fam.role, family_label=fam.label, comparison=c.comparison,
                                 estimate=e.ratio, ci_low=e.ci_low, ci_high=e.ci_high,
                                 p=e.p if tested else np.nan, p_adj=np.nan, significance="", N=e.n, df=m.df,
                                 confounded=e.confounded, note="; ".join(notes)))
        tested_idx = [i for i, r in enumerate(fam_rows) if np.isfinite(r["p"])]
        for i, pa in zip(tested_idx, holm([fam_rows[i]["p"] for i in tested_idx])):
            fam_rows[i]["p_adj"] = pa
            if fam.role == "confirmatory" and readout.role == "confirmatory":
                fam_rows[i]["significance"] = stars(pa)
        rows += fam_rows
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def _ref_cell(model: CellModel, cell: Cell, focus: Factor) -> Cell:
    i = [f.name for f in model.factors].index(focus.name)
    return cell[:i] + (focus.reference,) + cell[i + 1:]


def normalised(table: pd.DataFrame, model: CellModel, focus: Factor) -> pd.Series:
    """Each tumour as % of the model-estimated mean of the focus factor's reference level in the same cell of every
    other factor (ADR 0003 amendment). Reference tumours scatter around 100 %; NaN where that reference is absent."""
    ok = usable(table, model.factors, model.column)
    out = pd.Series(np.nan, index=table.index, dtype=float)
    keys = cell_of(table[ok], model.factors)
    v = pd.to_numeric(table.loc[ok, model.column]).astype(float)
    for i in keys.index:
        ref = _ref_cell(model, keys[i], focus)
        if ref in model.mean:
            out[i] = 100.0 * v[i] / np.exp(model.mean[ref])
    return out


def cell_means(model: CellModel, level: float = 0.95) -> pd.DataFrame:
    """Model-estimated (geometric) mean per present cell with a CI from the pooled variance."""
    q = float(sps.t.ppf(0.5 + level / 2, model.df)) if model.df > 0 else float("nan")
    rows = []
    for cell in sorted(model.n, key=lambda c: tuple(f.levels.index(x) for f, x in zip(model.factors, c))):
        se = float(np.sqrt(model.s2 / model.n[cell])) if model.df > 0 else float("nan")
        mu = model.mean[cell]
        rows.append({**{f.name: x for f, x in zip(model.factors, cell)}, "n": model.n[cell],
                     "mean": float(np.exp(mu)), "ci_low": float(np.exp(mu - q * se)),
                     "ci_high": float(np.exp(mu + q * se))})
    return pd.DataFrame(rows, columns=[*(f.name for f in model.factors), "n", "mean", "ci_low", "ci_high"])


def normalised_groups(model: CellModel, focus: Factor) -> pd.DataFrame:
    """Per present cell: model mean as % of its reference cell, with the CI of that fold-change (×100) — the same
    numbers as the simple-effect contrast of the focus factor in the forest plot. The reference cell is 100 % (no CI)."""
    from .families import Contrast
    from .config import Family

    fam = Family("normalised", "simple", (focus.name,), "descriptive")
    rows = []
    for _, r in cell_means(model).iterrows():
        cell = tuple(str(r[f.name]) for f in model.factors)
        ref = _ref_cell(model, cell, focus)
        pct = lo = hi = float("nan")
        if cell == ref:
            pct = 100.0
        elif ref in model.mean:
            e = estimate(model, Contrast(fam, "", (), "", {cell: 1.0, ref: -1.0}))
            pct, lo, hi = 100 * e.ratio, 100 * e.ci_low, 100 * e.ci_high
        rows.append({**{f.name: x for f, x in zip(model.factors, cell)}, "n": int(r["n"]), "pct": pct,
                     "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows, columns=[*(f.name for f in model.factors), "n", "pct", "ci_low", "ci_high"])
