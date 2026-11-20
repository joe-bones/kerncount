"""Factorial statistics for tumour readouts (ADR 0003 + amendment): one N-way linear model on log(readout).

The model has one mean per **cell** (combination of all declared factor levels = full interaction) and one residual
variance pooled across all cells; it is fitted with whatever cells are present. A contrast (families.py) is a set of
weights on cell means: estimate = Σ w·mean, SE = s·√(wᵀVw) (V = diag(1/n) without block effects), t with the pooled
residual df; exp() gives a fold-change (ratio of geometric means) or a ratio of ratios. Nothing here knows the
readout's or the factors' names.

Batch (ADR 0003 amendment 2026-11-20, ticket 13): the batch of every tumour (`batch` column, `[tumour] batch`, default
the culture batch) enters as a **fixed block effect** (additive on the log scale = a multiplicative shift of all its
tumours). Batches that share a cell are *connected*; within a connected set (here: the WT culture batches 190626 and
100726, which share cntrl and HERV −clodro) the block effects are estimated, so contrasts inside it are within-batch
estimates. Between connected sets a shift cannot be estimated: a contrast whose weights do not cancel within every
connected set is **confounded** (never tested); ratios of within-set ratios (F4/F5) cancel. With one batch per set the
model is the plain cell-means model (identical numbers). Block effects are centred within each set (tumour-weighted
mean 0), so a cell mean is the mean at the set's average batch.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats as sps

from .config import Factor
from .families import Cell, Contrast

CONFOUNDED = "confounded with imaging session"
CONFOUNDED_DAY = "confounded with acquisition day"
CONFOUNDED_BATCH = "confounded with culture batch"
SESSION_COLUMN = "imaging_session"
BATCH_COLUMN = "batch"  # tumours.csv: the imaging batch per `[tumour] batch` (table.build); else imaging_session


def batch_of(table: pd.DataFrame) -> tuple[str, str]:
    """(column holding each tumour's batch, label of contrasts confounded with it). `batch` (tickets 12, 13) when
    present: it equals culture_batch (batch = "culture_batch"), imaging_session ("session_folder") or acquisition_day
    ("acquisition_day")."""
    if BATCH_COLUMN in table.columns:
        b = table[BATCH_COLUMN].astype(str)
        for col, label in (("culture_batch", CONFOUNDED_BATCH), (SESSION_COLUMN, CONFOUNDED),
                           ("acquisition_day", CONFOUNDED_DAY)):
            if col in table.columns:
                c = table[col].fillna("").astype(str)
                if b.equals(c) or (b[c != ""].equals(c[c != ""]) and b[c == ""].str.startswith("(").all()):
                    return BATCH_COLUMN, label
        return BATCH_COLUMN, "confounded with imaging batch"
    return SESSION_COLUMN, CONFOUNDED


@dataclass
class CellModel:
    factors: list[Factor]
    column: str
    n: dict[Cell, int]  # tumours per present cell
    mean: dict[Cell, float]  # mean log(readout) per present cell
    s2: float  # pooled residual variance (log scale); nan when df = 0
    df: int  # residual degrees of freedom = tumours − present cells
    sessions: dict[Cell, dict[str, float]] = field(default_factory=dict)  # cell -> {batch: share of its tumours}
    dropped: int = 0  # rows left out (missing, non-positive or outside the declared levels)
    confounded_label: str = CONFOUNDED  # what a batch is (culture batch / imaging session / day), for notes
    # ticket 13: fixed block effects (log scale, centred within each connected set; 0 for a set of one batch)
    blocks: dict[str, float] = field(default_factory=dict)
    component: dict[Cell, int] = field(default_factory=dict)  # cell -> connected set of batches it belongs to
    cov: dict[tuple[Cell, Cell], float] = field(default_factory=dict)  # unscaled covariance of cell means (≠ diag)
    batch_column: str = ""
    batch_sets: list[list[str]] = field(default_factory=list)  # connected sets of batches (sorted)

    def v(self, a: Cell, b: Cell) -> float:
        """Unscaled covariance of two cell means (× s² = covariance); diag(1/n) without block effects."""
        if (a, b) in self.cov:
            return self.cov[(a, b)]
        return 1.0 / self.n[a] if a == b and not self.cov else 0.0

    @property
    def n_blocks(self) -> int:
        """Block-effect parameters estimated (batches beyond the first of each connected set)."""
        return len(self.blocks) - len(set(self.component.values())) if self.blocks else 0

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


def fit(table: pd.DataFrame, facs: list[Factor], column: str, session: str | None = None) -> CellModel:
    """`session` = the batch column (default: `batch_of(table)`); confounding is judged per batch."""
    label = CONFOUNDED
    if session is None:
        session, label = batch_of(table)
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
    batches = (t[session].astype(str).to_numpy() if session in t.columns else np.array([""] * len(t), dtype=object))
    component, sets = _connected(list(keys), list(batches))
    blocks = {b: 0.0 for bs in sets for b in bs}
    m = CellModel(list(facs), column, n, mean, float("nan"), len(y) - len(n), sessions, int((~ok).sum()), label,
                  blocks, component, {}, session if session in t.columns else "", sets)
    extra = [b for bs in sets for b in bs[1:]]
    if extra:
        ss = _fit_blocks(m, y, list(keys), list(batches), sets, extra)
    m.s2 = ss / m.df if m.df > 0 else float("nan")
    return m


def _connected(cells: list[Cell], batches: list[str]) -> tuple[dict[Cell, int], list[list[str]]]:
    """Connected sets of batches (two batches are connected when they share a cell): {cell: set index}, and the
    sorted batches of every set."""
    parent: dict[tuple, tuple] = {}

    def find(x: tuple) -> tuple:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for c, b in zip(cells, batches):
        ra, rb = find(("c", c)), find(("b", b))
        if ra != rb:
            parent[ra] = rb
    roots: dict[tuple, int] = {}
    sets: list[list[str]] = []
    for b in sorted(set(batches)):
        r = find(("b", b))
        if r not in roots:
            roots[r] = len(sets)
            sets.append([])
        sets[roots[r]].append(b)
    return {c: roots[find(("c", c))] for c in set(cells)}, sets


def _fit_blocks(m: CellModel, y: np.ndarray, cells: list[Cell], batches: list[str], sets: list[list[str]],
                extra: list[str]) -> float:
    """Least squares of cell means + block effects (the first batch of every set is the reference); then centred:
    within each set the tumour-weighted mean block effect is 0. Fills m.mean, m.blocks, m.cov, m.df; returns RSS."""
    order = list(m.n)
    ci = {c: i for i, c in enumerate(order)}
    bi = {b: len(order) + j for j, b in enumerate(extra)}
    p = len(order) + len(extra)
    X = np.zeros((len(y), p))
    for r, (c, b) in enumerate(zip(cells, batches)):
        X[r, ci[c]] = 1.0
        if b in bi:
            X[r, bi[b]] = 1.0
    xtx_inv = np.linalg.pinv(X.T @ X)
    theta = xtx_inv @ X.T @ y
    rss = float(((y - X @ theta) ** 2).sum())
    rank = int(np.linalg.matrix_rank(X))
    # centring: mean block effect of set k (tumour-weighted, reference = 0) as a row on theta
    counts = pd.Series(batches).value_counts().to_dict()
    bar = np.zeros((len(sets), p))
    for k, bs in enumerate(sets):
        tot = sum(counts.get(b, 0) for b in bs)
        for b in bs[1:]:
            bar[k, bi[b]] = counts.get(b, 0) / tot
    A = np.zeros((len(order) + len(m.blocks), p))
    for c, i in ci.items():
        A[i, i] = 1.0
        A[i] += bar[m.component[c]]
    names = [b for bs in sets for b in bs]
    for j, b in enumerate(names):
        k = next(k for k, bs in enumerate(sets) if b in bs)
        row = len(order) + j
        if b in bi:
            A[row, bi[b]] = 1.0
        A[row] -= bar[k]
    est = A @ theta
    cov = A[: len(order)] @ xtx_inv @ A[: len(order)].T
    m.mean = {c: float(est[ci[c]]) for c in order}
    m.blocks = {b: float(est[len(order) + j]) for j, b in enumerate(names)}
    m.cov = {(a, b): float(cov[ci[a], ci[b]]) for a in order for b in order}
    m.df = len(y) - rank
    return rss


def batch_adjusted(table: pd.DataFrame, model: CellModel) -> pd.DataFrame:
    """`table` with the model's readout divided by exp(its batch's block effect): the values the sensitivity models
    (Welch-type, permutation) and residual plots use, so that they compare the same within-batch contrasts. Their cell
    means are the primary model's cell means. Unchanged when there are no block effects."""
    if not model.n_blocks or not model.batch_column:
        return table
    out = table.copy()
    shift = np.exp(table[model.batch_column].astype(str).map(model.blocks).fillna(0.0).astype(float))
    out[model.column] = pd.to_numeric(table[model.column], errors="coerce") / shift
    return out


def confounded(model: CellModel, weights: dict[Cell, float]) -> bool:
    """True when the contrast's weight does not cancel within every connected set of batches (`[tumour] batch`):
    a multiplicative shift of one set (e.g. a genotype prepared in its own culture batches) would then move the
    estimate. Shifts between batches of one set are estimated (block effects), so contrasts inside a set are not
    confounded; ratios of within-set ratios cancel."""
    if model.component:
        per_set: dict[int, float] = {}
        for cell, w in weights.items():
            if cell in model.component:
                per_set[model.component[cell]] = per_set.get(model.component[cell], 0.0) + w
        return any(abs(v) > 1e-9 for v in per_set.values())
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
    se = float(np.sqrt(model.s2 * sum(wa * wb * model.v(a, b) for a, wa in c.weights.items()
                                      for b, wb in c.weights.items())))
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
                    notes.append(m.confounded_label + " (not tested)")
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
    other factor (ADR 0003 amendment) and in the tumour's own batch (block effect, ticket 13). Reference tumours
    scatter around 100 %; NaN where that reference is absent."""
    ok = usable(table, model.factors, model.column)
    out = pd.Series(np.nan, index=table.index, dtype=float)
    keys = cell_of(table[ok], model.factors)
    v = pd.to_numeric(batch_adjusted(table, model).loc[ok, model.column]).astype(float)
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
        se = float(np.sqrt(model.s2 * model.v(cell, cell))) if model.df > 0 else float("nan")
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
