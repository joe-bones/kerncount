"""Hierarchical (mixed) model across all n, per metric; comparison families; equivalence; robustness.

Unit of observation: one field (one row of the long table). Nesting: field in well in treatment in n (one plate).

Response scale (multiplicative effects -> log scale, results reported as ratios):
    count     y = log(count + 0.5)    0.5 keeps empty fields finite and is negligible for counts of ~10 or more
    positive  y = log(value)          (if a column has zeros: log(value + half its smallest positive value), noted)
    fraction  y = empirical logit log((k + 0.5) / (N - k + 0.5)) with k = fraction x N, N = nuclei in that field;
              without a usable N the fraction is clipped to [0.01, 0.99]. Results are odds ratios.

Model (REML; groups = n, everything else nested in n; design built by statsmodels MixedLM with vc_formula,
likelihood maximised over variances >= 0 by reml.py, see there why):
    y = treatment (fixed) [+ timepoint (fixed, only with >= 2 known timepoints)]
        + n                                    random intercept per experiment
        + n x treatment                        random, with a separate variance per treatment class
        + well within n x treatment            random (only when some condition has >= 2 wells)
        + field position within n              random n:position (only when [design] field_positions = "fixed")
        + field residual
    Position is a random n:position term rather than fixed position-in-n effects: every treatment is imaged at the
    same positions, so position is (nearly) orthogonal to treatment and only needs to soak up variance; the random
    term spends no degrees of freedom on 6 x n position means, shrinks noisy position means, and gives the position
    variance that the variance-components table reports. Its pattern may differ between days (n:position).

Fallbacks (first that fits wins, recorded in `model`): full -> common n x treatment variance -> without position.
A model "fails" when it raises, does not converge, or gives a non-finite / singular fixed-effect covariance.
A variance component estimated at exactly zero is a legitimate (boundary) estimate, not a failure: it is noted.
If every model fails the caller falls back to the RM ANOVA (and says so).

Small-sample inference (a Wald z is far too liberal with 3-6 n's): t tests with Satterthwaite df (from the
expected REML information), capped by a containment ceiling = the df of the n x treatment variance(s) the
comparison relies on (per class: cells - treatments - n's + 1; a one-treatment class: n's - 1), never above the RM
ANOVA error df, so a class variance estimated at zero cannot buy the comparison extra df. A comparison needs >= 2
n's containing both groups. Means: Satterthwaite capped at (n's with the treatment) - 1. Fixed-effect covariance =
GLS (X' V^-1 X)^-1 at the REML estimates. Null simulation: docs/statistics.md and tests/test_models.py.

Multiplicity: within each family, single-step multivariate t (Dunnett-type; correlation from the model covariance,
each comparison with its own df), giving adjusted p and simultaneous 95% CIs. The equivalence family uses Holm.
"""
from __future__ import annotations

import fnmatch
import re
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from .metrics import Metric

COUNT_OFFSET = 0.5
FRACTION_CLIP = 0.01
ALPHA = 0.05


def stars(p) -> str:
    if p is None or not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


# ---------------------------------------------------------------- configuration

@dataclass
class Family:
    name: str  # e.g. "vs LyoVec"
    kind: str  # vs_control | vehicle_check | positive_control | among
    pairs: list[tuple[str, str]]  # (treatment, comparator)
    reference: str | None = None
    adjust: str = "single-step multivariate t"

    @property
    def id(self) -> str:
        """File-name friendly id: vs_<control> | vehicle_check | vs_positive_control | within_<class>."""
        if self.kind == "vs_control":
            return "vs_" + _safe(self.reference)
        if self.kind == "positive_control":
            return "vs_positive_control"
        if self.kind == "among":
            return "within_" + _safe(self.name.split("within class ", 1)[-1].split(" (")[0])
        return self.kind


@dataclass
class Config:
    controls: list[str]
    primary_ref: str
    primary_test: str
    classes: dict[str, str]  # treatment -> class name
    class_order: list[str]
    vehicle_pairs: list[tuple[str, str]]
    positive_control: str
    among_classes: list[str]
    margin: float
    order: list[str]

    def families(self, present: list[str]) -> list[Family]:
        """Comparison families for the treatments present (in `order`)."""
        pres = [t for t in self.order if t in set(present)]
        fams = []
        for c in self.controls:
            if c in pres:
                fams.append(Family(f"vs {c}", "vs_control", [(t, c) for t in pres if t != c], c))
        vp = [(v, c) for v, c in self.vehicle_pairs if v in pres and c in pres]
        if vp:
            fams.append(Family("vehicle check (equivalence)", "vehicle_check", vp, vp[0][1], adjust="Holm"))
        pc = self.positive_control
        if pc and pc in pres:
            others = [t for t in pres if t != pc and t not in self.controls]
            if others:
                fams.append(Family(f"vs positive control {pc}", "positive_control", [(t, pc) for t in others], pc))
        for cl in self.among_classes:
            members = [t for t in pres if self.classes.get(t) == cl]
            pairs = [(members[i], members[j]) for i in range(len(members)) for j in range(i + 1, len(members))]
            if pairs:
                fams.append(Family(f"within class {cl} (all pairs)", "among", pairs, None))
        return fams


def _as_list(v) -> list[str]:
    if v is None or v == "":
        return []
    if isinstance(v, str):
        return [s.strip() for s in v.split(",") if s.strip()]
    return [str(s).strip() for s in v if str(s).strip()]


def parse_classes(v) -> dict[str, list[str]]:
    """[analysis] classes: a TOML table {name = ["pattern", ...]} (or its text form)."""
    if isinstance(v, str):
        v = v.strip()
        if not v:
            return {}
        import tomllib
        v = tomllib.loads(f"x = {v}")["x"]
    if not isinstance(v, dict):
        return {}
    return {str(k): _as_list(p) for k, p in v.items()}


def config(settings: dict, treatments) -> Config:
    from .analysis import control_groups, treatment_order

    a = settings["analysis"]
    order = treatment_order(list(treatments), a.get("treatment_order", []))
    controls = control_groups(a)
    primary = str(a.get("primary_family") or "").strip() or a["reference"]
    user = parse_classes(a.get("classes", {}))
    cls, corder = {}, []
    for t in order:
        hit = next((name for name, pats in user.items() if any(fnmatch.fnmatchcase(t, p) for p in pats)), None)
        if hit is None:
            hit = "controls" if t in controls else ("other" if user else "treatments")
        cls[t] = hit
        if hit not in corder:
            corder.append(hit)
    vc = str(a.get("vehicle_check", "auto")).strip()
    if vc.lower() == "auto":
        present_controls = [c for c in controls]
        pairs = [(c, present_controls[0]) for c in present_controls[1:]] if present_controls else []
    elif vc.lower() in ("", "none", "off"):
        pairs = []
    else:
        pairs = []
        for part in vc.split(";"):
            bits = re.split(r"\s+vs\.?\s+", part.strip())
            if len(bits) == 2 and all(bits):
                pairs.append((bits[0].strip(), bits[1].strip()))
    pt = str(a.get("primary_test", "mixed")).strip().lower()
    return Config(controls=controls, primary_ref=primary, primary_test=pt if pt in ("mixed", "rm_anova") else "mixed",
                  classes=cls, class_order=corder, vehicle_pairs=pairs,
                  positive_control=str(a.get("positive_control") or "").strip(),
                  among_classes=_as_list(a.get("among_classes", [])),
                  margin=float(a.get("equivalence_margin", 0.2) or 0.2), order=order)


# ---------------------------------------------------------------- response scale

@dataclass
class Scale:
    kind: str
    shift: float = 0.0  # added before the log (count offset, or zero fix for positive metrics)
    note: str = ""

    def back(self, mu):
        mu = np.asarray(mu, float)
        if self.kind == "fraction":
            return 1 / (1 + np.exp(-mu))
        return np.exp(mu) - self.shift

    @property
    def estimate_type(self) -> str:
        return "odds_ratio" if self.kind == "fraction" else "ratio"

    @property
    def describe(self) -> str:
        if self.kind == "count":
            return f"log(count + {self.shift:g})"
        if self.kind == "fraction":
            return "logit (empirical)"
        return "log" if not self.shift else f"log(value + {self.shift:.3g})"


def transform(values: pd.Series, kind: str, counts: pd.Series | None = None) -> tuple[pd.Series, Scale]:
    """Response on the analysis scale; rows that cannot be transformed become NaN."""
    v = pd.to_numeric(values, errors="coerce").astype(float)
    if kind == "count":
        v = v.where(v >= 0)
        return np.log(v + COUNT_OFFSET), Scale("count", COUNT_OFFSET, f"y = log(count + {COUNT_OFFSET:g})")
    if kind == "fraction":
        p = v.where((v >= 0) & (v <= 1))
        N = pd.to_numeric(counts, errors="coerce") if counts is not None else None
        if N is not None and (N > 0).any():
            k = p * N
            y = np.log((k + 0.5) / (N - k + 0.5))
            fb = ~(N > 0)
            pc = p.clip(FRACTION_CLIP, 1 - FRACTION_CLIP)
            y = y.where(~fb, np.log(pc / (1 - pc)))
            return y, Scale("fraction", note="y = empirical logit log((k + 0.5) / (N - k + 0.5)), N = nuclei in the field")
        pc = p.clip(FRACTION_CLIP, 1 - FRACTION_CLIP)
        return np.log(pc / (1 - pc)), Scale("fraction", note=f"y = logit(fraction clipped to [{FRACTION_CLIP}, {1 - FRACTION_CLIP}])")
    v = v.where(v >= 0)
    pos = v[v > 0]
    if (v == 0).any() and len(pos):
        shift = float(pos.min()) / 2
        return np.log(v + shift), Scale("positive", shift, f"y = log(value + {shift:.3g}); {int((v == 0).sum())} zero value(s)")
    neg = int((pd.to_numeric(values, errors="coerce") < 0).sum())
    return np.log(v.where(v > 0)), Scale("positive", 0.0, "y = log(value)" + (f"; {neg} negative value(s) left out" if neg else ""))


# ---------------------------------------------------------------- data preparation

def _safe(s: str) -> str:
    return re.sub(r"\W", "_", str(s))


@dataclass
class Prepared:
    d: pd.DataFrame
    scale: Scale
    treatments: list[str]  # present, in order; first = baseline of the fixed effects
    ns: list[str]
    timepoint_levels: list[str]
    can_well: bool
    can_position: bool
    notes: list[str] = field(default_factory=list)
    nt_df: float = 1.0  # df of the n x treatment interaction (the RM ANOVA error df): ceiling for contrast df


def prepare(df: pd.DataFrame, metric: Metric, cfg: Config, positions_ok: bool) -> Prepared:
    d = df.copy()
    d["n"] = d["n"].astype(str)
    d["treatment"] = d["treatment"].astype(str)
    y, scale = transform(d["value"], metric.kind, d["count"] if "count" in d and metric.kind == "fraction" else None)
    d["y"] = y
    notes = [scale.note] if scale.note else []
    bad = d["y"].isna() | ~np.isfinite(d["y"])
    if bad.any():
        notes.append(f"{int(bad.sum())} field(s) without a usable value left out")
    d = d[~bad].copy()
    present = [t for t in cfg.order if t in set(d["treatment"])]
    base = cfg.primary_ref if cfg.primary_ref in present else next((c for c in cfg.controls if c in present),
                                                                    present[0] if present else None)
    treatments = ([base] + [t for t in present if t != base]) if base else []
    ns = sorted(d["n"].unique(), key=str)
    if "well" not in d:
        d["well"] = 1
    if "position" not in d:
        d["position"] = 1
    if "timepoint" not in d:
        d["timepoint"] = ""
    d["well"] = pd.to_numeric(d["well"], errors="coerce").fillna(1).astype(int)
    d["position"] = pd.to_numeric(d["position"], errors="coerce").fillna(1).astype(int)
    d["timepoint"] = d["timepoint"].fillna("").astype(str)
    d["cell"] = d["treatment"].map(_safe) + "_w" + d["well"].astype(str)
    d["pos"] = "p" + d["position"].astype(str)
    for i, t in enumerate(treatments):
        d[f"I{i}"] = (d["treatment"] == t).astype(float)
        if i:
            d[f"T{i}"] = d[f"I{i}"]
    wells_per_cell = d.groupby(["n", "treatment"])["well"].nunique()
    can_well = bool((wells_per_cell >= 2).any())
    can_pos = bool(positions_ok and d["position"].nunique() >= 2)
    if not positions_ok:
        notes.append("field positions are random ([design] field_positions): position not modelled")
    tps = sorted(set(d["timepoint"]) - {""})
    tp_by_n = d.groupby("n")["timepoint"].first()
    timepoint_levels = []
    if len(tps) >= 2:
        timepoint_levels = sorted(set(tp_by_n.replace("", "unknown")))
        d["tp"] = d["timepoint"].replace("", "unknown")
        for j, lv in enumerate(timepoint_levels[1:], 1):
            d[f"P{j}"] = (d["tp"] == lv).astype(float)
        notes.append(f"timepoint ({', '.join(timepoint_levels)}) is a fixed effect; it varies only between n's, so it "
                     "is confounded with n: its effect is a between-experiment comparison")
    elif tps:
        unknown = int((tp_by_n == "").sum())
        notes.append(f"timepoint {tps[0]}" + (f" (unknown for {unknown} n)" if unknown else "") +
                     ": only one known level, not modelled")
    cells = d.groupby(["n", "treatment"]).ngroups
    nt_df = max(1.0, float(cells - len(treatments) - len(ns) + 1))
    return Prepared(d, scale, treatments, ns, timepoint_levels, can_well, can_pos, notes, nt_df)


# ---------------------------------------------------------------- fitting

SPECS = [("full", "class-specific n×treatment variances + well + position"),
         ("common", "common n×treatment variance + well + position"),
         ("no_position", "common n×treatment variance + well, no position")]


@dataclass
class Fit:
    spec: str
    description: str
    formula: str
    reml: object  # reml.REMLFit
    beta: np.ndarray
    cov: np.ndarray
    names: list[str]  # fixed-effect names (Intercept, T1, ..., P1, ...)
    vc_names: list[str]  # random terms in theta order: n, then the variance components
    vc_labels: dict[str, str]
    group_labels: list[str]  # n of each group, in the order of the REML groups
    rows: list[np.ndarray]  # positions (0..len(d)-1) of each group's rows in Prepared.d
    converged: bool
    warnings: list[str]
    method: str
    vc_colnames: dict = field(default_factory=dict)  # term -> column names of its design (e.g. I0, I3, ...)

    def df(self, L: np.ndarray, cap: float) -> float:
        """Satterthwaite df for L beta, capped at `cap` (containment ceiling), at least 1."""
        try:
            s = self.reml.satterthwaite(L)
        except Exception:  # noqa: BLE001
            s = np.inf
        if not np.isfinite(s):
            s = cap
        return float(max(1.0, min(s, cap)))

    @property
    def theta(self) -> dict[str, float]:
        return dict(zip(self.vc_names + ["residual"], map(float, self.reml.theta)))

    def fitted_resid(self, n_rows: int) -> tuple[np.ndarray, np.ndarray]:
        """Conditional fitted values and residuals in the row order of Prepared.d."""
        f, r = self.reml.fitted_resid()
        order = np.concatenate(self.rows)
        fo, ro = np.empty(n_rows), np.empty(n_rows)
        fo[order], ro[order] = f, r
        return fo, ro

    def random_effects(self) -> dict[str, dict[str, np.ndarray]]:
        """{n: {term: BLUP vector}} (term 'n' = the n intercept)."""
        return {lab: dict(zip(self.vc_names, b)) for lab, b in zip(self.group_labels, self.reml.blups())}


def _class_groups(p: Prepared, cfg: Config) -> dict[str, list[int]]:
    """class -> treatment indices (in p.treatments) usable for a separate n x treatment variance."""
    out: dict[str, list[int]] = {}
    for i, t in enumerate(p.treatments):
        out.setdefault(cfg.classes.get(t, "treatments"), []).append(i)
    return out


def _vc(p: Prepared, cfg: Config, spec: str) -> tuple[dict[str, str], dict[str, str], list[str]]:
    vc, labels, notes = {}, {}, []
    groups = _class_groups(p, cfg)
    cells = p.d.groupby(["n", "treatment"]).size()

    def ncells(idx):
        return int(sum(1 for (n, t) in cells.index if t in {p.treatments[i] for i in idx}))

    split = spec == "full" and len(groups) >= 2 and all(ncells(ix) >= 3 for ix in groups.values())
    if spec == "full" and len(groups) >= 2 and not split:
        notes.append("a treatment class has fewer than 3 n×treatment cells: one common n×treatment variance used")
    if split:
        for cl, ix in groups.items():
            key = f"nt_{_safe(cl)}"
            vc[key] = "0 + " + " + ".join(f"I{i}" for i in ix)
            labels[key] = f"n × treatment [{cl}]"
    else:
        vc["nt_all"] = "0 + " + " + ".join(f"I{i}" for i in range(len(p.treatments)))
        labels["nt_all"] = "n × treatment [all treatments]"
    if p.can_well:
        vc["well"] = "0 + C(cell)"
        labels["well"] = "well (within n × treatment)"
    if p.can_position and spec != "no_position":
        vc["position"] = "0 + C(pos)"
        labels["position"] = "field position (within n)"
    return vc, labels, notes


def _statsmodels_start(model) -> np.ndarray | None:
    """One statsmodels REML fit, used only as an extra starting point (it may fail at a zero variance)."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = model.fit(reml=True, method=["lbfgs"], maxiter=300)
        vc = np.asarray(r.vcomp, float) if getattr(model, "k_vc", 0) else np.array([])
        return np.r_[float(np.asarray(r.cov_re)[0, 0]), vc, float(r.scale)]
    except Exception:  # noqa: BLE001
        return None


def _fit_spec(p: Prepared, cfg: Config, spec: str, starts: list | None = None,
              quick: bool = False) -> tuple[Fit | None, str, list[str]]:
    """Fit one model of the ladder. statsmodels MixedLM builds the design (vc_formula) and gives one starting
    point; reml.REML maximises the REML likelihood over theta >= 0 from several starts (best likelihood wins)."""
    from statsmodels.regression.mixed_linear_model import MixedLM

    from . import reml as reml_mod

    vc, labels, notes = _vc(p, cfg, spec)
    fixed = [f"T{i}" for i in range(1, len(p.treatments))] + [f"P{j}" for j in range(1, len(p.timepoint_levels))]
    formula = "y ~ " + (" + ".join(fixed) if fixed else "1")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = MixedLM.from_formula(formula, p.d.reset_index(drop=True), groups="n", re_formula="1", vc_formula=vc)
    except Exception as e:  # noqa: BLE001
        return None, f"design: {type(e).__name__}: {e}", notes
    vc_names = list(model.exog_vc.names) if getattr(model, "k_vc", 0) else []
    groups, rows = [], []
    for g, lab in enumerate(model.group_labels):
        M = [np.asarray(model.exog_re_li[g], float)] + [np.asarray(model.exog_vc.mats[j][g], float)
                                                         for j in range(len(vc_names))]
        groups.append(reml_mod.Group(np.asarray(model.endog_li[g], float), np.asarray(model.exog_li[g], float), M))
        rows.append(np.asarray(model.row_indices[lab]))
    X = np.vstack([g.X for g in groups])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None, "fixed effects not estimable (rank-deficient design)", notes
    eng = reml_mod.REML(groups, ["n"] + vc_names)
    extra = list(starts or [])
    method = "REML, bounded L-BFGS-B (multi-start)"
    if not quick:
        s0 = _statsmodels_start(model)
        if s0 is not None:
            extra.append(s0)
            method += " incl. statsmodels MixedLM start"
    try:
        rf = eng.fit(extra_starts=[s for s in extra if s is not None and len(s) == eng.K + 1],
                     n_default=1 if quick else None)
    except Exception as e:  # noqa: BLE001
        return None, f"{type(e).__name__}: {e}", notes
    if not rf.converged:
        return None, f"did not converge ({rf.message})", notes
    se = np.sqrt(np.clip(np.diag(rf.C), 0, None))
    if not (np.all(np.isfinite(rf.beta)) and np.all(np.isfinite(se)) and np.all(se > 0)):
        return None, "non-finite fixed-effect estimates", notes
    desc = dict(SPECS)[spec]
    if spec == "full" and not any(k.startswith("nt_") and k != "nt_all" for k in vc):
        desc = dict(SPECS)["common"]
    if not p.can_position or spec == "no_position":
        desc = desc.replace(" + position", "").replace(", no position", "")
    if not p.can_well:
        desc = desc.replace(" + well", "")
    labels = {"n": "n (experiment)"} | labels
    return Fit(spec, desc, formula + "  |  random: n" + "".join(f" + {labels[k]}" for k in vc), rf, rf.beta, rf.C,
               list(model.exog_names), ["n"] + vc_names, labels, [str(x) for x in model.group_labels], rows,
               True, notes, method,
               {k: list(model.exog_vc.colnames[j][0]) for j, k in enumerate(vc_names)}), "", notes


def _same_terms(p: Prepared, cfg: Config, a: str, b: str) -> bool:
    return _vc(p, cfg, a)[0].keys() == _vc(p, cfg, b)[0].keys()


def fit(p: Prepared, cfg: Config, start: str = "full", starts: list | None = None,
        quick: bool = False) -> tuple[Fit | None, list[str]]:
    """Walk the fallback ladder from `start`; returns (fit or None, log of failed attempts)."""
    log = []
    names = [s for s, _ in SPECS]
    tried = []
    for spec in names[names.index(start):]:
        if any(_same_terms(p, cfg, spec, t) for t in tried):
            continue
        tried.append(spec)
        f, err, _ = _fit_spec(p, cfg, spec, starts, quick)
        if f is not None:
            return f, log
        log.append(f"{dict(SPECS)[spec]}: failed ({err})")
    return None, log


# ---------------------------------------------------------------- contrasts

def df_cap(p: Prepared, f: Fit, cfg: Config, a: str, b: str) -> float:
    """Containment ceiling for the df of contrast a vs b: the degrees of freedom of the n x treatment variance(s)
    the comparison relies on. With class-specific variances that is the sum over the classes of a and b of
    (cells - treatments - n's + 1) of that class; with one common variance it is the whole n x treatment df.
    Never more than the RM ANOVA error df. Keeps the df honest when a class variance is estimated as zero."""
    split = [k for k in f.vc_names if k.startswith("nt_") and k != "nt_all"]
    if not split:
        return p.nt_df
    cap = 0.0
    for cl in dict.fromkeys([cfg.classes.get(a, "treatments"), cfg.classes.get(b, "treatments")]):
        members = {t for t in p.treatments if cfg.classes.get(t, "treatments") == cl}
        sub = p.d[p.d["treatment"].isin(members)]
        cells = sub.groupby(["n", "treatment"]).ngroups
        if len(members) >= 2:
            cap += max(1.0, float(cells - len(members) - sub["n"].nunique() + 1))
        else:  # a class of one treatment: its n x treatment variance rests on (n's with it) - 1 df
            cap += max(1.0, float(sub["n"].nunique() - 1))
    return float(min(cap, p.nt_df))


def _n_with(d: pd.DataFrame, *ts: str) -> int:
    sets = [set(d.loc[d["treatment"] == t, "n"]) for t in ts]
    return len(set.intersection(*sets)) if sets else 0


def _L(p: Prepared, f: Fit, t: str) -> np.ndarray:
    """Linear predictor of treatment t, averaged over timepoint levels."""
    L = np.zeros(len(f.names))
    L[f.names.index("Intercept")] = 1.0
    i = p.treatments.index(t)
    if i:
        L[f.names.index(f"T{i}")] = 1.0
    k = len(p.timepoint_levels)
    for j in range(1, k):
        L[f.names.index(f"P{j}")] = 1.0 / k
    return L


def _mvt_adjust(tvals: np.ndarray, dfs: np.ndarray, R: np.ndarray, simultaneous: bool = True):
    """Single-step multivariate-t adjusted p-values and 95% critical values (per comparison, using its own df)."""
    m = len(tvals)
    padj = np.full(m, np.nan)
    crit = np.full(m, np.nan)
    ok = np.isfinite(tvals) & (dfs >= 1)
    if not ok.any():
        return padj, crit
    idx = np.where(ok)[0]
    Rs = R[np.ix_(idx, idx)]
    k = len(idx)
    cache: dict[float, object] = {}

    def cover(c, df):
        if k == 1:
            return 1 - 2 * stats.t.sf(c, df)
        mv = cache.get(df)
        if mv is None:
            mv = cache[df] = stats.multivariate_t(shape=Rs, df=df, allow_singular=True)
        return float(mv.cdf(np.full(k, c), lower_limit=np.full(k, -c), random_state=0))

    crit_by_df: dict[float, float] = {}
    for j, i in enumerate(idx):
        df = round(float(dfs[i]), 1)
        padj[i] = min(1.0, max(0.0, 1 - cover(abs(tvals[i]), df)))
        if simultaneous:
            if df not in crit_by_df:
                from scipy.optimize import brentq
                lo, hi = stats.t.ppf(0.975, df), stats.t.ppf(1 - 0.025 / k, df) + 1e-6
                try:
                    crit_by_df[df] = lo if k == 1 else brentq(lambda c: cover(c, df) - 0.95, lo, hi, xtol=1e-4)
                except ValueError:
                    crit_by_df[df] = hi
            crit[i] = crit_by_df[df]
    return padj, crit


def holm(p) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    run = 0.0
    for rank, i in enumerate(ok[np.argsort(p[ok])]):
        run = max(run, min(1.0, (len(ok) - rank) * p[i]))
        out[i] = run
    return out


def contrasts(p: Prepared, f: Fit, cfg: Config, simultaneous: bool = True) -> pd.DataFrame:
    """One row per comparison in every family (model-based)."""
    rows = []
    et = p.scale.estimate_type
    lm = np.log(1 + cfg.margin)
    for fam in cfg.families(p.treatments):
        Ls, info = [], []
        for a, b in fam.pairs:
            L = _L(p, f, a) - _L(p, f, b)
            est = float(L @ f.beta)
            se = float(np.sqrt(L @ f.cov @ L))
            N = _n_with(p.d, a, b)
            df = f.df(L, df_cap(p, f, cfg, a, b)) if N >= 2 else 0.0
            Ls.append(L)
            info.append((a, b, est, se, N, df))
        if not Ls:
            continue
        C = np.array(Ls) @ f.cov @ np.array(Ls).T
        sd = np.sqrt(np.diag(C))
        R = C / np.outer(sd, sd)
        tv = np.array([e / s if s > 0 and df >= 1 else np.nan for _, _, e, s, _, df in info])
        dfs = np.array([df for *_, df in info], float)
        if fam.kind == "vehicle_check":
            for (a, b, est, se, N, df) in info:
                r = dict(family=fam.name, family_kind=fam.kind, family_id=fam.id, comparison=f"{a} vs {b}", group=a, comparator=b,
                         estimate_type=et, log_estimate=est, se=se, df=df, N=N,
                         test=f"TOST equivalence (mixed model; margin ±{cfg.margin:.0%}: ratio within "
                              f"[{1 / (1 + cfg.margin):.3g}, {1 + cfg.margin:.3g}]; 90% CI)")
                if df >= 1:
                    q = stats.t.ppf(0.95, df)
                    lo, hi = est - q * se, est + q * se
                    p1 = stats.t.sf((est + lm) / se, df)
                    p2 = stats.t.cdf((est - lm) / se, df)
                    ptost = float(max(p1, p2))
                    verdict = ("equivalent" if lo > -lm and hi < lm else
                               "not equivalent" if lo >= lm or hi <= -lm else "inconclusive")
                    r.update(estimate=np.exp(est), ci_low=np.exp(lo), ci_high=np.exp(hi), t=np.nan, p=ptost,
                             conclusion=verdict,
                             note=f"90% CI {np.exp(lo):.3g}-{np.exp(hi):.3g} vs margin {1 / (1 + cfg.margin):.3g}-"
                                  f"{1 + cfg.margin:.3g}: {verdict}")
                else:
                    r.update(estimate=np.exp(est), ci_low=np.nan, ci_high=np.nan, t=np.nan, p=np.nan,
                             conclusion="not estimable", note="needs at least 2 n's containing both groups")
                rows.append(r)
            ps = holm([r["p"] for r in rows if r["family"] == fam.name])
            for r, pa in zip([r for r in rows if r["family"] == fam.name], ps):
                r["p_adj"] = pa
            continue
        padj, crit = _mvt_adjust(tv, dfs, R, simultaneous)
        for j, (a, b, est, se, N, df) in enumerate(info):
            r = dict(family=fam.name, family_kind=fam.kind, family_id=fam.id, comparison=f"{a} vs {b}", group=a, comparator=b,
                     estimate_type=et, log_estimate=est, se=se, df=df, N=N, estimate=np.exp(est),
                     test=f"mixed model contrast (t, df = matched n's - 1; {fam.adjust} within family)")
            if df >= 1:
                q = stats.t.ppf(0.975, df)
                r.update(ci_low=np.exp(est - q * se), ci_high=np.exp(est + q * se), t=tv[j],
                         p=float(2 * stats.t.sf(abs(tv[j]), df)), p_adj=padj[j],
                         ci_low_adj=np.exp(est - crit[j] * se) if np.isfinite(crit[j]) else np.nan,
                         ci_high_adj=np.exp(est + crit[j] * se) if np.isfinite(crit[j]) else np.nan, note="")
            else:
                r.update(ci_low=np.nan, ci_high=np.nan, t=np.nan, p=np.nan, p_adj=np.nan, ci_low_adj=np.nan,
                         ci_high_adj=np.nan,
                         note="not estimable: fewer than 2 n's contain both groups (estimate shown, no test)")
            rows.append(r)
    cols = ["family", "family_kind", "family_id", "comparison", "group", "comparator", "estimate_type", "estimate", "ci_low", "ci_high",
            "ci_low_adj", "ci_high_adj", "log_estimate", "se", "t", "df", "N", "p", "p_adj", "significance",
            "conclusion", "test", "note"]
    out = pd.DataFrame(rows).reindex(columns=cols)
    if len(out):
        out["significance"] = [r.conclusion if r.family_kind == "vehicle_check" else stars(r.p_adj)
                               for r in out.itertuples()]
        if et == "odds_ratio":
            out["test"] = out["test"].str.replace("ratio within", "odds ratio within")
    return out


def emm(p: Prepared, f: Fit) -> pd.DataFrame:
    """Estimated marginal means per treatment, back-transformed to absolute units, with 95% CI."""
    rows = []
    for t in p.treatments:
        L = _L(p, f, t)
        mu = float(L @ f.beta)
        se = float(np.sqrt(L @ f.cov @ L))
        N = _n_with(p.d, t)
        df = f.df(L, N - 1) if N >= 2 else 0.0
        q = stats.t.ppf(0.975, df) if df >= 1 else np.nan
        rows.append(dict(treatment=t, N=N, fields=int((p.d["treatment"] == t).sum()),
                         estimate=float(p.scale.back(mu)),
                         ci_low=float(p.scale.back(mu - q * se)) if df >= 1 else np.nan,
                         ci_high=float(p.scale.back(mu + q * se)) if df >= 1 else np.nan,
                         log_mean=mu, se=se, df=df,
                         estimate_type="mean",
                         mean_type="proportion (back-transformed logit)" if p.scale.kind == "fraction" else
                         "geometric mean" + (f" (minus offset {p.scale.shift:g})" if p.scale.shift else "")))
    return pd.DataFrame(rows)


def variance_components(p: Prepared, f: Fit | None, cfg: Config) -> pd.DataFrame:
    """Variance components on the analysis scale, with approximate Wald CIs on the log-variance scale."""
    allc = [("n", "n (experiment: day, plate, passage)")]
    groups = _class_groups(p, cfg)
    allc += [(f"nt_{_safe(c)}", f"n × treatment [{c}]") for c in groups]
    allc += [("nt_all", "n × treatment [all treatments]"), ("well", "well (within n × treatment)"),
             ("position", "field position (within n)"), ("residual", "field (residual)")]
    rows = []
    if f is None:
        return pd.DataFrame([dict(component=lab, variance=np.nan, note="not estimable (no model fitted)") for _, lab in allc])
    vals = f.theta
    ses = {}
    try:  # SEs from the inverse expected REML information (components at zero excluded)
        keys = f.vc_names + ["residual"]
        keep, cov = f.reml.theta_cov()
        for k, sd_ in zip([k for k, kp in zip(keys, keep) if kp], np.sqrt(np.clip(np.diag(cov), 0, None))):
            ses[k] = float(sd_)
    except Exception:  # noqa: BLE001
        pass
    total = sum(v for v in vals.values() if np.isfinite(v))
    for key, lab in allc:
        if key not in vals:
            if key in ("nt_all",) and any(k.startswith("nt_") for k in vals):
                continue
            if key.startswith("nt_") and key != "nt_all" and "nt_all" in vals:
                continue
            why = {"well": "not estimable (one well per condition)" if not p.can_well else "not in the fitted model",
                   "position": ("not modelled (positions random)" if not p.can_position else "dropped (fallback model)")}
            rows.append(dict(component=lab, variance=np.nan, sd=np.nan, cv_pct=np.nan, pct_of_total=np.nan,
                             ci_low=np.nan, ci_high=np.nan, note=why.get(key, "not estimable")))
            continue
        v = vals[key]
        se = ses.get(key, np.nan)
        if key == "residual":
            df_res = max(len(p.d) - len(f.beta), 1)
            lo, hi = v * df_res / stats.chi2.ppf(0.975, df_res), v * df_res / stats.chi2.ppf(0.025, df_res)
            note = "chi-square CI (residual df)"
        elif v > 1e-8 and np.isfinite(se) and se > 1.5 * v:
            lo, hi = np.nan, np.nan
            note = f"too uncertain for a CI (SE {se:.2g} > 1.5 x estimate; few levels)"
        elif v > 1e-8 and np.isfinite(se) and se > 0:
            z = 1.96 * se / v
            lo, hi = v * np.exp(-z), v * np.exp(z)
            note = "approximate Wald CI on log-variance scale (rough with few levels)"
        else:
            lo, hi = np.nan, np.nan
            note = "estimate at or near zero: CI not available" if v <= 1e-8 else "CI not available"
        rows.append(dict(component=lab, variance=v, sd=np.sqrt(max(v, 0)),
                         cv_pct=100 * np.sqrt(np.expm1(max(v, 0))) if p.scale.kind != "fraction" else np.nan,
                         pct_of_total=100 * v / total if total > 0 else np.nan, ci_low=lo, ci_high=hi, note=note))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- per-n effects and RM ANOVA paths

def per_n_ratios(p: Prepared, cfg: Config) -> pd.DataFrame:
    """Ratio of each treatment to each control within each n: model scale (geometric) and ratio of arithmetic means."""
    d = p.d
    g = d.groupby(["n", "treatment"])
    mu, am = g["y"].mean(), g["value"].mean()
    rows = []
    for n in p.ns:
        for c in cfg.controls:
            if (n, c) not in mu.index:
                continue
            for t in p.treatments:
                if t == c or (n, t) not in mu.index:
                    continue
                rows.append(dict(n=n, treatment=t, control=c, ratio=float(np.exp(mu[(n, t)] - mu[(n, c)])),
                                 ratio_of_means=float(am[(n, t)] / am[(n, c)]) if am[(n, c)] else np.nan,
                                 fields=int(g.size()[(n, t)]),
                                 estimate_type=p.scale.estimate_type))
    return pd.DataFrame(rows, columns=["n", "treatment", "control", "ratio", "ratio_of_means", "fields", "estimate_type"])


def rm_anova_contrasts(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """The Prism-standard RM one-way ANOVA + Dunnett (raw scale, per-n means) as model-like rows (vs each control)."""
    from .analysis import rm_anova_dunnett

    means = df.groupby(["n", "treatment"])["value"].mean().unstack("treatment")
    order = [t for t in cfg.order if t in means]
    rows = []
    for c in cfg.controls:
        if c not in means:
            continue
        _, dun = rm_anova_dunnett(means, order, c)
        for r in dun:
            rows.append(dict(family=f"vs {c}", family_kind="vs_control", comparison=r["comparison"],
                             group=r["comparison"].split(" vs ")[0], comparator=c, estimate_type="difference",
                             estimate=r["effect"], ci_low=r["ci_low"], ci_high=r["ci_high"], N=r["N"], p=r["p"],
                             p_adj=r["p_adj"], significance=stars(r["p_adj"]), test=r["test"], note=""))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- top level

@dataclass
class ModelResult:
    metric: Metric
    status: str  # ok | fallback_rm_anova | not_estimable
    model: str
    spec: str
    formula: str
    converged: bool
    df_rule: str
    notes: list[str]
    attempts: list[str]
    scale: Scale | None = None
    prepared: Prepared | None = None
    fit: Fit | None = None
    emm: pd.DataFrame = field(default_factory=pd.DataFrame)
    contrasts: pd.DataFrame = field(default_factory=pd.DataFrame)
    varcomp: pd.DataFrame = field(default_factory=pd.DataFrame)
    per_n: pd.DataFrame = field(default_factory=pd.DataFrame)
    loo: pd.DataFrame = field(default_factory=pd.DataFrame)
    loo_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    primary_test: str = "mixed"
    order: list = field(default_factory=list)  # treatment display order (settings)

    @property
    def ok(self) -> bool:
        return self.fit is not None

    def status_dict(self) -> dict:
        return dict(status=self.status, model=self.model, spec=self.spec, formula=self.formula,
                    converged=bool(self.converged), df_rule=self.df_rule,
                    response=self.scale.describe if self.scale else "", primary_test_used=self.primary_test,
                    notes=list(self.notes), attempts=list(self.attempts),
                    warnings=list(self.fit.warnings) if self.fit else [],
                    optimiser=self.fit.method if self.fit else "")


DF_RULE = ("t with Satterthwaite df, capped at the df of the n×treatment variance(s) the comparison relies on "
           "(containment: cells - treatments - n's + 1 per class involved); means: capped at (n's with the treatment) - 1; "
           "a comparison needs >= 2 n's containing both groups")


def analyse(df: pd.DataFrame, metric: Metric, settings: dict, positions_ok: bool | None = None,
            leave_one_out: bool = True) -> ModelResult:
    """Fit the hierarchical model for one metric and derive every model-based table. Never raises."""
    try:
        res = _analyse(df, metric, settings, positions_ok, leave_one_out)
        res.order = config(settings, df["treatment"].unique()).order
        return res
    except Exception as e:  # noqa: BLE001 - a generation build must never crash on the statistics
        cfg = config(settings, df["treatment"].unique()) if "treatment" in df else None
        return ModelResult(metric, "fallback_rm_anova", "none (error)", "", "", False, DF_RULE,
                           [f"model analysis failed ({type(e).__name__}: {e}); the RM ANOVA drives the stars"], [],
                           primary_test="rm_anova" if cfg else "rm_anova")


def _analyse(df, metric, settings, positions_ok, leave_one_out) -> ModelResult:
    from . import design

    if positions_ok is None:
        try:
            positions_ok = design.positions_meaningful(settings)
        except ValueError:
            positions_ok = False
    cfg = config(settings, df["treatment"].unique())
    p = prepare(df, metric, cfg, positions_ok)
    notes = list(p.notes)
    if len(p.ns) < 2 or len(p.treatments) < 2:
        why = "only one n" if len(p.ns) < 2 else "fewer than two treatments"
        res = ModelResult(metric, "not_estimable", f"not estimable ({why})", "", "", False, DF_RULE,
                          notes + [f"mixed model not estimable: {why}"], [], p.scale, p,
                          primary_test="none")
        res.per_n = per_n_ratios(p, cfg)
        res.varcomp = variance_components(p, None, cfg)
        return res
    f, attempts = fit(p, cfg)
    if f is None:
        res = ModelResult(metric, "fallback_rm_anova", "none (all mixed models failed)", "", "", False, DF_RULE,
                          notes + ["every mixed model failed; the RM ANOVA (Prism standard) drives the stars"],
                          attempts, p.scale, p, primary_test="rm_anova")
        res.per_n = per_n_ratios(p, cfg)
        res.varcomp = variance_components(p, None, cfg)
        if leave_one_out:
            res.loo, res.loo_summary = loo(df, metric, settings, positions_ok, use_model=False)
        return res
    vcs = variance_components(p, f, cfg)
    zero = vcs[(vcs["variance"] <= 1e-8) & vcs["variance"].notna()]["component"].tolist()
    if zero:
        notes.append("estimated as zero (boundary): " + ", ".join(zero))
    res = ModelResult(metric, "ok", f.description, f.spec, f.formula, f.converged, DF_RULE, notes, attempts, p.scale, p, f,
                      emm(p, f), contrasts(p, f, cfg), vcs, per_n_ratios(p, cfg),
                      primary_test=cfg.primary_test)
    if leave_one_out:
        res.loo, res.loo_summary = loo(df, metric, settings, positions_ok, use_model=cfg.primary_test == "mixed",
                                       start=f.spec, full=res.contrasts, theta0=f.reml.theta)
    return res


def loo(df: pd.DataFrame, metric: Metric, settings: dict, positions_ok: bool, use_model: bool = True,
        start: str = "full", full: pd.DataFrame | None = None, theta0=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Leave-one-n-out with the primary test: estimate, CI, p, p_adj with each n dropped, and a stability summary."""
    ns = sorted(df["n"].astype(str).unique())
    cfg = config(settings, df["treatment"].unique())
    rows = []
    if full is None or not use_model:
        full = rm_anova_contrasts(df, cfg) if not use_model else full
    for n in ns:
        sub = df[df["n"].astype(str) != n]
        tab = pd.DataFrame()
        note = ""
        if sub["n"].nunique() >= 2:
            if use_model:
                p = prepare(sub, metric, cfg, positions_ok)
                f, _ = fit(p, cfg, start=start, starts=[theta0] if theta0 is not None else None, quick=True)
                if f is not None:
                    tab = contrasts(p, f, cfg, simultaneous=False)
                    tab = tab[tab["family_kind"] != "vehicle_check"]
                else:
                    note = "model failed"
            else:
                tab = rm_anova_contrasts(sub, cfg)
        else:
            note = "fewer than 2 n's left"
        if tab.empty:
            rows.append(dict(dropped_n=n, family="", comparison="", note=note or "not estimable"))
            continue
        for r in tab.itertuples():
            rows.append(dict(dropped_n=n, family=r.family, comparison=r.comparison, estimate_type=r.estimate_type,
                             estimate=r.estimate, ci_low=r.ci_low, ci_high=r.ci_high, p=r.p, p_adj=r.p_adj,
                             N=r.N, note=getattr(r, "note", "") or ""))
    out = pd.DataFrame(rows, columns=["dropped_n", "family", "comparison", "estimate_type", "estimate", "ci_low", "ci_high",
                                      "p", "p_adj", "N", "note"])
    summ = []
    if full is not None and len(full):
        fl = full[full["family_kind"] != "vehicle_check"] if "family_kind" in full else full
        null = 0.0 if use_model is False else 1.0
        for r in fl.itertuples():
            drops = out[(out["family"] == r.family) & (out["comparison"] == r.comparison)]
            sig_full = bool(np.isfinite(r.p_adj) and r.p_adj < ALPHA)
            sign_full = np.sign((r.estimate - null) if np.isfinite(r.estimate) else np.nan)
            ests = drops["estimate"].to_numpy(float)
            pads = drops["p_adj"].to_numpy(float)
            same_dir = bool(len(ests) and np.all(np.sign(ests - null) == sign_full))
            same_sig = bool(len(pads) and np.all((pads < ALPHA) == sig_full) and np.all(np.isfinite(pads)))
            missing = len(ns) - len(drops)
            summ.append(dict(family=r.family, comparison=r.comparison, estimate=r.estimate, p_adj=r.p_adj,
                             significant=sig_full, loo_min=np.nanmin(ests) if len(ests) and np.isfinite(ests).any() else np.nan,
                             loo_max=np.nanmax(ests) if len(ests) and np.isfinite(ests).any() else np.nan,
                             loo_p_adj_max=np.nanmax(pads) if len(pads) and np.isfinite(pads).any() else np.nan,
                             drops=len(drops), same_direction=same_dir, same_significance=same_sig,
                             stable="yes" if (same_dir and same_sig and missing == 0) else "no",
                             note="" if missing == 0 else f"{missing} drop(s) not estimable (comparison needs that n)"))
    return out, pd.DataFrame(summ)
