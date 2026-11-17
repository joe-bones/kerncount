"""Hierarchical (mixed) model across all n, per metric; comparison families; equivalence; robustness.

Unit of observation ([analysis] model_unit):
    "well" (default)  one row per well. Biologists mean "ratio of mean nuclei per field", so each well contributes the
                      arithmetic mean of its fields (count, positive metrics) or its pooled fraction
                      sum(positives) / sum(nuclei) (fraction metrics). Every well is imaged at the same positions, so
                      position effects cancel in well means. A well with missing / excluded fields gets a
                      position-adjusted mean: sum of its fields / sum of the n-specific position pattern at the
                      positions it has (pattern = mean of field / well mean over the complete wells of that n,
                      normalised to mean 1), a ratio estimator that equals the plain mean for complete wells.
    "field"           one row per field (the earlier model): compares geometric means of fields, which sparse, uneven
                      fields drag down; kept as an option.

Response scale (multiplicative effects -> log scale, results reported as ratios):
    count     y = log(mean + 0.5)        (field unit: log(count + 0.5)); 0.5 is negligible for ~10 nuclei or more
    positive  y = log(mean)              (with zeros: log(value + half the smallest positive value), noted)
    fraction  y = empirical logit log((K + 0.5) / (N - K + 0.5)), K = positives, N = nuclei pooled over the well's
              fields. Binomial-type model with overdispersion: each well has a known binomial sampling variance
              1/(K + 0.5) + 1/(N - K + 0.5) (the variance of the empirical logit) added to the estimated well-to-well
              variance, so wells with few nuclei count less, but a well with many nuclei cannot dominate when the
              biology (well, n, n x treatment variances) is the main source of noise. (A plain nuclei weight would let
              big wells dominate; a pure binomial model would ignore overdispersion.)
              Headline effect: difference in percentage points (back-transformed model draws, see contrasts()),
              with odds ratios and % of reference as secondary columns.

Model (REML; groups = n, everything else nested in n; design built by statsmodels MixedLM with vc_formula,
likelihood maximised over variances >= 0 by reml.py):
    y = treatment (fixed) [+ timepoint (fixed, >= 2 known timepoints)] [+ density covariate (adjusted model only)]
        + n                                    random intercept per experiment
        + n x treatment                        random, a separate variance per treatment class
        + (field unit) well within n x treatment, field position within n (random n:position)
        + residual                             well-to-well (well unit) or field (field unit)
    Classes ([analysis] classes): "auto" (default) = one class per treatment, each control included.
    Pooling n x treatment variance across treatments gave a type I error of 0.16 for a treatment that varies more
    than the rest of its class, and pooling the controls gave 0.13 for "vs LyoVec" when the vehicle varies between
    n's but Control does not; one class per treatment is the conservative default (docs/statistics.md).
    A treatment present in only one n has no n x treatment term (its single cell is absorbed by its fixed effect).
    With one well per condition (well unit) the residual is not identifiable next to n x treatment; it is held at ~0.

Fallbacks (first that fits wins, recorded in `model`): full -> common n x treatment variance -> without position.
If every model fails the caller falls back to the RM ANOVA (and says so).

Small-sample inference: t tests with Satterthwaite df (expected REML information), capped by a containment ceiling =
the df of the n x treatment variance(s) the comparison relies on (multi-treatment class: cells - treatments - n's + 1;
one-treatment class: n's - 1), never above the RM ANOVA error df. A comparison needs >= 2 n's containing both groups.

One n (or a comparison with fewer than 2 n's containing both groups): estimates are filled in from the data
(descriptive), every across-n CI, p and star is empty with the note "needs ≥ 2 n". Field / well variation is never
used to make a single experiment look like a biological test.

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

from .metrics import Metric, estimand

COUNT_OFFSET = 0.5
FRACTION_CLIP = 0.01
ALPHA = 0.05
NEEDS_2N = "needs ≥ 2 n"
NEEDS_2N_PAIR = "needs ≥ 2 n containing both groups: not estimable (estimate shown, no test)"
OVERLAP_MIN = 0.5  # density overlap: each condition needs >= half of its wells inside the other's density range
DENSITY_FLAG = "conditions do not overlap in density; adjusted estimate not trustworthy"
RAW_SCALE_NOTE = ("RAW SCALE (difference): some well means are <= 0 (background-subtracted intensity), so the log is "
                  "impossible; analysed as a DIFFERENCE in intensity units (identity link), not a ratio")


def needs_raw_scale(df: pd.DataFrame, metric: Metric, unit: str = "well", positions_ok: bool = True) -> int:
    """Number of analysed values <= 0 of a positive-kind metric (well means at the well unit, field values < 0 at
    the field unit, where zero fields keep the log(value + shift) rule). > 0 -> the metric is analysed as a
    difference on the raw scale: dropping those wells (log of a non-positive value) biases the effect toward 0."""
    if metric.kind != "positive" or df.empty:
        return 0
    d = _normalise(df)
    if unit == "well":
        w, _ = well_means(d[d["value"].notna()], "value", positions_ok)
        return int((w["mean"] <= 0).sum())
    return int((d["value"] < 0).sum())


def raw_scale_metric(metric: Metric) -> Metric:
    from dataclasses import replace

    return replace(metric, kind="difference")


def stars(p) -> str:
    if p is None or not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def expit(x):
    from scipy.special import expit as _expit

    return _expit(np.asarray(x, float))


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
    unit: str = "well"  # well | field
    density: str = ""  # density covariate column for fraction metrics ("" = off)
    classes_mode: str = "auto"  # auto | two (controls vs the rest) | user
    margin_pp: float = 5.0  # equivalence margin of fraction metrics, in percentage points

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
    """[analysis] classes: a TOML table {name = ["pattern", ...]} (or its text form); "auto" -> {}."""
    if isinstance(v, str):
        v = v.strip()
        if not v or v.lower() == "auto":
            return {}
        import tomllib
        v = tomllib.loads(f"x = {v}")["x"]
    if not isinstance(v, dict):
        return {}
    return {str(k): _as_list(p) for k, p in v.items()}


def classes_mode(v) -> str:
    if v is None or (isinstance(v, str) and v.strip().lower() == "auto"):
        return "auto"
    return "user" if parse_classes(v) else "two"


def config(settings: dict, treatments) -> Config:
    from .analysis import control_groups, treatment_order

    a = settings["analysis"]
    order = treatment_order(list(treatments), a.get("treatment_order", []))
    controls = control_groups(a)
    primary = str(a.get("primary_family") or "").strip() or a["reference"]
    raw = a.get("classes", "auto")
    mode = classes_mode(raw)
    user = parse_classes(raw) if mode == "user" else {}
    cls, corder = {}, []
    for t in order:
        if mode == "auto":
            hit = t  # every treatment, each control included, carries its own n x treatment variance
        else:
            hit = next((name for name, pats in user.items() if any(fnmatch.fnmatchcase(t, p) for p in pats)), None)
            if hit is None:
                hit = "controls" if t in controls else ("other" if user else "treatments")
        cls[t] = hit
        if hit not in corder:
            corder.append(hit)
    vc = str(a.get("vehicle_check", "auto")).strip()
    if vc.lower() == "auto":
        pairs = [(c, controls[0]) for c in controls[1:]] if controls else []
    elif vc.lower() in ("", "none", "off"):
        pairs = []
    else:
        pairs = []
        for part in vc.split(";"):
            bits = re.split(r"\s+vs\.?\s+", part.strip())
            if len(bits) == 2 and all(bits):
                pairs.append((bits[0].strip(), bits[1].strip()))
    pt = str(a.get("primary_test", "mixed")).strip().lower()
    unit = str(a.get("model_unit", "well") or "well").strip().lower()
    return Config(controls=controls, primary_ref=primary, primary_test=pt if pt in ("mixed", "rm_anova") else "mixed",
                  classes=cls, class_order=corder, vehicle_pairs=pairs,
                  positive_control=str(a.get("positive_control") or "").strip(),
                  among_classes=_as_list(a.get("among_classes", [])),
                  margin=float(a.get("equivalence_margin", 0.2) or 0.2), order=order,
                  unit=unit if unit in ("well", "field") else "well",
                  density=str(a.get("density_covariate", "") or "").strip(), classes_mode=mode,
                  margin_pp=float(a.get("equivalence_margin_pp", 5.0) or 5.0))


# ---------------------------------------------------------------- response scale

@dataclass
class Scale:
    kind: str
    shift: float = 0.0  # added before the log (count offset, or zero fix for positive metrics)
    note: str = ""
    unit: str = "field"  # well | field: what one analysed value is

    def back(self, mu):
        mu = np.asarray(mu, float)
        if self.kind == "fraction":
            return expit(mu)
        if self.kind == "difference":
            return mu
        return np.maximum(np.exp(mu) - self.shift, 0.0)  # an offset never makes a mean (or its CI) negative

    @property
    def estimate_type(self) -> str:
        return {"fraction": "odds_ratio", "difference": "difference"}.get(self.kind, "ratio")

    @property
    def describe(self) -> str:
        what = "well mean" if self.unit == "well" else ("count" if self.kind == "count" else "value")
        if self.kind == "count":
            return f"log({what} + {self.shift:g})"
        if self.kind == "fraction":
            return "logit of the pooled well fraction (empirical)" if self.unit == "well" else "logit (empirical)"
        if self.kind == "difference":
            return f"{what} on the raw scale (identity link; some {self.unit} means <= 0)"
        return f"log({what})" if not self.shift else f"log({what} + {self.shift:.3g})"


def transform(values: pd.Series, kind: str, counts: pd.Series | None = None) -> tuple[pd.Series, Scale]:
    """Response on the analysis scale; rows that cannot be transformed become NaN. For fractions `counts` is the
    number of nuclei N the fraction is taken of (k = fraction x N positives)."""
    v = pd.to_numeric(values, errors="coerce").astype(float)
    if kind == "difference":
        return v, Scale("difference", note=RAW_SCALE_NOTE)
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
            return y, Scale("fraction", note="y = empirical logit log((k + 0.5) / (N - k + 0.5)), k positives of N nuclei")
        pc = p.clip(FRACTION_CLIP, 1 - FRACTION_CLIP)
        return np.log(pc / (1 - pc)), Scale("fraction", note=f"y = logit(fraction clipped to [{FRACTION_CLIP}, {1 - FRACTION_CLIP}])")
    v = v.where(v >= 0)
    pos = v[v > 0]
    if (v == 0).any() and len(pos):
        shift = float(pos.min()) / 2
        return np.log(v + shift), Scale("positive", shift, f"y = log(value + {shift:.3g}); {int((v == 0).sum())} zero value(s)")
    neg = int((pd.to_numeric(values, errors="coerce") < 0).sum())
    return np.log(v.where(v > 0)), Scale("positive", 0.0, "y = log(value)" + (f"; {neg} negative value(s) left out" if neg else ""))


def nuclei(fields: pd.DataFrame, metric: Metric) -> pd.Series | None:
    """Nuclei N each field's fraction is taken of (None when not available): the metric's own N column
    (Metric.n_column, e.g. n_measured_Ki67) when set and present, else from the count (Metric.denominator)."""
    ncol = getattr(metric, "n_column", "")
    if ncol and ncol in fields:
        return pd.to_numeric(fields[ncol], errors="coerce").astype(float)
    if "count" not in fields:
        return None
    c = pd.to_numeric(fields["count"], errors="coerce").astype(float)
    if metric.denominator == "count_plus_fraction":
        f = pd.to_numeric(fields["value"], errors="coerce").astype(float)
        return (c / (1 - f)).where(f < 1)
    return c


# ---------------------------------------------------------------- well aggregation

def _safe(s: str) -> str:
    return re.sub(r"\W", "_", str(s))


def well_means(f: pd.DataFrame, col: str, positions_ok: bool) -> tuple[pd.DataFrame, list[str]]:
    """Mean per field of each well (n, treatment, well) of column `col`, position-adjusted where fields are missing.

    A well missing some of its n's positions gets sum(fields) / sum(r_p over its positions), where r_p is the n's
    position pattern (mean of field / well mean over that n's complete wells, normalised to mean 1). For a complete
    well this is the plain mean. Returns columns n, treatment, well, mean, fields, adjusted."""
    f = f[["n", "treatment", "well", "position", col]].copy()
    f[col] = pd.to_numeric(f[col], errors="coerce")
    f = f[f[col].notna()]
    keys = ["n", "treatment", "well"]
    out = f.groupby(keys)[col].agg(mean="mean", fields="size", var_mean="var").reset_index()
    out["var_mean"] = out["var_mean"] / out["fields"]  # variance of the well mean (from its fields)
    out["adjusted"] = False
    notes = []
    if not positions_ok or f.empty:
        return out, notes
    n_adj, n_noref = 0, 0
    for n, fn in f.groupby("n"):
        P = sorted(fn["position"].unique())
        if len(P) < 2:
            continue
        by_well = {k: w.groupby("position")[col].mean() for k, w in fn.groupby(["treatment", "well"])}
        complete = {k: w for k, w in by_well.items() if len(w) == len(P) and w.mean() > 0}
        incomplete = [k for k, w in by_well.items() if len(w) < len(P)]
        if not incomplete:
            continue
        if not complete:
            n_noref += len(incomplete)
            continue
        r = pd.concat([w / w.mean() for w in complete.values()], axis=1).mean(axis=1).reindex(P)
        r = r / r.mean()
        for (t, w) in incomplete:
            x = by_well[(t, w)]
            den = float(r.reindex(x.index).sum())
            if den > 0:
                m = (out["n"] == n) & (out["treatment"] == t) & (out["well"] == w)
                out.loc[m, "mean"] = float(x.sum()) / den
                out.loc[m, "adjusted"] = True
                n_adj += 1
    if n_adj:
        notes.append(f"{n_adj} well(s) with missing fields: position-adjusted well mean (ratio estimator with the "
                     "n's position pattern from its complete wells)")
    if n_noref:
        notes.append(f"{n_noref} incomplete well(s) in n's without a complete well: plain mean of the fields present")
    return out, notes


def aggregate_wells(d: pd.DataFrame, metric: Metric, positions_ok: bool, covariate: bool = False
                    ) -> tuple[pd.DataFrame, list[str], bool]:
    """Field rows -> one row per well: value (mean per field, or pooled fraction), fields, adjusted, K and nuclei
    (fractions), timepoint, cell_line and (covariate=True) zraw = well mean of column `cov`.
    Returns (wells, notes, has_nuclei)."""
    keys = ["n", "treatment", "well"]
    d = d[d["value"].notna()].copy()
    notes: list[str] = []
    has_n = False
    meta = d.groupby(keys).agg(timepoint=("timepoint", "first"), fields_total=("value", "size")).reset_index()
    if metric.kind == "fraction":
        N = nuclei(d, metric)
        if N is not None and (N > 0).any():
            has_n = True
            d["_N"] = N.where(N > 0)
            d["_K"] = d["value"] * d["_N"]
            g = d[d["_N"].notna()].groupby(keys)
            w = g.agg(K=("_K", "sum"), nuclei=("_N", "sum"), fields=("value", "size")).reset_index()
            w["value"] = (w["K"] / w["nuclei"]).where(w["nuclei"] > 0)
            notes.append("fraction pooled per well: sum of positives / sum of nuclei over its fields "
                         + (f"(nuclei = {metric.n_column})" if metric.n_column and metric.n_column in d
                            else "(nuclei = counted + the fraction's own nuclei)" if metric.denominator == "count_plus_fraction"
                            else "(nuclei = counted nuclei)") + "; fields are not position-adjusted for fractions")
        else:
            w = d.groupby(keys)["value"].agg(value="mean", fields="size").reset_index()
            notes.append("fraction: nuclei per field unknown, well value = mean of the field fractions (no weights)")
        w["adjusted"] = False
    else:
        w, wn = well_means(d, "value", positions_ok)
        w = w.rename(columns={"mean": "value"})
        notes += wn
    w = w.merge(meta, on=keys, how="left")
    if covariate:
        z, _ = well_means(d.assign(cov=pd.to_numeric(d.get("cov"), errors="coerce")), "cov", positions_ok)
        w = w.merge(z[keys + ["mean"]].rename(columns={"mean": "zraw"}), on=keys, how="left")
    for c in ("cell_line",):
        if c in d:
            w = w.merge(d.groupby(keys)[c].first().reset_index(), on=keys, how="left")
    return w, notes, has_n


# ---------------------------------------------------------------- data preparation

@dataclass
class Prepared:
    d: pd.DataFrame  # analysis rows (wells or fields) with y
    scale: Scale
    treatments: list[str]  # present, in order; first = baseline of the fixed effects
    ns: list[str]
    timepoint_levels: list[str]
    can_well: bool
    can_position: bool
    notes: list[str] = field(default_factory=list)
    nt_df: float = 1.0  # df of the n x treatment interaction (the RM ANOVA error df): ceiling for contrast df
    unit: str = "well"
    fields: pd.DataFrame | None = None  # field-level rows (value, count, ...) for plots and per-n descriptions
    replicated: bool = True  # some n x treatment cell has >= 2 rows (else the residual is not identifiable)
    known: bool = False  # d["known"] = known sampling variance per row (fractions, well unit)
    covariate: str = ""  # density covariate column (adjusted model): d["Z"] = centred, d["zraw"] = raw
    z_center: float = 0.0
    metric: Metric | None = None
    has_nuclei: bool = False


def _normalise(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["n"] = d["n"].astype(str)
    d["treatment"] = d["treatment"].astype(str)
    for c, v in (("well", 1), ("position", 1)):
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(v).astype(int) if c in d else v
    d["timepoint"] = d["timepoint"].fillna("").astype(str) if "timepoint" in d else ""
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    return d


def prepare(df: pd.DataFrame, metric: Metric, cfg: Config, positions_ok: bool, unit: str | None = None,
            covariate: str = "") -> Prepared:
    unit = unit or cfg.unit
    d = _normalise(df)
    fields = d.copy()
    notes: list[str] = []
    if covariate:
        d["cov"] = pd.to_numeric(d[covariate], errors="coerce") if covariate in d else np.nan
    has_n = False
    if unit == "well":
        d, wn, has_n = aggregate_wells(d, metric, positions_ok, bool(covariate))
        notes.append("unit: well (mean per field of each well" + ("; pooled fraction" if metric.kind == "fraction" else "")
                     + ")")
        notes += wn
        counts = d["nuclei"] if (metric.kind == "fraction" and has_n) else None
    else:
        counts = nuclei(d, metric) if metric.kind == "fraction" else None
        has_n = counts is not None and bool((counts > 0).any())
        if covariate:
            d["zraw"] = d["cov"]
    y, scale = transform(d["value"], metric.kind, counts)
    scale.unit = unit
    d["y"] = y
    if scale.note:
        notes.insert(0, scale.note)
    known = False
    if unit == "well" and metric.kind == "fraction" and has_n:
        K, N = d["value"] * d["nuclei"], d["nuclei"]
        d["known"] = 1 / (K + 0.5) + 1 / (N - K + 0.5)
        known = True
    bad = d["y"].isna() | ~np.isfinite(d["y"])
    if covariate:
        bad |= d["zraw"].isna() | ~np.isfinite(d["zraw"])
    if bad.any():
        notes.append(f"{int(bad.sum())} {unit}(s) without a usable value" + (" or density" if covariate else "") + " left out")
    d = d[~bad].copy().reset_index(drop=True)
    zc = 0.0
    if covariate:
        zc = float(d["zraw"].mean()) if len(d) else 0.0
        d["Z"] = d["zraw"] - zc
    present = [t for t in cfg.order if t in set(d["treatment"])]
    base = cfg.primary_ref if cfg.primary_ref in present else next((c for c in cfg.controls if c in present),
                                                                    present[0] if present else None)
    treatments = ([base] + [t for t in present if t != base]) if base else []
    ns = sorted(d["n"].unique(), key=str)
    d["cell"] = d["treatment"].map(_safe) + "_w" + d["well"].astype(str)
    if "position" not in d:
        d["position"] = 1
    d["pos"] = "p" + d["position"].astype(str)
    for i, t in enumerate(treatments):
        d[f"I{i}"] = (d["treatment"] == t).astype(float)
        if i:
            d[f"T{i}"] = d[f"I{i}"]
    rows_per_cell = d.groupby(["n", "treatment"]).size()
    replicated = bool((rows_per_cell >= 2).any())
    if unit == "field":
        wells_per_cell = d.groupby(["n", "treatment"])["well"].nunique()
        can_well = bool((wells_per_cell >= 2).any())
        can_pos = bool(positions_ok and d["position"].nunique() >= 2)
        if not positions_ok:
            notes.append("field positions are random ([design] field_positions): position not modelled")
    else:
        can_well = can_pos = False
        if not replicated and len(d):
            notes.append("one well per condition: the well-to-well residual is not identifiable and is absorbed by "
                         "n × treatment")
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
    return Prepared(d, scale, treatments, ns, timepoint_levels, can_well, can_pos, notes, nt_df, unit, fields,
                    replicated, known, covariate, zc, metric, has_n)


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
    names: list[str]  # fixed-effect names (Intercept, T1, ..., P1, ..., Z)
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


def _cells_by_treatment(p: Prepared) -> pd.Series:
    return p.d.groupby("treatment")["n"].nunique()


def _class_groups(p: Prepared, cfg: Config, drop_single: bool = True) -> dict[str, list[int]]:
    """class -> treatment indices (in p.treatments) with an n x treatment term. Treatments present in fewer than 2
    n's are left out (drop_single): their only cell is absorbed by their fixed effect."""
    cells = _cells_by_treatment(p)
    out: dict[str, list[int]] = {}
    for i, t in enumerate(p.treatments):
        if drop_single and cells.get(t, 0) < 2:
            continue
        out.setdefault(cfg.classes.get(t, "treatments"), []).append(i)
    return out


def _vc(p: Prepared, cfg: Config, spec: str) -> tuple[dict[str, str], dict[str, str], list[str]]:
    vc, labels, notes = {}, {}, []
    groups = _class_groups(p, cfg)
    cells = _cells_by_treatment(p)

    def enough(ix):  # a multi-treatment class needs >= 3 cells, a one-treatment class >= 2 (i.e. 2 n's)
        return sum(int(cells.get(p.treatments[i], 0)) for i in ix) >= (3 if len(ix) >= 2 else 2)

    split = spec == "full" and len(groups) >= 2 and all(enough(ix) for ix in groups.values())
    if spec == "full" and len(groups) >= 2 and not split:
        notes.append("a treatment class has too few n×treatment cells: one common n×treatment variance used")
    if split:
        for cl, ix in groups.items():
            key = f"nt_{_safe(cl)}"
            vc[key] = "0 + " + " + ".join(f"I{i}" for i in ix)
            labels[key] = f"n × treatment [{cl}]"
    elif groups:
        ix = [i for v in groups.values() for i in v]
        vc["nt_all"] = "0 + " + " + ".join(f"I{i}" for i in sorted(ix))
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


def _fixed_names(p: Prepared) -> list[str]:
    return ([f"T{i}" for i in range(1, len(p.treatments))] + [f"P{j}" for j in range(1, len(p.timepoint_levels))]
            + (["Z"] if p.covariate else []))


def _fit_spec(p: Prepared, cfg: Config, spec: str, starts: list | None = None,
              quick: bool = False) -> tuple[Fit | None, str, list[str]]:
    """Fit one model of the ladder. statsmodels MixedLM builds the design (vc_formula) and gives one starting
    point; reml.REML maximises the REML likelihood over theta >= 0 from several starts (best likelihood wins)."""
    from statsmodels.regression.mixed_linear_model import MixedLM

    from . import reml as reml_mod

    vc, labels, notes = _vc(p, cfg, spec)
    fixed = _fixed_names(p)
    formula = "y ~ " + (" + ".join(fixed) if fixed else "1")
    data = p.d.reset_index(drop=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = MixedLM.from_formula(formula, data, groups="n", re_formula="1", vc_formula=vc or None)
    except Exception as e:  # noqa: BLE001
        return None, f"design: {type(e).__name__}: {e}", notes
    vc_names = list(model.exog_vc.names) if getattr(model, "k_vc", 0) else []
    known = data["known"].to_numpy(float) if p.known else None
    groups, rows = [], []
    for g, lab in enumerate(model.group_labels):
        M = [np.asarray(model.exog_re_li[g], float)] + [np.asarray(model.exog_vc.mats[j][g], float)
                                                         for j in range(len(vc_names))]
        idx = np.asarray(model.row_indices[lab])
        groups.append(reml_mod.Group(np.asarray(model.endog_li[g], float), np.asarray(model.exog_li[g], float), M,
                                     known[idx] if known is not None else None))
        rows.append(idx)
    X = np.vstack([g.X for g in groups])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None, "fixed effects not estimable (rank-deficient design)", notes
    eng = reml_mod.REML(groups, ["n"] + vc_names)
    if p.unit == "well" and not p.replicated:
        eng.fix_residual = 1e-6 * eng.var_y
    extra = list(starts or [])
    method = "REML, bounded L-BFGS-B (multi-start)"
    if not quick:  # quick = leave-one-n-out refits: still multi-start, without the statsmodels start
        s0 = _statsmodels_start(model)
        if s0 is not None:
            extra.append(s0)
            method += " incl. statsmodels MixedLM start"
    try:
        # quick (leave-one-n-out refits): every default start plus the full fit's variances, only the (slow)
        # statsmodels start is skipped, so a refit is a multi-start fit like the main one
        rf = eng.fit(extra_starts=[s for s in extra if s is not None and len(s) == eng.K + 1])
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
    if p.unit == "well":
        desc += " (unit: well" + ("; + known binomial sampling variance" if p.known else "") + ")"
    if p.covariate:
        desc += f" + density covariate {p.covariate}"
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
    (cells - treatments - n's + 1) of that class (one-treatment class: n's - 1); with one common variance it is the
    whole n x treatment df. Never more than the RM ANOVA error df."""
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


def _L(p: Prepared, f, t: str) -> np.ndarray:
    """Linear predictor of treatment t, averaged over timepoint levels, at the mean density (Z = 0)."""
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
    """Single-step multivariate-t adjusted p-values and 95% critical values (per comparison, using its own df).
    Adjusted p-values are clamped into [raw p, min(1, k * raw p)] (k = comparisons tested in the family)."""
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
        # numerical integration loses precision in the far tail: keep p_adj within its exact bounds
        # [raw p, min(1, k * raw p)] (single-step adjustment is never below the raw p nor above Bonferroni)
        praw = float(2 * stats.t.sf(abs(tvals[i]), df))
        padj[i] = float(np.clip(1 - cover(abs(tvals[i]), df), praw, min(1.0, k * praw)))
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


CONTRAST_COLUMNS = ["family", "family_kind", "family_id", "comparison", "group", "comparator", "estimand",
                    "estimate_type", "estimate", "ci_low", "ci_high", "ci_low_adj", "ci_high_adj", "pct_of_comparator",
                    "pp_difference", "pp_ci_low", "pp_ci_high", "pp_ci_low_adj", "pp_ci_high_adj", "pct_group",
                    "pct_comparator", "log_estimate", "se", "t", "df", "N", "p", "p_adj", "significance",
                    "conclusion", "test", "note", "margin_low", "margin_high"]


PP_DRAWS = 4000


def _pp_interval(mu: np.ndarray, S: np.ndarray, df: float, levels: list[float]) -> list[tuple[float, float]]:
    """Intervals for 100 * (expit(mu_a) - expit(mu_b)) from draws of (mu_a, mu_b) ~ bivariate t(df) with the model
    covariance S: respects the possible range, is asymmetric near 0 % and carries the uncertainty of both levels."""
    rng = np.random.default_rng(12345)
    try:
        z = rng.multivariate_normal(np.zeros(2), S, size=PP_DRAWS, method="eigh")
    except (np.linalg.LinAlgError, ValueError):
        return [(np.nan, np.nan) for _ in levels]
    w = np.sqrt(df / rng.chisquare(df, PP_DRAWS)) if np.isfinite(df) else np.ones(PP_DRAWS)
    m = mu[None, :] + z * w[:, None]
    d = 100 * (expit(m[:, 0]) - expit(m[:, 1]))
    out = []
    for a in levels:
        if not np.isfinite(a):
            out.append((np.nan, np.nan))
        else:
            out.append((float(np.quantile(d, a / 2)), float(np.quantile(d, 1 - a / 2))))
    return out


def _effect(p: Prepared, mu_a: float, mu_b: float, est: float, se: float, q: float, qadj: float,
            g: tuple | None = None, cov: np.ndarray | None = None) -> dict:
    """Effect columns for one comparison on the analysis scale (est = mu_a - mu_b, se) with critical values q (95%)
    and qadj (simultaneous); g = (La, Lb, df) for the percentage-point interval of fractions."""
    ok = np.isfinite(q) and np.isfinite(se)
    oka = ok and np.isfinite(qadj)
    if p.scale.kind == "difference":  # raw scale: est is already the difference in the metric's units
        return dict(log_estimate=np.nan, se=se, estimate=est,
                    ci_low=est - q * se if ok else np.nan, ci_high=est + q * se if ok else np.nan,
                    ci_low_adj=est - qadj * se if oka else np.nan, ci_high_adj=est + qadj * se if oka else np.nan,
                    pp_difference=np.nan, pp_ci_low=np.nan, pp_ci_high=np.nan, pp_ci_low_adj=np.nan,
                    pp_ci_high_adj=np.nan, pct_group=np.nan, pct_comparator=np.nan,
                    pct_of_comparator=100 * mu_a / mu_b if mu_b > 0 else np.nan)
    r = dict(log_estimate=est, se=se, estimate=float(np.exp(est)),
             ci_low=float(np.exp(est - q * se)) if ok else np.nan, ci_high=float(np.exp(est + q * se)) if ok else np.nan,
             ci_low_adj=float(np.exp(est - qadj * se)) if oka else np.nan,
             ci_high_adj=float(np.exp(est + qadj * se)) if oka else np.nan,
             pp_difference=np.nan, pp_ci_low=np.nan, pp_ci_high=np.nan, pp_ci_low_adj=np.nan, pp_ci_high_adj=np.nan,
             pct_group=np.nan, pct_comparator=np.nan)
    if p.scale.kind == "fraction":
        pa, pb = float(expit(mu_a)), float(expit(mu_b))
        (l1, h1), (l2, h2) = (np.nan, np.nan), (np.nan, np.nan)
        if ok and g is not None and cov is not None:
            La, Lb, df = g
            Lm = np.vstack([La, Lb])
            # coverage levels matching the t critical values of the logit-scale CI (95%) and the simultaneous CI
            lev = [2 * stats.t.sf(c, df) if np.isfinite(c) else np.nan for c in (q, qadj)]
            (l1, h1), (l2, h2) = _pp_interval(np.array([mu_a, mu_b]), Lm @ cov @ Lm.T, df, lev)
        r.update(pp_difference=100 * (pa - pb), pp_ci_low=l1, pp_ci_high=h1, pp_ci_low_adj=l2, pp_ci_high_adj=h2,
                 pct_group=100 * pa, pct_comparator=100 * pb, pct_of_comparator=100 * pa / pb if pb > 0 else np.nan)
    else:
        r["pct_of_comparator"] = 100 * float(np.exp(est))
    return r


def _pp_draws(mu: np.ndarray, S: np.ndarray, df: float) -> np.ndarray | None:
    """Draws of 100 * (expit(mu_a) - expit(mu_b)) (same generator as _pp_interval)."""
    rng = np.random.default_rng(12345)
    try:
        z = rng.multivariate_normal(np.zeros(2), S, size=PP_DRAWS, method="eigh")
    except (np.linalg.LinAlgError, ValueError):
        return None
    w = np.sqrt(df / rng.chisquare(df, PP_DRAWS)) if np.isfinite(df) else np.ones(PP_DRAWS)
    m = mu[None, :] + z * w[:, None]
    return 100 * (expit(m[:, 0]) - expit(m[:, 1]))


def _equivalence_rows(p: Prepared, f, cfg: Config, info: list, base: dict) -> list[dict]:
    """TOST rows of the vehicle-check family. Holm across the pairs; the verdict uses the same adjusted level:
    equivalent = Holm-adjusted TOST p < 0.05; not equivalent = the (1 - 2 * 0.05 / m) CI (m pairs; 90% for one pair)
    lies wholly outside the margin; else inconclusive. Margins: ratio within [1/(1+m), 1+m] (count / positive;
    odds ratio is not used for fractions), fractions +- equivalence_margin_pp percentage points, raw-scale
    differences +- margin x the comparator mean."""
    m = max(1, sum(1 for x in info if x[5] >= 1))
    a2 = 2 * ALPHA / m  # two-sided level of the reported CI (0.10 -> 90% CI for one pair)
    cl = f"{100 * (1 - a2):.0f}% CI" if m == 1 else f"{100 * (1 - a2):.1f}% CI (Bonferroni-consistent with Holm, {m} pairs)"
    kind = p.scale.kind
    frows = []
    for (a, b, est, se, N, df, mu_a, mu_b, g) in info:
        if kind == "fraction":
            margin_txt = f"±{cfg.margin_pp:g} percentage points"
        elif kind == "difference":
            margin_txt = f"±{cfg.margin:.0%} of the comparator mean (raw-scale difference)"
        else:
            margin_txt = f"ratio within [{1 / (1 + cfg.margin):.3g}, {1 + cfg.margin:.3g}]"
        r = dict(base, comparison=f"{a} vs {b}", group=a, comparator=b, df=df, N=N,
                 test=f"TOST equivalence (mixed model; margin {margin_txt}; {cl})")
        if df < 1:
            r.update(_effect(p, mu_a, mu_b, est, np.nan, np.nan, np.nan), t=np.nan, p=np.nan,
                     conclusion=NEEDS_2N, note=NEEDS_2N_PAIR, margin_low=np.nan, margin_high=np.nan)
            frows.append(r)
            continue
        q = stats.t.ppf(1 - a2 / 2, df)
        r.update(_effect(p, mu_a, mu_b, est, se, q, np.nan, g, f.cov), t=np.nan)
        if kind == "fraction":
            La, Lb, _ = g
            Lm = np.vstack([La, Lb])
            d = _pp_draws(np.array([mu_a, mu_b]), Lm @ f.cov @ Lm.T, df)
            M = cfg.margin_pp
            if d is None:
                pv, lo, hi = np.nan, np.nan, np.nan
            else:
                pv = float(max(np.mean(d <= -M), np.mean(d >= M), 1 / PP_DRAWS))
                lo, hi = float(np.quantile(d, a2 / 2)), float(np.quantile(d, 1 - a2 / 2))
            r.update(pp_ci_low=lo, pp_ci_high=hi, margin_low=-M, margin_high=M)
            fmt = lambda v: f"{v:+.2f} pp"  # noqa: E731
        else:
            if kind == "difference":
                M = cfg.margin * mu_b if mu_b > 0 else np.nan
                lo_m, hi_m = -M, M
                fmt = lambda v: f"{v:.4g}"  # noqa: E731
            else:
                M = np.log(1 + cfg.margin)
                lo_m, hi_m = np.exp(-M), np.exp(M)
                fmt = lambda v: f"{np.exp(v):.3g}"  # noqa: E731
            lo, hi = est - q * se, est + q * se
            pv = (float(max(stats.t.sf((est + M) / se, df), stats.t.cdf((est - M) / se, df)))
                  if np.isfinite(M) else np.nan)
            r.update(margin_low=lo_m, margin_high=hi_m)
        r.update(p=pv, _lo=lo, _hi=hi, _M=M, _fmt=fmt)
        frows.append(r)
    for r, pa_ in zip(frows, holm([r["p"] for r in frows])):
        r["p_adj"] = pa_
        if "_M" not in r:
            continue
        lo, hi, M, fmt = r.pop("_lo"), r.pop("_hi"), r.pop("_M"), r.pop("_fmt")
        if not np.isfinite(M) or not np.isfinite(r["p"]):
            verdict = "not testable (comparator mean <= 0)" if kind == "difference" else "not testable"
        elif pa_ < ALPHA:
            verdict = "equivalent"
        elif lo >= M or hi <= -M:
            verdict = "not equivalent"
        else:
            verdict = "inconclusive"
        if kind == "fraction":
            note = f"{cl} {fmt(lo)} to {fmt(hi)} vs margin ±{M:g} pp: {verdict}"
        else:
            note = f"{cl} {fmt(lo)}-{fmt(hi)} vs margin {fmt(-M)}-{fmt(M)}: {verdict}"
        r.update(conclusion=verdict, note=note + (f" (Holm-adjusted TOST p {pa_:.3g})" if np.isfinite(pa_) else ""))
    return frows


def contrasts(p: Prepared, f: Fit, cfg: Config, simultaneous: bool = True, kinds: tuple | None = None) -> pd.DataFrame:
    """One row per comparison in every family (model-based). Fractions also get the difference in percentage points
    at the model means, expit(mu_a) - expit(mu_b), with intervals from draws of (mu_a, mu_b) from the model's
    bivariate t (same df and coverage as the logit CI; see _pp_interval); p-value = the logit-scale test, which tests
    the same null hypothesis pa = pb."""
    rows = []
    et = p.scale.estimate_type
    est_name = estimand(p.metric, p.unit) if p.metric is not None else ""
    if p.scale.kind == "fraction":
        est_name = estimand(p.metric, p.unit) + "; odds ratio and % of comparator as secondary"
    for fam in cfg.families(p.treatments):
        if kinds and fam.kind not in kinds:
            continue
        Ls, info = [], []
        for a, b in fam.pairs:
            La, Lb = _L(p, f, a), _L(p, f, b)
            L = La - Lb
            est = float(L @ f.beta)
            se = float(np.sqrt(L @ f.cov @ L))
            N = _n_with(p.d, a, b)
            df = f.df(L, df_cap(p, f, cfg, a, b)) if N >= 2 else 0.0
            mu_a, mu_b = float(La @ f.beta), float(Lb @ f.beta)
            g = (La, Lb, df)
            Ls.append(L)
            info.append((a, b, est, se, N, df, mu_a, mu_b, g))
        if not Ls:
            continue
        C = np.array(Ls) @ f.cov @ np.array(Ls).T
        sd = np.sqrt(np.diag(C))
        R = C / np.outer(sd, sd)
        tv = np.array([e / s if s > 0 and df >= 1 else np.nan for _, _, e, s, _, df, *_ in info])
        dfs = np.array([x[5] for x in info], float)
        base = dict(family=fam.name, family_kind=fam.kind, family_id=fam.id, estimate_type=et, estimand=est_name)
        if fam.kind == "vehicle_check":
            rows += _equivalence_rows(p, f, cfg, info, base)
            continue
        padj, crit = _mvt_adjust(tv, dfs, R, simultaneous)
        for j, (a, b, est, se, N, df, mu_a, mu_b, g) in enumerate(info):
            r = dict(base, comparison=f"{a} vs {b}", group=a, comparator=b, df=df, N=N,
                     test=f"mixed model contrast (t, Satterthwaite df; {fam.adjust} within family)")
            if df >= 1:
                q = stats.t.ppf(0.975, df)
                r.update(_effect(p, mu_a, mu_b, est, se, q, crit[j], g, f.cov), t=tv[j],
                         p=float(2 * stats.t.sf(abs(tv[j]), df)), p_adj=padj[j], note="")
            else:
                r.update(_effect(p, mu_a, mu_b, est, np.nan, np.nan, np.nan), t=np.nan, p=np.nan, p_adj=np.nan,
                         note=NEEDS_2N_PAIR)
            rows.append(r)
    out = pd.DataFrame(rows).reindex(columns=CONTRAST_COLUMNS)
    if len(out):
        out["significance"] = [(r.conclusion if r.conclusion != NEEDS_2N else "") if r.family_kind == "vehicle_check"
                               else stars(r.p_adj) for r in out.itertuples()]
    return out


def _observed(p: Prepared) -> pd.Series:
    """Per (n, treatment): the observed value on the natural scale - mean per field (count / positive) or pooled
    fraction sum(K) / sum(N) over the fields (fraction)."""
    f = p.fields if p.fields is not None else p.d
    f = f[f["value"].notna()]
    if p.scale.kind == "fraction" and p.metric is not None:
        N = nuclei(f, p.metric)
        if N is not None and (N > 0).any():
            g = f.assign(_N=N.where(N > 0), _K=f["value"] * N.where(N > 0)).groupby(["n", "treatment"])
            return (g["_K"].sum() / g["_N"].sum()).astype(float)
    return f.groupby(["n", "treatment"])["value"].mean()


def emm(p: Prepared, f: Fit) -> pd.DataFrame:
    """Estimated marginal means per treatment, back-transformed to absolute units, with 95% CI."""
    rows = []
    obs = _observed(p)
    for t in p.treatments:
        L = _L(p, f, t)
        mu = float(L @ f.beta)
        se = float(np.sqrt(L @ f.cov @ L))
        N = _n_with(p.d, t)
        df = f.df(L, N - 1) if N >= 2 else 0.0
        q = stats.t.ppf(0.975, df) if df >= 1 else np.nan
        rows.append(_emm_row(p, t, mu, se, q, df, N, obs, "" if N >= 2 else NEEDS_2N))
    return pd.DataFrame(rows)


EMM_COLUMNS = ["treatment", "N", "wells", "fields", "estimand", "estimate", "ci_low", "ci_high", "percent",
               "percent_ci_low", "percent_ci_high", "observed_mean", "log_mean", "se", "df", "estimate_type",
               "mean_type", "note"]


def _emm_row(p, t, mu, se, q, df, N, obs, note) -> dict:
    ok = np.isfinite(q) and np.isfinite(se)
    est = float(p.scale.back(mu))
    lo = float(p.scale.back(mu - q * se)) if ok else np.nan
    hi = float(p.scale.back(mu + q * se)) if ok else np.nan
    frac = p.scale.kind == "fraction"
    fl = p.fields if p.fields is not None else p.d
    o = obs[obs.index.get_level_values("treatment") == t] if len(obs) else obs
    mean_type = ("proportion of nuclei (typical well; back-transformed logit)" if frac else
                 ("mean per field" if p.unit == "well" else "geometric mean of fields")
                 + (f" (minus offset {p.scale.shift:g})" if p.scale.shift else ""))
    return dict(treatment=t, N=N, wells=int((p.d["treatment"] == t).sum()) if p.unit == "well" else np.nan,
                fields=int(((fl["treatment"] == t) & fl["value"].notna()).sum()),
                estimand=estimand(p.metric, p.unit, "mean") if p.metric is not None else "",
                estimate=est, ci_low=lo, ci_high=hi,
                percent=100 * est if frac else np.nan, percent_ci_low=100 * lo if frac else np.nan,
                percent_ci_high=100 * hi if frac else np.nan,
                observed_mean=float(o.mean()) if len(o) else np.nan,
                log_mean=mu, se=se, df=df, estimate_type="mean", mean_type=mean_type, note=note)


# ---------------------------------------------------------------- descriptive estimates (one n)

class _OLS:
    """Fixed-effect least-squares fit used only for point estimates when there is one n (no CIs)."""

    def __init__(self, p: Prepared):
        fixed = _fixed_names(p)
        self.names = ["Intercept"] + fixed
        X = np.column_stack([np.ones(len(p.d))] + [p.d[c].to_numpy(float) for c in fixed])
        w = 1 / (p.d["known"].to_numpy(float) + p.d["y"].var()) if p.known and len(p.d) > 1 else np.ones(len(p.d))
        sw = np.sqrt(w)
        self.beta, *_ = np.linalg.lstsq(X * sw[:, None], p.d["y"].to_numpy(float) * sw, rcond=None)
        self.cov = np.full((len(self.names), len(self.names)), np.nan)


def descriptive(p: Prepared, cfg: Config, kinds: tuple | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Estimates without any across-n inference: means and comparisons filled in from the data, every CI, p and
    star empty with the note 'needs ≥ 2 n'."""
    f = _OLS(p)
    obs = _observed(p)
    em = pd.DataFrame([_emm_row(p, t, float(_L(p, f, t) @ f.beta), np.nan, np.nan, 0.0, _n_with(p.d, t), obs, NEEDS_2N)
                       for t in p.treatments], columns=EMM_COLUMNS)
    rows = []
    et = p.scale.estimate_type
    est_name = estimand(p.metric, p.unit) if p.metric is not None else ""
    for fam in cfg.families(p.treatments):
        if kinds and fam.kind not in kinds:
            continue
        for a, b in fam.pairs:
            mu_a, mu_b = float(_L(p, f, a) @ f.beta), float(_L(p, f, b) @ f.beta)
            r = dict(family=fam.name, family_kind=fam.kind, family_id=fam.id, comparison=f"{a} vs {b}", group=a,
                     comparator=b, estimate_type=et, estimand=est_name, df=0.0, N=_n_with(p.d, a, b), t=np.nan,
                     p=np.nan, p_adj=np.nan, test="descriptive (one n: no across-n test)",
                     conclusion=NEEDS_2N if fam.kind == "vehicle_check" else "", note=NEEDS_2N)
            r.update(_effect(p, mu_a, mu_b, mu_a - mu_b, np.nan, np.nan, np.nan))
            rows.append(r)
    out = pd.DataFrame(rows).reindex(columns=CONTRAST_COLUMNS)
    if len(out):
        out["significance"] = ""
    return em, out


def variance_components(p: Prepared, f: Fit | None, cfg: Config) -> pd.DataFrame:
    """Variance components on the analysis scale, with approximate Wald CIs on the log-variance scale."""
    allc = [("n", "n (experiment: day, plate, passage)")]
    groups = _class_groups(p, cfg)
    allc += [(f"nt_{_safe(c)}", f"n × treatment [{c}]") for c in groups]
    allc += [("nt_all", "n × treatment [all treatments]")]
    if p.unit == "well":
        res_lab = "well (residual: well-to-well within n × treatment)"
        allc += [("residual", res_lab)]
    else:
        allc += [("well", "well (within n × treatment)"), ("position", "field position (within n)"),
                 ("residual", "field (residual)")]
    rows = []
    cols = ["component", "variance", "sd", "cv_pct", "pct_of_total", "ci_low", "ci_high", "note"]
    if f is None:
        for key, lab in allc:
            if key == "nt_all" and len(groups) > 1:
                continue
            v, note = np.nan, ("not estimable (no model fitted)" if len(p.ns) >= 2 else f"{NEEDS_2N}: not estimable")
            if key == "residual" and len(p.ns) < 2:
                cell_var = p.d.groupby(["n", "treatment"])["y"].var(ddof=1).dropna()
                if len(cell_var):
                    v = float(cell_var.mean())
                    note = "descriptive: technical replicates (pooled within-condition variance of this n)"
            rows.append(dict(component=lab, variance=v, sd=np.sqrt(v) if np.isfinite(v) else np.nan,
                             cv_pct=100 * np.sqrt(np.expm1(v)) if np.isfinite(v) and p.scale.kind not in ("fraction", "difference") else np.nan,
                             pct_of_total=np.nan, ci_low=np.nan, ci_high=np.nan, note=note))
        return pd.DataFrame(rows, columns=cols)
    vals = f.theta
    ses = {}
    try:  # SEs from the inverse expected REML information (components at zero excluded)
        keys = f.vc_names + ["residual"]
        keep, cov = f.reml.theta_cov()
        for k, sd_ in zip([k for k, kp in zip(keys, keep) if kp], np.sqrt(np.clip(np.diag(cov), 0, None))):
            ses[k] = float(sd_)
    except Exception:  # noqa: BLE001
        pass
    fixed_res = p.unit == "well" and not p.replicated
    total = sum(v for k, v in vals.items() if np.isfinite(v) and not (k == "residual" and fixed_res))
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
        if key == "residual" and fixed_res:
            rows.append(dict(component=lab, variance=np.nan, sd=np.nan, cv_pct=np.nan, pct_of_total=np.nan,
                             ci_low=np.nan, ci_high=np.nan,
                             note="not identifiable (one well per condition): absorbed by n × treatment"))
            continue
        if key == "residual":
            df_res = max(len(p.d) - len(f.beta), 1)
            lo, hi = v * df_res / stats.chi2.ppf(0.975, df_res), v * df_res / stats.chi2.ppf(0.025, df_res)
            note = "chi-square CI (residual df)" + ("; each well also has its known binomial sampling variance"
                                                    if p.known else "")
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
                         cv_pct=100 * np.sqrt(np.expm1(max(v, 0))) if p.scale.kind not in ("fraction", "difference") else np.nan,
                         pct_of_total=100 * v / total if total > 0 else np.nan, ci_low=lo, ci_high=hi, note=note))
    return pd.DataFrame(rows, columns=cols)


# ---------------------------------------------------------------- per-n effects and RM ANOVA paths

TECH = "descriptive: technical replicates (fields)"
PER_N_COLUMNS = ["n", "treatment", "control", "ratio", "ratio_of_means", "pp_difference", "difference",
                 "tech_ci_low", "tech_ci_high", "interval", "fields", "estimate_type"]


def _field_summary(p: Prepared) -> pd.DataFrame:
    """Per (n, treatment) from the fields: natural-scale value, its technical SE (log scale for count/positive,
    proportion scale for fractions) and the number of fields."""
    f = p.fields if p.fields is not None else p.d
    f = f[f["value"].notna()]
    rows = []
    N = nuclei(f, p.metric) if (p.scale.kind == "fraction" and p.metric is not None) else None
    for (n, t), g in f.groupby(["n", "treatment"]):
        k = len(g)
        if p.scale.kind == "fraction":
            Ng = N.loc[g.index].where(N.loc[g.index] > 0) if N is not None else None
            if Ng is not None and Ng.notna().any():
                Kg = g["value"] * Ng
                ok = Ng.notna()
                prop = float(Kg[ok].sum() / Ng[ok].sum())
                m = int(ok.sum())
                se = (float(np.sqrt(((Kg[ok] - prop * Ng[ok]) ** 2).sum() / (m * (m - 1)))) / float(Ng[ok].mean())
                      if m > 1 else np.nan)
            else:
                prop = float(g["value"].mean())
                se = float(g["value"].std(ddof=1) / np.sqrt(k)) if k > 1 else np.nan
            rows.append(dict(n=n, treatment=t, value=prop, se=se, fields=k))
        else:
            m = float(g["value"].mean())
            if p.scale.kind == "difference":  # raw scale: absolute SE of the field mean
                se = float(g["value"].std(ddof=1) / np.sqrt(k)) if k > 1 else np.nan
            else:
                se = float(g["value"].std(ddof=1) / np.sqrt(k) / m) if (k > 1 and m > 0) else np.nan
            rows.append(dict(n=n, treatment=t, value=m, se=se, fields=k))
    return pd.DataFrame(rows, columns=["n", "treatment", "value", "se", "fields"]).set_index(["n", "treatment"])


def per_n_ratios(p: Prepared, cfg: Config) -> pd.DataFrame:
    """Effect of each treatment vs each control within each n: ratio on the model scale, ratio of the arithmetic
    means (count / positive) or difference in percentage points of the pooled fractions, with a technical 95%
    interval from the fields (descriptive: it describes the measurement of that one experiment)."""
    mu = p.d.groupby(["n", "treatment"])["y"].mean()
    fs = _field_summary(p)
    frac = p.scale.kind == "fraction"
    raw = p.scale.kind == "difference"
    rows = []
    for n in p.ns:
        for c in cfg.controls:
            if (n, c) not in mu.index or (n, c) not in fs.index:
                continue
            for t in p.treatments:
                if t == c or (n, t) not in mu.index or (n, t) not in fs.index:
                    continue
                vt, vc = fs.loc[(n, t)], fs.loc[(n, c)]
                se = float(np.sqrt(vt["se"] ** 2 + vc["se"] ** 2))
                rom = float(vt["value"] / vc["value"]) if vc["value"] > 0 else np.nan
                diff = np.nan
                if raw:  # raw-scale metric: the per-n effect is the difference of the means (intensity units)
                    pp = np.nan
                    diff = float(vt["value"] - vc["value"])
                    lo, hi = (diff - 1.96 * se, diff + 1.96 * se) if np.isfinite(se) else (np.nan, np.nan)
                elif frac:
                    pp = 100 * float(vt["value"] - vc["value"])
                    lo, hi = (pp - 196 * se, pp + 196 * se) if np.isfinite(se) else (np.nan, np.nan)
                else:
                    pp = np.nan
                    ok = np.isfinite(se) and np.isfinite(rom) and rom > 0
                    lo, hi = (rom * np.exp(-1.96 * se), rom * np.exp(1.96 * se)) if ok else (np.nan, np.nan)
                rows.append(dict(n=n, treatment=t, control=c,
                                 ratio=np.nan if raw else float(np.exp(mu[(n, t)] - mu[(n, c)])),
                                 ratio_of_means=rom, pp_difference=pp, difference=diff, tech_ci_low=lo, tech_ci_high=hi,
                                 interval=TECH, fields=int(vt["fields"]), estimate_type=p.scale.estimate_type))
    return pd.DataFrame(rows, columns=PER_N_COLUMNS)


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


# ---------------------------------------------------------------- density adjustment (fraction metrics)

@dataclass
class DensityResult:
    covariate: str
    prepared: Prepared
    fit: Fit | None
    emm: pd.DataFrame
    contrasts: pd.DataFrame
    slope: dict
    model: str
    notes: list[str]


def density_overlap(p: Prepared, a: str, b: str) -> float:
    """Overlap of the density ranges of two conditions: the smaller of (share of a's wells inside b's min-max range)
    and (share of b's wells inside a's range). 1 = fully interleaved, 0 = separate."""
    xa = p.d.loc[p.d["treatment"] == a, "zraw"].to_numpy(float)
    xb = p.d.loc[p.d["treatment"] == b, "zraw"].to_numpy(float)
    if not len(xa) or not len(xb):
        return np.nan
    ina = np.mean((xa >= xb.min()) & (xa <= xb.max()))
    inb = np.mean((xb >= xa.min()) & (xb <= xa.max()))
    return float(min(ina, inb))


def density_adjusted(df: pd.DataFrame, metric: Metric, settings: dict, positions_ok: bool,
                     covariate: str | None = None, start: str = "full") -> DensityResult | None:
    """Second model for a fraction metric with local crowding as a covariate (well mean of the covariate, centred):
    estimates at the mean density, slope, and per comparison the density overlap and its flag."""
    cfg = config(settings, df["treatment"].unique())
    cov = covariate if covariate is not None else cfg.density
    if not cov or cov not in df or metric.column == cov or pd.to_numeric(df[cov], errors="coerce").notna().sum() == 0:
        return None
    p = prepare(df, metric, cfg, positions_ok, covariate=cov)
    if len(p.treatments) < 2 or len(p.d) < 3:
        return None
    notes = list(p.notes)
    f = None
    kinds = ("vs_control", "positive_control")
    if len(p.ns) >= 2:
        f, attempts = fit(p, cfg, start=start)
        notes += attempts
    if f is not None:
        em, ct = emm(p, f), contrasts(p, f, cfg, kinds=kinds)
        j = f.names.index("Z")
        b, se = float(f.beta[j]), float(np.sqrt(f.cov[j, j]))
        L = np.zeros(len(f.names))
        L[j] = 1
        dfz = f.df(L, p.nt_df)
        q = stats.t.ppf(0.975, dfz)
        slope = dict(covariate=cov, slope=b, ci_low=b - q * se, ci_high=b + q * se, df=dfz,
                     p=float(2 * stats.t.sf(abs(b / se), dfz)) if se > 0 else np.nan,
                     scale="logit per unit of the covariate")
        model = f.description
    else:
        em, ct = descriptive(p, cfg, kinds=kinds)
        o = _OLS(p)
        slope = dict(covariate=cov, slope=float(o.beta[o.names.index("Z")]), ci_low=np.nan, ci_high=np.nan, df=0.0,
                     p=np.nan, scale="logit per unit of the covariate")
        model = "descriptive (one n; least squares)" if len(p.ns) < 2 else "none (all mixed models failed)"
    if len(ct):
        ov = [density_overlap(p, r.group, r.comparator) for r in ct.itertuples()]
        ct.insert(ct.columns.get_loc("note"), "density_overlap", ov)
        ct.insert(ct.columns.get_loc("note"), "density_flag",
                  [DENSITY_FLAG if (np.isfinite(o) and o < OVERLAP_MIN) else "" for o in ov])
        ct["family"] = ct["family"] + " (density-adjusted)"
    em.insert(1, "density_at", p.z_center)
    return DensityResult(cov, p, f, em, ct, slope, model, notes)


def density_curves(dens: DensityResult, points: int = 40) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All n: (points: one per well with n, treatment, x = density, y = %), (curves: treatment, x, y = model % at
    that density for a typical n, drawn over the treatment's own density range)."""
    p = dens.prepared
    pts = p.d[["n", "treatment", "zraw", "value"]].rename(columns={"zraw": "x"}).assign(y=lambda d: 100 * d["value"])
    f = dens.fit if dens.fit is not None else _OLS(p)
    j = f.names.index("Z")
    rows = []
    for t in p.treatments:
        xs = p.d.loc[p.d["treatment"] == t, "zraw"]
        if xs.empty:
            continue
        grid = np.linspace(xs.min(), xs.max(), points)
        mu = float(_L(p, f, t) @ f.beta)
        rows += [dict(treatment=t, x=x, y=100 * float(expit(mu + f.beta[j] * (x - p.z_center)))) for x in grid]
    return pts[["n", "treatment", "x", "y"]], pd.DataFrame(rows, columns=["treatment", "x", "y"])


def density_fields(fields: pd.DataFrame, metric: Metric, covariate: str, points: int = 40
                   ) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """One n: (points: one per field), (curves: descriptive least-squares fit of the empirical logit on treatment +
    a common density slope, per treatment over its own density range). None if there is nothing to plot."""
    d = _normalise(fields)
    if covariate not in d:
        return None
    d["x"] = pd.to_numeric(d[covariate], errors="coerce")
    N = nuclei(d, metric)
    y, _ = transform(d["value"], "fraction", N)
    d = d.assign(ly=y)
    d = d[d["x"].notna() & np.isfinite(d["ly"]) & d["value"].notna()]
    trs = list(dict.fromkeys(d["treatment"]))
    if len(d) < 3 or not trs:
        return None
    X = np.column_stack([(d["treatment"] == t).to_numpy(float) for t in trs] + [d["x"].to_numpy(float)])
    beta, *_ = np.linalg.lstsq(X, d["ly"].to_numpy(float), rcond=None)
    rows = []
    for i, t in enumerate(trs):
        xs = d.loc[d["treatment"] == t, "x"]
        for x in np.linspace(xs.min(), xs.max(), points):
            rows.append(dict(treatment=t, x=x, y=100 * float(expit(beta[i] + beta[-1] * x))))
    pts = d.assign(y=100 * d["value"])[["n", "treatment", "x", "y"]]
    return pts, pd.DataFrame(rows, columns=["treatment", "x", "y"])


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
    descriptive: bool = False  # estimates filled in without across-n inference (one n)
    unit: str = "well"
    density: DensityResult | None = None

    @property
    def ok(self) -> bool:
        return self.fit is not None

    @property
    def has_estimates(self) -> bool:
        return self.ok or (self.descriptive and len(self.emm) > 0)

    @property
    def estimand(self) -> str:
        return estimand(self.metric, self.unit) if self.metric is not None else ""

    def status_dict(self) -> dict:
        return dict(status=self.status, model=self.model, spec=self.spec, formula=self.formula,
                    converged=bool(self.converged), df_rule=self.df_rule,
                    response=self.scale.describe if self.scale else "", primary_test_used=self.primary_test,
                    unit=self.unit, estimand=self.estimand,
                    mean_estimand=estimand(self.metric, self.unit, "mean") if self.metric is not None else "",
                    n_count=len(self.prepared.ns) if self.prepared else 0,
                    across_n_inference=bool(self.ok),
                    notes=list(self.notes), attempts=list(self.attempts),
                    warnings=list(self.fit.warnings) if self.fit else [],
                    optimiser=self.fit.method if self.fit else "",
                    density_model=self.density.model if self.density else "",
                    effect_type=self.scale.estimate_type if self.scale else "",
                    raw_scale=bool(self.scale is not None and self.scale.kind == "difference"),
                    raw_scale_warning=RAW_SCALE_NOTE if (self.scale is not None and self.scale.kind == "difference")
                    else "")


DF_RULE = ("t with Satterthwaite df, capped at the df of the n×treatment variance(s) the comparison relies on "
           "(containment: cells - treatments - n's + 1 per class involved, n's - 1 for a one-treatment class); "
           "means: capped at (n's with the treatment) - 1; a comparison needs >= 2 n's containing both groups")


def analyse(df: pd.DataFrame, metric: Metric, settings: dict, positions_ok: bool | None = None,
            leave_one_out: bool = True, density: bool = True) -> ModelResult:
    """Fit the hierarchical model for one metric and derive every model-based table. Never raises."""
    try:
        res = _analyse(df, metric, settings, positions_ok, leave_one_out, density)
        res.order = config(settings, df["treatment"].unique()).order
        return res
    except Exception as e:  # noqa: BLE001 - a generation build must never crash on the statistics
        cfg = config(settings, df["treatment"].unique()) if "treatment" in df else None
        return ModelResult(metric, "fallback_rm_anova", "none (error)", "", "", False, DF_RULE,
                           [f"model analysis failed ({type(e).__name__}: {e}); the RM ANOVA drives the stars"], [],
                           primary_test="rm_anova", unit=cfg.unit if cfg else "well")


def _analyse(df, metric, settings, positions_ok, leave_one_out, density=True) -> ModelResult:
    from . import design

    if positions_ok is None:
        try:
            positions_ok = design.positions_meaningful(settings)
        except ValueError:
            positions_ok = False
    cfg = config(settings, df["treatment"].unique())
    n_bad = needs_raw_scale(df, metric, cfg.unit, positions_ok)
    if n_bad:  # log impossible for some wells: analyse the whole metric as a raw-scale difference (I3)
        metric = raw_scale_metric(metric)
    p = prepare(df, metric, cfg, positions_ok)
    notes = list(p.notes) + boundary_notes(df, metric)
    if n_bad:
        notes.insert(0, f"{RAW_SCALE_NOTE} ({n_bad} {cfg.unit} value(s) <= 0)")
    dens = None
    if density and metric.kind == "fraction" and cfg.density:
        try:
            dens = density_adjusted(df, metric, settings, positions_ok)
        except Exception as e:  # noqa: BLE001
            notes.append(f"density-adjusted model failed ({type(e).__name__}: {e})")
        if dens is None and cfg.density not in df:
            notes.append(f"density covariate {cfg.density} not available: no density-adjusted estimate")
    if len(p.ns) < 2 and len(p.treatments) >= 1:
        res = ModelResult(metric, "not_estimable", f"descriptive: one n ({NEEDS_2N} for CIs, p-values and stars)",
                          "", "", False, DF_RULE, notes + [f"one n: estimates from that n, no across-n inference "
                                                           f"({NEEDS_2N})"], [], p.scale, p,
                          primary_test="none", descriptive=True, unit=p.unit, density=dens)
        res.emm, res.contrasts = descriptive(p, cfg)
        res.per_n = per_n_ratios(p, cfg)
        res.varcomp = variance_components(p, None, cfg)
        if leave_one_out:
            res.loo, res.loo_summary = loo(df, metric, settings, positions_ok, use_model=True, full=res.contrasts)
        return res
    if len(p.ns) < 2 or len(p.treatments) < 2:
        why = "only one n" if len(p.ns) < 2 else "fewer than two treatments"
        res = ModelResult(metric, "not_estimable", f"not estimable ({why})", "", "", False, DF_RULE,
                          notes + [f"mixed model not estimable: {why}"], [], p.scale, p,
                          primary_test="none", unit=p.unit, density=dens)
        res.per_n = per_n_ratios(p, cfg)
        res.varcomp = variance_components(p, None, cfg)
        return res
    f, attempts = fit(p, cfg)
    if f is None:
        res = ModelResult(metric, "fallback_rm_anova", "none (all mixed models failed)", "", "", False, DF_RULE,
                          notes + ["every mixed model failed; the RM ANOVA (Prism standard) drives the stars"],
                          attempts, p.scale, p, primary_test="rm_anova", unit=p.unit, density=dens)
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
                      primary_test=cfg.primary_test, unit=p.unit, density=dens)
    if leave_one_out:
        res.loo, res.loo_summary = loo(df, metric, settings, positions_ok, use_model=cfg.primary_test == "mixed",
                                       start=f.spec, full=res.contrasts, theta0=f.reml.theta)
    return res


def boundary_notes(df: pd.DataFrame, metric: Metric) -> list[str]:
    """Conditions sitting at the edge of the scale in every field (docs/design.md, "Unusual designs"): no nuclei at
    all (count 0: ratios rest on the +0.5 offset), or 0% / 100% positive (the empirical logit keeps them finite, so
    the odds ratio is arbitrary; the percentage-point difference is the readout)."""
    if metric.kind not in ("count", "fraction") or "value" not in df:
        return []
    v = pd.to_numeric(df["value"], errors="coerce")
    out = []
    for t, g in v.groupby(df["treatment"], sort=True):
        g = g.dropna()
        if g.empty:
            continue
        if metric.kind == "count" and (g == 0).all():
            out.append(f"{t}: 0 nuclei in every field; its ratios rest on the +{COUNT_OFFSET:g} offset (read them as "
                       "'close to 0', an upper bound)")
        elif metric.kind == "fraction" and ((g == 0).all() or (g == 1).all()):
            pct = "0%" if (g == 0).all() else "100%"
            out.append(f"{t}: {pct} in every field; the empirical logit keeps it finite, so its odds ratios are "
                       "arbitrary: read the percentage-point difference")
    return out


LOO_COLUMNS = ["dropped_n", "family", "comparison", "estimate_type", "estimate", "ci_low", "ci_high",
               "pp_difference", "p", "p_adj", "N", "note"]


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
            note = f"fewer than 2 n's left ({NEEDS_2N})"
        if tab.empty:
            rows.append(dict(dropped_n=n, family="", comparison="", note=note or "not estimable"))
            continue
        for r in tab.itertuples():
            rows.append(dict(dropped_n=n, family=r.family, comparison=r.comparison, estimate_type=r.estimate_type,
                             estimate=r.estimate, ci_low=r.ci_low, ci_high=r.ci_high,
                             pp_difference=getattr(r, "pp_difference", np.nan), p=r.p, p_adj=r.p_adj,
                             N=r.N, note=getattr(r, "note", "") or ""))
    out = pd.DataFrame(rows, columns=LOO_COLUMNS)
    summ = []
    if full is not None and len(full):
        fl = full[full["family_kind"] != "vehicle_check"] if "family_kind" in full else full
        for r in fl.itertuples():
            # null value of the estimate: 0 for differences (RM ANOVA, raw-scale difference model), 1 for ratios
            null = 0.0 if (not use_model or getattr(r, "estimate_type", "") == "difference") else 1.0
            drops = out[(out["family"] == r.family) & (out["comparison"] == r.comparison)]
            if len(ns) < 3 and np.isfinite(r.p_adj):
                # dropping one of 2 n's leaves one n: no test is possible, so stability cannot be judged
                summ.append(dict(family=r.family, comparison=r.comparison, estimate=r.estimate, p_adj=r.p_adj,
                                 significant=bool(r.p_adj < ALPHA), loo_min=np.nan, loo_max=np.nan,
                                 loo_p_adj_max=np.nan, drops=len(drops), same_direction=False,
                                 same_significance=False, stable="not testable",
                                 note=f"leave-one-n-out needs ≥ 3 n (dropping one of {len(ns)} leaves too few)"))
                continue
            if not np.isfinite(r.p_adj):
                summ.append(dict(family=r.family, comparison=r.comparison, estimate=r.estimate, p_adj=r.p_adj,
                                 significant=False, loo_min=np.nan, loo_max=np.nan, loo_p_adj_max=np.nan,
                                 drops=len(drops), same_direction=False, same_significance=False, stable="not tested",
                                 note=NEEDS_2N))
                continue
            sig_full = bool(r.p_adj < ALPHA)
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
    return out, pd.DataFrame(summ, columns=["family", "comparison", "estimate", "p_adj", "significant", "loo_min",
                                            "loo_max", "loo_p_adj_max", "drops", "same_direction",
                                            "same_significance", "stable", "note"])
