"""Forward-looking experiment planner: how many n's, wells and fields does the next experiment need?

The sampling analysis (sampling.py) looks backwards ("could we have used fewer images?"). This module fits a
variance model for planning, with honest uncertainty, and predicts power, confidence-interval width, minimum
detectable effect and cost for any design. See docs/planner.md for the method in plain language.

Variance model (relative / log scale: a CV of 30% ~ a log-scale SD of 0.3, delta method)
------------------------------------------------------------------------------------------
One image (field) value y in condition (n, t), well j, field position p:

    y / condition mean = 1 + well_j + pos_{n,p} + e

    field (residual) e   : sigma2_field   field-to-field noise within a well, after position
    position pos_{n,p}   : sigma2_pos     a field-position pattern shared by every well of one n (only when
                                          design.positions_meaningful); it cancels in a treatment/reference ratio
                                          when every well is imaged at the same positions ("standardised")
    well                 : sigma2_well    well-to-well within a condition
    biological           : sigma2_bio     n x treatment variation of the log ratio treatment / reference (biology +
                                          day/plate effects on the effect size), pooled within a treatment class

The main analysis is a paired comparison across n of d_n = log(mean_t / mean_ref) (ratio paired t-test,
df = N - 1). Its per-n variance for a design with w wells x k fields (reference: w_r wells) is

    V = sigma2_bio + sigma2_well (1/w + 1/w_r) + sigma2_f (1/(w k) + 1/(w_r k))
    sigma2_f = sigma2_field (standardised positions) or sigma2_field + sigma2_pos (random positions)

Power = noncentral t (df = N - 1, noncentrality |log effect| / sqrt(V / N)); expected 95% CI half-width =
t_{0.975, N-1} sqrt(V / N). Every component has a generalized-pivotal-quantity (GPQ) confidence interval;
joint GPQ draws give the optimistic / expected / pessimistic band of every prediction.

Public API (used by the browser planner; pure functions, no file writes):
    load(project_dir, generation=None) -> plan
    evaluate(plan, n, wells, fields, effect, reference_wells=None, alpha=0.05, target_power=0.8, ...) -> dict
    curve(plan, vary, fixed, effect, ...) -> list[dict]
run(df, settings, out, ...) builds the plan and writes every planner output into sampling/gN/.
"""
from __future__ import annotations

import fnmatch
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from . import analysis
from . import design as design_mod
from .settings import DEFAULTS as _SETTINGS_DEFAULTS

PLAN_VERSION = 1
BAND = (0.10, 0.90)  # prediction band: 80% of the component uncertainty (optimistic .. pessimistic)
CI = (0.025, 0.975)  # component confidence intervals: 95%
N_DRAWS = 1000
PLATE_FORMATS = (6, 12, 24, 48, 96)
CONTROLS, TREATMENTS, OTHER = "controls", "treatments", "other"
PLANNER_DEFAULTS = _SETTINGS_DEFAULTS["planner"]


# ================================================================ settings

def options(settings: dict) -> dict:
    """[planner] settings, typed and validated (strings from the settings form are accepted)."""
    p = {**PLANNER_DEFAULTS, **dict(settings.get("planner", {}) or {})}

    def num(key, lo=None, hi=None, cast=float):
        try:
            v = cast(float(str(p[key]).replace(",", ".")))
        except (TypeError, ValueError):
            raise ValueError(f"settings.toml [planner] {key} = {p[key]!r}: expected a number") from None
        if (lo is not None and v < lo) or (hi is not None and v > hi):
            raise ValueError(f"settings.toml [planner] {key} = {p[key]!r}: must be between {lo} and {hi}")
        return v

    eff = p["effects"] if isinstance(p["effects"], (list, tuple)) else str(p["effects"]).split(",")
    effects = []
    for e in eff:
        try:
            v = float(str(e).strip().replace(",", "."))
        except ValueError:
            raise ValueError(f"settings.toml [planner] effects: {e!r} is not a ratio (e.g. 0.7)") from None
        if v <= 0 or v == 1:
            raise ValueError(f"settings.toml [planner] effects: {e!r} must be a positive ratio other than 1")
        effects.append(v)
    rw = p["reference_wells"]
    if not (isinstance(rw, str) and rw.strip().lower() == "auto"):
        rw = int(num("reference_wells", 1, 96))
    else:
        rw = "auto"
    fmt = int(num("plate_format", cast=int))
    if fmt not in PLATE_FORMATS:
        raise ValueError(f"settings.toml [planner] plate_format = {fmt}: expected one of {PLATE_FORMATS}")
    return dict(effects=effects or list(PLANNER_DEFAULTS["effects"]), target_power=num("target_power", 0.05, 0.999),
                alpha=num("alpha", 1e-4, 0.5), cost_per_n=num("cost_per_n", 0), cost_per_well=num("cost_per_well", 0),
                cost_per_field=num("cost_per_field", 0), max_n=int(num("max_n", 3, 100)), plate_format=fmt,
                avoid_edge=bool(p["avoid_edge"]), reference_wells=rw)


def references(settings: dict, present) -> list[str]:
    """[reference] + normalise_to, de-duplicated, only those present in the data (same order as sampling.run)."""
    a = settings["analysis"]
    return [r for r in dict.fromkeys([a["reference"]] + list(a["normalise_to"])) if r in set(present)]


def treatment_classes(settings: dict, treatments, refs) -> tuple[dict[str, str], str, list[str]]:
    """{treatment: class}, source, notes. analysis.classes (class -> names or glob patterns such as "miR-*") wins;
    otherwise references (reference + normalise_to) are "controls" and everything else "treatments"."""
    spec = settings.get("analysis", {}).get("classes")
    notes = []
    if isinstance(spec, dict) and spec:
        out = {}
        for t in treatments:
            for cls, pats in spec.items():
                pats = [pats] if isinstance(pats, str) else list(pats)
                if any(t == p or fnmatch.fnmatchcase(t, p) for p in pats):
                    out[t] = str(cls)
                    break
            else:
                out[t] = CONTROLS if t in refs else OTHER
        rest = [t for t in treatments if out[t] == OTHER]
        if rest:
            notes.append(f"not in any analysis.classes entry, grouped as '{OTHER}': {', '.join(rest)}")
        return out, "analysis.classes", notes
    return {t: CONTROLS if t in refs else TREATMENTS for t in treatments}, "default (references vs others)", notes


# ================================================================ variance components

def _assigned(df: pd.DataFrame, settings: dict, project_dir) -> pd.DataFrame:
    if "well" in df and "position" in df:
        return df
    a = design_mod.assign(df, settings, project_dir)
    return df.assign(well=a["well"], position=a["position"])


def within_components(df: pd.DataFrame, positions: bool) -> dict:
    """Nested ANOVA on values scaled by their condition (n x treatment) mean, so components are CV^2 (~ log-scale
    variances). Returns mean squares with df (for the pivots) and the point estimates. Unbalanced data use the
    standard k0 = (N - sum n_i^2 / N) / (a - 1) coefficients; single wells / single fields leave a component
    'not estimable' (None) with a note."""
    d = df[np.isfinite(df["value"])].copy()
    cm = d.groupby(["n", "treatment"])["value"].transform("mean")
    d = d[cm > 0]
    d["r"] = d["value"] / cm[cm > 0]
    keys = ["n", "treatment", "well"]
    d["wm"] = d.groupby(keys)["r"].transform("mean")
    d["dev"] = d["r"] - d["wm"]
    n_wells = d.groupby(keys).ngroups
    ss_within, df_within = float((d["dev"] ** 2).sum()), len(d) - n_wells
    notes = []

    # position: a pattern shared by all wells of one n (balanced: two-way well x position within n)
    pos = None
    if positions and df_within > 0:
        ss_p = df_p = 0.0
        m_num = 0.0
        for _, g in d.groupby("n"):
            pm = g.groupby("position")["dev"].agg(["mean", "size"])
            P = len(pm)
            if P < 2 or g.groupby(["treatment", "well"]).ngroups < 2:
                continue
            ss_p += float((pm["size"] * pm["mean"] ** 2).sum())
            df_p += P - 1
            Nn = pm["size"].sum()
            m_num += (Nn - (pm["size"] ** 2).sum() / Nn)  # = (P - 1) * m0
        if df_p > 0 and df_within - df_p > 0:
            pos = dict(ss=ss_p, df=df_p, m0=m_num / df_p)
        else:
            notes.append("field position effect not estimable (too few positions or wells per n)")
    elif not positions:
        notes.append("field positions are random ([design] field_positions): no position component; "
                     "all field-to-field differences count as field noise")

    if df_within <= 0:
        notes.append("every well has a single field: field-to-field noise cannot be separated from well-to-well")
        field = None
        ms_e = 0.0
        df_e = 0
    else:
        ss_e = ss_within - (pos["ss"] if pos else 0.0)
        df_e = df_within - (pos["df"] if pos else 0)
        ms_e = ss_e / df_e
        field = dict(ms=ms_e, df=df_e)
    if pos:
        pos["ms"] = pos["ss"] / pos["df"]

    # well: between wells within a condition
    ss_w = df_w = kw = 0.0
    for _, g in d.groupby(["n", "treatment"]):
        wm = g.groupby("well")["r"].agg(["mean", "size"])
        a = len(wm)
        if a < 2:
            continue
        N = wm["size"].sum()
        ss_w += float((wm["size"] * (wm["mean"] - g["r"].mean()) ** 2).sum())
        df_w += a - 1
        kw += N - (wm["size"] ** 2).sum() / N  # = (a - 1) k0
    well = dict(ms=ss_w / df_w, df=df_w, k0=kw / df_w) if df_w else None
    if well is None:
        notes.append("every condition sits in a single well: well-to-well variance is not estimable; it is "
                     "contained in the biological component, so predictions for more wells are conservative")

    # overdispersion: within-well variance of the raw values / their mean (Poisson = 1)
    raw = d.groupby(keys)["value"].agg(["var", "mean", "size"])
    raw = raw[(raw["size"] > 1) & (raw["mean"] > 0)]
    od = float((raw["var"] / raw["mean"]).median()) if len(raw) else float("nan")
    return dict(field=field, position=pos, well=well, overdispersion=od, notes=notes, data=d)


def _points(w: dict) -> dict:
    """Point estimates (clipped at 0) from within_components."""
    f = w["field"]["ms"] if w["field"] else 0.0
    if w["field"] is None and w["well"]:
        # single field per well: the well mean square contains field noise too; all of it is "well"
        return dict(field=0.0, position=0.0, well=w["well"]["ms"] / max(w["well"]["k0"], 1))
    pos = max((w["position"]["ms"] - f) / w["position"]["m0"], 0.0) if w["position"] else 0.0
    well = max((w["well"]["ms"] - f) / w["well"]["k0"], 0.0) if w["well"] else 0.0
    return dict(field=f, position=pos, well=well)


def _pivot(ms: float, dfree: float, chi: np.ndarray) -> np.ndarray:
    return ms * dfree / chi


def within_draws(w: dict, rng: np.random.Generator, n: int = N_DRAWS) -> dict[str, np.ndarray]:
    """GPQ draws of field, position and well variances (independent chi-square pivots of the mean squares)."""
    z = np.zeros(n)
    if w["field"] is None:
        if w["well"]:
            sw = _pivot(w["well"]["ms"], w["well"]["df"], rng.chisquare(w["well"]["df"], n))
            return dict(field=z, position=z, well=sw / max(w["well"]["k0"], 1))
        return dict(field=z, position=z, well=z)
    se = _pivot(w["field"]["ms"], w["field"]["df"], rng.chisquare(w["field"]["df"], n))
    pos = z
    if w["position"]:
        sp = _pivot(w["position"]["ms"], w["position"]["df"], rng.chisquare(w["position"]["df"], n))
        pos = np.maximum((sp - se) / w["position"]["m0"], 0)
    well = z
    if w["well"]:
        sw = _pivot(w["well"]["ms"], w["well"]["df"], rng.chisquare(w["well"]["df"], n))
        well = np.maximum((sw - se) / w["well"]["k0"], 0)
    return dict(field=se, position=pos, well=well)


def _meas_coef(d: pd.DataFrame) -> pd.Series:
    """Per condition: (a_well, a_field) so that Var(condition mean, relative) = a_well s2_well + a_field s2_field
    for its actual (possibly unbalanced) wells: a_well = sum n_j^2 / N^2, a_field = 1 / N."""
    def f(g):
        sz = g.groupby("well").size()
        N = sz.sum()
        return pd.Series(dict(a_well=float((sz ** 2).sum() / N ** 2), a_field=1.0 / N))
    return d.groupby(["n", "treatment"])[["well"]].apply(f)


def log_ratios(d: pd.DataFrame, ref: str) -> pd.DataFrame:
    """n x treatment log(mean / reference mean), n's without the reference dropped (as sampling._log_ratios)."""
    means = d.groupby(["n", "treatment"])["value"].mean().unstack("treatment")
    if ref not in means:
        return pd.DataFrame()
    m = means[means[ref].notna() & (means[ref] > 0)]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(m.drop(columns=ref).div(m[ref], axis=0)).replace([np.inf, -np.inf], np.nan)


def bio_fit(lr: pd.DataFrame, coef: pd.DataFrame, ref: str) -> dict | None:
    """Pooled per-n variance V of the log ratio for the treatments of one class vs one reference, at the actual
    design, and the measurement part it contains.

    d_{n,t} = delta_t + u_n + eps_{n,t}: u_n is shared by every treatment of that n (the reference's own noise and
    day effects on the response), eps is treatment-specific. V = Var(u) + Var(eps) is estimated without bias by
    the pooled within-treatment variance; its uncertainty comes from the n mean square (df N - 1) and the
    interaction mean square (df (cells - N - T + 1)), combined with Satterthwaite / GPQ:
        V = MS_n / t_bar + (1 - 1 / t_bar) MS_e   (t_bar = cells per n)."""
    lr = lr.loc[:, lr.notna().sum() >= 2]
    cells = lr.stack().dropna()
    if len(cells) < 2 or lr.shape[1] == 0:
        return None
    T = lr.shape[1]
    ns = sorted(set(cells.index.get_level_values(0)))
    N = len(ns)
    dev = lr - lr.mean()
    ss_pooled = float((dev ** 2).sum().sum())
    df_pooled = int((lr.notna().sum() - 1).sum())
    V = ss_pooled / df_pooled
    # measurement part at the actual design: treatment cells + their reference cells
    a_w = a_f = 0.0
    for (n, t) in cells.index:
        a_w += coef.loc[(n, t), "a_well"] + coef.loc[(n, ref), "a_well"]
        a_f += coef.loc[(n, t), "a_field"] + coef.loc[(n, ref), "a_field"]
    a_w, a_f = a_w / len(cells), a_f / len(cells)
    out = dict(V=V, df_pooled=df_pooled, cells=len(cells), n_count=N, treatments=list(lr.columns), a_well=a_w,
               a_field=a_f, kind="single")
    df_e = len(cells) - N - T + 1
    if T >= 2 and df_e >= 1 and N >= 2:
        # additive two-way fit (handles missing cells) for the interaction mean square
        X = pd.get_dummies(pd.DataFrame({"n": cells.index.get_level_values(0),
                                         "t": cells.index.get_level_values(1)}), drop_first=True, dtype=float)
        X.insert(0, "c", 1.0)
        beta, *_ = np.linalg.lstsq(X.to_numpy(), cells.to_numpy(), rcond=None)
        ms_e = float(((cells.to_numpy() - X.to_numpy() @ beta) ** 2).sum() / df_e)
        tbar = len(cells) / N
        ms_n = max(tbar * (V - (1 - 1 / tbar) * ms_e), 1e-12)
        a, b = ms_n / tbar, (1 - 1 / tbar) * ms_e
        out.update(kind="two-way", ms_n=ms_n, df_n=N - 1, ms_e=ms_e, df_e=df_e, tbar=tbar, V=a + b,
                   df=(a + b) ** 2 / (a ** 2 / (N - 1) + b ** 2 / df_e))
    else:
        out.update(df=float(df_pooled))
    return out


def bio_draws(b: dict, wd: dict, rng: np.random.Generator, n: int = N_DRAWS) -> tuple[np.ndarray, np.ndarray]:
    """(V draws at the actual design, biological variance draws) for one bio_fit."""
    if b["kind"] == "two-way":
        mn = _pivot(b["ms_n"], b["df_n"], rng.chisquare(b["df_n"], n))
        me = _pivot(b["ms_e"], b["df_e"], rng.chisquare(b["df_e"], n))
        V = mn / b["tbar"] + (1 - 1 / b["tbar"]) * me
    else:
        V = _pivot(b["V"], b["df"], rng.chisquare(b["df"], n))
    return V, np.maximum(V - b["a_well"] * wd["well"] - b["a_field"] * wd["field"], 0)


# ================================================================ prediction core (pure, vectorised)

def power(log_effect, sd, n, alpha=0.05):
    """Two-sided paired t-test power: df = n - 1, noncentrality |log_effect| / (sd / sqrt(n)). Vectorised."""
    log_effect, sd, n = np.broadcast_arrays(np.abs(np.asarray(log_effect, float)), np.asarray(sd, float),
                                            np.asarray(n, float))
    out = np.full(log_effect.shape, np.nan)
    ok = (n >= 2) & (sd > 0) & np.isfinite(sd) & np.isfinite(log_effect)
    if ok.any():
        dfree = n[ok] - 1
        nc = log_effect[ok] / (sd[ok] / np.sqrt(n[ok]))
        out[ok] = _two_sided(stats.t.ppf(1 - alpha / 2, dfree), dfree, nc)
    zero = (n >= 2) & (sd == 0) & (log_effect > 0)
    out[zero] = 1.0
    return out if out.ndim else float(out)


def _two_sided(tc, dfree, nc):
    """P(|T| > tc) for a noncentral t (nc >= 0). scipy's nct occasionally returns NaN deep in a tail: the lower
    tail is then ~0, and the upper tail is taken as 1 - cdf."""
    up = stats.nct.sf(tc, dfree, nc)
    up = np.where(np.isfinite(up), up, 1 - stats.nct.cdf(tc, dfree, nc))
    lo = np.nan_to_num(stats.nct.cdf(-tc, dfree, nc), nan=0.0)
    return np.clip(up + lo, 0, 1)


_NC_CACHE: dict = {}


def nc_for_power(dfree: int, target: float = 0.8, alpha: float = 0.05) -> float:
    """Noncentrality giving `target` power for a two-sided t-test with dfree df (cached)."""
    key = (int(dfree), round(target, 6), round(alpha, 6))
    if key not in _NC_CACHE:
        tc = stats.t.ppf(1 - alpha / 2, dfree)
        f = lambda nc: float(_two_sided(tc, dfree, nc)) - target
        hi = 4.0
        while f(hi) < 0:
            hi *= 1.5
        _NC_CACHE[key] = float(optimize.brentq(f, 0.0, hi, xtol=1e-6))
    return _NC_CACHE[key]


def design_variance(comp: dict, bio, wells, fields, reference_wells, standardised: bool):
    """V = bio + well (1/w + 1/w_r) + f (1/(w k) + 1/(w_r k)); comp holds field / position / well (scalars or
    arrays); f = field (standardised positions) or field + position (random positions). Broadcasts."""
    f = comp["field"] + (0 if standardised else comp["position"])
    w, k, wr = (np.asarray(x, float) for x in (wells, fields, reference_wells))
    return bio + comp["well"] * (1 / w + 1 / wr) + f * (1 / (w * k) + 1 / (wr * k))


def _pct(log_hw):
    return 100 * (np.exp(log_hw) - 1)


def _cost(plan: dict, n, wells, fields, reference_wells):
    c = plan["costs"]
    nref = len(plan["design"]["references"])
    ntreat = max(plan["design"]["conditions"] - nref, 0)
    per_well = c["per_well"] + np.asarray(fields) * c["per_field"]
    return np.asarray(n) * (c["per_n"] + (ntreat * np.asarray(wells) + nref * np.asarray(reference_wells)) * per_well)


# ================================================================ plan access

def _bio_key(ref: str, cls: str) -> str:
    return f"{ref}|{cls}"


def _np_draws(plan: dict) -> dict:
    cache = plan.get("_np")  # filled by load(); plans straight from JSON convert on the fly
    if cache is None:
        dr = plan["draws"]
        cache = dict(field=np.asarray(dr["field"]), position=np.asarray(dr["position"]), well=np.asarray(dr["well"]),
                     bio={k: np.asarray(v) for k, v in dr["bio"].items()})
    return cache


def _resolve(plan: dict, effect, reference, cls):
    """effect: a ratio (float) or a treatment name (its observed ratio vs the reference). Returns
    (ratio, reference, class, label)."""
    ref = reference or plan["design"]["main_reference"]
    if ref not in plan["design"]["references"]:
        raise ValueError(f"unknown reference {ref!r}; expected one of {plan['design']['references']}")
    label = None
    if isinstance(effect, str):
        hit = [e for e in plan["effects"] if e["reference"] == ref and e["treatment"] == effect]
        if not hit:
            raise ValueError(f"no observed effect for {effect!r} vs {ref}")
        label, ratio = effect, hit[0]["ratio"]
        cls = cls or plan["classes"].get(effect)
    else:
        ratio = float(effect)
        if not ratio > 0:
            raise ValueError("effect must be a positive ratio, e.g. 0.7")
    cls = cls or plan["default_class"]
    if _bio_key(ref, cls) not in plan["bio"]:
        avail = [k.split("|", 1)[1] for k in plan["bio"] if k.startswith(ref + "|")]
        if not avail:
            raise ValueError(f"no biological variance estimate vs {ref}")
        cls = plan["default_class"] if plan["default_class"] in avail else avail[0]
    return ratio, ref, cls, label


def _auto_ref_wells(plan, wells):
    k = max(plan["design"]["conditions"] - 1, 1)
    return max(int(round(wells * math.sqrt(k))), wells)


def evaluate(plan: dict, n: int, wells: int, fields: int, effect, reference_wells=None, alpha: float = 0.05,
             target_power: float = 0.8, reference: str | None = None, cls: str | None = None,
             positions: str | None = None) -> dict:
    """Predict one design. Pure and fast (< 5 ms); safe to call on every slider move.

    plan            the dict from load() (or from run())
    n               number of independent experiments (n's, one plate each), >= 2
    wells, fields   wells per treatment condition per n, fields (images) per well
    effect          true ratio treatment / reference to plan for (e.g. 0.7 = 30% fewer nuclei), or a treatment
                    name (uses its observed ratio and its class)
    reference_wells wells per reference condition (None = same as wells; "auto" = Dunnett sqrt(k) rule)
    alpha, target_power  test level and the power that defines the minimum detectable effect (MDE)
    reference, cls  which reference (default: main) and treatment class (default: the largest non-control class)
    positions       "standardised" (same field spots in every well; default when positions are meaningful) or
                    "random"

    Returns {power, ci_halfwidth_pct, mde_ratio, sd_log_ratio}: each {low, mid, high} in numeric order
    (mid = point estimate; low/high = the 10th/90th percentile of the component uncertainty: power low and
    CI-width high are pessimistic; mde_ratio low = pessimistic for a decrease, it is the ratio below 1 that is
    detectable with target_power), plus cost (planner cost units per experiment series), images, wells_per_n,
    measurement_share (share of V that more wells/fields can reduce), effect ratio, reference, class, and
    `notes` (short plain-language strings)."""
    if n < 2:
        raise ValueError("n must be at least 2 (a paired comparison needs 2 n's)")
    if wells < 1 or fields < 1:
        raise ValueError("wells and fields must be at least 1")
    ratio, ref, cls, label = _resolve(plan, effect, reference, cls)
    if reference_wells is None:
        reference_wells = wells
    elif reference_wells == "auto":
        reference_wells = _auto_ref_wells(plan, wells)
    std = (positions or plan["design"]["positions_default"]) == "standardised"
    pt = plan["components_point"]
    bio_pt = plan["bio"][_bio_key(ref, cls)]["var"]
    dr = _np_draws(plan)
    V_mid = float(design_variance(pt, bio_pt, wells, fields, reference_wells, std))
    V_d = design_variance(dr, dr["bio"][_bio_key(ref, cls)], wells, fields, reference_wells, std)
    q = np.quantile(V_d, BAND)
    V = np.array([min(q[0], V_mid), V_mid, max(q[1], V_mid)])  # optimistic, expected, pessimistic
    sd = np.sqrt(V)
    le = abs(math.log(ratio))
    pw = power(le, sd, n, alpha)
    tc = stats.t.ppf(1 - alpha / 2, n - 1)
    hw = _pct(tc * sd / math.sqrt(n))
    nc = nc_for_power(n - 1, target_power, alpha)
    mde = np.exp(-nc * sd / math.sqrt(n))  # decrease detectable with target power
    meas = V_mid - bio_pt
    wells_per_n = (plan["design"]["conditions"] - len(plan["design"]["references"])) * wells + \
        len(plan["design"]["references"]) * reference_wells
    res = dict(n=int(n), wells=int(wells), fields=int(fields), reference_wells=int(reference_wells),
               effect_ratio=ratio, effect_label=label or f"ratio {ratio:g}", reference=ref, cls=cls,
               positions="standardised" if std else "random",
               power=dict(low=float(pw[2]), mid=float(pw[1]), high=float(pw[0])),
               ci_halfwidth_pct=dict(low=float(hw[0]), mid=float(hw[1]), high=float(hw[2])),
               mde_ratio=dict(low=float(mde[2]), mid=float(mde[1]), high=float(mde[0])),
               sd_log_ratio=dict(low=float(sd[0]), mid=float(sd[1]), high=float(sd[2])),
               measurement_share=float(meas / V_mid) if V_mid > 0 else float("nan"),
               cost=float(_cost(plan, n, wells, fields, reference_wells)), images=int(wells * fields),
               wells_per_n=int(wells_per_n), t_critical=float(tc))
    res["notes"] = _notes(plan, res, target_power)
    return res


def _notes(plan, r, target):
    notes = []
    p = r["power"]
    if p["mid"] >= target and p["low"] >= target:
        notes.append(f"Reaches {target:.0%} power even in the pessimistic case.")
    elif p["mid"] >= target:
        notes.append(f"Reaches {target:.0%} power on the expected variance, but not in the pessimistic case "
                     f"({p['low']:.0%}).")
    else:
        notes.append(f"Below {target:.0%} power ({p['mid']:.0%} expected).")
    if r["n"] <= 3:
        notes.append(f"With {r['n']} n's the t critical value is {r['t_critical']:.2f}: very few n's make any test weak.")
    ms = r["measurement_share"]
    if np.isfinite(ms):
        if ms < 0.25:
            notes.append(f"Only {ms:.0%} of the noise is measurement: more wells or fields barely help; add n's.")
        elif ms > 0.6:
            notes.append(f"{ms:.0%} of the noise is measurement: more fields or wells per n help a lot.")
    fmt = plan["options"]["plate_format"]
    if r["wells_per_n"] > fmt:
        notes.append(f"Needs {r['wells_per_n']} wells per n: more than one {fmt}-well plate.")
    if plan["components"]["well"]["var"] is None:
        notes.append("Well-to-well variance could not be estimated (single wells); more wells may help more than shown.")
    if r["effect_label"] in {e["treatment"] for e in plan["effects"]}:
        notes.append("Observed effects from few n's tend to be overestimated: plan with a smaller effect too.")
    return notes


def curve(plan: dict, vary: str = "n", fixed: dict | None = None, effect=0.7, values=None, **kw) -> list[dict]:
    """Points for a chart: evaluate() while one of n / wells / fields / reference_wells varies.

    vary    "n" | "wells" | "fields" | "reference_wells"
    fixed   the other design values, e.g. {"n": 4, "wells": 2, "fields": 6}; missing ones default to the current
            design (n defaults to the current number of n's, at least 2)
    effect  ratio or treatment name (as evaluate)
    values  the x values (default: n 2..max_n, wells 1..6, fields 1..max(3 x current, 12), reference_wells 1..12)
    kw      passed to evaluate (alpha, target_power, reference, cls, positions)
    Returns [{"x": value, **evaluate(...)}, ...]."""
    d = plan["design"]
    base = dict(n=max(d["n_count"], 2), wells=d["wells"], fields=d["fields_per_well"], reference_wells=None)
    base.update(fixed or {})
    if values is None:
        values = dict(n=range(2, plan["options"]["max_n"] + 1), wells=range(1, 7),
                      fields=range(1, max(3 * d["fields_per_well"], 12) + 1), reference_wells=range(1, 13))[vary]
    if vary not in ("n", "wells", "fields", "reference_wells"):
        raise ValueError("vary must be n, wells, fields or reference_wells")
    out = []
    for v in values:
        args = {**base, vary: int(v)}
        out.append(dict(x=int(v), **evaluate(plan, args["n"], args["wells"], args["fields"], effect,
                                             reference_wells=args["reference_wells"], **kw)))
    return out


def load(project_dir, generation: str | None = None) -> dict:
    """Read sampling/<generation>/planner.json (default: the latest generation that has one) and return the plan,
    with numpy copies of the draws cached under "_np" for fast evaluate() calls. Raises FileNotFoundError."""
    root = Path(project_dir) / "sampling"
    if generation:
        path = root / generation / "planner.json"
    else:
        cands = sorted((p for p in root.glob("g*") if p.name[1:].isdigit() and (p / "planner.json").exists()),
                       key=lambda p: int(p.name[1:]))
        if not cands:
            raise FileNotFoundError(f"no planner.json under {root} - run `kerncount sampling` first")
        path = cands[-1] / "planner.json"
    plan = json.loads(path.read_text(encoding="utf-8"))
    plan["_np"] = _np_draws(plan)
    return plan


# ================================================================ building the plan

def _ci(x: np.ndarray, lv=CI) -> list[float]:
    return [float(v) for v in np.quantile(x, lv)]


def _class_of_ref_comparisons(classes, lr_cols):
    by = {}
    for t in lr_cols:
        by.setdefault(classes.get(t, OTHER), []).append(t)
    return by


def build(df: pd.DataFrame, settings: dict, project_dir=None, seed: int = 0, n_draws: int = N_DRAWS,
          metric_label: str = "nuclei per field") -> dict:
    """Fit the planning variance model on df (n, treatment, image, value; well/position optional) and return the
    plan dict (JSON-serialisable apart from nothing: numpy is converted)."""
    rng = np.random.default_rng(seed)
    opt = options(settings)
    d0 = _assigned(df, settings, project_dir)
    positions = design_mod.positions_meaningful(settings)
    desc = design_mod.describe(design_mod.assign(df[["n", "treatment", "image"]], settings, project_dir)
                               if "image" in df else d0)
    W, K = design_mod.nominal(desc)
    wc = within_components(d0, positions)
    d = wc.pop("data")
    pt = _points(wc)
    wd = within_draws(wc, rng, n_draws)
    notes = list(wc["notes"])

    treatments = analysis.treatment_order(sorted(d["treatment"].unique()), settings["analysis"]["treatment_order"])
    refs = references(settings, treatments)
    classes, class_source, cnotes = treatment_classes(settings, treatments, refs)
    notes += cnotes
    coef = _meas_coef(d)
    ns = sorted(d["n"].unique())

    bio, bio_draws_out, effects = {}, {}, []
    for ref in refs:
        lr = log_ratios(d, ref)
        for cls, ts in _class_of_ref_comparisons(classes, list(lr.columns)).items():
            b = bio_fit(lr[ts], coef, ref)
            if b is None:
                continue
            Vd, bd = bio_draws(b, wd, rng, n_draws)
            meas_act = b["a_well"] * pt["well"] + b["a_field"] * pt["field"]
            b.update(var=max(b["V"] - meas_act, 0.0), ci95=_ci(bd), ci80=_ci(bd, BAND), V_ci95=_ci(Vd),
                     meas_actual=meas_act, reference=ref, cls=cls)
            if b["var"] == 0:
                notes.append(f"biological variance vs {ref} for class '{cls}' is estimated as 0 (upper 95% limit "
                             f"{b['ci95'][1]:.3f}); predictions use the band to stay honest")
            bio[_bio_key(ref, cls)] = b
            bio_draws_out[_bio_key(ref, cls)] = bd
        for t in lr.columns:
            x = lr[t].dropna()
            if len(x) < 1:
                continue
            e = dict(reference=ref, treatment=t, cls=classes.get(t, OTHER), n_pairs=int(len(x)),
                     log_ratio=float(x.mean()), ratio=float(np.exp(x.mean())), n_list=list(x.index),
                     values=[float(v) for v in x])
            if len(x) >= 2:
                sd = float(x.std(ddof=1))
                tc = stats.t.ppf(0.975, len(x) - 1)
                e.update(sd_own=sd, p=float(stats.ttest_1samp(x, 0).pvalue) if sd > 0 else float("nan"),
                         ci95_ratio=[float(np.exp(x.mean() - tc * sd / np.sqrt(len(x)))),
                                     float(np.exp(x.mean() + tc * sd / np.sqrt(len(x))))])
            effects.append(e)

    # default class: the non-control class with most treatments
    counts = pd.Series([c for t, c in classes.items() if t not in refs]).value_counts()
    default_class = str(counts.index[0]) if len(counts) else CONTROLS

    def comp_entry(name, var, draws, dfree, note=""):
        if var is None:
            return dict(var=None, sd=None, ci95=None, ci80=None, df=None, note=note)
        return dict(var=float(var), sd=float(np.sqrt(var)), ci95=_ci(draws), ci80=_ci(draws, BAND),
                    df=None if dfree is None else float(dfree), note=note)

    comps = dict(
        field=comp_entry("field", pt["field"] if wc["field"] else None, wd["field"],
                         wc["field"]["df"] if wc["field"] else None,
                         "within-well field-to-field noise" + (" after removing the field-position pattern"
                                                                if wc["position"] else "")),
        position=comp_entry("position", pt["position"] if wc["position"] else None, wd["position"],
                            wc["position"]["df"] if wc["position"] else None,
                            "field-position pattern shared by the wells of one n; cancels with standardised positions"),
        well=comp_entry("well", pt["well"] if wc["well"] else None, wd["well"], wc["well"]["df"] if wc["well"] else None,
                        "well-to-well within a condition"),
    )
    pos_share = (pt["position"] / (pt["position"] + pt["field"])) if (pt["position"] + pt["field"]) > 0 else 0.0
    missing = {}
    for ref in refs:
        have = set(d.loc[d["treatment"] == ref, "n"])
        missing[ref] = [n for n in ns if n not in have]
    images = d.groupby(["n", "treatment"]).size()
    design = dict(n_count=len(ns), ns=ns, conditions=len(treatments), condition_names=treatments, references=refs,
                  main_reference=refs[0] if refs else None, wells=int(W), fields_per_well=int(K),
                  images_per_condition=float(images.median()), images_min=int(images.min()),
                  images_max=int(images.max()), positions_meaningful=bool(positions),
                  positions_default="standardised" if positions else "random",
                  unbalanced=int(desc["unbalanced"].sum()) if len(desc) else 0,
                  single_well=int(desc["single_well"].sum()) if len(desc) else 0, missing_reference=missing)
    plan = dict(version=PLAN_VERSION, created=time.strftime("%Y-%m-%d %H:%M:%S"), metric=metric_label,
                scale="relative (CV^2 ~ log-scale variance)", band=list(BAND), ci=list(CI),
                options=opt, costs=dict(per_n=opt["cost_per_n"], per_well=opt["cost_per_well"],
                                        per_field=opt["cost_per_field"], unit="minutes (or any unit you choose)"),
                design=design, components=comps, components_point=pt,
                position_share_of_field_variance=float(pos_share), overdispersion=wc["overdispersion"],
                classes=classes, class_source=class_source, default_class=default_class,
                bio={k: {kk: vv for kk, vv in v.items()} for k, v in bio.items()},
                effects=effects, notes=notes,
                draws=dict(field=wd["field"], position=wd["position"], well=wd["well"], bio=bio_draws_out))
    return _jsonable(plan)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return None if not math.isfinite(v) else float(f"{v:.6g}")
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def _py(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    return v


# ================================================================ derived analyses (grids, cost, allocation, more n)

def focus_effect(plan: dict) -> float:
    """The planning effect used for single-effect summaries: the middle one by size (0.7 of [0.8, 0.7, 0.5])."""
    effs = sorted(plan["options"]["effects"], key=lambda e: abs(math.log(e)))
    return effs[len(effs) // 2]


def _grid(plan: dict):
    d = plan["design"]
    ns = np.arange(2, plan["options"]["max_n"] + 1)
    ws = np.arange(1, 7)
    ks = np.arange(1, max(3 * d["fields_per_well"], 12) + 1)
    N, Wg, Kg = (a.ravel() for a in np.meshgrid(ns, ws, ks, indexing="ij"))
    return N, Wg, Kg


def _bands(plan, key, W, K, WR, std=None):
    """V (optimistic, expected, pessimistic) arrays for designs W x K (reference wells WR)."""
    std = plan["design"]["positions_default"] == "standardised" if std is None else std
    pt, dr = plan["components_point"], _np_draws(plan)
    W, K, WR = (np.asarray(x, float) for x in (W, K, WR))
    mid = design_variance(pt, plan["bio"][key]["var"], W, K, WR, std)
    uw, inv = np.unique(np.stack([W, K, WR]), axis=1, return_inverse=True)
    inv = np.asarray(inv).ravel()
    Vd = design_variance({k: dr[k][:, None] for k in ("field", "position", "well")}, dr["bio"][key][:, None],
                         uw[0][None, :], uw[1][None, :], uw[2][None, :], std)
    q = np.quantile(Vd, BAND, axis=0)
    return np.minimum(q[0][inv], mid), mid, np.maximum(q[1][inv], mid)


def n_needed(log_eff, V, target=0.8, alpha=0.05, cap=60):
    """Smallest n (2..cap) with power >= target, per V (NaN if not reached by cap)."""
    V = np.atleast_1d(np.asarray(V, float))
    out = np.full(V.shape, np.nan)
    ns = np.arange(2, cap + 1)
    pw = power(log_eff, np.sqrt(V)[:, None], ns[None, :], alpha)
    ok = pw >= target
    hit = ok.any(axis=1)
    out[hit] = ns[ok[hit].argmax(axis=1)]
    return out


def forward_tables(plan: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """design_forward (power / CI width / cost for planning and observed effects over the design grid) and mde
    (minimum detectable ratio at target power) tables. Reference wells = treatment wells in the grid."""
    o = plan["options"]
    alpha, target = o["alpha"], o["target_power"]
    N, Wg, Kg = _grid(plan)
    cost = _cost(plan, N, Wg, Kg, Wg)
    tc = stats.t.ppf(1 - alpha / 2, N - 1)
    nc = np.array([nc_for_power(n - 1, target, alpha) for n in N])
    fwd, mde = [], []
    for key in plan["bio"]:
        ref, cls = key.split("|", 1)
        lo, mid, hi = _bands(plan, key, Wg, Kg, Wg)
        sd = np.sqrt(np.stack([lo, mid, hi]))
        hw = _pct(tc * sd / np.sqrt(N))
        md = np.exp(-nc * sd / np.sqrt(N))
        base = pd.DataFrame(dict(reference=ref, cls=cls, n=N, wells=Wg, fields_per_well=Kg, reference_wells=Wg,
                                 images_per_condition=Wg * Kg, cost=cost, sd_log_ratio_low=sd[0],
                                 sd_log_ratio_mid=sd[1], sd_log_ratio_high=sd[2],
                                 ci_halfwidth_pct_low=hw[0], ci_halfwidth_pct_mid=hw[1], ci_halfwidth_pct_high=hw[2]))
        mde.append(base.drop(columns=["cost"]).assign(mde_ratio_low=md[2], mde_ratio_mid=md[1], mde_ratio_high=md[0]))
        efx = [(f"planning {e:g}", e) for e in o["effects"]] if cls == plan["default_class"] else []
        efx += [(e["treatment"], e["ratio"]) for e in plan["effects"]
                if e["reference"] == ref and e["cls"] == cls and e["n_pairs"] >= 2]
        for label, ratio in efx:
            p = power(abs(math.log(ratio)), sd, N[None, :], alpha)
            fwd.append(base.assign(effect=label, effect_ratio=ratio, power_low=p[2], power_mid=p[1], power_high=p[0]))
    cols = ["reference", "cls", "effect", "effect_ratio", "n", "wells", "fields_per_well", "reference_wells",
            "images_per_condition", "power_low", "power_mid", "power_high", "ci_halfwidth_pct_low",
            "ci_halfwidth_pct_mid", "ci_halfwidth_pct_high", "sd_log_ratio_low", "sd_log_ratio_mid",
            "sd_log_ratio_high", "cost"]
    fwd = pd.concat(fwd, ignore_index=True)[cols] if fwd else pd.DataFrame(columns=cols)
    mde = pd.concat(mde, ignore_index=True) if mde else pd.DataFrame()
    return fwd.rename(columns={"cls": "class"}), mde.rename(columns={"cls": "class"})


def cost_analysis(plan: dict, fwd: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Pareto frontier of power (expected) vs cost per planning effect (main reference, default class), the
    cheapest design reaching target power (expected and pessimistic), and the classic optimum
    fields per well k* = sqrt((c_well / c_field) (s2_field / s2_well)) and wells per n
    w* = sqrt(B c_n / (A D)) for V(w) = A + B / w, A = bio, B = 2 (s2_well + s2_field / k*), D = conditions
    x (c_well + k* c_field)."""
    o, c, pt = plan["options"], plan["costs"], plan["components_point"]
    ref, cls = plan["design"]["main_reference"], plan["default_class"]
    target = o["target_power"]
    rows, cheapest = [], {}
    sub = fwd[(fwd["reference"] == ref) & (fwd["class"] == cls) & fwd["effect"].str.startswith("planning ")]
    for label, g in sub.groupby("effect", sort=False):
        g = g.sort_values(["cost", "power_mid"], ascending=[True, False])
        best = -1.0
        for r in g.itertuples():
            if r.power_mid > best + 1e-9:
                best = r.power_mid
                rows.append(dict(effect=label, effect_ratio=r.effect_ratio, n=r.n, wells=r.wells,
                                 fields_per_well=r.fields_per_well, reference_wells=r.reference_wells, cost=r.cost,
                                 power_low=r.power_low, power_mid=r.power_mid, power_high=r.power_high,
                                 ci_halfwidth_pct_mid=r.ci_halfwidth_pct_mid))
        ch = {}
        for band, col in (("expected", "power_mid"), ("pessimistic", "power_low")):
            ok = g[g[col] >= target]
            ch[band] = None if ok.empty else {k: _py(v) for k, v in ok.iloc[0][
                ["n", "wells", "fields_per_well", "cost", "power_low", "power_mid", "power_high"]].items()}
        cheapest[label] = ch
    front = pd.DataFrame(rows)
    if len(front):
        def is_cheapest(r):
            e = cheapest[r.effect]["expected"]
            return bool(e) and (e["n"], e["wells"], e["fields_per_well"]) == (r.n, r.wells, r.fields_per_well)
        front["cheapest_reaching_target"] = [is_cheapest(r) for r in front.itertuples()]
    f, wv = pt["field"], pt["well"]
    if plan["design"]["positions_default"] == "random":
        f += pt["position"]
    kmax = max(3 * plan["design"]["fields_per_well"], 12)
    note = ""
    if c["per_field"] <= 0:
        k_opt, note = float(kmax), "fields cost nothing: take as many as practical"
    elif wv <= 0:
        k_opt, note = float(kmax), "well-to-well variance estimated as 0: extra fields are as good as extra wells"
    else:
        k_opt = math.sqrt((c["per_well"] / c["per_field"]) * (f / wv))
    k_use = int(min(max(round(k_opt), 1), kmax))
    A = plan["bio"][_bio_key(ref, cls)]["var"] if _bio_key(ref, cls) in plan["bio"] else 0.0
    B = 2 * (wv + f / k_use)
    D = plan["design"]["conditions"] * (c["per_well"] + k_use * c["per_field"])
    w_opt = math.sqrt(B * c["per_n"] / (A * D)) if A > 0 and D > 0 else float("inf")
    optimum = dict(fields_per_well=k_opt, fields_per_well_rounded=k_use,
                   wells_per_n=w_opt if math.isfinite(w_opt) else None,
                   wells_per_n_rounded=int(min(max(round(w_opt), 1), 6)) if math.isfinite(w_opt) else 6, note=note,
                   formula_fields="sqrt((cost_per_well / cost_per_field) * (s2_field / s2_well))",
                   formula_wells="sqrt(B * cost_per_n / (A * D)); A = bio, B = 2 (s2_well + s2_field / k), "
                                 "D = conditions x (cost_per_well + k cost_per_field)")
    return front, dict(cheapest=cheapest, optimum=optimum)


def allocation(plan: dict) -> tuple[pd.DataFrame, dict]:
    """Reference wells: Dunnett's sqrt(k) rule, its honest gain from the variance components (only the reference's
    measurement part shrinks), and the best split of one plate between treatment and reference wells
    (minimum mean V over all treatment x reference comparisons, at the current fields per well)."""
    from . import platemap

    d, o = plan["design"], plan["options"]
    refs, C = d["references"], d["conditions"]
    R, T = len(refs), max(C - len(refs), 0)
    K, W = d["fields_per_well"], d["wells"]
    cls = plan["default_class"]
    eff = focus_effect(plan)
    cap_all = o["plate_format"]
    cap_inner = len(platemap.inner_wells(cap_all))
    keys = [_bio_key(r, cls) for r in refs if _bio_key(r, cls) in plan["bio"]]
    if not keys or T == 0 or R == 0:
        return pd.DataFrame(), dict(note="no treatment x reference comparisons to allocate")
    pt = plan["components_point"]
    std = d["positions_default"] == "standardised"

    def meanV(wt, wr):
        return float(np.mean([design_variance(pt, plan["bio"][k]["var"], wt, K, wr, std) for k in keys]))

    kcomp = C - 1
    dunnett = max(int(round(W * math.sqrt(kcomp))), W)
    V0 = meanV(W, W)
    rows = []
    for wt in range(1, 7):
        for wr in range(wt, 4 * wt + 1):
            used = T * wt + R * wr
            if used > cap_all and not (wt == W and wr in (W, dunnett)):
                continue
            V = meanV(wt, wr)
            rows.append(dict(treatment_wells=wt, reference_wells=wr, wells_used=used, fits_plate=used <= cap_all,
                             fits_inner_wells=used <= cap_inner, mean_V=V, sd_log_ratio=math.sqrt(V),
                             variance_vs_current=V / V0,
                             n_needed=_py(n_needed(abs(math.log(eff)), V, o["target_power"], o["alpha"])[0]),
                             rule="current" if (wt == W and wr == W) else
                             ("dunnett_sqrt_k" if (wt == W and wr == dunnett) else "")))
    tab = pd.DataFrame(rows)
    fit = tab[tab["fits_plate"]]
    if fit.empty:
        return tab, dict(note=f"not even 1 well per condition fits a {cap_all}-well plate")
    best = fit.sort_values(["mean_V", "wells_used", "reference_wells"]).iloc[0]
    tab["recommended"] = (tab["treatment_wells"] == best["treatment_wells"]) & \
                         (tab["reference_wells"] == best["reference_wells"])
    Vd = meanV(W, dunnett)
    summary = dict(comparisons_per_reference=kcomp, sqrt_k=math.sqrt(kcomp), dunnett_reference_wells=dunnett,
                   current_wells=W, gain_variance_pct=100 * (1 - Vd / V0), dunnett_wells_used=T * W + R * dunnett,
                   dunnett_fits=T * W + R * dunnett <= cap_all,
                   recommended=dict(treatment_wells=int(best["treatment_wells"]),
                                    reference_wells=int(best["reference_wells"]), wells_used=int(best["wells_used"]),
                                    gain_variance_pct=100 * (1 - best["mean_V"] / V0),
                                    n_needed=best["n_needed"], fits_inner_wells=bool(best["fits_inner_wells"])),
                   n_needed_current=_py(n_needed(abs(math.log(eff)), V0, o["target_power"], o["alpha"])[0]),
                   plate_format=cap_all, inner_wells=cap_inner, treatments=T, references=R, focus_effect=eff)
    return tab, summary


def more_n(plan: dict, k_max: int = 8, sims: int = 4000, seed: int = 0) -> tuple[pd.DataFrame, dict]:
    """Predictive probability that the final paired test at N_total = current + k reaches p < alpha, per
    treatment x reference. Approximate posterior: sigma2 ~ scaled-inverse-chi-square(df, V) with V and df the
    class-pooled per-n variance at the current design; mean ~ Normal(observed mean, sigma2 / n_now). New n's are
    simulated at the current design and added to the observed log ratios.
    Planning information only: testing after every new n and stopping when p < alpha (optional stopping) inflates
    false positives; the null simulation in the summary shows by how much."""
    rng = np.random.default_rng(seed)
    alpha = plan["options"]["alpha"]
    rows = []
    for e in plan["effects"]:
        key = _bio_key(e["reference"], e["cls"])
        if e["n_pairs"] < 2 or key not in plan["bio"]:
            continue
        b = plan["bio"][key]
        x = np.asarray(e["values"], float)
        n0 = len(x)
        s2 = b["V"] * b["df"] / rng.chisquare(b["df"], sims)
        mu = rng.normal(x.mean(), np.sqrt(s2 / n0))
        new = rng.normal(mu[:, None], np.sqrt(s2)[:, None], (sims, k_max))
        cs = x.sum() + np.cumsum(new, axis=1)
        cq = (x ** 2).sum() + np.cumsum(new ** 2, axis=1)
        Ns = n0 + np.arange(1, k_max + 1)
        m = cs / Ns
        var = np.maximum((cq - Ns * m ** 2) / (Ns - 1), 1e-300)
        p = 2 * stats.t.sf(np.abs(m) / np.sqrt(var / Ns), Ns - 1)
        prob = (p < alpha).mean(axis=0)
        for k in range(k_max):
            rows.append(dict(reference=e["reference"], treatment=e["treatment"], cls=e["cls"], n_now=n0,
                             ratio_now=e["ratio"], p_now=e.get("p"), k_added=k + 1, n_total=int(Ns[k]),
                             prob_p_below_alpha=float(prob[k])))
    tab = pd.DataFrame(rows).rename(columns={"cls": "class"})
    # optional stopping under the null: start at n0 n's, test after every new n up to n0 + k_max
    n0 = max(2, min((e["n_pairs"] for e in plan["effects"] if e["n_pairs"] >= 2), default=3))
    z = rng.normal(0, 1, (sims, n0 + k_max))
    cs, cq = np.cumsum(z, axis=1), np.cumsum(z ** 2, axis=1)
    Ns = np.arange(1, n0 + k_max + 1)
    m = cs / Ns
    with np.errstate(divide="ignore", invalid="ignore"):
        var = (cq - Ns * m ** 2) / (Ns - 1)
        p = 2 * stats.t.sf(np.abs(m) / np.sqrt(var / Ns), Ns - 1)
    seq = float((p[:, n0 - 1:] < alpha).any(axis=1).mean())
    summary = dict(start_n=int(n0), k_max=k_max, sims=sims, false_positive_fixed_n=alpha,
                   false_positive_optional_stopping=seq)
    return tab, summary


def position_effect(plan: dict) -> dict:
    """How much standardised field positions (or a model adjusting for position) reduce the noise."""
    d, pt = plan["design"], plan["components_point"]
    out = dict(meaningful=d["positions_meaningful"], estimated=plan["components"]["position"]["var"] is not None,
               share_of_field_variance=plan["position_share_of_field_variance"])
    key = _bio_key(d["main_reference"], plan["default_class"])
    if not out["estimated"] or key not in plan["bio"]:
        return out
    W, K, bio = d["wells"], d["fields_per_well"], plan["bio"][key]["var"]
    Vs = float(design_variance(pt, bio, W, K, W, True))
    Vr = float(design_variance(pt, bio, W, K, W, False))
    k_eq = next((k for k in range(1, 500) if design_variance(pt, bio, W, k, W, False) <= Vs), None)
    eff, o = focus_effect(plan), plan["options"]
    out.update(field_sd_random=math.sqrt(pt["field"] + pt["position"]), field_sd_standardised=math.sqrt(pt["field"]),
               V_random=Vr, V_standardised=Vs, V_reduction_pct=100 * (1 - Vs / Vr),
               measurement_reduction_pct=100 * (1 - (Vs - bio) / max(Vr - bio, 1e-12)),
               fields_per_well_random_to_match=k_eq,
               n_needed_random=_py(n_needed(abs(math.log(eff)), Vr, o["target_power"], o["alpha"])[0]),
               n_needed_standardised=_py(n_needed(abs(math.log(eff)), Vs, o["target_power"], o["alpha"])[0]))
    return out


def missing_reference_cost(plan: dict) -> list[dict]:
    """What n's lacking a reference cost: CI half-width and power with the n's that have it vs all n's."""
    d, out = plan["design"], []
    eff, o = focus_effect(plan), plan["options"]
    for ref, miss in d["missing_reference"].items():
        if not miss or _bio_key(ref, plan["default_class"]) not in plan["bio"]:
            continue
        have, full = d["n_count"] - len(miss), d["n_count"]
        if have < 2:
            out.append(dict(reference=ref, missing=miss, n_have=have, n_all=full, note="fewer than 2 n's left"))
            continue
        a, b = (evaluate(plan, n, d["wells"], d["fields_per_well"], eff, reference=ref, alpha=o["alpha"],
                         target_power=o["target_power"]) for n in (have, full))
        out.append(dict(reference=ref, missing=miss, n_have=have, n_all=full, effect=eff,
                        t_critical_have=a["t_critical"], t_critical_all=b["t_critical"],
                        ci_halfwidth_pct_have=a["ci_halfwidth_pct"]["mid"], ci_halfwidth_pct_all=b["ci_halfwidth_pct"]["mid"],
                        power_have=a["power"]["mid"], power_all=b["power"]["mid"]))
    return out


# ================================================================ run: write everything into sampling/gN/

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

# descriptions of every file in sampling/gN/ for manifest.json: name -> (kind, title, description)
FILES = {
    "REPORT.txt": ("report", "Report", "Plain-language report: recommendations for the next experiment, then the sampling analysis."),
    "planner.json": ("data", "Planner model", "Variance components with confidence intervals, current design, observed effects and costs; input for the interactive planner."),
    "design_forward.csv": ("table", "Forward design table", "Predicted power and 95% CI width (with bands) for n, wells and fields, per effect."),
    "mde.csv": ("table", "Minimum detectable effect", "Smallest ratio detectable with the target power for every design, with bands."),
    "cost_frontier.csv": ("table", "Cost frontier", "Designs that give the most power for their cost, per planning effect."),
    "allocation.csv": ("table", "Reference well allocation", "Variance and n's needed for splits of the plate between treatment and reference wells."),
    "more_n.csv": ("table", "How many more n's", "Predicted chance that the final test reaches p < alpha after adding k more n's (planning only)."),
    "planner_components.csv": ("table", "Planning variance components", "Field, position, well and biological variance with 95% and 80% intervals."),
    "variance_components.png": ("plot", "Where the noise comes from", "Each variance component with its 95% interval, and how they add up for one n's ratio."),
    "power_vs_n.png": ("plot", "Power vs number of n's", "Chance of detecting each planning effect as n's are added, with the uncertainty band."),
    "ci_width_vs_fields.png": ("plot", "CI width vs fields per well", "How wide the 95% CI of the ratio is for more fields per well, for several numbers of n's."),
    "cost_frontier.png": ("plot", "Power for the effort", "The most powerful design for each total effort, per planning effect."),
    "more_n.png": ("plot", "Adding n's", "Predicted chance of p < alpha per treatment if k more n's are added (not a stopping rule)."),
    "variance_components.csv": ("table", "Variance components (sampling)", "Field and well CV^2 from the sampling analysis."),
    "effects.csv": ("table", "Observed effects", "Observed ratios, p values and biological variance per treatment and reference."),
    "design_power.csv": ("table", "Power for smaller designs", "Predicted power at the current n for fewer wells or fields, per treatment."),
    "subsampling.csv": ("table", "Subsampling real images", "How often a subsample of the real images gives the same conclusion."),
    "summary_by_design.csv": ("table", "Summary by design", "Median noise, n's needed and subsampling agreement per smaller design."),
    "batch_by_condition.csv": ("table", "Batch checks per condition", "Per-well means, ratios and flags for wells that disagree."),
    "batch_by_n.csv": ("table", "Batch checks per n", "Reference levels, missing references and effect-size shifts per n."),
    "images_needed.png": ("plot", "Images needed", "Noise, subsampling agreement and per-n error for fewer images per condition."),
}


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _save(fig, out: Path, stem: str):
    fig.savefig(out / f"{stem}.png", dpi=170)
    fig.savefig(out / f"{stem}.pdf")
    import matplotlib.pyplot as plt
    plt.close(fig)


def _plots(plan, out, mde, fwd, front, mn):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d, o = plan["design"], plan["options"]
    W, K, n_now = d["wells"], d["fields_per_well"], d["n_count"]
    main, cls = d["main_reference"], plan["default_class"]
    target = o["target_power"]
    effs = sorted(o["effects"], key=lambda e: -abs(math.log(e)))

    # 1. variance components with 95% intervals
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.3, 1]))
    items = []
    for name, lab in (("field", "field-to-field"), ("position", "field position"), ("well", "well-to-well")):
        c = plan["components"][name]
        if c["var"] is not None:
            items.append((lab, c["var"], c["ci95"]))
    for key, b in plan["bio"].items():
        ref, cl = key.split("|", 1)
        items.append((f"biological: {cl} vs {ref}", b["var"], b["ci95"]))
    y = np.arange(len(items))[::-1]
    cap = max(1.0, 1.4 * max(np.sqrt(v) for _, v, _ in items)) if items else 1.0
    for yi, (lab, v, ci) in zip(y, items):
        col = PALETTE[0] if lab.startswith("biological") else PALETTE[1]
        lo, hi = np.sqrt(ci)
        a1.plot([lo, min(hi, cap)], [yi, yi], color=col, lw=2.2, solid_capstyle="round")
        a1.plot(np.sqrt(v), yi, "o", color=col, ms=8, mec="white", mew=1.5)
        txt = f"{np.sqrt(v):.2f}  [{lo:.2f}-{hi:.2f}]"
        if hi > cap:
            a1.annotate("", (cap * 1.04, yi), (cap * 0.97, yi), arrowprops=dict(arrowstyle="->", color=col, lw=1.5))
        a1.text(min(hi, cap) + cap * 0.06, yi + 0.28, txt, va="center", ha="right" if hi > cap * 0.6 else "left",
                fontsize=8, color=INK2)
    a1.set_yticks(y, [i[0] for i in items], fontsize=9)
    a1.set_xlim(0, cap * 1.1)
    a1.set_ylim(-0.6, len(items) - 0.3)
    a1.set_xlabel("SD on the log scale (~ CV); line = 95% interval", color=INK2)
    a1.set_title("Variance components", loc="left", fontsize=11, color=INK)
    _style(a1)
    pt = plan["components_point"]
    key = _bio_key(main, cls)
    if key in plan["bio"]:
        bio = plan["bio"][key]["var"]
        parts = lambda std: [bio, pt["well"] * 2 / W, pt["field"] * 2 / (W * K), 0 if std else pt["position"] * 2 / (W * K)]
        labels = ["biological (n-to-n)", "wells", "fields", "field position"]
        for yi, (std, lab) in enumerate([(False, "random positions"), (True, "standardised positions")]):
            left = 0
            for j, v in enumerate(parts(std)):
                if v > 0:
                    a2.barh(yi, v, left=left, color=PALETTE[[0, 2, 1, 4][j]], edgecolor="white", lw=2, height=0.55,
                            label=labels[j] if yi == 0 or j < 3 else None)
                left += v
            a2.text(left * 1.02, yi, f"SD {np.sqrt(left):.2f}", va="center", fontsize=9, color=INK2)
        a2.set_yticks([0, 1], ["random\npositions", "standardised\npositions"], fontsize=9)
        a2.set_xlabel(f"variance of one n's log ratio vs {main}\n({W} wells x {K} fields)", color=INK2)
        h, l = a2.get_legend_handles_labels()
        uniq = dict(zip(l, h))
        a2.legend(uniq.values(), uniq.keys(), fontsize=8, frameon=False, loc="upper center", ncol=4,
                  bbox_to_anchor=(0.5, -0.25))
        a2.set_title("How the noise adds up", loc="left", fontsize=11, color=INK)
        a2.set_xlim(0, a2.get_xlim()[1] * 1.2)
        a2.set_ylim(-0.6, 1.6)
        _style(a2)
    fig.tight_layout()
    _save(fig, out, "variance_components")

    # 2. power vs n per reference, bands
    refs = [r for r in d["references"] if _bio_key(r, cls) in plan["bio"]]
    if refs:
        fig, axes = plt.subplots(1, len(refs), figsize=(5.6 * len(refs), 4.3), squeeze=False, sharey=True)
        for ax, ref in zip(axes[0], refs):
            for i, e in enumerate(effs):
                g = fwd[(fwd["reference"] == ref) & (fwd["class"] == cls) & (fwd["effect"] == f"planning {e:g}") &
                        (fwd["wells"] == W) & (fwd["fields_per_well"] == K)].sort_values("n")
                col = PALETTE[i % 8]
                ax.fill_between(g["n"], g["power_low"], g["power_high"], color=col, alpha=0.15, lw=0)
                ax.plot(g["n"], g["power_mid"], "-o", color=col, lw=2, ms=4, label=f"ratio {e:g}")
            ax.axhline(target, color=INK2, ls="--", lw=0.9)
            ax.axvline(n_now, color=INK2, ls=":", lw=0.9)
            ax.text(n_now + 0.1, 0.02, f"now: {n_now} n's", fontsize=8, color=INK2)
            ax.set(ylim=(0, 1.02), xlabel="number of n's (independent experiments)")
            ax.set_title(f"vs {ref}  ({W} wells x {K} fields)", loc="left", fontsize=11, color=INK)
            _style(ax)
        axes[0][0].set_ylabel(f"power (alpha {o['alpha']:g}); band = 80% uncertainty", color=INK2)
        axes[0][0].legend(fontsize=9, frameon=False, loc="upper left")
        fig.tight_layout()
        _save(fig, out, "power_vs_n")

    # 3. CI half-width vs fields per well
    g0 = mde[(mde["reference"] == main) & (mde["class"] == cls) & (mde["wells"] == W)]
    if len(g0):
        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        ns = [n for n in sorted({3, 4, n_now, 6, 8, 12}) if 2 <= n <= o["max_n"]]
        for i, n in enumerate(ns):
            g = g0[g0["n"] == n].sort_values("fields_per_well")
            col = PALETTE[i % 8]
            ax.fill_between(g["fields_per_well"], g["ci_halfwidth_pct_low"], np.minimum(g["ci_halfwidth_pct_high"], 2000),
                            color=col, alpha=0.10, lw=0)
            ax.plot(g["fields_per_well"], g["ci_halfwidth_pct_mid"], "-", color=col, lw=2, label=f"{n} n's")
            last = g.iloc[-1]
            ax.text(last["fields_per_well"] + 0.3, last["ci_halfwidth_pct_mid"], f"{n}", fontsize=8, color=INK2,
                    va="center")
        ax.axvline(K, color=INK2, ls=":", lw=0.9)
        ax.set_yscale("log")
        ticks = [t for t in (10, 20, 50, 100, 200, 500, 1000, 2000)
                 if g0["ci_halfwidth_pct_low"].min() * 0.8 <= t <= 2000]
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"±{v:g}%"))
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.set(xlabel=f"fields per well ({W} wells per condition)",
               ylabel="expected 95% CI half-width of the ratio")
        ax.set_title(f"CI width vs fields (vs {main}); band = 80% uncertainty", loc="left", fontsize=11, color=INK)
        ax.legend(fontsize=9, frameon=False, ncol=2)
        _style(ax)
        fig.tight_layout()
        _save(fig, out, "ci_width_vs_fields")

    # 4. cost frontier
    if len(front):
        fig, ax = plt.subplots(figsize=(8, 4.8))
        unit = plan["costs"]["unit"].split(" ")[0]
        for i, e in enumerate(effs):
            g = front[front["effect"] == f"planning {e:g}"].sort_values("cost")
            col = PALETTE[i % 8]
            ax.step(g["cost"], g["power_mid"], where="post", color=col, lw=2, label=f"ratio {e:g}")
            ch = g[g["cheapest_reaching_target"]]
            for r in ch.itertuples():
                ax.plot(r.cost, r.power_mid, "o", color=col, ms=9, mec="white", mew=1.5)
                ax.annotate(f"{r.n} n x {r.wells} w x {r.fields_per_well} f", (r.cost, r.power_mid),
                            xytext=(6, -14), textcoords="offset points", fontsize=8, color=INK)
            cur = fwd[(fwd["reference"] == main) & (fwd["class"] == cls) & (fwd["effect"] == f"planning {e:g}") &
                      (fwd["n"] == max(n_now, 2)) & (fwd["wells"] == W) & (fwd["fields_per_well"] == K)]
            if len(cur):
                ax.plot(cur["cost"], cur["power_mid"], "D", color=col, ms=7, mec=INK, mew=0.8)
        ax.axhline(target, color=INK2, ls="--", lw=0.9)
        ax.plot([], [], "D", color="white", mec=INK, label="current design")
        ax.plot([], [], "o", color=INK2, mec="white", label=f"cheapest reaching {target:.0%}")
        ax.set_xscale("log")
        ax.set(ylim=(0, 1.02), xlabel=f"total effort ({unit}; cost per n {plan['costs']['per_n']:g}, per well "
                                        f"{plan['costs']['per_well']:g}, per field {plan['costs']['per_field']:g})",
               ylabel="power (expected)")
        ax.set_title(f"Most power for the effort (vs {main}; expected noise, bands in cost_frontier.csv)", loc="left",
                     fontsize=11, color=INK)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        _style(ax)
        fig.tight_layout()
        _save(fig, out, "cost_frontier")

    # 5. more n
    if len(mn):
        refs = list(dict.fromkeys(mn["reference"]))
        names = [t for t in d["condition_names"] if t in set(mn["treatment"])]
        cmap = {t: PALETTE[i % 8] for i, t in enumerate(names)}  # colour follows the treatment in every panel
        ls = {t: ("-" if i < 8 else "--") for i, t in enumerate(names)}
        fig, axes = plt.subplots(1, len(refs), figsize=(5.4 * len(refs) + 1.8, 4.5), squeeze=False, sharey=True)
        for ax, ref in zip(axes[0], refs):
            g = mn[mn["reference"] == ref]
            for t, gt in g.groupby("treatment", sort=False):
                ax.plot(gt["n_total"], gt["prob_p_below_alpha"], ls[t], marker="o", color=cmap[t], lw=1.8, ms=3.5,
                        label=t)
            ax.axhline(o["target_power"], color=INK2, ls="--", lw=0.8)
            ax.set(ylim=(0, 1.02), xlabel="total n's in the final analysis")
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            ax.set_title(f"vs {ref} (now {int(g['n_now'].max())} n's)", loc="left", fontsize=11, color=INK)
            _style(ax)
        axes[0][0].set_ylabel(f"predicted chance of p < {o['alpha']:g}", color=INK2)
        hs = {}
        for ax in axes[0]:
            for h, l in zip(*ax.get_legend_handles_labels()):
                hs.setdefault(l, h)
        fig.legend([hs[t] for t in names if t in hs], [t for t in names if t in hs], fontsize=9, frameon=False,
                   loc="center right")
        fig.text(0.01, 0.005, "Planning information, not a stopping rule: fix the final number of n's in advance "
                 "and test once.", fontsize=8, color=INK2)
        fig.tight_layout(rect=(0, 0.03, 1 - 1.6 / (5.4 * len(refs) + 1.8), 1))
        _save(fig, out, "more_n")


def _manifest(out: Path):
    items = []
    for p in sorted(out.iterdir()):
        if p.suffix == ".pdf" or p.name == "manifest.json":
            continue
        name = p.name
        if name in FILES:
            kind, title, desc = FILES[name]
        elif name.startswith("images_needed_vs_"):
            ref = p.stem.removeprefix("images_needed_vs_")
            kind, title, desc = "plot", f"Images needed (vs {ref})", f"Noise and subsampling for fewer images, vs {ref}."
        else:
            kind = "plot" if p.suffix == ".png" else "table" if p.suffix == ".csv" else "file"
            title, desc = p.stem.replace("_", " "), ""
        item = dict(path=name, kind=kind, title=title, description=desc)
        if p.with_suffix(".pdf").exists() and p.suffix == ".png":
            item["pdf"] = p.with_suffix(".pdf").name
        items.append(item)
    (out / "manifest.json").write_text(json.dumps(dict(files=items), indent=1), encoding="utf-8")


def _r(x, nd=2):
    return "-" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}"


def _n_txt(x):
    return f">{60}" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{int(x)}"


def recommendations(plan: dict) -> list[str]:
    """The 'Recommendations for the next experiment' section of REPORT.txt (plain language)."""
    d, o, s = plan["design"], plan["options"], plan["summaries"]
    main, cls = d["main_reference"], plan["default_class"]
    W, K, N = d["wells"], d["fields_per_well"], d["n_count"]
    key = _bio_key(main, cls)
    L = ["Recommendations for the next experiment", "=" * 70]
    if key not in plan["bio"]:
        return L + ["  Not enough data to plan (need at least 2 n's with the reference and treatments).", ""]
    alpha, target = o["alpha"], o["target_power"]
    b = plan["bio"][key]
    cur = evaluate(plan, max(N, 2), W, K, focus_effect(plan), alpha=alpha, target_power=target)
    L += [f"Main comparison: each treatment vs {main}, ratio paired t-test across n's (df = n's - 1).",
          f"Now: {N} n's x {W} wells x {K} fields. Noise in one n's log ratio: SD {cur['sd_log_ratio']['mid']:.2f} "
          f"(biological {math.sqrt(b['var']):.2f}, measurement {math.sqrt(max(cur['sd_log_ratio']['mid'] ** 2 - b['var'], 0)):.2f});",
          f"plausible range {cur['sd_log_ratio']['low']:.2f}-{cur['sd_log_ratio']['high']:.2f} (80% band; few n's make it wide).",
          "",
          f"1. How many n's? ({target:.0%} power, alpha {alpha:g}, {W} wells x {K} fields; expected [optimistic-pessimistic])"]
    pt = plan["components_point"]
    std = d["positions_default"] == "standardised"
    for e in sorted(o["effects"], key=lambda e: -abs(math.log(e))):
        lo, mid, hi = _bands(plan, key, np.array([W]), np.array([K]), np.array([W]))
        nn = [n_needed(abs(math.log(e)), v, target, alpha)[0] for v in (lo, mid, hi)]
        L.append(f"   ratio {e:g} ({abs(1 - e):.0%} {'fewer' if e < 1 else 'more'}): {_n_txt(nn[1])} n's "
                 f"[{_n_txt(nn[0])}-{_n_txt(nn[2])}]")
    eff_obs = [e for e in plan["effects"] if e["reference"] == main and e["cls"] == cls and e["n_pairs"] >= 2]
    if eff_obs:
        med = float(np.exp(np.median([math.log(e["ratio"]) for e in eff_obs])))
        Vm = float(design_variance(pt, b["var"], W, K, W, std))
        L.append(f"   observed effects (median ratio {med:.2f}): {_n_txt(n_needed(abs(math.log(med)), Vm, target, alpha)[0])} "
                 "n's - but effects seen in few n's are usually overestimated")
    e0 = focus_effect(plan)
    n_inf = n_needed(abs(math.log(e0)), max(b["var"], 1e-12), target, alpha)[0]
    L += [f"   Even with unlimited images per n, ratio {e0:g} needs {_n_txt(n_inf)} n's: the number of n's is the",
          f"   bottleneck (t critical value {stats.t.ppf(1 - alpha / 2, 2):.2f} at 3 n's, {stats.t.ppf(1 - alpha / 2, 4):.2f} at 5).",
          ""]
    ms = cur["measurement_share"]
    opt = s["cost"]["optimum"]
    L.append("2. Fields or wells?")
    L.append(f"   Measurement (wells + fields) is {ms:.0%} of the noise in one ratio; the rest is n-to-n biology.")
    for (w2, k2, lab) in ((W, 2 * K, f"{2 * K} fields per well"), (W + 1, K, f"{W + 1} wells x {K} fields"),
                          (W, max(K // 2, 1), f"{max(K // 2, 1)} fields per well")):
        V2 = float(design_variance(pt, b["var"], w2, k2, w2, std))
        V0 = float(design_variance(pt, b["var"], W, K, W, std))
        L.append(f"   {lab:<22}: noise {100 * (math.sqrt(V2 / V0) - 1):+.0f}%, "
                 f"n's for ratio {e0:g}: {_n_txt(n_needed(abs(math.log(e0)), V2, target, alpha)[0])} "
                 f"(now {_n_txt(n_needed(abs(math.log(e0)), V0, target, alpha)[0])})")
    if opt["note"]:
        L.append(f"   Cost-optimal fields per well: {opt['note']}.")
    else:
        L.append(f"   Cost-optimal fields per well = sqrt((cost per well / cost per field) x (field var / well var)) = "
                 f"sqrt(({plan['costs']['per_well']:g}/{plan['costs']['per_field']:g}) x "
                 f"({pt['field']:.3f}/{pt['well']:.3f})) = {opt['fields_per_well']:.1f}"
                 + (f"; wells per n ~ {opt['wells_per_n']:.1f}" if opt["wells_per_n"] else ""))
    ch = s["cost"]["cheapest"].get(f"planning {e0:g}", {})
    for band in ("expected", "pessimistic"):
        c = ch.get(band)
        if c:
            L.append(f"   Cheapest design for {target:.0%} power at ratio {e0:g} ({band} noise): {c['n']} n's x "
                     f"{c['wells']} wells x {c['fields_per_well']} fields (effort {c['cost']:.0f})")
        else:
            L.append(f"   No design up to {o['max_n']} n's reaches {target:.0%} power at ratio {e0:g} ({band} noise).")
    L.append("")
    al = s["allocation"]
    L.append("3. Reference wells")
    if "recommended" in al:
        rec = al["recommended"]
        L += [f"   Each reference is compared with {al['comparisons_per_reference']} conditions: Dunnett's rule gives "
              f"{W} x sqrt({al['comparisons_per_reference']}) = {W * al['sqrt_k']:.1f} -> {al['dunnett_reference_wells']} "
              "wells per reference.",
              f"   Honest gain: it only shrinks the reference's measurement noise, so the variance of a ratio drops by "
              f"{al['gain_variance_pct']:.0f}%" + ("" if al["dunnett_fits"] else
                                                   f" (but {al['dunnett_wells_used']} wells do not fit a {al['plate_format']}-well plate)") + ".",
              f"   Best use of a {al['plate_format']}-well plate ({al['treatments']} treatments, {al['references']} references): "
              f"{rec['treatment_wells']} wells per treatment and {rec['reference_wells']} per reference "
              f"({rec['wells_used']} wells; variance {-rec['gain_variance_pct']:+.0f}% vs now; n's for ratio {e0:g}: "
              f"{_n_txt(rec['n_needed'])} vs {_n_txt(al['n_needed_current'])})."]
        if not rec["fits_inner_wells"]:
            L.append(f"   That uses edge wells (only {al['inner_wells']} inner wells): randomise and balance them "
                     "(kerncount platemap).")
    else:
        L.append(f"   {al.get('note', '')}")
    L.append("")
    pe = s["positions"]
    L.append("4. Field positions")
    if pe.get("estimated"):
        L += [f"   Field position explains {pe['share_of_field_variance']:.0%} of field-to-field variance (the pattern "
              "differs from day to day).",
              f"   Image the same positions in every well (standardised), or adjust for position: field SD "
              f"{pe['field_sd_random']:.2f} -> {pe['field_sd_standardised']:.2f}, ratio noise {-pe['V_reduction_pct']:+.0f}% "
              f"in variance; with random positions you would need {pe['fields_per_well_random_to_match'] or 'many'} fields "
              f"per well for the same precision."]
    elif not pe.get("meaningful"):
        L.append("   Positions are random ([design] field_positions): fixing the stage positions would let the "
                 "analysis remove the position pattern (on comparable data it was ~35% of field noise).")
    else:
        L.append("   The position effect could not be estimated from this design.")
    L.append("")
    L.append("5. Missing references")
    mr = s["missing_reference"]
    if not mr:
        L.append("   Every n has every reference.")
    for m in mr:
        if "note" in m:
            L.append(f"   {', '.join(m['missing'])} lack {m['reference']}: {m['note']}.")
            continue
        L += [f"   {', '.join(x.split('_')[0] for x in m['missing'])} lack{'s' if len(m['missing']) == 1 else ''} "
              f"{m['reference']}, so comparisons vs "
              f"{m['reference']} use {m['n_have']} instead of {m['n_all']} n's:",
              f"   t critical {m['t_critical_have']:.2f} instead of {m['t_critical_all']:.2f}; CI half-width "
              f"±{m['ci_halfwidth_pct_have']:.0f}% instead of ±{m['ci_halfwidth_pct_all']:.0f}%; power at ratio "
              f"{m['effect']:g}: {m['power_have']:.0%} instead of {m['power_all']:.0%}. Put every reference on every plate."]
    L.append("")
    mn = s["more_n"]
    L.append("6. Adding n's (planning information, NOT a stopping rule)")
    tab = pd.DataFrame(plan.get("_more_n_rows", []))
    if len(tab):
        for ref, g in tab.groupby("reference", sort=False):
            parts = []
            for t, gt in g.groupby("treatment", sort=False):
                p = {int(r.k_added): r.prob_p_below_alpha for r in gt.itertuples()}
                parts.append(f"{t} {p.get(1, float('nan')):.0%}/{p.get(4, float('nan')):.0%}")
            L.append(f"   vs {ref}, chance of p < {alpha:g} after +1 / +4 n's: " + "; ".join(parts))
    L += [f"   Optional stopping: if there were no effect at all, testing after every new n (from {mn['start_n']} up to "
          f"{mn['start_n'] + mn['k_max']} n's) and stopping at p < {alpha:g} gives a false positive in "
          f"{mn['false_positive_optional_stopping']:.0%} of series instead of {alpha:.0%}.",
          "   Decide the final number of n's now (e.g. from section 1), run them all, then test once.", ""]
    L += ["7. Caution",
          f"   All of this rests on {N} n's. The biological variance has only ~{b['df']:.0f} degrees of freedom, so its",
          "   95% interval is wide (see planner.json / variance_components.png); use the pessimistic numbers when",
          "   the experiment matters. Classes pool treatments to stabilise the estimate"
          f" ({plan['class_source']}).",
          "   The planner separates the field-position pattern from field noise, so its field / well SDs"
          f" ({math.sqrt(pt['field']):.2f} / {math.sqrt(pt['well']):.2f})",
          "   differ from section 1 of the sampling report below, which counts position as field noise.",
          "   CI half-widths are given as +x% on the ratio (upper limit = ratio x (1 + x/100); lower = ratio / the same).",
          ""]
    return L


def run(df: pd.DataFrame, settings: dict, out: Path, project_dir=None, seed: int = 0,
        metric_label: str = "nuclei per field", sims: int = 4000) -> dict:
    """Build the plan from df (n, treatment, image, value; well/position optional) and write planner.json,
    planner_components.csv, design_forward.csv, mde.csv, cost_frontier.csv, allocation.csv, more_n.csv, the plots,
    manifest.json, and prepend 'Recommendations for the next experiment' to REPORT.txt (kept otherwise)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    plan = build(df, settings, project_dir, seed=seed, metric_label=metric_label)
    fwd, mde = forward_tables(plan)
    front, cost = cost_analysis(plan, fwd) if len(fwd) else (pd.DataFrame(), dict(cheapest={}, optimum=dict(note="no data")))
    alloc, alloc_sum = allocation(plan)
    mn, mn_sum = more_n(plan, sims=sims, seed=seed)
    plan["summaries"] = _jsonable(dict(cost=cost, allocation=alloc_sum, more_n=mn_sum, positions=position_effect(plan),
                                       missing_reference=missing_reference_cost(plan), focus_effect=focus_effect(plan)))
    comp_rows = []
    for name in ("field", "position", "well"):
        c = plan["components"][name]
        comp_rows.append(dict(component=name, reference="", cls="", variance=c["var"], sd=c["sd"], df=c["df"],
                              ci95_low=(c["ci95"] or [None])[0], ci95_high=(c["ci95"] or [None, None])[1],
                              ci80_low=(c["ci80"] or [None])[0], ci80_high=(c["ci80"] or [None, None])[1], note=c["note"]))
    for key, b in plan["bio"].items():
        ref, cl = key.split("|", 1)
        comp_rows.append(dict(component="biological", reference=ref, cls=cl, variance=b["var"], sd=math.sqrt(b["var"]),
                              df=b["df"], ci95_low=b["ci95"][0], ci95_high=b["ci95"][1], ci80_low=b["ci80"][0],
                              ci80_high=b["ci80"][1], note=f"pooled over {', '.join(b['treatments'])}"))
    pd.DataFrame(comp_rows).rename(columns={"cls": "class"}).to_csv(out / "planner_components.csv", index=False)
    f4 = "%.4g"
    fwd.to_csv(out / "design_forward.csv", index=False, float_format=f4)
    mde.to_csv(out / "mde.csv", index=False, float_format=f4)
    front.to_csv(out / "cost_frontier.csv", index=False, float_format=f4)
    alloc.to_csv(out / "allocation.csv", index=False, float_format=f4)
    mn.to_csv(out / "more_n.csv", index=False, float_format=f4)
    if len(fwd):
        _plots(plan, out, mde, fwd, front, mn)
    plan["_more_n_rows"] = mn.to_dict("records")
    rec = recommendations(plan)
    plan.pop("_more_n_rows")
    (out / "planner.json").write_text(json.dumps(_jsonable(plan), indent=None, separators=(",", ":")),
                                      encoding="utf-8")
    rp = out / "REPORT.txt"
    old = rp.read_text(encoding="utf-8") if rp.exists() else ""
    marker = "Recommendations for the next experiment"
    if old.startswith(marker):  # re-run: replace the previous recommendations
        cut = old.find("\n\n\n")
        old = old[cut + 3:] if cut >= 0 else ""
    rp.write_text("\n".join(rec) + "\n\n" + old, encoding="utf-8")
    _manifest(out)
    return plan
