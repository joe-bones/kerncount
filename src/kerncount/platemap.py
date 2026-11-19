"""Randomised plate maps for the next experiment(s): which condition goes into which well of each plate.

One plate = one n (CONTEXT.md). Every plate is randomised separately; with `balance`, the plates are chosen so
that each condition visits different plate regions (edge / inner wells, rows, columns) across plates, and the
replicate wells of one condition are kept apart. Unused wells stay empty (medium only).

    plates = generate(["Control", "LyoVec", "miR-1"], plate_format=24, wells_per_condition={"Control": 4, ...},
                      n_plates=3, avoid_edge=True, seed=1)
    render(plates[0], "plate_1")            # plate_1.png + plate_1.pdf, printable on A4
    edge_usage(plates)                      # edge-well use per condition

CLI: kerncount platemap <experiment> --plates 3 --format 24 [--wells 2] [--reference-wells 4] [--avoid-edge]
[--seed 1] writes analysis/<experiment>/plans/platemap_<timestamp>/.
"""
from __future__ import annotations

import math
import string
import time
from pathlib import Path

import numpy as np
import pandas as pd

FORMATS = {6: (2, 3), 12: (3, 4), 24: (4, 6), 48: (6, 8), 96: (8, 12)}  # wells -> (rows, columns)
COLUMNS = ["plate", "plate_well", "row", "col", "condition", "replicate", "edge"]
EMPTY = ""  # condition of an unused well (medium only)
# categorical palette (validated, colour-blind safe on adjacent pairs); extra conditions get lighter tints
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def layout(plate_format: int) -> tuple[int, int]:
    if plate_format not in FORMATS:
        raise ValueError(f"plate format {plate_format}: expected one of {', '.join(map(str, FORMATS))}")
    return FORMATS[plate_format]


def well_name(r: int, c: int) -> str:
    """0-based row/column -> "B3"."""
    return f"{string.ascii_uppercase[r]}{c + 1}"


def all_wells(plate_format: int) -> list[tuple[int, int]]:
    R, C = layout(plate_format)
    return [(r, c) for r in range(R) for c in range(C)]


def is_edge(plate_format: int, r: int, c: int) -> bool:
    R, C = layout(plate_format)
    return r in (0, R - 1) or c in (0, C - 1)


def inner_wells(plate_format: int) -> list[tuple[int, int]]:
    return [(r, c) for r, c in all_wells(plate_format) if not is_edge(plate_format, r, c)]


def _counts(conditions, wells_per_condition) -> dict[str, int]:
    if not conditions:
        raise ValueError("no conditions given")
    if len(set(conditions)) != len(conditions):
        raise ValueError("conditions must be unique")
    if isinstance(wells_per_condition, dict):
        missing = [c for c in conditions if c not in wells_per_condition]
        if missing:
            raise ValueError(f"wells_per_condition has no entry for {', '.join(missing)}")
        out = {c: int(wells_per_condition[c]) for c in conditions}
    else:
        out = {c: int(wells_per_condition) for c in conditions}
    bad = [c for c, n in out.items() if n < 1]
    if bad:
        raise ValueError(f"every condition needs at least 1 well ({', '.join(bad)})")
    return out


def check_fits(conditions, plate_format: int, wells_per_condition) -> dict:
    """Raise ValueError with a clear message if the wells requested do not fit one plate; else return a summary."""
    R, C = layout(plate_format)
    counts = _counts(list(conditions), wells_per_condition)
    need = sum(counts.values())
    if need > plate_format:
        detail = ", ".join(f"{c} {n}" for c, n in counts.items())
        raise ValueError(f"{need} wells requested ({detail}) but a {plate_format}-well plate has only "
                         f"{plate_format}: use fewer wells per condition or a larger plate format")
    return dict(wells_needed=need, wells_total=plate_format, inner_wells=len(inner_wells(plate_format)), counts=counts)


def _score(assign, plate_format, cum, conds, edge_frac):
    """Lower = better. assign: {(r, c): condition}. cum: cumulative per-condition exposure from earlier plates."""
    R, C = layout(plate_format)
    ci = {c: i for i, c in enumerate(conds)}
    edge = cum["edge"].copy()
    rows = cum["rows"].copy()
    cols = cum["cols"].copy()
    tot = cum["total"].copy()
    adj = 0
    for (r, c), cond in assign.items():
        i = ci[cond]
        edge[i] += is_edge(plate_format, r, c)
        rows[i, r] += 1
        cols[i, c] += 1
        tot[i] += 1
        for dr, dc in ((0, 1), (1, 0)):
            if assign.get((r + dr, c + dc)) == cond:
                adj += 1
    s_edge = float(((edge - edge_frac * tot) ** 2).sum())
    s_rows = float(((rows - tot[:, None] / R) ** 2).sum()) / R
    s_cols = float(((cols - tot[:, None] / C) ** 2).sum()) / C
    return 4 * s_edge + s_rows + s_cols + 2 * adj


def generate(conditions, plate_format: int = 24, wells_per_condition=2, n_plates: int = 1, avoid_edge: bool = False,
             seed: int | None = None, balance: bool = True, candidates: int = 400,
             negative_control: str | None = None, negative_control_wells: int = 2) -> list[pd.DataFrame]:
    """One randomised plate map per plate (n).

    conditions           condition names, in the order used for colours and the legend
    plate_format         6, 12, 24, 48 or 96 wells
    wells_per_condition  an int (same for all) or {condition: wells} (e.g. more wells for references)
    n_plates             number of plates (= n's) to plan
    avoid_edge           use inner wells first; edge wells are used only when the inner ones are not enough
                         (then the edge wells are shared out evenly, see edge_usage)
    seed                 random seed (None = a new random seed; it is recorded in attrs["seed"])
    balance              choose each plate (from `candidates` random layouts) so that conditions rotate through
                         edge / inner wells, rows and columns across plates and replicates are not neighbours;
                         False = every plate is an independent random permutation
    negative_control     name of a negative-control condition (secondary-only / isotype / unstained) that gets
                         `negative_control_wells` (1 or 2) wells on every plate; it is added to `conditions` if missing.
                         Its wells set that plate's positivity cut-off (docs/intensity.md).
    Returns one DataFrame per plate with columns plate, plate_well ("B3"), row ("B"), col (1-based), condition
    ("" = empty, medium only), replicate (1.. per condition, 0 for empty), edge (bool); every well of the plate
    is listed. attrs: conditions, seed, plate_format, n_plates, avoid_edge, balance, notes.
    Raises ValueError if the wells do not fit."""
    conditions = list(conditions)
    if negative_control:
        nw = int(negative_control_wells)
        if not 1 <= nw <= 2:
            raise ValueError("negative_control_wells must be 1 or 2")
        if not isinstance(wells_per_condition, dict):
            wells_per_condition = {c: int(wells_per_condition) for c in conditions}
        wells_per_condition = {**wells_per_condition, negative_control: nw}
        if negative_control not in conditions:
            conditions.append(negative_control)
    info = check_fits(conditions, plate_format, wells_per_condition)
    counts = info["counts"]
    need = info["wells_needed"]
    if seed is None:
        seed = int(np.random.default_rng().integers(0, 2 ** 31 - 1))
    rng = np.random.default_rng(seed)
    R, C = layout(plate_format)
    allw = all_wells(plate_format)
    inner = inner_wells(plate_format)
    edges = [w for w in allw if w not in set(inner)]
    notes = []
    if avoid_edge and need > len(inner):
        notes.append(f"{need} wells needed but the plate has only {len(inner)} inner wells: "
                     f"{need - len(inner)} edge wells are used on every plate (shared out evenly between conditions)")
    labels = [c for c in conditions for _ in range(counts[c])]
    n_edge_used = (max(need - len(inner), 0) if avoid_edge else None)
    edge_frac = (n_edge_used / need) if avoid_edge else len(edges) / len(allw)
    cum = dict(edge=np.zeros(len(conditions)), rows=np.zeros((len(conditions), R)),
               cols=np.zeros((len(conditions), C)), total=np.zeros(len(conditions)))
    if not avoid_edge:
        edge_frac = None  # computed per candidate from the wells actually used

    def candidate():
        if avoid_edge:
            used = inner + [edges[i] for i in rng.choice(len(edges), n_edge_used, replace=False)] if n_edge_used \
                else [inner[i] for i in rng.choice(len(inner), need, replace=False)]
        else:
            used = [allw[i] for i in rng.choice(len(allw), need, replace=False)]
        perm = rng.permutation(len(labels))
        return {w: labels[perm[i]] for i, w in enumerate(used)}

    plates = []
    for p in range(n_plates):
        if balance:
            best, best_s = None, math.inf
            for _ in range(candidates):
                a = candidate()
                ef = edge_frac if edge_frac is not None else sum(is_edge(plate_format, *w) for w in a) / need
                s = _score(a, plate_format, cum, conditions, ef)
                if s < best_s:
                    best, best_s = a, s
            a = best
        else:
            a = candidate()
        ci = {c: i for i, c in enumerate(conditions)}
        for (r, c), cond in a.items():
            i = ci[cond]
            cum["edge"][i] += is_edge(plate_format, r, c)
            cum["rows"][i, r] += 1
            cum["cols"][i, c] += 1
            cum["total"][i] += 1
        rows, rep = [], {c: 0 for c in conditions}
        for r, c in allw:  # reading order: replicate numbers run left to right, top to bottom
            cond = a.get((r, c), EMPTY)
            if cond != EMPTY:
                rep[cond] += 1
            rows.append(dict(plate=p + 1, plate_well=well_name(r, c), row=string.ascii_uppercase[r], col=c + 1,
                             condition=cond, replicate=rep[cond] if cond != EMPTY else 0,
                             edge=is_edge(plate_format, r, c)))
        df = pd.DataFrame(rows, columns=COLUMNS)
        df.attrs.update(conditions=conditions, seed=seed, plate_format=plate_format, n_plates=n_plates,
                        avoid_edge=avoid_edge, balance=balance, notes=list(notes))
        plates.append(df)
    return plates


def edge_usage(plates: list[pd.DataFrame]) -> pd.DataFrame:
    """Per condition: wells over all plates, edge wells, edge share, and the number of distinct rows / columns /
    well positions it occupied (more = better spread across plate regions)."""
    allp = pd.concat(plates, ignore_index=True)
    allp = allp[allp["condition"] != EMPTY]
    conds = plates[0].attrs.get("conditions", sorted(allp["condition"].unique()))
    out = []
    for c in conds:
        g = allp[allp["condition"] == c]
        row = dict(condition=c, wells=len(g), edge_wells=int(g["edge"].sum()),
                   edge_pct=100 * g["edge"].mean() if len(g) else np.nan, rows_used=g["row"].nunique(),
                   cols_used=g["col"].nunique(), distinct_positions=g["plate_well"].nunique())
        for p, gp in g.groupby("plate"):
            row[f"plate{p}_edge"] = int(gp["edge"].sum())
        out.append(row)
    return pd.DataFrame(out)


# ---------------------------------------------------------------- rendering

def colors_for(conditions) -> dict[str, str]:
    """Fixed colour per condition (palette order; beyond 8, lighter tints of the same hues)."""
    out = {}
    for i, c in enumerate(conditions):
        base = PALETTE[i % len(PALETTE)]
        k = i // len(PALETTE)
        if k:
            rgb = np.array([int(base[j:j + 2], 16) for j in (1, 3, 5)], float)
            rgb = rgb + (255 - rgb) * min(0.45 * k, 0.8)
            base = "#" + "".join(f"{int(v):02x}" for v in rgb)
        out[c] = base
    return out


def _text_color(hex_color: str) -> str:
    r, g, b = (int(hex_color[j:j + 2], 16) / 255 for j in (1, 3, 5))
    lin = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in (r, g, b)]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#111111" if lum > 0.33 else "#ffffff"


def _wrap(name: str, width: int) -> str:
    if len(name) <= width:
        return name
    parts, line, lines = name.replace("-", "-\u0000").split("\u0000"), "", []
    for p in parts:
        if len(line) + len(p) > width and line:
            lines.append(line)
            line = p
        else:
            line += p
    lines.append(line)
    return "\n".join(lines[:3])


def render(plate_df: pd.DataFrame, out_path, title: str | None = None, subtitle: str | None = None,
           colors: dict | None = None) -> list[Path]:
    """Draw one plate map as a printable A4-landscape diagram (rows A.., columns 1.., one colour per condition,
    condition name / code and replicate number in each well, legend, title with plate number and seed).
    out_path: file stem or a .png/.pdf path; both <stem>.png and <stem>.pdf are written. Returns the paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyBboxPatch

    at = plate_df.attrs
    fmt = int(at.get("plate_format", len(plate_df)))
    R, C = layout(fmt)
    conds = list(at.get("conditions") or [c for c in plate_df["condition"].unique() if c != EMPTY])
    colors = colors or colors_for(conds)
    codes = fmt >= 48  # small wells: a short code in the well, names in the legend
    plate_no = int(plate_df["plate"].iloc[0])
    if title is None:
        title = f"Plate {plate_no} of {at.get('n_plates', '?')}  ·  {fmt}-well  ·  seed {at.get('seed', '?')}"

    fig = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
    ax = fig.add_axes([0.04, 0.10, 0.66, 0.74])
    lg = fig.add_axes([0.72, 0.10, 0.26, 0.74])
    lg.axis("off")
    ax.set_xlim(-0.9, C - 0.1)
    ax.set_ylim(R - 0.1, -0.95)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((-0.6, -0.65), C + 0.2, R + 0.3, boxstyle="round,pad=0.02,rounding_size=0.25",
                                fc="#f6f6f4", ec="#9a9994", lw=1.5))
    rad = 0.43
    fs_lab = {6: 30, 12: 24, 24: 17, 48: 12, 96: 9}[fmt]
    fs_in = {6: 15, 12: 13, 24: 10, 48: 9, 96: 7}[fmt]
    width = {6: 12, 12: 10, 24: 9, 48: 6, 96: 5}[fmt]
    for c in range(C):
        ax.text(c, -0.82, str(c + 1), ha="center", va="center", fontsize=max(fs_lab * 0.7, 11), color="#52514e")
    for r in range(R):
        ax.text(-0.78, r, string.ascii_uppercase[r], ha="center", va="center", fontsize=max(fs_lab * 0.7, 11), color="#52514e")
    for w in plate_df.itertuples():
        r, c = string.ascii_uppercase.index(w.row), w.col - 1
        if w.condition == EMPTY:
            ax.add_patch(Circle((c, r), rad, fc="white", ec="#b8b7b1", lw=1.0, ls=(0, (3, 2))))
            ax.text(c, r, "empty", ha="center", va="center", fontsize=fs_in * 0.8, color="#9a9994")
            continue
        col = colors[w.condition]
        ax.add_patch(Circle((c, r), rad, fc=col, ec="white", lw=2))
        name = f"{conds.index(w.condition) + 1}" if codes else _wrap(w.condition, width)
        ax.text(c, r, f"{name}\n#{w.replicate}", ha="center", va="center", fontsize=fs_in,
                color=_text_color(col), fontweight="bold", linespacing=1.1)
    fig.text(0.04, 0.93, title, fontsize=18, fontweight="bold", color="#0b0b0b")
    sub = subtitle or ("Randomised layout. Empty wells: medium only."
                       + (" Inner wells used first (avoid edge)." if at.get("avoid_edge") else ""))
    fig.text(0.04, 0.89, sub, fontsize=11, color="#52514e")
    # legend
    used = plate_df[plate_df["condition"] != EMPTY].groupby("condition").size()
    y = 0.98
    lg.text(0, y, "Conditions", fontsize=13, fontweight="bold", va="top", transform=lg.transAxes)
    step = min(0.075, 0.9 / max(len(conds) + 1, 1))
    for i, cname in enumerate(conds):
        y = 0.98 - (i + 1) * step - 0.02
        lg.scatter([0.06], [y], s=260, c=[colors[cname]], edgecolors="white", transform=lg.transAxes, clip_on=False)
        lab = f"{i + 1}  {cname}" if codes else cname
        lg.text(0.14, y, f"{lab}  ({int(used.get(cname, 0))} wells)", fontsize=11, va="center",
                transform=lg.transAxes, color="#0b0b0b")
    n_empty = int((plate_df["condition"] == EMPTY).sum())
    if n_empty:
        y -= step
        lg.scatter([0.06], [y], s=260, c="white", edgecolors="#9a9994", linestyles="--", transform=lg.transAxes, clip_on=False)
        lg.text(0.14, y, f"empty, medium only ({n_empty} wells)", fontsize=11, va="center", transform=lg.transAxes,
                color="#52514e")
    lg.text(0, 0.0, "#k = replicate well k of that condition", fontsize=9, color="#52514e", transform=lg.transAxes)
    fig.text(0.04, 0.04, f"kerncount platemap  ·  {time.strftime('%Y-%m-%d')}  ·  write the date and n here: "
             "________________", fontsize=9, color="#52514e")
    out_path = Path(out_path)
    stem = out_path.with_suffix("") if out_path.suffix.lower() in (".png", ".pdf") else out_path
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = [stem.with_suffix(".png"), stem.with_suffix(".pdf")]
    fig.savefig(paths[0], dpi=200)
    fig.savefig(paths[1])
    plt.close(fig)
    return paths


# ---------------------------------------------------------------- project command

def project_conditions(project_dir: Path, settings: dict) -> list[str]:
    """Treatments of a project (counts.csv, else layout.csv), in [analysis] treatment_order."""
    from . import analysis
    from .discover import read_table

    project_dir = Path(project_dir)
    names: list[str] = []
    for f in ("counts.csv", "layout.csv"):
        p = project_dir / f
        if p.exists():
            t = read_table(p)
            if "treatment" in t:
                if "include" in t:
                    t = t[t["include"].astype(str).str.strip().str.lower().isin(["yes", "true", "1", ""])]
                names = sorted({str(x).strip() for x in t["treatment"] if str(x).strip() and str(x) != "nan"})
                if names:
                    break
    return analysis.treatment_order(names, settings["analysis"]["treatment_order"]) if names else []


def reference_wells_auto(wells: int, n_conditions: int) -> int:
    """Dunnett's rule: reference wells = wells x sqrt(comparisons), comparisons = conditions - 1."""
    return max(int(round(wells * math.sqrt(max(n_conditions - 1, 1)))), wells)


def plan_project(project_dir: Path, plates: int = 3, plate_format: int | None = None, wells: int = 2,
                 reference_wells=None, avoid_edge: bool | None = None, seed: int | None = None,
                 balance: bool = True, conditions: list[str] | None = None, negative_control: str | None = None,
                 negative_control_wells: int = 2) -> Path:
    """Write analysis/<experiment>/plans/platemap_<timestamp>/ with plate_K.csv/.png/.pdf, all_plates.csv,
    edge_usage.csv and README.txt. reference_wells: int, "auto" (Dunnett sqrt(k), reduced until the plate fits)
    or None ([planner] reference_wells). negative_control: name of the negative-control condition to reserve
    negative_control_wells wells for on every plate; None = [experiment] negative_control ("" = none).
    Returns the folder."""
    from . import settings as settings_mod
    from .planner import options

    project_dir = Path(project_dir)
    s = settings_mod.load(project_dir)
    opt = options(s)
    fmt = int(plate_format or opt["plate_format"])
    avoid = opt["avoid_edge"] if avoid_edge is None else avoid_edge
    conds = conditions or project_conditions(project_dir, s)
    if not conds:
        raise ValueError("no conditions found (counts.csv / layout.csv); pass them explicitly")
    nc = (s.get("experiment", {}).get("negative_control") or "") if negative_control is None else negative_control
    nc = str(nc).strip()
    conds = [c for c in conds if c != nc]
    a = s["analysis"]
    refs = [r for r in dict.fromkeys([a["reference"]] + list(a["normalise_to"])) if r in conds]
    rw = opt["reference_wells"] if reference_wells is None else reference_wells
    notes = []
    auto = isinstance(rw, str) and rw.lower() == "auto"
    if auto:
        rw = reference_wells_auto(wells, len(conds))
        while rw > wells and (len(conds) - len(refs)) * wells + len(refs) * rw + (negative_control_wells if nc else 0) > fmt:
            rw -= 1
        notes.append(f"reference wells: Dunnett's sqrt(k) rule ({wells} x sqrt({len(conds) - 1}) = "
                     f"{wells * math.sqrt(max(len(conds) - 1, 1)):.1f}), {rw} after fitting the {fmt}-well plate")
    counts = {c: (int(rw) if c in refs else wells) for c in conds}
    if nc:
        notes.append(f"negative control '{nc}': {negative_control_wells} well(s) per plate; it sets each plate's "
                     "positivity cut-off and is left out of the statistics")
    maps = generate(conds, fmt, counts, plates, avoid, seed, balance, negative_control=nc or None,
                    negative_control_wells=negative_control_wells)
    counts = {**counts, **({nc: int(negative_control_wells)} if nc else {})}
    conds = maps[0].attrs["conditions"]
    seed = maps[0].attrs["seed"]
    out = project_dir / "plans" / f"platemap_{time.strftime('%Y%m%d-%H%M%S')}"
    out.mkdir(parents=True, exist_ok=False)
    delim = s["project"].get("csv_delimiter", ",")
    colors = colors_for(conds)
    for m in maps:
        p = int(m["plate"].iloc[0])
        m.to_csv(out / f"plate_{p}.csv", index=False, sep=delim)
        render(m, out / f"plate_{p}", colors=colors)
    pd.concat(maps, ignore_index=True).to_csv(out / "all_plates.csv", index=False, sep=delim)
    eu = edge_usage(maps)
    eu.to_csv(out / "edge_usage.csv", index=False, sep=delim)
    notes += maps[0].attrs["notes"]
    L = [f"Plate maps for {project_dir.name}", "=" * 60,
         f"{plates} plate(s) (one per n), {fmt}-well, seed {seed}"
         f"{', balanced across plates' if balance else ''}{', inner wells first' if avoid else ''}.",
         "Wells per condition: " + ", ".join(f"{c} {n}" for c, n in counts.items()),
         f"Wells used per plate: {sum(counts.values())} of {fmt}; the rest stay empty (medium only).", ""]
    L += [f"NOTE: {x}" for x in notes]
    L += ["", "Edge-well use per condition over all plates:"]
    L += ["  " + x for x in eu.round(1).to_string(index=False).splitlines()]
    L += ["", "Files: plate_K.png / .pdf (print on A4), plate_K.csv (plate, plate_well, row, col, condition, replicate,",
          "edge), all_plates.csv, edge_usage.csv.",
          "Re-create exactly the same maps with the same seed:",
          f"  kerncount platemap {project_dir.name} --plates {plates} --format {fmt} --wells {wells} "
          f"--reference-wells {rw}{' --avoid-edge' if avoid else ''} --seed {seed}"
          f"{f' --negative-control {chr(34)}{nc}{chr(34)} --control-wells {negative_control_wells}' if nc else ''}"
          f"{'' if balance else ' --no-balance'}",
          "",
          "Why randomise: plate position effects (edge evaporation, temperature, imaging order) then spread over",
          "the conditions instead of biasing one of them. Use the map for seeding, treating and imaging.",
          "Note: the image numbering (which images belong to which well) still has to match [design] in",
          "settings.toml (or wells.csv) when the images are analysed."]
    (out / "README.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    return out
