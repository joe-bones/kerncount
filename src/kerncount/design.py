"""Plate design: which well and which field position every image belongs to, and the timepoint of every n.

This is the single place that knows the design. Everything else calls `assign` (adds well / position / timepoint)
and `describe` (one row per n x treatment with flags). See docs/design.md.

Nesting (CONTEXT.md): field (= image) in well in condition (= treatment) in n (= one plate).

Settings, all in [design]:
    wells            int (default 2) or "auto" (= wells.csv if present, else 2)
    fields_per_well  "auto" (default; per condition ceil(max image number / wells)) or an int (fixed)
    well_layout      "consecutive" (1..k = well 1, k+1..2k = well 2, ...) | "alternating" (odd/even, round robin)
                     | "file" (wells.csv in the project folder)
    field_positions  "fixed" (default; field k is the same physical spot in every well) | "random"
Old names plates / fields_per_plate / plate_layout are still read (with a note).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

WELLS_FILE = "wells.csv"
DEFAULT_WELLS = 2
LAYOUTS = ("consecutive", "alternating", "file")
POSITIONS = ("fixed", "random")
_OLD_NAMES = {"plates": "wells", "fields_per_plate": "fields_per_well", "plate_layout": "well_layout"}
KEY = ["n", "treatment", "image"]
COLUMNS = ["well", "position", "timepoint"]  # added by assign(), in this order


@dataclass
class Spec:
    """Resolved [design] settings. `wells` / `fields_per_well` are an int or "auto"."""
    wells: int | str = DEFAULT_WELLS
    fields_per_well: int | str = "auto"
    well_layout: str = "consecutive"
    field_positions: str = "fixed"
    notes: list[str] = field(default_factory=list)


def _auto(v) -> bool:
    return isinstance(v, str) and v.strip().lower() == "auto"


def spec(settings: dict) -> Spec:
    """Read [design] (old key names accepted) and validate it. Raises ValueError on an invalid value."""
    d = dict(settings.get("design", {}))
    notes = []
    for old, new in _OLD_NAMES.items():
        if old in d:  # old names are not in the defaults, so their presence means the user's file uses them
            d[new] = d[old]
            notes.append(f"settings.toml [design] uses the old name '{old}'; please rename it to '{new}'")

    def count(key, default):
        v = d.get(key, default)
        if _auto(v):
            return "auto"
        try:
            n = int(v)
        except (TypeError, ValueError):
            raise ValueError(f"settings.toml [design] {key} = {v!r}: expected a whole number or \"auto\"") from None
        if n < 1:
            raise ValueError(f"settings.toml [design] {key} = {v!r}: must be at least 1")
        return n

    layout = str(d.get("well_layout", "consecutive")).strip().lower()
    if layout not in LAYOUTS:
        raise ValueError(f"settings.toml [design] well_layout = {layout!r}: expected one of {', '.join(LAYOUTS)}")
    pos = str(d.get("field_positions", "fixed")).strip().lower()
    if pos not in POSITIONS:
        raise ValueError(f"settings.toml [design] field_positions = {pos!r}: expected one of {', '.join(POSITIONS)}")
    return Spec(count("wells", DEFAULT_WELLS), count("fields_per_well", "auto"), layout, pos, notes)


def positions_meaningful(settings: dict) -> bool:
    """True when field position k is the same physical spot in every well ([design] field_positions = "fixed").
    Models that use the `position` column (e.g. as a fixed effect) must check this first and ignore position when
    False: the column is still filled in, but then it is only the order in which fields were taken."""
    return spec(settings).field_positions == "fixed"


# ---------------------------------------------------------------- wells.csv

def wells_file(project_dir: Path | None) -> Path | None:
    return Path(project_dir) / WELLS_FILE if project_dir is not None else None


def read_wells_file(path: Path) -> pd.DataFrame:
    """wells.csv: columns image, well, position (position may be blank) and optional n, treatment (blank = all)."""
    from .discover import read_table

    t = read_table(path)
    t.columns = [c.strip().lower() for c in t.columns]
    missing = {"image", "well"} - set(t.columns)
    if missing:
        raise ValueError(f"{path.name}: missing column(s) {', '.join(sorted(missing))} (need image, well, position)")
    for c in ("n", "treatment", "position"):
        if c not in t:
            t[c] = ""
    t = t[[s.strip() != "" for s in t["image"]]].copy()
    for c in ("n", "treatment", "position"):
        t[c] = t[c].astype(str).str.strip()
    try:
        t["image"] = t["image"].astype(int)
        t["well"] = t["well"].astype(int)
        t["position"] = [int(p) if p else 0 for p in t["position"]]  # 0 = work it out from the image order
    except ValueError as e:
        raise ValueError(f"{path.name}: image, well and position must be whole numbers ({e})") from None
    if (t["well"] < 1).any() or (t["position"] < 0).any():
        raise ValueError(f"{path.name}: well and position must be at least 1")
    return t.reset_index(drop=True)


def _n_matches(rule: str, n: str) -> bool:
    """A wells.csv n of "n1" also matches the full label "n1_2026-06-05"."""
    return rule == n or n.startswith(rule + "_")


def _from_file(rules: pd.DataFrame, n: str, treatment: str, image: int) -> tuple[int, int] | None:
    """Most specific rule wins: n + treatment > n only > treatment only > neither. Later rows win ties."""
    best, best_score = None, -1
    for r in rules[rules["image"] == image].itertuples():
        if r.n and not _n_matches(r.n, n):
            continue
        if r.treatment and r.treatment != treatment:
            continue
        score = 2 * bool(r.n) + bool(r.treatment) + (1 if r.n and r.treatment else 0)
        if score >= best_score:
            best, best_score = (int(r.well), int(r.position)), score
    return best


# ---------------------------------------------------------------- timepoint

def timepoints(layout: pd.DataFrame | None) -> dict[str, str]:
    """{n: timepoint} from layout.csv. Uses its timepoint column, or (old files without it) parses the n folder name."""
    if layout is None or len(layout) == 0 or "n" not in layout:
        return {}
    from .discover import parse_timepoint

    out: dict[str, str] = {}
    for _, r in layout.iterrows():
        n = str(r["n"])
        tp = str(r["timepoint"]).strip() if "timepoint" in layout else parse_timepoint(str(r.get("n_folder", "")))
        if tp and tp.lower() != "nan" and not out.get(n):
            out[n] = tp
        out.setdefault(n, "")
    return out


def _read_layout(project_dir: Path | None) -> pd.DataFrame | None:
    if project_dir is None or not (Path(project_dir) / "layout.csv").exists():
        return None
    from .discover import read_table

    return read_table(Path(project_dir) / "layout.csv")


# ---------------------------------------------------------------- assign

def _consecutive(image: int, k: int, w: int) -> tuple[int, int]:
    well = min((image - 1) // k + 1, w)  # images beyond w x k stay in the last well (flagged by describe)
    return well, image - (well - 1) * k


def _alternating(image: int, w: int) -> tuple[int, int]:
    return (image - 1) % w + 1, (image - 1) // w + 1


def assign(table: pd.DataFrame, settings: dict, project_dir: Path | None = None,
           layout: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return a copy of `table` (needs n, treatment, image) with columns well (int, 1-based), position (int, 1-based
    field position within its well) and timepoint (str, "" if unknown) added or replaced. Rows are neither dropped
    nor reordered and the index is kept.

    project_dir: where wells.csv (well_layout = "file", or wells = "auto") and layout.csv (timepoint) are read.
    layout: a layout table to take timepoints from instead of project_dir/layout.csv. If neither gives a timepoint
    and `table` already has a timepoint column, that column is kept.

    The result carries the resolved design in .attrs["design"] (a per-condition dict used by describe)."""
    sp = spec(settings)
    out = table.copy()
    if out.empty:
        for c, dt in zip(COLUMNS, ("int64", "int64", "object")):
            out[c] = pd.Series(dtype=dt)
        for c in metadata_columns():
            out[c] = pd.Series(dtype="object")
        out.attrs["design"] = dict(spec=sp, conditions={})
        return out
    images = out["image"].astype("int64")
    wf = wells_file(project_dir)
    use_file = sp.well_layout == "file" or (_auto(sp.wells) and wf is not None and wf.exists())
    rules = None
    if use_file:
        if wf is None or not wf.exists():
            raise ValueError(f"[design] well_layout = \"file\" but there is no {WELLS_FILE} in the project folder "
                             "(kerncount wells-template <experiment> writes one)")
        rules = read_wells_file(wf)
    W = DEFAULT_WELLS if _auto(sp.wells) else sp.wells
    fallback = "consecutive" if sp.well_layout == "file" else sp.well_layout

    img = images.to_numpy()
    wells = np.zeros(len(out), dtype="int64")
    pos = np.zeros(len(out), dtype="int64")
    conds = {}
    keys = [out["n"].astype(str).to_numpy(), out["treatment"].astype(str).to_numpy()]
    for (n, t), ix in pd.Series(np.arange(len(out))).groupby(keys, sort=False).indices.items():
        mx = int(img[ix].max())
        K = math.ceil(mx / W) if _auto(sp.fields_per_well) else sp.fields_per_well
        unmapped, blank = [], []
        for i in ix:
            hit = _from_file(rules, n, t, int(img[i])) if rules is not None else None
            if hit is None:
                if rules is not None:
                    unmapped.append(int(img[i]))
                hit = _alternating(int(img[i]), W) if fallback == "alternating" else _consecutive(int(img[i]), K, W)
            elif hit[1] == 0:
                blank.append(i)
            wells[i], pos[i] = hit
        if blank:  # wells.csv left position blank: number the fields of each well in image order
            b = np.asarray(blank)
            ranks = pd.DataFrame({"well": wells[b], "image": img[b]}).groupby("well")["image"].rank(method="first")
            pos[b] = ranks.to_numpy().astype("int64")
        c = dict(wells=W, fields_per_well=K, layout=fallback, fields_source="auto" if _auto(sp.fields_per_well) else "fixed",
                 unmapped=sorted(unmapped))
        if rules is not None:  # the file defines the design: its wells and positions are the expected ones
            c.update(wells=int(wells[ix].max()), fields_per_well=int(pos[ix].max()), layout="file", fields_source="file")
        conds[(n, t)] = c
    out["well"] = pd.Series(wells, index=out.index)
    out["position"] = pd.Series(pos, index=out.index)

    lay = layout if layout is not None else _read_layout(project_dir)
    tps = timepoints(lay)
    if tps or "timepoint" not in out:
        prev = out["timepoint"].fillna("").astype(str) if "timepoint" in out else pd.Series("", index=out.index)
        out["timepoint"] = [tps.get(str(n)) or p for n, p in zip(out["n"], prev)]
    else:
        out["timepoint"] = out["timepoint"].fillna("").astype(str)
    _add_metadata(out, settings, lay)
    out.attrs["design"] = dict(spec=sp, conditions=conds)
    return out


# ---------------------------------------------------------------- experiment metadata

META_KEYS = ["magnification", "camera_pixel_um", "seeding_density", "seeding_unit", "stain", "negative_control",
             "negative_control_type"]  # each becomes a meta_<key> column; cell_line has its own column


def metadata_columns() -> list[str]:
    """Metadata columns added by `assign` after well / position / timepoint, in this order."""
    return ["cell_line"] + [f"meta_{k}" for k in META_KEYS]


def output_columns(assigned: pd.DataFrame) -> list[str]:
    """Every design + metadata column `assign` added (for merging into per-field tables such as counts_used.csv)."""
    return [c for c in COLUMNS + metadata_columns() if c in assigned]


def _text(v) -> str:
    s = "" if v is None else str(v).strip()
    if s.lower() == "nan":
        return ""
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    return s


def metadata(settings: dict, layout: pd.DataFrame | None) -> tuple[dict, dict[tuple[str, str], dict]]:
    """(experiment defaults, {(n, treatment): overrides}) for cell_line + META_KEYS. Overrides come from
    same-named columns in layout.csv; a blank cell keeps the experiment value. A layout row is one n x treatment
    folder, so put a value in every row of an n to override it for the whole n."""
    e = settings.get("experiment", {}) or {}
    base = {k: _text(e.get(k, "")) for k in ["cell_line"] + META_KEYS}
    over: dict[tuple[str, str], dict] = {}
    if layout is not None and len(layout) and {"n", "treatment"} <= set(layout.columns):
        cols = [k for k in base if k in layout.columns]
        for _, r in layout.iterrows():
            vals = {k: _text(r[k]) for k in cols if _text(r[k])}
            if vals:
                over.setdefault((str(r["n"]), str(r["treatment"])), {}).update(vals)
    return base, over


def _add_metadata(out: pd.DataFrame, settings: dict, layout: pd.DataFrame | None) -> None:
    base, over = metadata(settings, layout)
    keys = list(zip(out["n"].astype(str), out["treatment"].astype(str)))
    for k in ["cell_line"] + META_KEYS:
        col = "cell_line" if k == "cell_line" else f"meta_{k}"
        out[col] = [over.get(key, {}).get(k, base[k]) for key in keys]


# ---------------------------------------------------------------- describe

DESCRIBE_COLUMNS = ["n", "treatment", "timepoint", "images", "wells_expected", "wells_present", "fields_per_well",
                    "fields_source", "fields_min", "fields_max", "max_image", "missing_images", "beyond_design",
                    "unbalanced", "single_well", "flags"]


def _numbers(xs) -> str:
    return "; ".join(str(x) for x in xs)


def describe(table: pd.DataFrame) -> pd.DataFrame:
    """One row per n x treatment of an `assign`ed table:
    images, wells_expected / wells_present, fields_per_well (the design size; fields_source auto | fixed | file),
    fields_min / fields_max (images per present well), max_image, missing_images (image numbers absent from
    1..max(max_image, wells x fields_per_well), "; "-separated), beyond_design (image numbers above wells x
    fields_per_well), the booleans unbalanced (wells have different numbers of fields, or an expected well is
    empty) and single_well, and `flags` (plain-language, "; "-separated; "" = nothing unusual)."""
    info = table.attrs.get("design", {}).get("conditions", {})
    rows = []
    for (n, t), g in table.groupby([table["n"].astype(str), table["treatment"].astype(str)], sort=False):
        c = info.get((n, t), {})
        W = int(c.get("wells", g["well"].max()))
        K = int(c.get("fields_per_well", g["position"].max()))
        imgs = sorted(set(int(i) for i in g["image"]))
        per_well = g.groupby("well").size()
        file = c.get("fields_source") == "file"
        expected = range(1, (imgs[-1] if file else max(imgs[-1], W * K)) + 1)
        missing = [i for i in expected if i not in set(imgs)]
        beyond = [] if file else [i for i in imgs if i > W * K]
        present = sorted(int(w) for w in per_well.index)
        unbalanced = per_well.nunique() > 1 or len(present) < W
        flags = []
        if len(imgs) == 1:
            flags.append("single image")
        if len(present) == 1:
            flags.append("single well")
        if unbalanced:
            flags.append("unbalanced: " + ", ".join(f"well {w} {int(per_well.get(w, 0))}" for w in range(1, max(W, present[-1]) + 1)))
        if missing:
            flags.append(f"missing image(s) {_numbers(missing)}")
        if beyond:
            where = f"kept in well {W}" if c.get("layout") == "consecutive" else f"positions above {K}"
            flags.append(f"image(s) {_numbers(beyond)} beyond the design ({W} wells x {K} fields): {where}")
        if c.get("unmapped"):
            flags.append(f"image(s) {_numbers(c['unmapped'])} not in {WELLS_FILE}: {c.get('layout', 'consecutive')} rule used")
        if len(g) != len(imgs):
            flags.append("duplicate image numbers")
        tp = g["timepoint"].iloc[0] if "timepoint" in g else ""
        rows.append(dict(n=n, treatment=t, timepoint=tp, images=len(imgs), wells_expected=W, wells_present=len(present),
                         fields_per_well=K, fields_source=c.get("fields_source", ""), fields_min=int(per_well.min()),
                         fields_max=int(per_well.max()), max_image=imgs[-1], missing_images=_numbers(missing),
                         beyond_design=_numbers(beyond), unbalanced=bool(unbalanced), single_well=len(present) == 1,
                         flags="; ".join(flags)))
    out = pd.DataFrame(rows, columns=DESCRIBE_COLUMNS)
    # auto fields per well read from the highest image number: a condition imaged in one well only (1-6) or missing
    # its last images looks like a smaller design. Say so when the rest of its n has a different size.
    for n, g in out[out["fields_source"] == "auto"].groupby("n"):
        if len(g) < 2:
            continue
        mode = g["fields_per_well"].mode()
        k = int(mode.max())
        if (g["fields_per_well"] == k).sum() * 2 <= len(g):
            continue  # no clear majority
        for i in g.index[g["fields_per_well"] != k]:
            r = out.loc[i]
            note = (f"fields per well read as {r.fields_per_well} (highest image {r.max_image} / {r.wells_expected} "
                    f"wells) but {k} in the other conditions of this n: if only one well was imaged or the last "
                    f"images are missing, the well split is wrong; set it in {WELLS_FILE} or fix fields_per_well")
            out.at[i, "flags"] = "; ".join(x for x in (r.flags, note) if x)
    return out


def absent(assigned: pd.DataFrame, settings: dict, layout: pd.DataFrame | None = None) -> pd.DataFrame:
    """Conditions an n does not have, as design_check rows with images = 0 and a flag (docs/design.md, "Unusual
    designs"): a condition present in another n of the same cell line but not in this one (a missing reference is
    named as such), and an included layout.csv folder that gave no images at all. Cell lines are never compared."""
    cols = DESCRIBE_COLUMNS
    if assigned.empty and (layout is None or layout.empty):
        return pd.DataFrame(columns=cols)
    a = settings.get("analysis", {})
    refs = [r for r in dict.fromkeys([a.get("reference", "")] + list(a.get("normalise_to") or [])) if r]
    neg = str((settings.get("experiment") or {}).get("negative_control", "") or "").strip()
    W = DEFAULT_WELLS if _auto(spec(settings).wells) else spec(settings).wells
    have = {(str(n), str(t)) for n, t in zip(assigned["n"], assigned["treatment"])} if len(assigned) else set()
    line = {}
    if "cell_line" in assigned and len(assigned):
        for n, g in assigned.groupby(assigned["n"].astype(str)):
            line[n] = "; ".join(sorted(set(g["cell_line"].fillna("").astype(str)) - {""}))
    tps = {}
    if "timepoint" in assigned and len(assigned):
        tps = assigned.groupby(assigned["n"].astype(str))["timepoint"].first().to_dict()
    ns = sorted({n for n, _ in have})
    rows, seen = [], set()

    def row(n, t, flag):
        seen.add((n, t))
        rows.append(dict(n=n, treatment=t, timepoint=tps.get(n, ""), images=0, wells_expected=W, wells_present=0,
                         fields_per_well=np.nan, fields_source="", fields_min=0, fields_max=0, max_image=np.nan,
                         missing_images="", beyond_design="", unbalanced=False, single_well=False, flags=flag))

    if layout is not None and len(layout) and {"n", "treatment"} <= set(layout.columns):
        inc = layout["include"].astype(str).str.strip().str.lower().isin(("yes", "y", "true", "1")) \
            if "include" in layout else pd.Series(True, index=layout.index)
        for r in layout[inc].itertuples():
            n, t = str(r.n), str(r.treatment)
            if (n, t) not in have and (n, t) not in seen:
                folder = f"{getattr(r, 'n_folder', '')}/{getattr(r, 'treatment_folder', '')}".strip("/")
                row(n, t, f"no images: folder {folder or t} is empty (or holds no numbered .tif files); "
                          f"{t} is missing from {n.split('_')[0]}")
    for t in sorted({t for _, t in have}):
        if t == neg:
            continue
        with_t = [n for n in ns if (n, t) in have]
        for n in ns:
            if (n, t) in have or (n, t) in seen:
                continue
            peers = [m for m in with_t if line.get(m, "") == line.get(n, "")]
            if not peers:
                continue  # only other cell lines have it: never compared
            short = ", ".join(m.split("_")[0] for m in peers)
            if t in refs:
                flag = (f"reference {t} absent from this n (present in {short}): comparisons vs {t} use only the "
                        f"n's that have it; this n enters them only through the other reference(s)")
            else:
                flag = (f"condition absent from this n (present in {short}): its comparisons use only the n's that "
                        "have it")
            row(n, t, flag)
    for r in rows:  # an empty folder of a condition the other n's have: say both
        if r["flags"].startswith("no images") and r["treatment"] in refs:
            r["flags"] += f" (reference {r['treatment']}: comparisons vs {r['treatment']} use only the n's that have it)"
    return pd.DataFrame(rows, columns=cols)


def summary_line(desc: pd.DataFrame, settings: dict) -> str:
    """One log line, e.g. 'Design: 4 n, 8 conditions, 2 wells x 6 fields (auto); 1 condition unbalanced
    (n1 / CL075: 11 images)'. Rows with images = 0 (`absent`) are listed as missing conditions."""
    gone = desc[desc["images"] == 0] if len(desc) else desc
    desc = desc[desc["images"] != 0] if len(desc) else desc
    if desc.empty:
        return "Design: no images" + (f"; {len(gone)} condition(s) without images, see design_check.csv" if len(gone) else "")
    sp = spec(settings)

    def rng(col):
        lo, hi = int(desc[col].min()), int(desc[col].max())
        return f"{lo}" if lo == hi else f"{lo}-{hi}"

    src = "/".join(sorted(set(desc["fields_source"]) - {""})) or "fixed"
    ns, conds = desc["n"].nunique(), desc["treatment"].nunique()
    line = f"Design: {ns} n, {conds} condition{'s' * (conds != 1)}, {rng('wells_expected')} wells x {rng('fields_per_well')} fields ({src})"
    if sp.field_positions == "random":
        line += ", field positions random"
    tps = sorted(set(desc["timepoint"]) - {""})
    if tps:
        line += f", timepoint {', '.join(tps)}" + (" (some n unknown)" if (desc["timepoint"] == "").any() else "")

    def listing(d):
        bits = [f"{r.n.split('_')[0]} / {r.treatment}: {r.images} image{'s' * (r.images != 1)}" for r in d.itertuples()]
        return ", ".join(bits[:3]) + (f", +{len(bits) - 3} more" if len(bits) > 3 else "")

    unb = desc[desc["unbalanced"]]
    other = desc[~desc["unbalanced"] & (desc["flags"] != "")]
    if len(unb):
        line += f"; {len(unb)} condition{'s' * (len(unb) != 1)} unbalanced ({listing(unb)})"
    if len(other):
        line += f"; {len(other)} other condition{'s' * (len(other) != 1)} flagged ({listing(other)})"
    if len(gone):
        bits = [f"{r.n.split('_')[0]} / {r.treatment}" for r in gone.itertuples()]
        line += (f"; {len(gone)} condition{'s' * (len(gone) != 1)} missing from an n ("
                 + ", ".join(bits[:4]) + (f", +{len(bits) - 4} more" if len(bits) > 4 else "") + ")")
    if len(unb) or len(other) or len(gone):
        line += "; see design_check.csv"
    return line


def nominal(desc: pd.DataFrame) -> tuple[int, int]:
    """(wells, fields_per_well) of the largest condition design: the full design used by the sampling analysis."""
    if desc.empty:
        return DEFAULT_WELLS, 1
    return int(desc["wells_expected"].max()), int(desc["fields_per_well"].max())


# ---------------------------------------------------------------- wells.csv template

def wells_template(table: pd.DataFrame) -> pd.DataFrame:
    """A wells.csv for an `assign`ed table that reproduces its current assignment: one generic row (n, treatment
    blank) per image number with the most common well/position, plus n + treatment rows only where a condition
    differs. Edit it, then set [design] well_layout = "file"."""
    rows = []
    t = table[KEY + ["well", "position"]].astype({"image": int})
    for img, g in t.groupby("image"):
        common = g.groupby(["well", "position"]).size().sort_values(ascending=False).index[0]
        rows.append(dict(n="", treatment="", image=img, well=int(common[0]), position=int(common[1])))
        for r in g.itertuples():
            if (r.well, r.position) != common:
                rows.append(dict(n=r.n, treatment=r.treatment, image=img, well=int(r.well), position=int(r.position)))
    return pd.DataFrame(rows, columns=["n", "treatment", "image", "well", "position"])


def write_wells_template(project_dir: Path, table: pd.DataFrame, settings: dict, overwrite: bool = False) -> Path:
    """Write wells.csv (or wells.new.csv if wells.csv exists and overwrite is False) for the images in `table`
    (n, treatment, image), assigned with the current settings. Returns the path written."""
    project_dir = Path(project_dir)
    s = {**settings, "design": {**settings.get("design", {}), "well_layout": _no_file(settings)}}
    s["design"].pop("plate_layout", None)
    tmpl = wells_template(assign(table, s))
    path = project_dir / WELLS_FILE
    if path.exists() and not overwrite:
        path = project_dir / "wells.new.csv"
    elif path.exists():
        from .settings import backup
        backup(project_dir, WELLS_FILE)  # overwrite: keep the hand-edited version in _work/backups/
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    tmpl.to_csv(path, index=False, sep=delim)
    return path


def _no_file(settings: dict) -> str:
    """The current layout, or consecutive if it is "file" (a template must not depend on the file it replaces)."""
    sp = spec(settings)
    return "consecutive" if sp.well_layout == "file" else sp.well_layout
