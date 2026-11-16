"""Settings as a form: field descriptions for the web interface, and a comment-preserving settings.toml writer.

Standard library only (the web interface must start without the scientific stack).
"""
from __future__ import annotations

import json
import re
import tomllib

from .settings import DEFAULTS

# kind: bool | int | float | choice | multi | list | num_or (a number, or one of `special`) | measure | text
# reseg = changing it re-detects every nucleus (the slow step, ~2-3 s per image)
FIELDS: dict[str, dict[str, dict]] = {
    "project": {
        "raw_dir": dict(kind="text", label="Raw images folder",
                        help="Where the raw images are. Relative paths are relative to the project folder."),
        "date_order": dict(kind="choice", options=["DMY", "MDY", "YMD"], label="Date order in n-folder names",
                           help="DMY: 05.06.26 = 5 June 2026."),
        "csv_delimiter": dict(kind="choice", options=[",", ";"], label="CSV delimiter",
                              help='Use ";" if Excel (German/European locale) shows CSVs in one column.'),
    },
    "channels": {
        "names": dict(kind="list", label="Channel names", help='e.g. DAPI, marker. Empty = C1, C2, ...', reseg=True),
        "segment": dict(kind="int", min=1, label="Channel used to find nuclei", help="1 = first channel.", reseg=True),
    },
    "segmentation": {
        "backend": dict(kind="choice", options=["stardist", "classical"], label="Detector", reseg=True,
                        help="stardist = deep-learning detector (recommended); classical = Fiji-style threshold + watershed."),
        "nucleus_diameter_px": dict(kind="num_or", special=["auto"], label="Nucleus diameter (px)", reseg=True,
                                    help="auto measures nucleus size per folder (adapts to magnification), or a fixed number."),
        "target_diameter_px": dict(kind="float", label="Detector working diameter (px)", reseg=True,
                                   help="Images are rescaled internally so nuclei are this size."),
        "prob_thresh": dict(kind="num_or", special=["default"], label="Detector confidence threshold", reseg=True,
                            help="Lower = more objects. default = the model's own value."),
        "pixel_size_um": dict(kind="num_or", special=[""], label="Pixel size (µm)",
                              help="If known: adds µm² areas and cells/mm². Leave empty if unknown. Does not re-detect."),
    },
    "counting": {
        "edge_rule": dict(kind="choice", options=["half", "exclude", "include"], label="Edge nuclei",
                          help="half = unbiased frame (top/left-touching count, bottom/right don't); exclude = drop all; include = count all."),
        "min_area_frac": dict(kind="float", label="Debris cut-off (× typical area)",
                              help="Smaller objects are debris / micronuclei and are not counted."),
        "exclude_out_of_focus": dict(kind="bool", label="Drop out-of-focus nuclei",
                                     help="Nuclei outside the focal plane (floating or dead cells) are not counted."),
        "focus_frac": dict(kind="float", label="Focus: edge sharpness (× typical)",
                           help="Out of focus = sharpness below this ..."),
        "texture_frac": dict(kind="float", label="Focus: chromatin texture (× typical)",
                             help="... AND texture below this. Both must be low."),
        "merge_area_frac": dict(kind="float", label="Possible merge (× typical area)",
                                help="Larger and non-convex = flagged as a possible merge (counted as 1)."),
    },
    "qc": {
        "high_max_uncertainty": dict(kind="float", label="High: max plausible range", help="Fraction of the count, e.g. 0.10 = 10%."),
        "high_max_disagreement": dict(kind="float", label="High: max disagreement with 2nd method", help="e.g. 0.15 = 15%."),
        "max_saturated_frac": dict(kind="float", label="High: max saturated nuclear pixels", help="e.g. 0.02 = 2%."),
        "low_min_uncertainty": dict(kind="float", label="Low: plausible range at least", help="e.g. 0.25 = 25%."),
        "low_min_disagreement": dict(kind="float", label="Low: disagreement at least", help="e.g. 0.35 = 35%."),
    },
    "intensity": {
        "enabled": dict(kind="bool", label="Measure intensities", reseg=True, help="Per-cell mean intensity (replaces the Mean-intensity macros)."),
        "measure": dict(kind="measure", label="Channels to measure", reseg=True, help='"all", or channel numbers, e.g. 2 or 2, 3.'),
        "regions": dict(kind="multi", options=["nucleus", "ring", "cell"], label="Regions", reseg=True,
                        help="ring = band around the nucleus (cytoplasm); cell = nucleus + ring (Fiji's 'enlarge')."),
        "ring_px": dict(kind="num_or", special=["auto"], label="Ring width (px)", reseg=True, help="auto = 0.25 × nucleus diameter."),
        "subtract_background": dict(kind="bool", label="Subtract background", reseg=True, help="Median of empty areas, per image."),
    },
    "analysis": {
        "reference": dict(kind="text", label="Main comparison group", help="Used for the headline plots."),
        "normalise_to": dict(kind="list", label="Control groups",
                             help="Each gets its own % plots and tests, e.g. Control, LyoVec."),
        "two_way_anova": dict(kind="bool", label="Two-way ANOVA (extra)",
                              help="Per n: treatment × well; all n: treatment × n; Tukey all pairs."),
        "treatment_order": dict(kind="list", label="Treatment order", help="Order in tables and plots; others follow alphabetically."),
        "exclude_low_qc": dict(kind="bool", label="Leave out QC-low images",
                               help="Images with QC low (and no manual value) are left out of the statistics."),
    },
    "design": {
        "wells": dict(kind="int", min=1, label="Wells per condition per n", help="Used by the sampling analysis."),
        "fields_per_well": dict(kind="int", min=1, label="Fields (images) per well"),
        "well_layout": dict(kind="choice", options=["consecutive", "alternating"], label="Well layout",
                            help="consecutive = images 1-6 well 1, 7-12 well 2; alternating = odd/even."),
    },
    "review": {
        "overlay_alpha": dict(kind="float", label="Overlay transparency", help="0 = invisible, 1 = opaque."),
        "jpeg_quality": dict(kind="int", min=1, max=100, label="Review image JPEG quality"),
    },
}

SECTION_TITLES = {
    "project": "Project", "channels": "Channels", "segmentation": "Finding nuclei", "counting": "What counts as a nucleus",
    "qc": "Quality levels", "intensity": "Intensity per cell", "analysis": "Statistics",
    "design": "Plate design (sampling analysis)", "review": "Review images",
}


def schema() -> list[dict]:
    return [dict(section=s, title=SECTION_TITLES.get(s, s),
                 fields=[dict(key=k, default=DEFAULTS[s][k], **meta) for k, meta in keys.items()])
            for s, keys in FIELDS.items()]


def _num(v, as_int=False):
    if isinstance(v, bool):
        raise ValueError
    f = float(str(v).strip().replace(",", "."))
    if as_int:
        if f != int(f):
            raise ValueError
        return int(f)
    return int(f) if f == int(f) and "." not in str(v) else f


def coerce(values: dict) -> dict:
    """Turn form values (strings, booleans, lists) into typed settings; raise ValueError naming the bad field."""
    out: dict[str, dict] = {}
    for sec, keys in values.items():
        if sec not in FIELDS:
            continue
        for key, v in keys.items():
            meta = FIELDS[sec].get(key)
            if meta is None:
                continue
            k, where = meta["kind"], f"{SECTION_TITLES[sec]} → {meta.get('label', key)}"
            try:
                if k == "bool":
                    val = bool(v)
                elif k in ("int", "float"):
                    val = _num(v, as_int=k == "int")
                    if "min" in meta and val < meta["min"] or "max" in meta and val > meta["max"]:
                        raise ValueError
                elif k == "choice":
                    if v not in meta["options"]:
                        raise ValueError
                    val = v
                elif k == "multi":
                    val = [o for o in meta["options"] if o in v]
                    if not val:
                        raise ValueError
                elif k == "list":
                    items = v if isinstance(v, list) else str(v).split(",")
                    val = [str(i).strip() for i in items if str(i).strip()]
                elif k == "num_or":
                    s = str(v).strip()
                    val = s if s in meta["special"] else _num(s)
                elif k == "measure":
                    s = str(v).strip().lower()
                    val = "all" if s in ("", "all") else [_num(c, as_int=True) for c in s.replace(";", ",").split(",") if c.strip()]
                else:
                    val = str(v).strip()
            except (ValueError, TypeError):
                raise ValueError(f"{where}: '{v}' is not a valid value") from None
            out.setdefault(sec, {})[key] = val
    return out


def toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(i) for i in v) + "]"
    return json.dumps(str(v), ensure_ascii=False)  # a JSON string is a valid TOML basic string


_SECTION = re.compile(r"^\s*\[([A-Za-z0-9_\-]+)\]\s*(#.*)?$")
_KEY = re.compile(r"^(\s*)([A-Za-z0-9_\-]+)(\s*=\s*)(.*)$")


def _split_comment(rest: str) -> tuple[str, str]:
    """'"a#b"   # note' -> ('"a#b"', '   # note'): find a # outside quotes."""
    quote = None
    i = 0
    while i < len(rest):
        c = rest[i]
        if quote:
            if c == "\\" and quote == '"':
                i += 1
            elif c == quote:
                quote = None
        elif c in "\"'":
            quote = c
        elif c == "#":
            value = rest[:i].rstrip()
            return value, rest[len(value):]
        i += 1
    return rest.rstrip(), ""


def update_text(text: str, typed: dict) -> str:
    """Write `typed` values into settings.toml text, keeping every comment and every key not being changed."""
    lines = text.splitlines()
    todo = {s: dict(kv) for s, kv in typed.items()}
    section, last_line_of = None, {}
    for i, line in enumerate(lines):
        m = _SECTION.match(line)
        if m:
            section = m.group(1)
            last_line_of[section] = i
            continue
        m = _KEY.match(line)
        if m and section:
            last_line_of[section] = i
            key = m.group(2)
            if key in todo.get(section, {}):
                value, comment = _split_comment(m.group(4))
                new = toml_value(todo[section].pop(key))
                if comment and len(value) > len(new):  # keep comments roughly aligned
                    new = new.ljust(len(value))
                lines[i] = f"{m.group(1)}{key}{m.group(3)}{new}{comment}"
    # keys or sections not in the file yet: append (insert from the bottom so line numbers stay valid)
    for sec in sorted((s for s in todo if todo[s] and s in last_line_of), key=lambda s: -last_line_of[s]):
        at = last_line_of[sec] + 1
        lines[at:at] = [f"{k} = {toml_value(v)}" for k, v in todo[sec].items()]
    for sec, kv in todo.items():
        if kv and sec not in last_line_of:
            lines += ["", f"[{sec}]"] + [f"{k} = {toml_value(v)}" for k, v in kv.items()]
    out = "\n".join(lines) + "\n"
    parsed = tomllib.loads(out)  # never write a file that cannot be read back
    for sec, kv in typed.items():
        for k, v in kv.items():
            if parsed.get(sec, {}).get(k) != v:
                raise ValueError(f"could not write {sec}.{k} (please edit settings.toml directly)")
    return out
