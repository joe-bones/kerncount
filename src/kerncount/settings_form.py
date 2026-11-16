"""Settings as a form: field descriptions for the web interface, and a comment-preserving settings.toml writer.

Standard library only (the web interface must start without the scientific stack).
"""
from __future__ import annotations

import json
import re
import tomllib

from .settings import DEFAULTS

# short names of classify.READOUTS (without "readout_"); kept here because this module must not import numpy/pandas
READOUT_NAMES = ["frac_out_of_focus", "micronuclei_per_nucleus", "nuclear_area_median", "nuclear_area_median_um2",
                 "nuclear_area_cv", "frac_irregular", "coverage", "nn_distance_median"]

# kind: bool | int | float | choice | multi | list | num_or (a number, or one of `special`) | measure | names ("all" or some of `options`) | text
# | table (TOML inline table of name = [strings], e.g. { mirna = ["miR-*"] }; shown and edited as text)
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
    "experiment": {
        "cell_line": dict(kind="text", label="Cell line",
                          help="e.g. LN229. Experiments (n's) of different cell lines are never pooled. Override per n "
                               "or condition with a cell_line column in layout.csv."),
        "magnification": dict(kind="num_or", special=[""], label="Objective magnification",
                              help="e.g. 20 for a 20x objective. With the camera pixel size it gives the pixel size."),
        "camera_pixel_um": dict(kind="num_or", special=[""], label="Camera pixel size (µm)",
                                help="Physical size of one camera pixel, e.g. 6.5. Pixel size in the image = this ÷ "
                                     "magnification (used when 'Pixel size (µm)' under Finding nuclei is empty)."),
        "seeding_density": dict(kind="num_or", special=[""], label="Seeding density",
                                help="Cells seeded, e.g. 30000. Override per n or condition in layout.csv."),
        "seeding_unit": dict(kind="text", label="Seeding unit", help='e.g. "cells/well" or "cells/cm²".'),
        "stain": dict(kind="text", label="Stain / antibodies",
                      help="e.g. Ki-67 (rabbit) + anti-rabbit AF488, DAPI. For the record only."),
        "notes": dict(kind="text", label="Notes", help="Anything worth knowing about this experiment."),
        "negative_control": dict(kind="text", label="Negative-control condition",
                                 help="Name of the condition with negative-control wells (stained without the specific "
                                      "signal). If present in an n, its nuclei set that n's positivity cut-off and it is "
                                      "left out of the statistics. Empty = no control wells."),
        "negative_control_type": dict(kind="choice", options=["secondary_only", "isotype", "unstained"],
                                      label="Negative-control type",
                                      help="secondary_only = secondary antibody without the primary (indirect "
                                           "staining); isotype = non-binding antibody of the same isotype; unstained "
                                           "= no antibody (only autofluorescence)."),
    },
    "channels": {
        "names": dict(kind="list", label="Channel names",
                      help='e.g. DAPI, marker. Empty = C1, C2, ... (separate files: empty = the file order below).',
                      reseg=True),
        "segment": dict(kind="text", label="Channel used to find nuclei",
                        help='A number (1 = first channel) or a channel name, e.g. DAPI.', reseg=True),
        "files": dict(kind="choice", options=["single", "interleaved"], label="How channels are saved", reseg=True,
                      help="single = all channels inside one file per image; interleaved = one file per channel, "
                           "alternating in file order (e.g. 1 = Ki67, 2 = DAPI, 3 = Ki67, 4 = DAPI, ...). Gaps in the "
                           "numbering are ignored. Run init again (--force-layout) after changing it."),
        "file_order": dict(kind="list", label="Channel order of the files",
                           help="Interleaved only: channel names in file order, e.g. Ki67, DAPI.", reseg=True),
        "align": dict(kind="choice", options=["check", "correct", "off"], label="Channel alignment",
                      help="Interleaved only: check = measure the shift between the files of a field and flag big "
                           "ones; correct = also shift the marker back (whole pixels, up to the limit below); off."),
        "align_warn_px": dict(kind="float", min=0, label="Alignment warning above (px)",
                              help="Fields whose channels are shifted by more than this are flagged."),
        "align_max_correct_px": dict(kind="float", min=0, label="Correct shifts up to (px)",
                                     help="Bigger shifts are not corrected (they suggest mis-paired files) but flagged."),
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
        "regions": dict(kind="multi", options=["nucleus", "ring", "cell", "fiji_cell"], label="Regions", reseg=True,
                        help="ring = band around the nucleus (cytoplasm); cell = nucleus + ring; fiji_cell = nucleus "
                             "enlarged by the Fiji ring width below (the old macro's 'Enlarge 6 px'). Changing only "
                             "fiji_cell does not re-detect."),
        "ring_px": dict(kind="num_or", special=["auto"], label="Ring width (px)", reseg=True, help="auto = 0.25 × nucleus diameter."),
        "subtract_background": dict(kind="bool", label="Subtract background", reseg=True, help="Median of empty areas, per image."),
        "per_nucleus": dict(kind="text", label="Marker channels measured per nucleus",
                            help='auto = every channel except the nuclear one; or names, e.g. Ki67; off = none. '
                                 "Per nucleus: mean, total, 90th percentile, local background. Does not re-detect."),
        "erode_px": dict(kind="num_or", special=["auto"], label="Shrink nucleus outline by (px)",
                         help="Measured inside the outline shrunk by this much, so the blurry rim and neighbours' "
                              "light are left out. auto = 1 px, 2 px for nuclei at least 16 px across."),
        "bg_gap_px": dict(kind="int", min=0, label="Background: gap around nuclei (px)",
                          help="The local background ring starts this far outside the nucleus and keeps this far "
                               "from every other nucleus."),
        "bg_ring_px": dict(kind="num_or", special=["auto"], label="Background ring width (px)",
                           help="auto = half a nucleus diameter. The median of the ring is the nucleus' background."),
        "bg_min_px": dict(kind="int", min=1, label="Background: minimum ring pixels",
                          help="Crowded nuclei with fewer ring pixels get a wider ring, then the image background."),
        "fiji_ring_px": dict(kind="int", min=1, label="Fiji ring width (px)",
                             help="Region fiji_cell = nucleus enlarged by this (the old macro used 6 px)."),
        "flatfield": dict(kind="text", label="Flat-field correction",
                          help="Uneven illumination (e.g. darker edges). auto = estimated per n and channel from all "
                               "of that n's images; off = none; or a folder of reference images (uniform slide; "
                               "files named with the channel; 'dark' in the name = dark frames). Intensities only, "
                               "never counts."),
        "flatfield_dark": dict(kind="bool", label="Flat field: remove background pattern first",
                               help="Estimate and remove the additive background (camera offset, autofluorescence) "
                                    "pattern before estimating the illumination."),
    },
    "readouts": {
        "enabled": dict(kind="bool", label="Extra readouts",
                        help="Per-image readouts besides the count (out-of-focus fraction, micronuclei, nuclear size and "
                             "shape, crowding), computed from the detected nuclei. Fast: does not re-detect."),
        "include": dict(kind="names", options=READOUT_NAMES, label="Readouts to include",
                        help='"all", or names separated by commas: ' + ", ".join(READOUT_NAMES) + "."),
    },
    "positivity": {
        "enabled": dict(kind="bool", label="% positive nuclei",
                        help="Classify each nucleus as positive or negative for the marker; per image the % positive "
                             "(primary readout) and the same with the 90th percentile."),
        "channels": dict(kind="text", label="Channels", help="auto = every per-nucleus channel; or names, e.g. Ki67."),
        "control_percentile": dict(kind="float", min=50, max=100, label="Negative control: percentile",
                                   help="With negative-control wells in an n: cut-off = this percentile of their nuclei "
                                        "(99 = 1% of control nuclei would be called positive)."),
        "manual_cutoff": dict(kind="text", label="Manual cut-off",
                              help="A fixed cut-off (corrected intensity units), used when there are no control wells. "
                                   "Empty = the value locked in the tuning view (cutoffs.json), else the automatic "
                                   "rule. Fix it before looking at the statistics."),
        "lod_k": dict(kind="float", min=0.5, max=20, label="Automatic rule: k x noise",
                      help="Limit of detection: a nucleus is positive when its signal exceeds k x the background "
                           "noise expected for its size. 3 = about 1 false positive in 700 unstained nuclei."),
        "min_nuclei": dict(kind="int", min=1, label="Minimum nuclei per image",
                           help="Images with fewer measurable nuclei get no % positive."),
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
        "primary_test": dict(kind="choice", options=["mixed", "rm_anova"], label="Test behind the stars",
                             help="mixed = hierarchical model on the log scale (fields in wells in n; each treatment "
                                  "class gets its own n-to-n variability; results as ratios). rm_anova = Prism-style "
                                  "repeated-measures ANOVA + Dunnett on the per-n means. Both are always computed; this "
                                  "only decides which one draws the stars. If the model cannot be fitted, rm_anova is "
                                  "used and the plot says so."),
        "primary_family": dict(kind="text", label="Main claim: compare with",
                               help="The control whose comparisons are the main result (headline stars), e.g. LyoVec. "
                                    "Empty = the main comparison group above."),
        "vehicle_check": dict(kind="text", label="Vehicle check (equivalence test)",
                              help='e.g. "LyoVec vs Control": tests whether the vehicle has NO meaningful effect '
                                   "(equivalence, not just 'not significant'). auto = every other control group vs the "
                                   "main comparison group; empty = off. Several: separate with ;"),
        "equivalence_margin": dict(kind="float", min=0.01, max=1, label="Equivalence margin",
                                   help="0.2 = the vehicle counts as equivalent if its ratio to the control (90% CI) "
                                        "lies entirely within 1/1.2 to 1.2 (about ±20%)."),
        "positive_control": dict(kind="text", label="Positive control",
                                 help="e.g. CL075. Every treatment that is not a control is also compared with it. "
                                      "Empty = none."),
        "classes": dict(kind="table", label="Treatment classes",
                        help='Groups of treatments whose effects may vary differently between n\'s, each with its own '
                             'n-to-n variance in the model. TOML inline table of name = [patterns], * = wildcard, e.g. '
                             '{ mirna = ["miR-*"] }. Empty {} = controls form one class, everything else another.'),
        "among_classes": dict(kind="list", label="Compare within classes",
                              help="Class names whose members are compared with each other (all pairs), e.g. mirna."),
    },
    "acceptance": {
        "vehicle_within": dict(kind="num_or", special=[""], label="Vehicle within ± (ratio)",
                               help="Per n, the vehicle / control ratio of nuclei per field must lie within 1 ± this "
                                    "(0.3 = 0.7 to 1.3). Uses the vehicle check pair. Empty = off."),
        "positive_control_max_ratio": dict(kind="num_or", special=[""], label="Positive control must reach (ratio)",
                                           help="Per n, positive control / main comparison group must be at most this "
                                                "ratio (e.g. 0.8 = at least 20% fewer nuclei); a value above 1 means "
                                                "at least this ratio. Needs a positive control. Empty = off."),
        "min_images_per_condition": dict(kind="int", min=0, label="Minimum images per condition",
                                         help="Every condition in that n needs at least this many usable images."),
        "max_low_qc_fraction": dict(kind="float", min=0, max=1, label="Maximum share of QC-low images",
                                    help="e.g. 0.5 = at most half of that n's images may have QC low."),
        "enforce": dict(kind="bool", label="Leave out n's that fail",
                        help="Off = only report (acceptance.csv). On = n's that fail any criterion are left out of "
                             "every all-n statistic; they are listed in the About sheet and the manifest."),
    },
    "design": {
        "wells": dict(kind="num_or", special=["auto"], label="Wells per condition per n",
                      help="Usually 2. auto = take the wells from wells.csv in the project folder if it exists, else 2."),
        "fields_per_well": dict(kind="num_or", special=["auto"], label="Fields (images) per well",
                                help="auto = per condition, highest image number ÷ wells, rounded up (11 images = 6 + 5, "
                                     "18 images = 9 + 9). A number fixes it; images beyond wells × fields are flagged."),
        "well_layout": dict(kind="choice", options=["consecutive", "alternating", "file"], label="Well layout",
                            help="consecutive = images 1-6 well 1, 7-12 well 2; alternating = odd/even; "
                                 "file = wells.csv in the project folder (columns image, well, position; optional n, treatment)."),
        "field_positions": dict(kind="choice", options=["fixed", "random"], label="Field positions",
                                help="fixed = field 1 of every well is the same spot (e.g. centre), so position can "
                                     "explain differences between fields; random = fields were taken anywhere, "
                                     "position is ignored."),
    },
    "review": {
        "overlay_alpha": dict(kind="float", label="Overlay transparency", help="0 = invisible, 1 = opaque."),
        "jpeg_quality": dict(kind="int", min=1, max=100, label="Review image JPEG quality"),
    },
    "planner": {
        "effects": dict(kind="list", label="Effects to plan for (ratios)",
                        help="Treatment ÷ reference ratios, e.g. 0.8, 0.7, 0.5 (0.7 = 30% fewer nuclei)."),
        "target_power": dict(kind="float", min=0.05, max=0.999, label="Target power",
                             help="Chance of finding a real effect (p < alpha); 0.8 is usual."),
        "alpha": dict(kind="float", min=0.0001, max=0.5, label="Significance level (alpha)", help="Usually 0.05."),
        "cost_per_n": dict(kind="float", min=0, label="Effort per n (plate / day)",
                           help="Seeding, transfection, staining for one independent experiment. Any unit, e.g. minutes."),
        "cost_per_well": dict(kind="float", min=0, label="Effort per well", help="Same unit; pipetting and handling one more well."),
        "cost_per_field": dict(kind="float", min=0, label="Effort per field (image)",
                               help="Same unit; taking and checking one more image."),
        "max_n": dict(kind="int", min=3, max=100, label="Largest number of n's in the tables"),
        "plate_format": dict(kind="int", min=6, max=96, label="Plate format (wells)", help="6, 12, 24, 48 or 96."),
        "avoid_edge": dict(kind="bool", label="Plate maps: avoid edge wells",
                           help="Use inner wells first; edge wells (evaporation) only when the inner ones are not enough."),
        "reference_wells": dict(kind="num_or", special=["auto"], label="Wells per reference condition",
                                help="auto = Dunnett's rule: wells × √(number of comparisons), as far as the plate allows."),
    },
}

SECTION_TITLES = {
    "project": "Project", "experiment": "Experiment (cell line, imaging, controls)", "channels": "Channels", "segmentation": "Finding nuclei", "counting": "What counts as a nucleus",
    "qc": "Quality levels", "intensity": "Intensity per cell / nucleus",
    "positivity": "% positive nuclei (e.g. Ki-67)", "readouts": "Extra readouts (beyond the count)",
    "analysis": "Statistics", "acceptance": "Acceptance criteria per n",
    "design": "Plate design (wells, fields, positions)", "review": "Review images",
    "planner": "Planning the next experiment",
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
                elif k == "table":
                    val = parse_table(v)
                elif k == "measure":
                    s = str(v).strip().lower()
                    val = "all" if s in ("", "all") else [_num(c, as_int=True) for c in s.replace(";", ",").split(",") if c.strip()]
                elif k == "names":  # "all", or a subset of `options` (the readout_ prefix is optional)
                    items = v if isinstance(v, list) else str(v).replace(";", ",").split(",")
                    items = [str(i).strip().removeprefix("readout_") for i in items if str(i).strip()]
                    if not items or "all" in items:
                        val = "all"
                    else:
                        if any(i not in meta["options"] for i in items):
                            raise ValueError
                        val = [o for o in meta["options"] if o in items]
                else:
                    val = str(v).strip()
            except (ValueError, TypeError):
                raise ValueError(f"{where}: '{v}' is not a valid value") from None
            out.setdefault(sec, {})[key] = val
    return out


def parse_table(v) -> dict:
    """kind "table": a TOML inline table of name = [strings] (dict or its text form, e.g. '{ mirna = ["miR-*"] }')."""
    if isinstance(v, dict):
        t = v
    else:
        s = str(v).strip()
        if not s:
            return {}
        if not s.startswith("{"):
            s = "{ " + s + " }"
        t = tomllib.loads(f"x = {s}")["x"]
    if not isinstance(t, dict):
        raise ValueError
    out = {}
    for name, pats in t.items():
        if isinstance(pats, str):
            pats = [pats]
        if not isinstance(pats, list) or not all(isinstance(p, str) for p in pats):
            raise ValueError
        out[str(name)] = [p.strip() for p in pats if p.strip()]
    return out


def table_text(v: dict) -> str:
    """Display form of a kind "table" value."""
    return toml_value(v or {})


def toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(i) for i in v) + "]"
    if isinstance(v, dict):
        if not v:
            return "{}"
        keys = [k if re.fullmatch(r"[A-Za-z0-9_\-]+", k) else json.dumps(k) for k in v]
        return "{ " + ", ".join(f"{k} = {toml_value(x)}" for k, x in zip(keys, v.values())) + " }"
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
