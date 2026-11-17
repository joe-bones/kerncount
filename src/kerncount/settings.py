"""Project settings: a commented settings.toml merged over built-in defaults."""
from __future__ import annotations

import copy
import hashlib
import json
import tomllib
from pathlib import Path

DEFAULTS: dict = {
    "project": {
        "raw_dir": "",
        "date_order": "DMY",  # how dates in n-folder names are written (05.06.26 = 5 June 2026)
        "csv_delimiter": ",",  # use ";" if your Excel is set to German/European locale
    },
    "experiment": {
        # descriptive metadata (docs/intensity.md); per-n / per-condition overrides: same-named columns in layout.csv
        "cell_line": "",  # e.g. "LN229"; n's of different cell lines are never pooled
        "magnification": "",  # objective, e.g. 20 (for 20x)
        "camera_pixel_um": "",  # camera pixel size in µm; pixel size = camera_pixel_um / magnification
        "seeding_density": "",  # e.g. 30000
        "seeding_unit": "cells/well",
        "stain": "",  # e.g. "Ki-67 (rabbit) + anti-rabbit AF488; DAPI"
        "notes": "",
        "negative_control": "",  # treatment name of the negative-control wells ("" = none on the plate)
        "negative_control_type": "secondary_only",  # secondary_only | isotype | unstained
    },
    "channels": {
        "names": [],  # e.g. ["DAPI", "marker"]; empty = C1, C2, ... (interleaved: empty = file_order)
        "segment": 1,  # channel used to find nuclei: 1-based number or a channel name (e.g. "DAPI")
        "files": "single",  # "single": all channels inside one file; "interleaved": one file per channel, alternating
        "file_order": [],  # interleaved: channel names in file order, e.g. ["Ki67", "DAPI"]
        "align": "check",  # interleaved: "check" (report shifts), "correct" (also shift small offsets back), "off"
        "align_warn_px": 2.0,  # warn when the channels of a field are shifted by more than this
        "align_max_correct_px": 10.0,  # "correct" only fixes shifts up to this (bigger = suspicious, flagged instead)
    },
    "segmentation": {
        "backend": "stardist",  # "stardist" (deep learning, recommended) or "classical"
        "nucleus_diameter_px": "auto",  # "auto" adapts to magnification; or a number
        "target_diameter_px": 30,  # size StarDist works best at; images are rescaled to this
        "prob_thresh": "default",
        "pixel_size_um": "",  # optional; if known, areas/densities are also given in µm
    },
    "counting": {
        "edge_rule": "half",  # "half" = unbiased frame (keep top/left-touching), "exclude", "include"
        "min_area_frac": 0.2,  # objects smaller than this x median nucleus area = debris/micronuclei
        "exclude_out_of_focus": True,
        "focus_frac": 0.5,  # nuclei sharper than this x typical sharpness are "in focus"
        "texture_frac": 0.8,  # ...or with chromatin texture above this x typical (rescues sharp nuclei in glare)
        "merge_area_frac": 2.0,  # objects bigger than this x median & non-convex = suspected merges
    },
    "qc": {
        "high_max_uncertainty": 0.10,
        "low_min_uncertainty": 0.25,
        "high_max_disagreement": 0.15,
        "low_min_disagreement": 0.35,
        "max_saturated_frac": 0.02,
    },
    "intensity": {
        "enabled": True,
        "measure": "all",  # "all" or list of 1-based channels
        "regions": ["nucleus", "cell"],  # nucleus, ring (cytoplasm band), cell (nucleus + ring)
        "ring_px": "auto",  # band width; "auto" = 0.25 x nucleus diameter (Fiji macro used 6 px)
        "subtract_background": True,
        # per-nucleus measurement of marker channels (docs/intensity.md); never changes detection or counts
        "per_nucleus": "auto",  # channels measured per nucleus: "auto" = every channel except the nuclear one; list; "off"
        "erode_px": "auto",  # shrink the nucleus outline before measuring: "auto" = 1 px, 2 px for nuclei >= 16 px across
        "bg_gap_px": 3,  # local background starts this far outside the nucleus and avoids this margin around every nucleus
        "bg_ring_px": "auto",  # width of the local background ring; "auto" = 0.5 x nucleus diameter
        "bg_min_px": 40,  # fewer ring pixels than this: wider ring, then the image-level background
        "fiji_ring_px": 6,  # legacy region "fiji_cell" = nucleus enlarged by this many px (the Fiji macro)
        "flatfield": "auto",  # "auto" (estimated per n and channel from that n's images), "off", or a folder of reference images
        "flatfield_dark": True,  # remove the additive background pattern (camera offset, autofluorescence) before estimating
    },
    "readouts": {
        "enabled": True,  # extra per-image readouts from the classified nuclei (docs/readouts.md); no re-detection
        "include": "all",  # "all" or a list of names, e.g. ["frac_out_of_focus", "coverage"]
    },
    "positivity": {
        "enabled": True,  # % positive nuclei per marker channel (docs/intensity.md)
        "channels": "auto",  # "auto" = every per-nucleus channel, or a list of names
        "control_percentile": 99.0,  # with negative-control wells: cut-off = this percentile of their nuclei
        "manual_cutoff": "",  # a number (all channels) or { Ki67 = 1200 }; "" = cutoffs.json lock, else the LoD rule
        "lod_k": 3.0,  # limit-of-detection rule: positive when signal > k x background noise SD for the nucleus area
        "min_nuclei": 5,  # fewer measurable nuclei in an image: the image's % positive is left empty
    },
    "analysis": {
        "reference": "Control",  # main comparison group (used for the headline plots)
        "normalise_to": ["Control", "LyoVec"],  # every one of these also gets its own % plots and tests
        "two_way_anova": True,  # extra two-way ANOVA + all-pairs Tukey (per n: treatment x well; all n: treatment x n)
        "treatment_order": ["Control", "LyoVec", "CL075"],  # remaining treatments follow alphabetically
        "exclude_low_qc": False,  # g1 is fully automatic; set true to drop QC=low images from stats
        # all-n statistics (docs/statistics.md)
        "primary_test": "mixed",  # drives the stars: "mixed" (hierarchical model, log scale) or "rm_anova" (Prism standard)
        "primary_family": "",  # control whose comparisons are the main claim; "" = `reference`
        "vehicle_check": "auto",  # equivalence test, e.g. "LyoVec vs Control"; "auto" = other controls vs reference; "" = off
        "equivalence_margin": 0.2,  # vehicle check: ratio must lie within [1/1.2, 1.2] to call it equivalent
        "positive_control": "",  # e.g. "CL075": every non-control treatment is also compared with it
        "model_unit": "well",  # all-n model unit: "well" (mean per field of each well; ratios of means) or "field"
        "classes": "auto",  # own n x treatment variance per class: "auto" = one per treatment, each control
                            # included; {} = controls vs everything else; or e.g. { mirna = ["miR-*"] }
        "density_covariate": "readout_coverage",  # fraction metrics: also a density-adjusted model; "" = off
        "among_classes": [],  # classes whose members are compared with each other (all pairs), e.g. ["mirna"]
    },
    "acceptance": {
        # per-n acceptance criteria (count metric); report only unless enforce = true
        "vehicle_within": 0.3,  # vehicle / control ratio in that n must lie within 1 +- this ("" = off)
        "positive_control_max_ratio": "",  # positive control / reference must reach this ratio (e.g. 0.8) ("" = off)
        "min_images_per_condition": 3,  # every condition in that n needs at least this many usable images
        "max_low_qc_fraction": 0.5,  # at most this fraction of that n's images may be QC low
        "no_duplicate_images": True,  # an n fails if one of its raw images is byte-identical to an image of another condition
        "enforce": False,  # true: n's that fail are left out of all-n statistics (listed in About and manifest)
    },
    "design": {
        # one plate per n; each condition sits in `wells` wells, each imaged `fields_per_well` times (docs/design.md)
        "wells": 2,  # wells imaged per condition per n; "auto" = from wells.csv if present, else 2
        "fields_per_well": "auto",  # "auto" = per condition ceil(highest image number / wells); or a fixed number
        "well_layout": "consecutive",  # "consecutive" (1-6 = well 1, 7-12 = well 2), "alternating" (odd/even), "file"
        "field_positions": "fixed",  # "fixed": field k is the same spot in every well; "random": position means nothing
    },
    "review": {
        "overlay_alpha": 0.35,
        "jpeg_quality": 88,
    },
    "planner": {
        # experiment planner (`kerncount sampling` -> planner.json, `kerncount platemap`); see docs/planner.md
        "effects": [0.8, 0.7, 0.5],  # treatment / reference ratios to plan for (0.7 = 30% fewer nuclei)
        "target_power": 0.8,  # chance of p < alpha if the effect is real
        "alpha": 0.05,
        "cost_per_n": 240.0,  # effort per independent experiment (one plate, one day): seeding, treating, staining
        "cost_per_well": 5.0,  # effort per extra well (pipetting, handling)
        "cost_per_field": 1.0,  # effort per image (acquisition + checking); any unit, default minutes
        "max_n": 12,  # largest number of n's in the design tables
        "plate_format": 24,  # 6, 12, 24, 48 or 96 wells
        "avoid_edge": True,  # plate maps: inner wells first, edge wells only when needed
        "reference_wells": "auto",  # wells per reference condition; "auto" = Dunnett's sqrt(k) rule
        # power / MDE for a star on the headline plot: the family's multivariate-t adjusted critical value
        # (false = unadjusted paired t-test of one comparison, as the ratio t-test)
        "adjust_for_family": True,
    },
}

TEMPLATE = '''# kerncount settings — edit freely; anything left out falls back to defaults.

[project]
raw_dir = "{raw_dir}"   # raw images; relative paths are relative to this folder
date_order = "DMY"        # n-folder dates: 05.06.26 = 5 June 2026
csv_delimiter = ","       # set to ";" if Excel (German/European locale) shows everything in one column

[experiment]              # what was imaged; per-n / per-condition overrides: same-named columns in layout.csv
cell_line = ""            # e.g. "LN229"; different cell lines are never pooled
magnification = ""        # objective, e.g. 20
camera_pixel_um = ""      # camera pixel size (µm); pixel size = camera_pixel_um / magnification
seeding_density = ""      # e.g. 30000
seeding_unit = "cells/well"
stain = ""                # e.g. "Ki-67 + secondary AF488; DAPI"
notes = ""
negative_control = ""     # treatment name of negative-control wells ("" = none): sets the positivity cut-off
negative_control_type = "secondary_only"   # secondary_only | isotype | unstained

[channels]
names = []                # e.g. ["DAPI", "marker"]; empty = C1, C2, ... (interleaved: empty = file_order)
segment = 1               # channel used to find nuclei: number (1-based) or name, e.g. "DAPI"
files = "single"          # "single": channels inside one file; "interleaved": one file per channel, alternating
file_order = []           # interleaved: channel names in file order, e.g. ["Ki67", "DAPI"]
align = "check"           # interleaved: "check" shifts between the files of a field, "correct" small ones, "off"
align_warn_px = 2.0       # warn above this shift (px)
align_max_correct_px = 10.0   # "correct" fixes shifts up to this; larger ones are only flagged

[segmentation]
backend = "stardist"      # "stardist" (recommended) or "classical" (Fiji-like threshold + watershed)
nucleus_diameter_px = "auto"   # "auto" adapts to magnification/camera; or a fixed number
target_diameter_px = 30
prob_thresh = "default"
pixel_size_um = ""        # optional, e.g. 0.65 — adds µm² areas and cells/mm²

[counting]
edge_rule = "half"        # "half" (unbiased: count nuclei touching top/left edges only), "exclude", "include"
min_area_frac = 0.2       # smaller than 0.2 x typical nucleus = debris / micronucleus (not counted)
exclude_out_of_focus = true
focus_frac = 0.5          # edge sharpness below 0.5 x typical ...
texture_frac = 0.8        # ... AND chromatin texture below 0.8 x typical = out of focus (not counted)
merge_area_frac = 2.0     # larger than 2 x typical and non-convex = suspected merged nuclei (flagged)

[qc]
high_max_uncertainty = 0.10
low_min_uncertainty = 0.25
high_max_disagreement = 0.15
low_min_disagreement = 0.35
max_saturated_frac = 0.02

[intensity]
enabled = true
measure = "all"           # "all" or a list of channels, e.g. [2]
regions = ["nucleus", "cell"]   # nucleus, ring, cell (= nucleus + ring), fiji_cell (= nucleus + fiji_ring_px, the old macro)
ring_px = "auto"          # "auto" = 0.25 x nucleus diameter
subtract_background = true
per_nucleus = "auto"      # marker channels measured per nucleus: "auto" = all except the nuclear channel; "off" = none
erode_px = "auto"         # shrink outlines before measuring: 1 px, 2 px for nuclei >= 16 px across
bg_gap_px = 3             # local background ring starts this far out and avoids this margin around all nuclei
bg_ring_px = "auto"       # ring width; "auto" = 0.5 x nucleus diameter
bg_min_px = 40            # too few ring pixels: wider ring, then the image background
fiji_ring_px = 6          # region "fiji_cell" = nucleus enlarged by 6 px (the old Fiji macro)
flatfield = "auto"        # uneven illumination: "auto" (estimated per n), "off", or a folder of reference images
flatfield_dark = true     # remove the additive background pattern before estimating the flat field

[readouts]                # extra per-image readouts from the detected nuclei (see docs/readouts.md); fast, no re-detection
enabled = true
include = "all"           # "all" or a list, e.g. ["frac_out_of_focus", "micronuclei_per_nucleus", "coverage"]

[positivity]              # % positive nuclei per marker channel (docs/intensity.md)
enabled = true
channels = "auto"         # "auto" = every per-nucleus channel
control_percentile = 99.0 # negative-control wells present: cut-off = this percentile of their nuclei
manual_cutoff = ""        # a number, or {{ Ki67 = 1200 }}; "" = locked value from the tuning view, else automatic
lod_k = 3.0               # automatic: signal > 3 x background noise SD for the nucleus area
min_nuclei = 5            # fewer measurable nuclei in an image: its % positive is left empty

[analysis]
reference = "Control"               # main comparison group (headline plots)
normalise_to = ["Control", "LyoVec"] # each also gets % plots and its own tests (vs Control AND vs LyoVec)
two_way_anova = true                 # extra: two-way ANOVA + Tukey all-pairs (per n: treatment x well; all n: treatment x n)
treatment_order = ["Control", "LyoVec", "CL075"]
exclude_low_qc = false
# all-n statistics (see docs/statistics.md)
primary_test = "mixed"               # what drives the stars: "mixed" (hierarchical model on the log scale) or "rm_anova"
primary_family = ""                  # control whose comparisons are the main claim, e.g. "LyoVec"; "" = reference
vehicle_check = "auto"               # equivalence test, e.g. "LyoVec vs Control"; "auto" = other controls vs reference; "" = off
equivalence_margin = 0.2             # equivalent = ratio (90% CI) within [1/1.2, 1.2]
positive_control = ""                # e.g. "CL075": every non-control treatment is also compared with it
model_unit = "well"                  # model unit: "well" (well mean per field; ratios of mean nuclei per field) or "field"
classes = "auto"                     # own n-to-n variance per class: "auto" = one per treatment, each control
                                     # included (conservative); {{}} = controls vs the rest; or e.g. {{ mirna = ["miR-*"] }}
density_covariate = "readout_coverage"   # fraction metrics: extra density-adjusted estimate; "" = off
among_classes = []                   # compare members of these classes with each other, e.g. ["mirna"]

[acceptance]                         # per-n acceptance criteria (on the counts); report only unless enforce = true
vehicle_within = 0.3                 # vehicle / control in that n within 1 +- 0.3; "" = off
positive_control_max_ratio = ""      # e.g. 0.8: positive control / reference must reach 0.8 (>= the value if above 1); "" = off
min_images_per_condition = 3         # every condition needs at least this many usable images
max_low_qc_fraction = 0.5            # at most this share of the n's images may be QC low
no_duplicate_images = true           # fail an n whose raw image is byte-identical to one of another condition
enforce = false                      # true: failing n's are left out of all-n statistics (listed in About)

[design]                  # which well / field each image is; one plate per n (see docs/design.md)
wells = 2                 # wells per condition per n; "auto" = take them from wells.csv if it exists, else 2
fields_per_well = "auto"  # "auto": per condition, highest image number / wells (rounded up), so 11 images = 6 + 5
                          # and 18 images = 9 + 9; or a fixed number (images above wells x fields are flagged)
well_layout = "consecutive"    # "consecutive": images 1-6 = well 1, 7-12 = well 2; "alternating" (odd/even);
                               # "file": read wells.csv (columns image, well, position, optional n, treatment)
field_positions = "fixed" # "fixed": field k is the same spot in every well (position is used); "random": it is not

[review]
overlay_alpha = 0.35
jpeg_quality = 88

[planner]                 # plan the next experiment (kerncount sampling -> planner; kerncount platemap); docs/planner.md
effects = [0.8, 0.7, 0.5] # treatment / reference ratios to plan for (0.7 = 30% fewer nuclei than the reference)
target_power = 0.8        # chance of p < alpha when the effect is real
alpha = 0.05
cost_per_n = 240.0        # effort per n (one plate, one day: seeding, transfection, staining); any unit, e.g. minutes
cost_per_well = 5.0       # effort per extra well
cost_per_field = 1.0      # effort per image (acquisition + checking)
max_n = 12                # largest number of n's in the design tables
plate_format = 24         # 6, 12, 24, 48 or 96 wells
avoid_edge = true         # plate maps: inner wells first, edge wells only when needed
reference_wells = "auto"  # wells per reference condition; "auto" = Dunnett's rule (wells x sqrt(comparisons))
adjust_for_family = true  # power = chance of a star on the headline plot (family-adjusted); false = one unadjusted test
'''


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load(project_dir: Path) -> dict:
    path = Path(project_dir) / "settings.toml"
    user = tomllib.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    return _merge(DEFAULTS, user)


BACKUPS = "_work/backups"


def backup(project_dir: Path, name: str) -> Path | None:
    """Copy <project>/<name> (settings.toml, layout.csv, wells.csv, ...) to _work/backups/<stem>_<timestamp><suffix>
    before it is overwritten. Returns the copy (None when there was nothing to back up)."""
    import shutil
    import time

    src = Path(project_dir) / name
    if not src.is_file():
        return None
    d = Path(project_dir) / BACKUPS
    d.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    dst = d / f"{src.stem}_{stamp}{src.suffix}"
    k = 1
    while dst.exists():
        k += 1
        dst = d / f"{src.stem}_{stamp}_{k}{src.suffix}"
    shutil.copy2(src, dst)
    return dst


def write_template(project_dir: Path, raw_dir: Path) -> Path | None:
    """raw_dir is stored relative to the project folder when possible, so a workspace can be moved. Symlinks are
    NOT resolved: a data/<experiment> link to a folder elsewhere is stored as ../../data/<experiment>, so the link
    (not its current target) stays the source of truth. An existing settings.toml is backed up first."""
    import os

    raw_abs, prj_abs = os.path.abspath(raw_dir), os.path.abspath(project_dir)
    try:
        raw = os.path.relpath(raw_abs, prj_abs)
    except ValueError:  # different drives on Windows
        raw = str(raw_abs)
    raw = raw.replace("\\", "/")
    bk = backup(project_dir, "settings.toml")
    from .fileio import atomic_write

    atomic_write(Path(project_dir) / "settings.toml", TEMPLATE.format(raw_dir=raw))
    return bk


def section_hash(settings: dict, *sections: str) -> str:
    blob = json.dumps({s: settings.get(s) for s in sections}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


def pixel_size_um(settings: dict) -> float | None:
    """[segmentation] pixel_size_um if set, else [experiment] camera_pixel_um / magnification, else None."""
    v = settings.get("segmentation", {}).get("pixel_size_um", "")
    if v not in ("", None):
        return float(v)
    e = settings.get("experiment", {})
    try:
        cam, mag = float(e.get("camera_pixel_um") or 0), float(str(e.get("magnification") or 0).lower().rstrip("x"))
    except ValueError:
        return None
    return cam / mag if cam > 0 and mag > 0 else None


# ---------------------------------------------------------------- channels

# Keys that existed before per-nucleus intensity / interleaved files. Only these enter the segmentation hash for
# single-file projects, so adding new settings never invalidates cached detections.
_LEGACY_SEG_KEYS = {"channels": ("names", "segment"),
                    "intensity": ("enabled", "measure", "regions", "ring_px", "subtract_background")}


def interleaved(settings: dict) -> bool:
    return str(settings.get("channels", {}).get("files", "single")).strip().lower() == "interleaved"


def channel_names(settings: dict, n_channels: int) -> list[str]:
    """Names of the channels of a field stack: channels.names, else (interleaved) file_order, else C1, C2, ..."""
    ch = settings.get("channels", {})
    names = [str(x) for x in (ch.get("names") or [])]
    if not names and interleaved(settings):
        names = [str(x) for x in (ch.get("file_order") or [])]
    return [names[i] if i < len(names) and names[i] else f"C{i + 1}" for i in range(n_channels)]


def channel_index(settings: dict, ref, n_channels: int) -> int:
    """0-based index of a channel given as a 1-based number ("2", 2) or a name ("DAPI", case-insensitive)."""
    names = channel_names(settings, n_channels)
    s = str(ref).strip()
    if s.lstrip("-").isdigit():
        i = int(s) - 1
    else:
        low = [n.lower() for n in names]
        if s.lower() not in low:
            raise ValueError(f"channel {s!r} not found; channels are {', '.join(names)}")
        i = low.index(s.lower())
    if not 0 <= i < n_channels:
        raise ValueError(f"channel {s} does not exist (the images have {n_channels} channel(s): {', '.join(names)})")
    return i


def n_channels_expected(settings: dict) -> int | None:
    """Number of channels of a field when it is known from the settings alone (interleaved), else None."""
    if interleaved(settings):
        return max(len(settings["channels"].get("file_order") or []), 1)
    return None


def seg_hash_input(settings: dict) -> dict:
    """The settings that decide the detection result (and cached per-object measurements). Single-file projects
    hash exactly the keys they always had; interleaved projects add their file layout."""
    out = {"segmentation": {**settings["segmentation"], "pixel_size_um": ""}}
    for sec, keys in _LEGACY_SEG_KEYS.items():
        d = settings.get(sec, {})
        out[sec] = {k: d[k] for k in keys if k in d}
    regions = out["intensity"].get("regions")
    if isinstance(regions, list):  # the legacy fiji_cell region is measured later, outside the detection cache
        out["intensity"]["regions"] = [r for r in regions if r != "fiji_cell"]
    if interleaved(settings):
        out["channels"].update(files="interleaved", file_order=list(settings["channels"].get("file_order") or []))
    return out
