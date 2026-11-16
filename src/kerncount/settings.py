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
    "channels": {
        "names": [],  # e.g. ["DAPI", "marker"]; empty = C1, C2, ...
        "segment": 1,  # 1-based channel used to find nuclei
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
        "classes": {},  # treatment classes with their own n x treatment variance, e.g. { mirna = ["miR-*"] }
        "among_classes": [],  # classes whose members are compared with each other (all pairs), e.g. ["mirna"]
    },
    "acceptance": {
        # per-n acceptance criteria (count metric); report only unless enforce = true
        "vehicle_within": 0.3,  # vehicle / control ratio in that n must lie within 1 +- this ("" = off)
        "positive_control_max_ratio": "",  # positive control / reference must reach this ratio (e.g. 0.8) ("" = off)
        "min_images_per_condition": 3,  # every condition in that n needs at least this many usable images
        "max_low_qc_fraction": 0.5,  # at most this fraction of that n's images may be QC low
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
}

TEMPLATE = '''# kerncount settings — edit freely; anything left out falls back to defaults.

[project]
raw_dir = "{raw_dir}"   # raw images; relative paths are relative to this folder
date_order = "DMY"        # n-folder dates: 05.06.26 = 5 June 2026
csv_delimiter = ","       # set to ";" if Excel (German/European locale) shows everything in one column

[channels]
names = []                # e.g. ["DAPI", "marker"]; empty = C1, C2, ...
segment = 1               # channel (1-based) used to find nuclei

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
regions = ["nucleus", "cell"]   # nucleus, ring, cell (= nucleus + ring, like the Fiji macro)
ring_px = "auto"          # "auto" = 0.25 x nucleus diameter
subtract_background = true

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
classes = {{}}                         # own n x treatment variance per class, e.g. {{ mirna = ["miR-*"] }}
among_classes = []                   # compare members of these classes with each other, e.g. ["mirna"]

[acceptance]                         # per-n acceptance criteria (on the counts); report only unless enforce = true
vehicle_within = 0.3                 # vehicle / control in that n within 1 +- 0.3; "" = off
positive_control_max_ratio = ""      # e.g. 0.8: positive control / reference must reach 0.8 (>= the value if above 1); "" = off
min_images_per_condition = 3         # every condition needs at least this many usable images
max_low_qc_fraction = 0.5            # at most this share of the n's images may be QC low
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


def write_template(project_dir: Path, raw_dir: Path) -> None:
    """raw_dir is stored relative to the project folder when possible, so a workspace can be moved."""
    import os

    raw_abs, prj_abs = Path(raw_dir).resolve(), Path(project_dir).resolve()
    try:
        raw = os.path.relpath(raw_abs, prj_abs)
    except ValueError:  # different drives on Windows
        raw = str(raw_abs)
    raw = raw.replace("\\", "/")
    (Path(project_dir) / "settings.toml").write_text(TEMPLATE.format(raw_dir=raw), encoding="utf-8")


def section_hash(settings: dict, *sections: str) -> str:
    blob = json.dumps({s: settings.get(s) for s in sections}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]
