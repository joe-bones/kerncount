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
        "names": [],  # e.g. ["DAPI", "pEIF2A"]; empty = C1, C2, ...
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
        "reference": "LyoVec",  # group every treatment is tested against
        "normalise_to": ["LyoVec", "Control"],
        "treatment_order": ["Control", "LyoVec", "CL075"],  # remaining treatments follow alphabetically
        "exclude_low_qc": False,  # g1 is fully automatic; set true to drop QC=low images from stats
    },
    "design": {
        "plates": 2,  # plates (wells) imaged per condition per n
        "fields_per_plate": 6,
        "plate_layout": "consecutive",  # "consecutive": images 1-6 = plate 1, 7-12 = plate 2; or "alternating"
    },
    "review": {
        "overlay_alpha": 0.35,
        "jpeg_quality": 88,
    },
}

TEMPLATE = '''# cellcount settings — edit freely; anything left out falls back to defaults.

[project]
raw_dir = "{raw_dir}"
date_order = "DMY"        # n-folder dates: 05.06.26 = 5 June 2026
csv_delimiter = ","       # set to ";" if Excel (German/European locale) shows everything in one column

[channels]
names = []                # e.g. ["DAPI", "pEIF2A"]; empty = C1, C2, ...
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
reference = "LyoVec"
normalise_to = ["LyoVec", "Control"]
treatment_order = ["Control", "LyoVec", "CL075"]
exclude_low_qc = false

[design]                  # used by `cellcount sampling` (how many images are needed?)
plates = 2
fields_per_plate = 6
plate_layout = "consecutive"   # images 1-6 = plate 1, 7-12 = plate 2; or "alternating" (odd/even)

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
    raw = str(Path(raw_dir).resolve()).replace("\\", "/")
    (Path(project_dir) / "settings.toml").write_text(TEMPLATE.format(raw_dir=raw), encoding="utf-8")


def section_hash(settings: dict, *sections: str) -> str:
    blob = json.dumps({s: settings.get(s) for s in sections}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]
