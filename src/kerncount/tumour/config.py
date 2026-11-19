"""[tumour] settings: defaults for the TLR7 tumour experiment and the declared factors (ADR 0004 amendment).

A factor is one crossed condition axis on every tumour. Nothing downstream assumes how many factors there are:
layout columns, tumour IDs, groups and plot facets are all built from this list.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

ROLES = ("row", "panel", "x", "colour")
SOURCES = ("filename", "session")

DEFAULTS: dict = {
    "smooth_sigma_um": 2.0,  # Gaussian smoothing before thresholding (physical scale, same for every tumour)
    "threshold": "intermedians",  # per-tumour automatic threshold, no rolling ball (docs/tumour.md)
    "min_object_um3": 1000.0,  # connected objects smaller than this are noise (not satellites)
    "factors": [
        {"name": "genotype", "levels": ["WT", "TLR7 KO"], "reference": "WT", "role": "row", "source": "session",
         "aliases": {"KO7": "TLR7 KO", "KO": "TLR7 KO"}},
        {"name": "tumour_line", "levels": ["cntrl", "HERV", "N4"], "reference": "cntrl", "role": "x",
         "source": "filename", "aliases": {}},
        {"name": "clodronate", "levels": ["-clodro", "+clodro"], "reference": "-clodro", "role": "panel",
         "source": "filename", "aliases": {"-clodron": "-clodro"}},
    ],
}

TEMPLATE = '''# kerncount settings — tumour volume experiment (docs/tumour.md); anything left out falls back to defaults.

[project]
raw_dir = "{raw_dir}"   # raw experiment folder: one sub-folder per imaging session, flat z-stack files
date_order = "DMY"        # dates in file names: 030826 = 3 Aug 2026
csv_delimiter = ","       # set to ";" if Excel (German/European locale) shows everything in one column

[experiment]
type = "tumour_volume"    # mCherry tumours in brain sections; one tumour (hemisphere A or B) = one n
notes = ""

[tumour]
smooth_sigma_um = 2.0     # Gaussian smoothing (µm) before the threshold; identical for every tumour
threshold = "intermedians"   # per-tumour automatic threshold without rolling ball: "intermedians" or "otsu"
min_object_um3 = 1000.0   # smaller connected objects are noise; larger ones besides the main mass are satellites

# Factors: one column per factor in layout.csv. levels = allowed values, reference = comparison baseline,
# role = where the factor goes in group plots (row, panel, x, colour), source = where the value is read
# ("filename" token or imaging-"session" folder name), aliases = other spellings mapped onto a level.
[[tumour.factors]]
name = "genotype"
levels = ["WT", "TLR7 KO"]
reference = "WT"
role = "row"
source = "session"
aliases = {{ KO7 = "TLR7 KO", KO = "TLR7 KO" }}

[[tumour.factors]]
name = "tumour_line"
levels = ["cntrl", "HERV", "N4"]
reference = "cntrl"
role = "x"
source = "filename"

[[tumour.factors]]
name = "clodronate"
levels = ["-clodro", "+clodro"]
reference = "-clodro"
role = "panel"
source = "filename"
aliases = {{ "-clodron" = "-clodro" }}
'''


@dataclass(frozen=True)
class Factor:
    name: str
    levels: tuple[str, ...]
    reference: str
    role: str = ""
    source: str = "filename"
    aliases: dict[str, str] = field(default_factory=dict)

    def match(self, token: str) -> tuple[str, bool] | None:
        """(level, via_alias) when `token` names a level of this factor (case-insensitive), else None."""
        t = token.strip().lower()
        for lv in self.levels:
            if t == lv.lower():
                return lv, False
        for k, v in self.aliases.items():
            if t == str(k).lower():
                return str(v), True
        return None


def tumour_settings(settings: dict) -> dict:
    """[tumour] merged over DEFAULTS (a user `factors` list replaces the default list as a whole)."""
    out = copy.deepcopy(DEFAULTS)
    out.update(copy.deepcopy(settings.get("tumour") or {}))
    return out


def factors(settings: dict) -> list[Factor]:
    """The declared factors, validated (levels non-empty, reference is a level, known role / source)."""
    out = []
    for f in tumour_settings(settings)["factors"]:
        name = str(f["name"])
        levels = tuple(str(x) for x in f["levels"])
        ref = str(f.get("reference") or levels[0])
        role, source = str(f.get("role", "")), str(f.get("source", "filename"))
        if not levels or ref not in levels:
            raise ValueError(f"[tumour.factors] {name}: reference {ref!r} is not one of its levels {list(levels)}")
        if role and role not in ROLES:
            raise ValueError(f"[tumour.factors] {name}: role {role!r} must be one of {', '.join(ROLES)}")
        if source not in SOURCES:
            raise ValueError(f"[tumour.factors] {name}: source {source!r} must be one of {', '.join(SOURCES)}")
        aliases = {str(k): str(v) for k, v in (f.get("aliases") or {}).items()}
        bad = [v for v in aliases.values() if v not in levels]
        if bad:
            raise ValueError(f"[tumour.factors] {name}: aliases point to unknown levels {bad}")
        out.append(Factor(name, levels, ref, role, source, aliases))
    return out
