"""[tumour] settings: defaults for the TLR7 tumour experiment and the declared factors (ADR 0004 amendment).

A factor is one crossed condition axis on every tumour. Nothing downstream assumes how many factors there are:
layout columns, tumour IDs, groups and plot facets are all built from this list.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

ROLES = ("row", "panel", "x", "colour")
SOURCES = ("filename", "session")
EFFECTS = ("simple", "interaction")  # simple effect of one factor / interaction of several, within cells of the rest
FAMILY_ROLES = ("confirmatory", "exploratory", "descriptive")
READOUT_ROLES = ("confirmatory", "exploratory")

DEFAULTS: dict = {
    "smooth_sigma_um": 2.0,  # Gaussian smoothing before thresholding (physical scale, same for every tumour)
    "threshold": "intermedians",  # per-tumour automatic threshold, no rolling ball (docs/tumour.md)
    "min_object_um3": 1000.0,  # connected objects smaller than this are noise (not satellites)
    # Sensitivity models (ticket 08): random permutations per contrast (when not enumerated exactly) and their seed.
    "sensitivity_permutations": 9999,
    "sensitivity_seed": 20261119,
    # Power planner (ticket 09, planning aid only): fold-change to plan for, target power, two-sided alpha per family.
    "planner_fold_change": 1.5,
    "planner_power": 0.8,
    "planner_alpha": 0.05,
    "factors": [
        {"name": "genotype", "levels": ["WT", "TLR7 KO"], "reference": "WT", "role": "row", "source": "session",
         "aliases": {"KO7": "TLR7 KO", "KO": "TLR7 KO"}},
        {"name": "tumour_line", "levels": ["cntrl", "HERV", "N4"], "reference": "cntrl", "role": "x",
         "source": "filename", "aliases": {}},
        {"name": "clodronate", "levels": ["-clodro", "+clodro"], "reference": "-clodro", "role": "panel",
         "source": "filename", "aliases": {"-clodron": "-clodro"}},
    ],
    # Readouts analysed by the factorial statistics (ADR 0003): column in tumours.csv + confirmatory/exploratory.
    # Only confirmatory readouts in confirmatory families carry significance stars.
    "readouts": [
        {"name": "main_mass_volume", "column": "main_mass_volume_mm3", "label": "Main-mass volume (mm³)",
         "role": "confirmatory"},
    ],
    # Factor whose reference level normalised values are expressed against ("" = the factor with role x).
    "normalise_by": "",
    # Contrast families (ADR 0003 amendment), Holm within each family and readout. effect = "simple": each level of
    # the factor vs its reference within every cell of the other factors; "interaction": ratio of ratios of the listed
    # factors (each non-reference vs reference) within every cell of the rest. role = confirmatory (stars),
    # exploratory (p-values, no stars) or descriptive (estimate + CI only, never tested).
    "families": [
        {"id": "F1", "effect": "simple", "factors": ["tumour_line"], "role": "confirmatory",
         "label": "tumour line vs cntrl within genotype × clodronate"},
        {"id": "F2", "effect": "simple", "factors": ["clodronate"], "role": "confirmatory",
         "label": "+clodro vs −clodro within genotype × tumour line"},
        {"id": "F3", "effect": "interaction", "factors": ["tumour_line", "clodronate"], "role": "exploratory",
         "label": "tumour line × clodronate interaction within genotype"},
        {"id": "F4", "effect": "interaction", "factors": ["genotype", "tumour_line"], "role": "confirmatory",
         "label": "genotype × tumour line (ratio of line effects) within clodronate"},
        {"id": "F5", "effect": "interaction", "factors": ["genotype", "clodronate"], "role": "exploratory",
         "label": "genotype × clodronate (ratio of clodronate effects) within tumour line"},
        {"id": "F6", "effect": "simple", "factors": ["genotype"], "role": "descriptive",
         "label": "direct TLR7 KO vs WT within tumour line × clodronate"},
        {"id": "F7", "effect": "interaction", "factors": ["genotype", "tumour_line", "clodronate"],
         "role": "exploratory", "label": "three-way interaction"},
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
sensitivity_permutations = 9999   # sensitivity: random permutations per contrast (exact enumeration when fewer splits)
sensitivity_seed = 20261119       # seed of the permutations (recorded in the report; same seed = same p-values)
planner_fold_change = 1.5         # power planner (planning aid): tumours per group needed to detect this fold-change

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

# Statistics (ADR 0003): LOCKED before results are viewed — the run records a hash of factors, readouts and families.
# Readouts: column in tumours.csv; role = "confirmatory" (stars) or "exploratory" (effect sizes, labelled).
[[tumour.readouts]]
name = "main_mass_volume"
column = "main_mass_volume_mm3"
label = "Main-mass volume (mm³)"
role = "confirmatory"

# Contrast families, Holm within each. effect = "simple" (each level of the factor vs its reference, within every
# combination of the other factors) or "interaction" (ratio of ratios of the listed factors, within the rest; the last
# factor is the inner ratio). role = confirmatory / exploratory / descriptive (never tested).
[[tumour.families]]
id = "F1"
effect = "simple"
factors = ["tumour_line"]
role = "confirmatory"
label = "tumour line vs cntrl within genotype × clodronate"

[[tumour.families]]
id = "F2"
effect = "simple"
factors = ["clodronate"]
role = "confirmatory"
label = "+clodro vs −clodro within genotype × tumour line"

[[tumour.families]]
id = "F3"
effect = "interaction"
factors = ["tumour_line", "clodronate"]
role = "exploratory"
label = "tumour line × clodronate interaction within genotype"

[[tumour.families]]
id = "F4"
effect = "interaction"
factors = ["genotype", "tumour_line"]
role = "confirmatory"
label = "genotype × tumour line (ratio of line effects) within clodronate"

[[tumour.families]]
id = "F5"
effect = "interaction"
factors = ["genotype", "clodronate"]
role = "exploratory"
label = "genotype × clodronate (ratio of clodronate effects) within tumour line"

[[tumour.families]]
id = "F6"
effect = "simple"
factors = ["genotype"]
role = "descriptive"
label = "direct TLR7 KO vs WT within tumour line × clodronate"

[[tumour.families]]
id = "F7"
effect = "interaction"
factors = ["genotype", "tumour_line", "clodronate"]
role = "exploratory"
label = "three-way interaction"
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


@dataclass(frozen=True)
class Readout:
    """One analysed readout: `column` in tumours.csv, a short `name` for file names, a plot `label` and its role."""
    name: str
    column: str
    label: str
    role: str = "exploratory"


@dataclass(frozen=True)
class Family:
    """One contrast family (ADR 0003 amendment): `effect` of `factors`, within every cell of the other factors."""
    id: str
    effect: str
    factors: tuple[str, ...]
    role: str
    label: str = ""


def readouts(settings: dict) -> list[Readout]:
    out = []
    for r in tumour_settings(settings)["readouts"]:
        col = str(r["column"])
        name = str(r.get("name") or col)
        role = str(r.get("role", "exploratory"))
        if role not in READOUT_ROLES:
            raise ValueError(f"[tumour.readouts] {name}: role {role!r} must be one of {', '.join(READOUT_ROLES)}")
        out.append(Readout(name, col, str(r.get("label") or name), role))
    return out


def families(settings: dict) -> list[Family]:
    """Declared contrast families, validated for shape. Factor names are checked against the factors by the caller
    (`families.usable`), so a family on an undeclared factor is reported, not silently dropped."""
    out = []
    for f in tumour_settings(settings)["families"]:
        fid = str(f["id"])
        effect, role = str(f.get("effect", "simple")), str(f.get("role", "exploratory"))
        fs = tuple(str(x) for x in f.get("factors") or ())
        if effect not in EFFECTS:
            raise ValueError(f"[tumour.families] {fid}: effect {effect!r} must be one of {', '.join(EFFECTS)}")
        if role not in FAMILY_ROLES:
            raise ValueError(f"[tumour.families] {fid}: role {role!r} must be one of {', '.join(FAMILY_ROLES)}")
        if not fs or (effect == "simple" and len(fs) != 1) or (effect == "interaction" and len(fs) < 2):
            raise ValueError(f"[tumour.families] {fid}: a simple effect names one factor, an interaction two or more")
        out.append(Family(fid, effect, fs, role, str(f.get("label", ""))))
    return out


def normalise_by(settings: dict, facs: list[Factor]) -> Factor | None:
    """The focus factor for normalised values: `[tumour] normalise_by`, else the factor plotted on x."""
    name = str(tumour_settings(settings).get("normalise_by") or "")
    if name:
        f = next((f for f in facs if f.name == name), None)
        if f is None:
            raise ValueError(f"[tumour] normalise_by = {name!r} is not a declared factor")
        return f
    return next((f for f in facs if f.role == "x"), facs[0] if facs else None)
