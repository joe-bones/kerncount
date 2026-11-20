"""Acquisition settings of every original stack (ticket 11, ADR 0002): recorded, checked for consistency, reported.

Read from each file's OME-XML (no pixels are read). Two groups of settings:

- `CHECKED` — settings that change the recorded brightness or the optical section (detector gain = PMT voltage,
  offset, laser line / power / attenuation, pinhole, objective, zoom, channel, pixel type, xy pixel size). They must
  be identical across tumours for grey values to be comparable between tumours; only then is one experiment-wide
  threshold (the shared-threshold robustness check, `shared.py`) meaningful.
- `REPORTED` — settings that may legitimately differ (z-step, frame size, number of layers, acquisition day). They are
  listed, not judged. The acquisition day matters for the design: groups imaged on different days can differ by a
  day effect even with identical settings (laser ageing, alignment). Whether the statistics treat a session folder
  or an acquisition day as one imaging batch is the `[tumour] batch` setting (ticket 12).

Outputs: `acquisition.csv` (one row per tumour) and a report section `## Acquisition settings`.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

FILENAME = "acquisition.csv"
CHECKED = ["microscope", "objective", "detector", "detector_voltage", "detector_offset", "detector_zoom",
           "laser_wavelength_nm", "laser_power", "laser_attenuation", "pinhole", "channel", "fluor",
           "excitation_nm", "emission_nm", "pixel_type", "voxel_x_um", "voxel_y_um"]
REPORTED = ["acquisition_day", "acquisition_time", "voxel_z_um", "size_x", "size_y", "size_z"]
COLUMNS = [*CHECKED, *REPORTED]


def _tag(el: ET.Element) -> str:
    return el.tag.rsplit("}", 1)[-1]


def _first(root: ET.Element, tag: str) -> ET.Element | None:
    return next((el for el in root.iter() if _tag(el) == tag), None)


def parse(xml: str) -> dict[str, str]:
    """Acquisition settings from an OME-XML string (first image); '' for anything not recorded."""
    out = {c: "" for c in COLUMNS}
    if not xml:
        return out
    root = ET.fromstring(xml)

    def get(tag: str, attr: str) -> str:
        el = _first(root, tag)
        return "" if el is None else (el.get(attr) or "")

    d = _first(root, "AcquisitionDate")
    date = (d.text or "").strip() if d is not None else ""
    out.update(microscope=get("Microscope", "Model"), objective=get("Objective", "Model"),
               detector=get("Detector", "Model"), detector_zoom=get("Detector", "Zoom"),
               detector_voltage=get("DetectorSettings", "Voltage") or get("Detector", "Voltage"),
               detector_offset=get("DetectorSettings", "Offset") or get("Detector", "Offset"),
               laser_wavelength_nm=get("LightSourceSettings", "Wavelength"),
               laser_power=get("Laser", "Power"), laser_attenuation=get("LightSourceSettings", "Attenuation"),
               pinhole=get("Channel", "PinholeSize"), channel=get("Channel", "Name"), fluor=get("Channel", "Fluor"),
               excitation_nm=get("Channel", "ExcitationWavelength"), emission_nm=get("Channel", "EmissionWavelength"),
               pixel_type=get("Pixels", "Type"), voxel_x_um=get("Pixels", "PhysicalSizeX"),
               voxel_y_um=get("Pixels", "PhysicalSizeY"), voxel_z_um=get("Pixels", "PhysicalSizeZ"),
               size_x=get("Pixels", "SizeX"), size_y=get("Pixels", "SizeY"), size_z=get("Pixels", "SizeZ"),
               acquisition_day=date[:10], acquisition_time=date[11:16])
    return out


def read(path: Path) -> dict[str, str]:
    """`parse` on the OME-XML of one original stack (header only)."""
    import tifffile

    with tifffile.TiffFile(path) as tf:
        return parse(tf.ome_metadata or "")


def table(raw: Path, tumours: pd.DataFrame) -> pd.DataFrame:
    """tumour_id + every acquisition setting, one row per tumour of `tumours` (needs its `file` column)."""
    rows = [{"tumour_id": r["tumour_id"], "file": r["file"], **read(Path(raw) / r["file"])}
            for _, r in tumours.iterrows()]
    return pd.DataFrame(rows, columns=["tumour_id", "file", *COLUMNS])


@dataclass
class Consistency:
    """identical: no CHECKED setting differs between tumours (settings not recorded anywhere are in `missing`)."""
    identical: bool
    differs: dict[str, list[str]]  # setting -> distinct recorded values
    missing: list[str]  # CHECKED settings recorded for no tumour
    partly: list[str]  # CHECKED settings recorded for some tumours only
    values: dict[str, str]  # setting -> the shared value (identical, recorded settings)

    def describe(self) -> str:
        if self.differs:
            return "acquisition settings differ between tumours: " + "; ".join(
                f"{k} ({', '.join(v)})" for k, v in self.differs.items())
        s = "acquisition settings identical across tumours"
        gaps = [*self.missing, *self.partly]
        return s + (f" (not recorded, so not checked: {', '.join(gaps)})" if gaps else "")


def consistency(acq: pd.DataFrame) -> Consistency:
    differs, missing, partly, values = {}, [], [], {}
    for c in CHECKED:
        v = acq[c].fillna("").astype(str).str.strip() if c in acq.columns else pd.Series([""] * len(acq))
        rec = v[v != ""]
        if not len(rec):
            missing.append(c)
            continue
        if len(rec) < len(v):
            partly.append(c)
        distinct = sorted(set(rec))
        if len(distinct) > 1:
            differs[c] = distinct
        else:
            values[c] = distinct[0]
    return Consistency(not differs, differs, missing, partly, values)


def report_lines(acq: pd.DataFrame, tumours: pd.DataFrame, facs: list, batch: str = "session_folder") -> list[str]:
    """Markdown: checked settings (shared values or what differs), settings that vary, groups x acquisition day."""
    c = consistency(acq)
    md = ["## Acquisition settings (from each original's OME metadata)", "",
          f"- Verdict: **{c.describe()}**."]
    if c.values:
        md.append("- Shared values: " + "; ".join(f"{k} = {v}" for k, v in c.values.items()) + ".")
    if "detector_voltage" in c.values:
        md.append("- Detector gain is recorded as the PMT voltage (`detector_voltage`).")
    for k in ("voxel_z_um", "size_z", "size_x", "size_y"):
        vals = sorted(set(acq[k].astype(str)) - {""}, key=lambda s: float(s) if s.replace(".", "", 1).isdigit() else 0)
        if len(vals) > 1:
            md.append(f"- {k} varies (expected; not a brightness setting): {', '.join(vals)}"
                      + (" — counts: " + ", ".join(f"{v} × {n}" for v, n in acq[k].astype(str).value_counts()
                                                    .sort_index().items()) if k == "voxel_z_um" else "") + ".")
    days = acq["acquisition_day"].astype(str)
    if (days != "").any():
        names = [f.name for f in facs if f.name in tumours.columns]
        t = tumours[["tumour_id", *names]].merge(acq[["tumour_id", "acquisition_day"]], on="tumour_id")
        t["group"] = t[names].astype(str).agg(" / ".join, axis=1) if names else "all"
        ct = pd.crosstab(t["group"], t["acquisition_day"])
        rule = ("`[tumour] batch = \"session_folder\"`: the confounding rule of the statistics treats each session "
                "folder as one batch, not each day" if batch == "session_folder" else
                "`[tumour] batch = \"acquisition_day\"`: the confounding rule of the statistics treats each day as "
                "one batch")
        md += ["", f"Tumours per group and acquisition day ({ct.shape[1]} day(s)); {rule}. Groups imaged on "
               "different days can differ by a day effect even with identical settings.", "",
               "| group | " + " | ".join(ct.columns) + " |", "|---|" + "---|" * ct.shape[1]]
        for g, r in ct.iterrows():
            md.append(f"| {g} | " + " | ".join(str(int(x)) if x else "" for x in r) + " |")
        single = [f"{g} ({ct.columns[r.to_numpy().nonzero()[0][0]]})" for g, r in ct.iterrows()
                  if (r > 0).sum() == 1]
        if ct.shape[1] > 1 and single:
            md += ["", "- Groups imaged on a single day: " + "; ".join(single) + "."]
    return md + [""]


def analyse(out: Path, tumours: pd.DataFrame, settings: dict, facs: list, report: Path | None = None
            ) -> tuple[list[str], Consistency | None]:
    """Write acquisition.csv and the report section; returns (written files, consistency) — (.., None) when the
    table has no file column or the raw folder is unavailable."""
    from .measure import raw_dir

    if "file" not in tumours.columns or not len(tumours) or not settings.get("project", {}).get("raw_dir"):
        return [], None
    raw = raw_dir(Path(out), settings)
    if not raw.is_dir():
        return [], None
    acq = table(raw, tumours)
    path = Path(out) / FILENAME
    acq.to_csv(path, index=False, sep=settings.get("project", {}).get("csv_delimiter", ","))
    if report is not None and report.exists():
        with report.open("a", encoding="utf-8") as fh:
            from .config import batch

            fh.write("\n" + "\n".join(report_lines(acq, tumours, facs, batch(settings))) + "\n")
    return [str(path)], consistency(acq)
