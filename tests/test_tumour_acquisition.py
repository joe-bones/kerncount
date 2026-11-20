"""Acquisition settings and the shared-threshold robustness check (ticket 11, ADR 0002).

Seams: OME-XML -> settings (acquisition.parse), settings table -> consistency verdict and per-day group table,
per-tumour thresholds -> one experiment-wide threshold -> re-segmentation (shared.measurements), gate + outputs
(shared.analyse), end to end through `kerncount run`."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tumour_synth import SynthTumour, write_ome

OME = """<?xml version="1.0"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
  <Instrument ID="Instrument:1">
    <Microscope Model="Axio Imager.Z2" Type="Upright" />
    <Laser ID="LightSource:1" Power="10" />
    <Detector Model="GaAsP-Pmt2" Offset="0" Zoom="1" ID="Detector:1" />
    <Objective Model="Plan-Apochromat 10x/0.45 M27" ID="Objective:1" />
  </Instrument>
  <Image ID="Image:1" Name="">
    <AcquisitionDate>{date}</AcquisitionDate>
    <Pixels ID="Pixels:1" DimensionOrder="XYZTC" Type="uint8" SizeX="1024" SizeY="1024" SizeZ="{nz}" SizeC="1"
            SizeT="1" PhysicalSizeX="0.6239258" PhysicalSizeY="0.6239258" PhysicalSizeZ="{dz}">
      <Channel ID="Channel:1:0" Name="mCher" PinholeSize="26.2007847" ExcitationWavelength="587"
               EmissionWavelength="610" Fluor="mCherry">
        <LightSourceSettings ID="LightSource:1" Attenuation="0.998" Wavelength="561" />
        <DetectorSettings ID="Detector:1" Offset="0" Voltage="{volt}" />
      </Channel>
    </Pixels>
  </Image>
</OME>"""


def _xml(date="2026-08-18T12:17:03.43Z", nz=53, dz=2.36, volt=596.9971) -> str:
    return OME.format(date=date, nz=nz, dz=dz, volt=volt)


def test_parse_reads_gain_laser_pinhole_channel_voxel_and_acquisition_day():
    from kerncount.tumour import acquisition

    a = acquisition.parse(_xml())
    assert a["detector_voltage"] == "596.9971" and a["detector_offset"] == "0"
    assert a["laser_power"] == "10" and a["laser_attenuation"] == "0.998" and a["laser_wavelength_nm"] == "561"
    assert a["pinhole"] == "26.2007847" and a["fluor"] == "mCherry" and a["channel"] == "mCher"
    assert a["objective"].startswith("Plan-Apochromat 10x") and a["pixel_type"] == "uint8"
    assert a["voxel_z_um"] == "2.36" and a["size_z"] == "53"
    assert a["acquisition_day"] == "2026-08-18" and a["acquisition_time"] == "12:17"
    assert acquisition.parse("")["detector_voltage"] == ""


def _acq(over: dict | None = None) -> pd.DataFrame:
    from kerncount.tumour import acquisition

    rows = []
    for i, (line, day) in enumerate([("cntrl", "2026-08-15T09:00Z"), ("cntrl", "2026-08-17T09:00Z"),
                                     ("N4", "2026-08-19T09:00Z"), ("N4", "2026-08-19T10:00Z")]):
        rows.append({"tumour_id": f"t{i}", "file": f"f{i}", **acquisition.parse(
            _xml(date=day, dz=2.36 if i else 1.75, volt=(over or {}).get(i, 596.9971)))})
    return pd.DataFrame(rows)


def test_consistency_identical_when_only_z_step_and_day_differ():
    from kerncount.tumour import acquisition

    c = acquisition.consistency(_acq())
    assert c.identical and not c.differs and c.values["detector_voltage"] == "596.9971"
    assert "identical" in c.describe()


def test_consistency_names_a_differing_gain_and_settings_never_recorded():
    from kerncount.tumour import acquisition

    a = _acq({2: 650.0})
    c = acquisition.consistency(a)
    assert not c.identical and c.differs["detector_voltage"] == ["596.9971", "650.0"]
    assert "detector_voltage" in c.describe()
    a["pinhole"] = ""
    assert "pinhole" in acquisition.consistency(a).missing


def test_report_lists_groups_per_acquisition_day_and_single_day_groups():
    from kerncount.tumour import acquisition
    from kerncount.tumour.config import Factor

    line = Factor("tumour_line", ("cntrl", "N4"), "cntrl", "x")
    tum = pd.DataFrame({"tumour_id": ["t0", "t1", "t2", "t3"], "tumour_line": ["cntrl", "cntrl", "N4", "N4"]})
    md = "\n".join(acquisition.report_lines(_acq(), tum, [line]))
    assert "| N4 |  |  | 2 |" in md and "| cntrl | 1 | 1 |  |" in md
    assert "single day: N4 (2026-08-19)" in md
    assert "voxel_z_um varies" in md and "Detector gain is recorded" in md


def test_shared_threshold_is_the_median_of_the_per_tumour_thresholds():
    from kerncount.tumour import shared

    assert shared.value(pd.DataFrame({"threshold": [80.0, 100.0, 130.0, np.nan]})) == 100.0
    with pytest.raises(ValueError):
        shared.value(pd.DataFrame({"threshold": [np.nan]}))


def _project(tmp_path, brightness=(1.0, 0.5)):
    """Two identical tumours, the second recorded at lower brightness; settings pointing at their folder."""
    from test_tumour_stats import _settings

    raw = tmp_path / "raw"
    files = []
    for i, b in enumerate(brightness):
        t = SynthTumour(semi_axes_um=(15, 40, 35), voxel_um=(2.36, 1.25, 1.25), shape=(24, 96, 96), plateau=90,
                        brightness=b, seed=7)
        write_ome(raw / f"t{i}.ome.tiff", t.stack(), t.voxel_um)
        files.append(f"t{i}.ome.tiff")
    s = _settings()
    s["project"]["raw_dir"] = str(raw)
    return s, files


def test_shared_threshold_shrinks_the_dim_tumour_while_its_own_threshold_does_not(tmp_path):
    """The primary per-tumour threshold measures both copies alike; at one grey value the dimmer copy shrinks —
    exactly the brightness effect the check exposes."""
    from kerncount.tumour import segment, shared
    from kerncount.tumour.config import tumour_settings

    s, files = _project(tmp_path)
    ts = tumour_settings(s)
    prim = [segment.measure_file(tmp_path / "raw" / f, ts).readouts for f in files]
    assert prim[1]["cross_sectional_area_um2"] == pytest.approx(prim[0]["cross_sectional_area_um2"], rel=0.05)
    table = pd.DataFrame({"tumour_id": ["a", "b"], "file": files, "threshold": [p["threshold"] for p in prim]})
    t = shared.value(table)
    m = shared.measurements(tmp_path, table, s, t).set_index("tumour_id")
    assert m["threshold"].tolist() == pytest.approx([t, t], abs=2e-3)
    assert m.loc["b", "cross_sectional_area_um2"] < 0.9 * prim[1]["cross_sectional_area_um2"]
    assert m.loc["a", "cross_sectional_area_um2"] > m.loc["b", "cross_sectional_area_um2"]
    cache = sorted(p.name for p in (tmp_path / "_work" / "tumour").glob("*.shared.json"))
    assert len(cache) == 2
    mtimes = [p.stat().st_mtime_ns for p in (tmp_path / "_work" / "tumour").glob("*.shared.json")]
    shared.measurements(tmp_path, table, s, t)  # cached
    assert [p.stat().st_mtime_ns for p in (tmp_path / "_work" / "tumour").glob("*.shared.json")] == mtimes


def test_shared_check_is_not_run_when_switched_off_or_acquisition_settings_differ(tmp_path):
    from kerncount.tumour import acquisition, shared
    from kerncount.tumour.config import factors, families

    s, files = _project(tmp_path)
    table = pd.DataFrame({"tumour_id": ["a", "b"], "file": files, "threshold": [50.0, 25.0],
                          "cross_sectional_area_mm2": [0.1, 0.1]})
    report = tmp_path / "stats_report.md"
    report.write_text("# r\n")
    assert shared.analyse(tmp_path, table, s, factors(s), families(s), report) == []
    assert "shared_threshold_check = false" in report.read_text()
    s["tumour"]["shared_threshold_check"] = True
    differ = acquisition.consistency(_acq({2: 650.0}))
    assert shared.analyse(tmp_path, table, s, factors(s), families(s), report, differ) == []
    assert "not run: acquisition settings differ" in report.read_text()
    assert not (tmp_path / shared.FILENAME).exists()


def test_a_check_that_is_not_run_removes_the_outputs_of_an_earlier_run(tmp_path):
    """Ticket 11 acceptance: the KO-only run wrote shared_threshold.csv + plots; after the WT set arrived the check
    was ruled out (settings differ) but the stale 48-tumour table stayed in the results (and the GUI lists it)."""
    from kerncount.tumour import acquisition, shared
    from kerncount.tumour.config import factors, families

    s, files = _project(tmp_path)
    s["tumour"]["shared_threshold_check"] = True
    table = pd.DataFrame({"tumour_id": ["a", "b"], "file": files, "threshold": [50.0, 25.0],
                          "cross_sectional_area_mm2": [0.1, 0.1]})
    stale = [tmp_path / shared.FILENAME, tmp_path / "plots" / "cross_sectional_area_vs_shared_threshold.png",
             tmp_path / "plots" / "perimeter_vs_shared_threshold.pdf"]
    (tmp_path / "plots").mkdir(exist_ok=True)
    for p in stale:
        p.write_text("from an earlier run")
    keep = tmp_path / "plots" / "cross_sectional_area_vs_fiji.png"
    keep.write_text("not ours")
    differ = acquisition.consistency(_acq({2: 650.0}))
    assert shared.analyse(tmp_path, table, s, factors(s), families(s), None, differ) == []
    assert [p for p in stale if p.exists()] == []
    assert keep.exists()


def test_run_writes_acquisition_table_and_the_shared_threshold_check(tmp_path, monkeypatch):
    from kerncount import cli
    from test_tumour import _experiment

    _experiment(tmp_path, monkeypatch)
    prj = tmp_path / "analysis" / "tumour_exp"
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    p = prj / "settings.toml"
    assert "shared_threshold_check = false" in p.read_text()
    p.write_text(p.read_text().replace("shared_threshold_check = false", "shared_threshold_check = true"))
    cli.main(["run", "tumour_exp", "--no-review"])
    acq = pd.read_csv(prj / "acquisition.csv")
    assert len(acq) == 4 and set(acq["pixel_type"]) == {"uint8"}
    cmp = pd.read_csv(prj / "shared_threshold.csv")
    assert set(cmp["readout"]) == {"cross_sectional_area", "perimeter"}
    assert {"shared_estimate", "direction_agree"} <= set(cmp.columns)
    assert (prj / "plots" / "cross_sectional_area_vs_shared_threshold.png").exists()
    md = (prj / "stats_report.md").read_text()
    assert "## Acquisition settings" in md and "## Shared-threshold robustness check" in md
    assert "Shared threshold = median" in md
    assert "Shared threshold" in pd.ExcelFile(prj / "results.xlsx").sheet_names
