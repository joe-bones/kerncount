"""Any number of imaging-session folders and genotypes (ticket 12).

Seams: session-folder name -> factor levels + session label (aliases, tokens, unknown), multi-folder scan -> unique
tumour IDs, file date vs OME acquisition day, layout migration across a folder rename / reorganisation, cache and
review-code stability across a rename, statistics / families / plots with 1, 2 and 3 genotype levels, and the
`[tumour] batch` setting -> which contrasts are confounded and how they are labelled."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import discover, experiment_types, settings as settings_mod

OME = ('<?xml version="1.0" encoding="UTF-8"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
       '<Image ID="Image:0"><AcquisitionDate>{date}</AcquisitionDate><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
       'Type="uint8" SizeX="4" SizeY="4" SizeZ="2" SizeC="1" SizeT="1" PhysicalSizeX="1.0" PhysicalSizeY="1.0" '
       'PhysicalSizeZ="2.0"><Channel ID="Channel:0:0" SamplesPerPixel="1"/><TiffData/></Pixels></Image></OME>')


def _settings(tmp_path: Path, extra: str = "") -> dict:
    (tmp_path / "settings.toml").write_text('[experiment]\ntype = "tumour_volume"\n' + extra)
    return settings_mod.load(tmp_path)


def _touch(folder: Path, *names: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for n in names:
        (folder / n).write_bytes(b"")


def _ome(path: Path, day: str) -> None:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, np.zeros((2, 4, 4), np.uint8), description=OME.format(date=f"{day}T10:00:00"),
                     metadata=None)


def _scan(tmp_path: Path, files: dict[str, list[str]], extra: str = "") -> tuple[Path, dict, pd.DataFrame]:
    raw = tmp_path / "raw"
    for session, names in files.items():
        _touch(raw / session if session else raw, *names)
    s = _settings(tmp_path, extra)
    return raw, s, experiment_types.get("tumour_volume").scan_layout(raw, s)


ONE = ["030826 HERV -clodro n1 A.ome.tiff"]

# ---------------------------------------------------------------- seam 1: folder name -> factor levels


@pytest.mark.parametrize("folder, genotype, label", [
    ("5d K07", "TLR7 KO", "5d"),
    ("5d WT", "WT", "5d"),
    ("ready 5d KO7 030826", "TLR7 KO", "ready 5d 030826"),
    ("TLR7 KO 10d", "TLR7 KO", "10d"),
    ("tlr7-ko_5d", "TLR7 KO", "5d"),
    ("5d TLR7KO", "TLR7 KO", "5d"),
    ("5d k07 repeat", "TLR7 KO", "5d repeat"),
    ("wt", "WT", ""),
])
def test_session_folder_names_are_read_token_wise_through_aliases(tmp_path, folder, genotype, label):
    _, _, layout = _scan(tmp_path, {folder: ONE})
    row = layout.iloc[0]
    assert row["genotype"] == genotype
    assert row["imaging_session"] == folder and row["session_label"] == label
    assert row["tumour_id"] == f"{genotype} / HERV / -clodro / n1A"


def test_a_new_genotype_and_a_timepoint_factor_are_a_settings_change_only(tmp_path):
    extra = '''
[[tumour.factors]]
name = "genotype"
levels = ["WT", "TLR7 KO", "TLR9 KO"]
reference = "WT"
role = "row"
source = "session"
aliases = { K07 = "TLR7 KO", K09 = "TLR9 KO" }

[[tumour.factors]]
name = "timepoint"
levels = ["5d", "10d"]
reference = "5d"
role = "colour"
source = "session"

[[tumour.factors]]
name = "tumour_line"
levels = ["cntrl", "HERV", "N4"]
reference = "cntrl"
role = "x"

[[tumour.factors]]
name = "clodronate"
levels = ["-clodro", "+clodro"]
reference = "-clodro"
role = "panel"
'''
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {"5d K07": ONE, "10d K09": ONE, "5d WT": ONE}, extra)
    got = layout.set_index("imaging_session")
    assert got.loc["10d K09", "genotype"] == "TLR9 KO" and got.loc["10d K09", "timepoint"] == "10d"
    assert got.loc["5d K07", "tumour_id"] == "TLR7 KO / 5d / HERV / -clodro / n1A"
    assert got.loc["5d WT", "session_label"] == ""  # every token was read as a factor level
    assert len(lay.check(layout, raw, s).included) == 3


def test_a_folder_without_a_readable_genotype_blocks_naming_the_folder_and_the_fixes(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {"5d K07": ONE, "5d TLR9-KO": ONE + ["030826 HERV -clodro n1 B.ome.tiff"]})
    with pytest.raises(lay.LayoutBlocked) as e:
        lay.check(layout, raw, s)
    assert len(e.value.problems) == 1
    msg = e.value.problems[0]
    assert "'5d TLR9-KO'" in msg and "genotype" in msg and "2 tumours" in msg
    assert "alias" in msg and "layout.csv" in msg and "TLR9-KO" in msg  # names the unread token to alias


def test_hidden_and_system_folders_are_ignored_and_flat_files_still_work(tmp_path):
    raw, s, layout = _scan(tmp_path, {"5d WT": ONE, ".snapshots": ONE, "_old": ONE, "@eaDir": ONE,
                                      "$RECYCLE.BIN": ONE})
    assert list(layout["imaging_session"]) == ["5d WT"]
    flat = tmp_path / "flat" / "5d K07"
    _touch(flat, *ONE)
    layout = experiment_types.get("tumour_volume").scan_layout(flat, s)
    assert list(layout["genotype"]) == ["TLR7 KO"] and layout["file"].iloc[0] == ONE[0]


# ---------------------------------------------------------------- seam 2: several folders -> unique tumour IDs


def test_repeat_sessions_of_one_genotype_get_the_session_in_the_id_only_where_needed(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {
        "5d WT": ["030826 HERV -clodro n1 A.ome.tiff", "030826 cntrl -clodro n1 A.ome.tiff"],
        "5d WT repeat": ["120926 HERV -clodro n1 A.ome.tiff", "120926 N4 -clodro n2 B.ome.tiff"],
        "5d K07": ["030826 HERV -clodro n1 A.ome.tiff"]})
    ids = dict(zip(layout["file"], layout["tumour_id"]))
    assert ids["5d WT/030826 HERV -clodro n1 A.ome.tiff"] == "WT / HERV / -clodro / n1A / 5d WT"
    assert ids["5d WT repeat/120926 HERV -clodro n1 A.ome.tiff"] == "WT / HERV / -clodro / n1A / 5d WT repeat"
    assert ids["5d WT/030826 cntrl -clodro n1 A.ome.tiff"] == "WT / cntrl / -clodro / n1A"  # no clash: plain ID
    assert ids["5d K07/030826 HERV -clodro n1 A.ome.tiff"] == "TLR7 KO / HERV / -clodro / n1A"
    v = lay.check(layout, raw, s)
    assert v.included["tumour_id"].is_unique and len(v.included) == 5


def test_duplicate_ids_within_one_session_still_block_naming_the_files(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {
        "5d WT": ["030826 HERV -clodro n1 A.ome.tiff", "030826 HERV -clodron n1 A.ome.tiff"],
        "5d WT repeat": ["030826 HERV -clodro n1 A.ome.tiff"]})
    with pytest.raises(lay.LayoutBlocked) as e:
        lay.check(layout, raw, s)
    (p,) = e.value.problems
    assert "5d WT/030826 HERV -clodro n1 A.ome.tiff" in p and "5d WT/030826 HERV -clodron n1 A.ome.tiff" in p
    assert "5d WT repeat/" not in p


# ---------------------------------------------------------------- seam 3: dates


def test_file_date_is_optional_and_the_acquisition_day_comes_from_the_ome_metadata(tmp_path):
    from kerncount.tumour import layout as lay

    raw = tmp_path / "raw"
    _ome(raw / "5d K07" / "030826 HERV -clodro n1 A.ome.tiff", "2026-08-15")
    _ome(raw / "5d K07" / "HERV -clodro n1 B.ome.tiff", "2026-08-17")  # no date in the name
    _ome(raw / "5d WT" / "120926 HERV -clodro n1 A.ome.tiff", "2026-09-12")
    _touch(raw / "5d WT", "2026-13-45 HERV -clodro n2 A.ome.tiff")  # no OME, nonsense date
    s = _settings(tmp_path)
    layout = lay.scan(raw, s).set_index("file")
    assert layout.loc["5d K07/030826 HERV -clodro n1 A.ome.tiff", "file_date"] == "2026-08-03"
    assert layout.loc["5d K07/HERV -clodro n1 B.ome.tiff", "file_date"] == ""
    assert layout.loc["5d K07/030826 HERV -clodro n1 A.ome.tiff", "acquisition_day"] == "2026-08-15"
    assert layout.loc["5d K07/HERV -clodro n1 B.ome.tiff", "acquisition_day"] == "2026-08-17"
    assert layout.loc["5d WT/2026-13-45 HERV -clodro n2 A.ome.tiff", "acquisition_day"] == ""
    v = lay.check(layout.reset_index(), raw, s)
    assert len(v.included) == 4  # neither date is required


# ---------------------------------------------------------------- seam 4: layout migration across a rename


def _workspace(tmp_path, monkeypatch, session: str, names: list[str]) -> tuple[Path, Path]:
    raw = tmp_path / "data" / "tumour_exp"
    _touch(raw / session, *names)
    monkeypatch.chdir(tmp_path)
    return raw, tmp_path / "analysis" / "tumour_exp"


def test_a_renamed_session_folder_keeps_layout_edits_and_reinit_carries_them_over(tmp_path, monkeypatch, capsys):
    from kerncount import cli
    from kerncount.tumour import layout as lay

    names = ["030826 HERV -clodro n1 A.ome.tiff", "030826 HERV -clodro n1 B.ome.tiff",
             "030826 cntrl +clodro n2 A.ome.tiff"]
    raw, prj = _workspace(tmp_path, monkeypatch, "ready 5d KO7 030826", names)
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    edited = discover.read_table(prj / "layout.csv")
    b = edited["file"].str.endswith("n1 B.ome.tiff")
    edited.loc[b, ["include", "note"]] = ["no", "out of focus"]
    c = edited["file"].str.contains("cntrl")
    edited.loc[c, "tumour_line"] = "N4"  # a user override
    edited = edited.drop(columns=["session_label", "acquisition_day"]).rename(columns={"file_date": "date"})  # old
    discover.write_layout(edited, prj / "layout.csv")  # layout.csv as written before ticket 12

    (raw / "ready 5d KO7 030826").rename(raw / "5d K07")
    s = settings_mod.load(prj)
    v = lay.check(discover.read_table(prj / "layout.csv"), raw, s)  # the run works before any re-init
    rows = v.layout.set_index("file")
    assert set(rows.index) == {f"5d K07/{n}" for n in names}
    assert rows.loc["5d K07/030826 HERV -clodro n1 B.ome.tiff", "include"] == "no"
    assert v.excluded["exclusion_reason"].tolist() == ["out of focus"]
    cn = rows.loc["5d K07/030826 cntrl +clodro n2 A.ome.tiff"]
    assert cn["tumour_line"] == "N4" and cn["tumour_id"] == "TLR7 KO / N4 / +clodro / n2A"
    assert "tumour_line set in layout.csv to 'N4'" in cn["note"] and "genotype set" not in cn["note"]
    assert (rows["imaging_session"] == "5d K07").all()
    assert any("moved" in w and "kerncount init" in w for w in v.warnings)

    capsys.readouterr()
    cli.main(["init", "tumour_exp"])
    out = capsys.readouterr().out
    assert "ready 5d KO7 030826/030826 HERV -clodro n1 A.ome.tiff -> 5d K07/" in out
    new = discover.read_table(prj / "layout.csv").set_index("file")
    assert set(new.index) == {f"5d K07/{n}" for n in names}
    assert new.loc["5d K07/030826 HERV -clodro n1 B.ome.tiff", "note"] == "out of focus"
    assert new.loc["5d K07/030826 cntrl +clodro n2 A.ome.tiff", "tumour_line"] == "N4"
    assert {"file_date", "session_label", "acquisition_day"} <= set(new.columns) and "date" not in new.columns
    assert (new["file_date"] == "2026-08-03").all()
    backups = list((prj / "_work" / "backups").glob("layout*.csv"))
    assert backups and "ready 5d KO7 030826" in backups[0].read_text()  # never overwritten silently

    capsys.readouterr()
    cli.main(["init", "tumour_exp"])  # nothing moved now: layout.csv kept, fresh scan beside it
    assert "layout.csv already exists (kept)" in capsys.readouterr().out
    assert discover.read_table(prj / "layout.csv").equals(new.reset_index()[discover.read_table(
        prj / "layout.csv").columns])


def test_a_file_moved_to_another_genotype_folder_takes_that_folders_genotype(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {"5d WT": ["030826 HERV -clodro n1 A.ome.tiff"],
                                      "5d K07": ["030826 HERV -clodro n2 A.ome.tiff"]})
    layout.loc[layout["imaging_session"] == "5d WT", "note"] = "checked"
    (raw / "5d WT" / "030826 HERV -clodro n1 A.ome.tiff").rename(raw / "5d K07" / "030826 HERV -clodro n1 A.ome.tiff")
    v = lay.check(layout, raw, s)
    row = v.layout.set_index("file").loc["5d K07/030826 HERV -clodro n1 A.ome.tiff"]
    assert row["genotype"] == "TLR7 KO" and row["note"] == "checked"  # misfiled: the folder decides, edit kept


def test_new_session_folders_are_added_at_reinit_without_losing_edits(tmp_path, monkeypatch):
    from kerncount import cli

    raw, prj = _workspace(tmp_path, monkeypatch, "5d K07", ONE)
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    lay_ = discover.read_table(prj / "layout.csv")
    lay_.loc[0, "note"] = "nice one"
    discover.write_layout(lay_, prj / "layout.csv")
    _touch(raw / "5d WT", *ONE)
    cli.main(["init", "tumour_exp"])
    new = discover.read_table(prj / "layout.csv").set_index("file")
    assert new.loc["5d K07/" + ONE[0], "note"] == "nice one"
    assert new.loc["5d WT/" + ONE[0], "genotype"] == "WT"


# ---------------------------------------------------------------- seam 5: caches and review codes across a rename


def test_caches_and_review_codes_survive_a_session_folder_rename(tmp_path, monkeypatch):
    from tumour_synth import SynthTumour, write_ome

    from kerncount import cli
    from kerncount.tumour import segment

    raw = tmp_path / "data" / "tumour_exp"
    for i, name in enumerate(["030826 HERV -clodro n1 A.ome.tiff", "030826 cntrl -clodro n1 B.ome.tiff"]):
        t = SynthTumour(semi_axes_um=(15, 35, 30), voxel_um=(2.36, 1.25, 1.25), shape=(24, 96, 96), plateau=90,
                        seed=i + 1)
        write_ome(raw / "ready 5d KO7 030826" / name, t.stack(), t.voxel_um)
    monkeypatch.chdir(tmp_path)
    prj = tmp_path / "analysis" / "tumour_exp"
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    cli.main(["run", "tumour_exp"])
    before = pd.read_csv(prj / "tumours.csv")
    key = (prj / "review_key.csv").read_text()
    work = {p.name: p.stat().st_mtime_ns for p in (prj / "_work" / "tumour").iterdir()}

    (raw / "ready 5d KO7 030826").rename(raw / "5d K07")

    def boom(*a, **k):
        raise AssertionError("re-segmented after a rename")
    monkeypatch.setattr(segment, "measure_file", boom)
    monkeypatch.setattr("kerncount.tumour.measure.measure_file", boom)
    monkeypatch.setattr(segment, "read_stack", boom)
    cli.main(["run", "tumour_exp"])  # before re-init: rows matched by file name
    cli.main(["init", "tumour_exp"])
    cli.main(["run", "tumour_exp"])  # after re-init
    after = pd.read_csv(prj / "tumours.csv")
    assert (prj / "review_key.csv").read_text() == key
    assert {p.name: p.stat().st_mtime_ns for p in (prj / "_work" / "tumour").iterdir()
            if p.name in work} == work  # no cache file rewritten
    assert list(after["tumour_id"]) == list(before["tumour_id"])
    assert (after["imaging_session"] == "5d K07").all()
    pd.testing.assert_series_equal(after["cross_sectional_area_mm2"], before["cross_sectional_area_mm2"])


# ---------------------------------------------------------------- seam 6: 1 / 2 / 3 genotype levels


GEN3 = {"name": "genotype", "levels": ["WT", "TLR7 KO", "TLR9 KO"], "reference": "WT", "role": "row",
        "source": "session"}


def _three_settings(**tumour) -> dict:
    from test_tumour_stats import _settings as stats_settings

    from kerncount.tumour.config import DEFAULTS

    facs = [GEN3, *DEFAULTS["factors"][1:]]
    return stats_settings(factors=facs, **tumour)


def _vol(genotypes, seed=3) -> pd.DataFrame:
    from test_tumour_stats import _vol_table

    return _vol_table(np.random.default_rng(seed), genotypes=genotypes)


def test_f4_compares_every_non_reference_genotype_with_the_reference(tmp_path):
    from kerncount.tumour import config, stats
    from kerncount.tumour.families import contrasts

    s = _three_settings()
    facs = config.factors(s)
    f4 = next(f for f in config.families(s) if f.id == "F4")
    cs = contrasts(facs, f4)
    assert len(cs) == 2 * 2 * 2  # (TLR7 KO, TLR9 KO) x (HERV, N4) x clodronate
    assert {c.level.split(", ")[0] for c in cs} == {"TLR7 KO", "TLR9 KO"}
    ro = config.Readout("v", "main_mass_volume_mm3", "v", "confirmatory")
    res = stats.readout_results(_vol(("WT", "TLR7 KO")), facs, [f4], ro)
    tlr9 = res[res["comparison"].str.contains("TLR9 KO")]
    tlr7 = res[~res["comparison"].str.contains("TLR9 KO")]
    assert len(tlr9) == 4 and tlr9["note"].str.contains("needs TLR9 KO data").all()
    assert tlr7["p"].notna().all()  # within-session ratio of ratios: tested


@pytest.mark.parametrize("genotypes", [("TLR7 KO",), ("WT", "TLR7 KO"), ("WT", "TLR7 KO", "TLR9 KO")])
def test_analysis_and_plots_run_with_one_two_and_three_genotypes(tmp_path, monkeypatch, genotypes):
    import matplotlib.pyplot as plt

    from kerncount.experiment_types import RunOptions, get

    shapes = []
    real = plt.subplots

    def spy(*a, **k):
        fig, ax = real(*a, **k)
        shapes.append(np.shape(ax))
        return fig, ax
    monkeypatch.setattr(plt, "subplots", spy)
    t = _vol(genotypes)
    get("tumour_volume").analyse(tmp_path, t, _three_settings(), RunOptions())
    st = pd.read_csv(tmp_path / "stats_long.csv")
    absent = [g for g in ("WT", "TLR7 KO", "TLR9 KO") if g not in genotypes]
    for g in absent:
        assert st["note"].str.contains(f"needs {g} data").any(), g
    assert (len(genotypes), 2) in shapes  # group plots: one row per present genotype, one panel per clodronate
    assert (tmp_path / "plots" / "cross_sectional_area_interaction.png").exists() or \
        (tmp_path / "plots" / "main_mass_volume_interaction.png").exists()


# ---------------------------------------------------------------- seam 7: batch setting -> confounding


def _batch_layout() -> tuple[pd.DataFrame, pd.DataFrame]:
    """KO only, one session folder; HERV imaged on another day than cntrl (as in the real KO set)."""
    rng = np.random.default_rng(5)
    rows, meas = [], []
    for line, day in (("cntrl", "2026-08-15"), ("HERV", "2026-08-17")):
        for c in ("-clodro", "+clodro"):
            for k in range(1, 5):
                tid = f"TLR7 KO / {line} / {c} / n{k}A"
                rows.append(dict(tumour_id=tid, imaging_session="5d K07", session_label="5d", file_date="2026-08-03",
                                 acquisition_day=day, genotype="TLR7 KO", tumour_line=line, clodronate=c,
                                 section=f"n{k}", hemisphere="A", tumour=f"n{k}A", file=f"5d K07/{tid}.ome.tiff",
                                 legacy_mask="", legacy_stats="", include="yes", note=""))
                meas.append(dict(tumour_id=tid, main_mass_volume_um3=float(np.exp(rng.normal(20, 0.3)))))
    return pd.DataFrame(rows), pd.DataFrame(meas)


@pytest.mark.parametrize("batch, f1_tested, label", [
    ("session_folder", True, ""),
    ("acquisition_day", False, "confounded with acquisition day"),
])
def test_batch_setting_decides_which_contrasts_are_confounded(tmp_path, batch, f1_tested, label):
    from test_tumour_stats import _settings as stats_settings

    from kerncount.tumour import config, planner, stats, table

    s = stats_settings(batch=batch)
    assert config.batch(s) == batch
    lay_, meas = _batch_layout()
    t = table.build(lay_, meas, s)
    expect = t["imaging_session"] if batch == "session_folder" else t["acquisition_day"]
    assert list(t["batch"]) == list(expect)
    facs = config.factors(s)
    fams = [f for f in config.families(s) if f.id in ("F1", "F2")]
    ro = config.Readout("v", "main_mass_volume_mm3", "v", "confirmatory")
    res = stats.readout_results(t, facs, fams, ro).set_index("comparison")
    herv = res.loc["HERV vs cntrl | TLR7 KO, -clodro"]
    assert bool(np.isfinite(herv["p"])) == f1_tested
    if label:
        assert label in herv["note"]
    assert np.isfinite(res.loc["+clodro vs -clodro | TLR7 KO, HERV", "p"])  # same day: tested either way
    plan = planner.plan_readout(stats.fit(t, facs, ro.column), facs, fams, ro, batch_kind=batch)
    assert any(plan["comparison"].str.startswith("HERV vs cntrl | TLR7 KO")) == f1_tested


def test_batch_setting_is_validated_and_outside_the_lock_hash_only_for_session_folder():
    """Ticket 13 (user decision 2026-11-20) changed the default to culture_batch; session_folder stays outside the
    hash so that locks made before ticket 12 remain valid when it is chosen explicitly."""
    from test_tumour_stats import _settings as stats_settings

    from kerncount.tumour import config, results

    assert config.batch(stats_settings()) == "culture_batch"
    with pytest.raises(ValueError, match="batch"):
        config.batch(stats_settings(batch="weekday"))
    session = results.lock(stats_settings(batch="session_folder"))
    assert "batch" not in session[1]
    assert results.lock(stats_settings())[0] != session[0]
    assert results.lock(stats_settings())[0] != results.lock(stats_settings(batch="acquisition_day"))[0]


def test_report_states_the_batch_definition(tmp_path):
    from test_tumour_stats import _settings as stats_settings

    from kerncount.experiment_types import RunOptions, get

    lay_, meas = _batch_layout()
    from kerncount.tumour import table

    s = stats_settings(batch="acquisition_day")
    t = table.build(lay_, meas, s)
    get("tumour_volume").analyse(tmp_path, t, s, RunOptions())
    md = (tmp_path / "stats_report.md").read_text()
    assert "Batch definition: acquisition_day" in md
    assert "Batches: 2026-08-15, 2026-08-17" in md


# ---------------------------------------------------------------- real data (slow)

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7"


@pytest.mark.slow
@pytest.mark.skipif(not (REAL / "5d K07").is_dir(), reason="real tumour data not present")
def test_real_5d_k07_folder_gives_48_ko_tumours_imaged_on_four_days(tmp_path):
    from kerncount.tumour import layout as lay

    s = _settings(tmp_path)
    layout = lay.scan(REAL, s)
    ko = layout[layout["imaging_session"] == "5d K07"]
    assert len(ko) == 48 and set(ko["genotype"]) == {"TLR7 KO"} and set(ko["session_label"]) == {"5d"}
    assert set(ko["acquisition_day"]) == {"2026-08-15", "2026-08-17", "2026-08-18", "2026-08-19"}
    assert set(ko["file_date"]) == {"2026-08-03"}
    assert len(lay.check(ko, None, s).included) == 48
