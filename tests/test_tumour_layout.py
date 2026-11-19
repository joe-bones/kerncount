"""Layout validation for tumour experiments (ticket 03): filename variants, unreadable names, duplicate tumour IDs,
undeterminable session factors, per-row overrides in layout.csv and exclusions with a reason."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from kerncount import experiment_types, settings as settings_mod

KO_SESSION = "ready 5d KO7 030826"

# The 144 files of data/tumour_TLR7/ready 5d KO7 030826/ as delivered (listed 2026-11-19).
REAL_NAMES = [
    '030826 HERV +clodro n1 A.ome.tiff',
    '030826 HERV +clodro n1 Af.tif',
    '030826 HERV +clodro n1 B.ome.tiff',
    '030826 HERV +clodro n1 Bf.tif',
    '030826 HERV +clodro n2 A.ome.tiff',
    '030826 HERV +clodro n2 Af.tif',
    '030826 HERV +clodro n2 B.ome.tiff',
    '030826 HERV +clodro n2 Bf.tif',
    '030826 HERV +clodro n3 A.ome.tiff',
    '030826 HERV +clodro n3 Af.tif',
    '030826 HERV +clodro n3 B.ome.tiff',
    '030826 HERV +clodro n3 Bf.tif',
    '030826 HERV -clodro n1 A.ome.tiff',
    '030826 HERV -clodro n1 Af.tif',
    '030826 HERV -clodro n1 B.ome.tiff',
    '030826 HERV -clodro n1 Bf.tif',
    '030826 HERV -clodro n2 A.ome.tiff',
    '030826 HERV -clodro n2 Af.tif',
    '030826 HERV -clodro n2 B.ome.tiff',
    '030826 HERV -clodro n2 Bf.tif',
    '030826 HERV -clodro n3 A.ome.tiff',
    '030826 HERV -clodro n3 Af.tif',
    '030826 HERV -clodro n3 B.ome.tiff',
    '030826 HERV -clodro n3 Bf.tif',
    '030826 HERV -clodro n4 A.ome.tiff',
    '030826 HERV -clodro n4 Af.tif',
    '030826 HERV -clodro n4 B.ome.tiff',
    '030826 HERV -clodro n4 B.txt',
    '030826 HERV -clodro n4 Bf.tif',
    '030826 HERV -clodro n5 A.ome.tiff',
    '030826 HERV -clodro n5 Af.tif',
    '030826 HERV -clodro n5 B.ome.tiff',
    '030826 HERV -clodro n5 Bf.tif',
    '030826 HERV -clodro n5 Bf_3Dresults.csv',
    '030826 HERV -clodro n6 A.ome.tiff',
    '030826 HERV -clodro n6 Af.tif',
    '030826 HERV -clodro n6 Af_3Dresults.csv',
    '030826 HERV -clodron n6 B.ome.tiff',
    '030826 HERV -clodron n6 Bf.tif',
    '030826 N4 +clodro n1 A.ome.tiff',
    '030826 N4 +clodro n1 Af.tif',
    '030826 N4 +clodro n1 B.ome.tiff',
    '030826 N4 +clodro n1 Bf.tif',
    '030826 N4 +clodro n2 A.ome.tiff',
    '030826 N4 +clodro n2 Af.tif',
    '030826 N4 +clodro n2 B.ome.tiff',
    '030826 N4 +clodro n2 Bf.tif',
    '030826 N4 +clodro n3 A.ome.tiff',
    '030826 N4 +clodro n3 Af.tif',
    '030826 N4 +clodro n3 B.ome.tiff',
    '030826 N4 +clodro n3 Bf.tif',
    '030826 N4 -clodro n1 A.ome.tiff',
    '030826 N4 -clodro n1 Af.tif',
    '030826 N4 -clodro n1 B.ome.tiff',
    '030826 N4 -clodro n1 Bf.tif',
    '030826 N4 -clodro n2 A.ome.tiff',
    '030826 N4 -clodro n2 Af.tif',
    '030826 N4 -clodro n2 B.ome.tiff',
    '030826 N4 -clodro n2 Bf.tif',
    '030826 N4 -clodro n3 A.ome.tiff',
    '030826 N4 -clodro n3 Af.tif',
    '030826 N4 -clodro n3 B.ome.tiff',
    '030826 N4 -clodro n3 Bf.tif',
    '030826 cntrl +clodro n1 A.ome.tiff',
    '030826 cntrl +clodro n1 Af.tif',
    '030826 cntrl +clodro n1 B.ome.tiff',
    '030826 cntrl +clodro n1 Bf.tif',
    '030826 cntrl +clodro n2 A.ome.tiff',
    '030826 cntrl +clodro n2 Af.tif',
    '030826 cntrl +clodro n2 B.ome.tiff',
    '030826 cntrl +clodro n2 Bf.tif',
    '030826 cntrl +clodro n3 A.ome.tiff',
    '030826 cntrl +clodro n3 Af.tif',
    '030826 cntrl +clodro n3 B.ome.tiff',
    '030826 cntrl +clodro n3 Bf.tif',
    '030826 cntrl -clodro n1 A.ome.tiff',
    '030826 cntrl -clodro n1 Af.tif',
    '030826 cntrl -clodro n1 B.ome.tiff',
    '030826 cntrl -clodro n1 Bf.tif',
    '030826 cntrl -clodro n2 A.ome.tiff',
    '030826 cntrl -clodro n2 Af.tif',
    '030826 cntrl -clodro n2 B.ome.tiff',
    '030826 cntrl -clodro n2 Bf.tif',
    '030826 cntrl -clodro n3 A.ome.tiff',
    '030826 cntrl -clodro n3 Af.tif',
    '030826 cntrl -clodro n3 B.ome.tiff',
    '030826 cntrl -clodro n3 Bf.tif',
    '030826 cntrl -clodro n4 A.ome.tiff',
    '030826 cntrl -clodro n4 Af.tif',
    '030826 cntrl -clodro n4 B.ome.tiff',
    '030826 cntrl -clodro n4 Bf.tif',
    '030826 cntrl -clodro n5 A.ome.tiff',
    '030826 cntrl -clodro n5 Af.tif',
    '030826 cntrl -clodro n5 B.ome.tiff',
    '030826 cntrl -clodro n5 Bf.tif',
    '030826 cntrl -clodro n6 A.ome.tiff',
    '030826 cntrl -clodro n6 Af.tif',
    '030826 cntrl -clodro n6 B.ome.tiff',
    '030826 cntrl -clodro n6 Bf.tif',
    'Statistics for 030826 HERV +clodro n1 Af.csv',
    'Statistics for 030826 HERV +clodro n1 Bf.csv',
    'Statistics for 030826 HERV +clodro n2 Af.csv',
    'Statistics for 030826 HERV +clodro n2 Bf.csv',
    'Statistics for 030826 HERV +clodro n3 Af.csv',
    'Statistics for 030826 HERV +clodro n3 Bf.csv',
    'Statistics for 030826 HERV -clodro n1 Af.csv',
    'Statistics for 030826 HERV -clodro n1 Bf.csv',
    'Statistics for 030826 HERV -clodro n2 Af.csv',
    'Statistics for 030826 HERV -clodro n2 Bf.csv',
    'Statistics for 030826 HERV -clodro n3 Af.csv',
    'Statistics for 030826 HERV -clodro n3 Bf.csv',
    'Statistics for 030826 HERV -clodro n4 Af.csv',
    'Statistics for 030826 HERV -clodro n5 Af.csv',
    'Statistics for 030826 HERV -clodron n6 Bf.csv',
    'Statistics for 030826 N4 +clodro n1 Af.csv',
    'Statistics for 030826 N4 +clodro n1 Bf.csv',
    'Statistics for 030826 N4 +clodro n2 Af.csv',
    'Statistics for 030826 N4 +clodro n2 Bf.csv',
    'Statistics for 030826 N4 +clodro n3 Af.csv',
    'Statistics for 030826 N4 +clodro n3 Bf.csv',
    'Statistics for 030826 N4 -clodro n1 Af.csv',
    'Statistics for 030826 N4 -clodro n1 Bf.csv',
    'Statistics for 030826 N4 -clodro n2 Af.csv',
    'Statistics for 030826 N4 -clodro n2 Bf.csv',
    'Statistics for 030826 N4 -clodro n3 Af.csv',
    'Statistics for 030826 N4 -clodro n3 Bf.csv',
    'Statistics for 030826 cntrl +clodro n1 Af.csv',
    'Statistics for 030826 cntrl +clodro n1 Bf.csv',
    'Statistics for 030826 cntrl +clodro n2 Af.csv',
    'Statistics for 030826 cntrl +clodro n2 Bf.csv',
    'Statistics for 030826 cntrl +clodro n3 Af.csv',
    'Statistics for 030826 cntrl +clodro n3 Bf.csv',
    'Statistics for 030826 cntrl -clodro n1 Af.csv',
    'Statistics for 030826 cntrl -clodro n1 Bf.csv',
    'Statistics for 030826 cntrl -clodro n2 Af.csv',
    'Statistics for 030826 cntrl -clodro n2 Bf.csv',
    'Statistics for 030826 cntrl -clodro n3 Af.csv',
    'Statistics for 030826 cntrl -clodro n3 Bf.csv',
    'Statistics for 030826 cntrl -clodro n4 Af.csv',
    'Statistics for 030826 cntrl -clodro n4 Bf.csv',
    'Statistics for 030826 cntrl -clodro n5 Af.csv',
    'Statistics for 030826 cntrl -clodro n5 Bf.csv',
    'Statistics for 030826 cntrl -clodro n6 Af.csv',
    'Statistics for 030826 cntrl -clodro n6 Bf.csv',
]


def _settings(tmp_path: Path) -> dict:
    (tmp_path / "settings.toml").write_text('[experiment]\ntype = "tumour_volume"\n')
    return settings_mod.load(tmp_path)


def _touch(folder: Path, *names: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for n in names:
        (folder / n).write_bytes(b"")


def _scan(tmp_path: Path, files: dict[str, list[str]]) -> tuple[Path, dict, pd.DataFrame]:
    raw = tmp_path / "raw"
    for session, names in files.items():
        _touch(raw / session, *names)
    s = _settings(tmp_path)
    return raw, s, experiment_types.get("tumour_volume").scan_layout(raw, s)


GROUPS = {("cntrl", "-clodro"): 12, ("cntrl", "+clodro"): 6, ("HERV", "-clodro"): 12, ("HERV", "+clodro"): 6,
          ("N4", "-clodro"): 6, ("N4", "+clodro"): 6}


def _assert_real_layout(layout: pd.DataFrame) -> None:
    assert len(layout) == 48
    assert (layout["include"] == "yes").all()
    assert layout["tumour_id"].is_unique
    assert set(layout["genotype"]) == {"TLR7 KO"}
    assert layout.groupby(["tumour_line", "clodronate"]).size().to_dict() == GROUPS
    noted = layout[layout["note"] != ""]
    assert list(noted["tumour_id"]) == ["TLR7 KO / HERV / -clodro / n6B"]  # the one -clodron typo
    assert "-clodron" in noted["note"].iloc[0]
    assert not layout["file"].str.contains(r"f\.tif$|Statistics|\.txt$|_3Dresults").any()
    assert (layout["legacy_mask"] != "").all()  # every tumour has its Fiji mask ...
    assert (layout["legacy_stats"] != "").all()  # ... and statistics, incl. the header-less .txt and _3Dresults.csv
    row = layout.set_index("tumour_id")
    assert row.loc["TLR7 KO / HERV / -clodro / n4B", "legacy_stats"].endswith("030826 HERV -clodro n4 B.txt")
    assert row.loc["TLR7 KO / HERV / -clodro / n6A", "legacy_stats"].endswith("n6 Af_3Dresults.csv")


def test_real_filename_list_gives_48_tumours_in_the_expected_groups(tmp_path):
    _, _, layout = _scan(tmp_path, {KO_SESSION: REAL_NAMES})
    _assert_real_layout(layout)


def test_known_naming_variants_are_normalised_and_noted(tmp_path):
    _, _, layout = _scan(tmp_path, {KO_SESSION: [
        "030826 HERV -clodro n1 A.ome.tiff.tif",  # original re-saved by Fiji
        "030826 herv -clodro n2 A.ome.tiff",  # case
        "030826 HERV -clodro n3A.ome.tiff",  # section + hemisphere glued
        "030826 HERV - clodro n4 A.ome.tiff",  # sign split from its word
        "030826_HERV_-clodro_n5_A.OME.TIFF",  # underscores, upper-case ending
        "030826 HERV −clodro n6  A.ome.tiff",  # typographic minus, double space
    ]})
    by_file = layout.set_index(layout["file"].str.split("/").str[-1])
    expect = {
        "030826 HERV -clodro n1 A.ome.tiff.tif": (1, ".ome.tiff.tif"),
        "030826 herv -clodro n2 A.ome.tiff": (2, "'herv' read as 'HERV'"),
        "030826 HERV -clodro n3A.ome.tiff": (3, "'n3A' read as section n3 + hemisphere A"),
        "030826 HERV - clodro n4 A.ome.tiff": (4, "sign joined"),
        "030826_HERV_-clodro_n5_A.OME.TIFF": (5, "spacing"),
        "030826 HERV −clodro n6  A.ome.tiff": (6, "dash"),
    }
    for name, (k, note) in expect.items():
        row = by_file.loc[name]
        assert row["tumour_id"] == f"TLR7 KO / HERV / -clodro / n{k}A", name
        assert row["include"] == "yes" and note in row["note"], (name, row["note"])
    assert layout["date"].eq("2026-08-03").all()


def test_unreadable_names_are_flagged_and_block_the_run_until_fixed_or_excluded(tmp_path):
    names = ["030826 HERV -clodro n1 A.ome.tiff", "030826 HREV -clodro n1 B.ome.tiff", "030826 HERV n2 A.ome.tiff",
             "030826 HERV -clodro A.ome.tiff"]
    raw, s, layout = _scan(tmp_path, {KO_SESSION: names})
    from kerncount.tumour import layout as lay

    flagged = layout[layout["note"].str.startswith("could not read")]
    assert len(flagged) == 3 and (flagged["tumour_id"] == "").all()
    note = dict(zip(flagged["file"].str.split("/").str[-1], flagged["note"]))
    assert "tumour_line" in note[names[1]] and "clodronate" in note[names[2]] and "section" in note[names[3]]

    with pytest.raises(lay.LayoutBlocked) as e:
        lay.check(layout, raw, s)
    assert len(e.value.problems) == 3  # every problem at once, not the first only
    for n in names[1:]:
        assert n in str(e.value.code)

    fixed = layout.copy()
    fixed.loc[fixed["file"].str.endswith("HREV -clodro n1 B.ome.tiff"), "tumour_line"] = "HERV"
    fixed.loc[fixed["file"].str.endswith("HERV n2 A.ome.tiff"), ["include", "note"]] = ["no", "no clodronate label"]
    fixed.loc[fixed["file"].str.endswith("-clodro A.ome.tiff"), "section"] = "n3"
    v = lay.check(fixed, raw, s)
    assert sorted(v.included["tumour_id"]) == ["TLR7 KO / HERV / -clodro / n1A", "TLR7 KO / HERV / -clodro / n1B",
                                            "TLR7 KO / HERV / -clodro / n3A"]
    hrev = v.included.set_index("tumour_id").loc["TLR7 KO / HERV / -clodro / n1B", "note"]
    assert "tumour_line set in layout.csv to 'HERV' (not readable from the file name)" in hrev
    assert list(v.excluded["exclusion_reason"]) == ["no clodronate label"]


def test_duplicate_tumour_ids_block_the_run_naming_the_conflicting_files(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {KO_SESSION: [
        "030826 HERV -clodro n6 B.ome.tiff", "030826 HERV -clodron n6 B.ome.tiff",  # typo copy -> same ID
        "030826 N4 +clodro n1 A.ome.tiff", "030826 N4 +clodro n1 A.ome.tiff.tif",
        "030826 cntrl -clodro n1 A.ome.tiff"]})
    with pytest.raises(lay.LayoutBlocked) as e:
        lay.check(layout, raw, s)
    assert len(e.value.problems) == 2
    msg = str(e.value.code)
    for tid, files in {"TLR7 KO / HERV / -clodro / n6B": ("HERV -clodro n6 B.ome.tiff", "HERV -clodron n6 B.ome.tiff"),
                       "TLR7 KO / N4 / +clodro / n1A": ("n1 A.ome.tiff", "n1 A.ome.tiff.tif")}.items():
        line = next(p for p in e.value.problems if tid in p)
        assert all(f in line for f in files), line
        assert tid in msg

    layout.loc[layout["file"].str.endswith("-clodron n6 B.ome.tiff"), ["include", "note"]] = ["no", "duplicate copy"]
    layout.loc[layout["file"].str.endswith(".tiff.tif"), "include"] = "no"
    v = lay.validate(layout, raw, s)
    assert v.problems == [] and len(v.included) == 3
    assert sorted(v.excluded["exclusion_reason"]) == ["duplicate copy", "excluded in layout.csv (no reason given)"]
    assert any("without a reason" in w for w in v.warnings)


def test_a_session_folder_without_a_readable_genotype_blocks_until_set_in_layout(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {KO_SESSION: ["030826 HERV -clodro n1 A.ome.tiff"],
                                      "ready 5d 120926": ["120926 HERV -clodro n1 A.ome.tiff",
                                                          "120926 HERV -clodro n1 B.ome.tiff"]})
    unknown = layout[layout["imaging_session"] == "ready 5d 120926"]
    assert (unknown["genotype"] == "").all() and unknown["note"].str.contains("genotype").all()
    with pytest.raises(lay.LayoutBlocked) as e:
        lay.check(layout, raw, s)
    assert len(e.value.problems) == 1  # one problem per session, not per tumour
    assert "ready 5d 120926" in e.value.problems[0] and "genotype" in e.value.problems[0]

    layout.loc[layout["imaging_session"] == "ready 5d 120926", "genotype"] = "wt"  # case-insensitive level
    v = lay.check(layout, raw, s)
    wt = v.included[v.included["imaging_session"] == "ready 5d 120926"]
    assert list(wt["tumour_id"]) == ["WT / HERV / -clodro / n1A", "WT / HERV / -clodro / n1B"]
    assert wt["note"].str.contains("genotype set in layout.csv to 'WT'").all()


def test_factor_values_can_be_overridden_per_row_with_the_override_noted(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _scan(tmp_path, {KO_SESSION: ["030826 HERV -clodro n1 A.ome.tiff",
                                                   "030826 HERV -clodro n1 B.ome.tiff"]})
    a = layout["file"].str.endswith("n1 A.ome.tiff")
    layout.loc[a, ["genotype", "tumour_line", "hemisphere"]] = ["WT", "N4", "b"]  # stale tumour_id left as is
    v = lay.check(layout, raw, s)
    row = v.included[v.included["file"].str.endswith("n1 A.ome.tiff")].iloc[0]
    assert row["tumour_id"] == "WT / N4 / -clodro / n1B" and row["tumour"] == "n1B"
    assert "genotype set in layout.csv to 'WT' (from the session folder: 'TLR7 KO')" in row["note"]
    assert "tumour_line set in layout.csv to 'N4' (from the file name: 'HERV')" in row["note"]
    assert "hemisphere set in layout.csv to 'B'" in row["note"]
    untouched = v.included[v.included["file"].str.endswith("n1 B.ome.tiff")].iloc[0]
    assert untouched["note"] == ""

    layout.loc[a, "tumour_line"] = "N5"
    with pytest.raises(lay.LayoutBlocked, match="'N5' is not one of cntrl, HERV, N4"):
        lay.check(layout, raw, s)


# ---------------------------------------------------------------- kerncount init / run


def _workspace(tmp_path, monkeypatch, names: list[str]) -> tuple[Path, Path]:
    raw = tmp_path / "data" / "tumour_exp"
    _touch(raw / KO_SESSION, *names)  # empty files: a run that got past validation would fail to read them
    monkeypatch.chdir(tmp_path)
    return raw, tmp_path / "analysis" / "tumour_exp"


def test_run_refuses_to_start_on_a_blocked_layout_and_lists_every_problem(tmp_path, monkeypatch, capsys):
    from kerncount import cli
    from kerncount.tumour import layout as lay

    _, prj = _workspace(tmp_path, monkeypatch, ["030826 HERV -clodro n6 B.ome.tiff",
                                                "030826 HERV -clodron n6 B.ome.tiff", "030826 HREV +clodro n1 A.ome.tiff"])
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    out = capsys.readouterr().out
    assert "2 problems" in out and "HREV" in out and "-clodron n6 B" in out  # shown at init already
    assert "include = no" in out and "n / treatment" not in out  # tumour wording of the final hint

    with pytest.raises(lay.LayoutBlocked) as e:
        cli.main(["run", "tumour_exp"])
    assert len(e.value.problems) == 2
    assert not (prj / "_work" / "tumour").exists() or not any((prj / "_work" / "tumour").iterdir())
    assert not (prj / "tumours.csv").exists()


def test_reinit_keeps_an_edited_layout_csv(tmp_path, monkeypatch):
    from kerncount import cli, discover

    _, prj = _workspace(tmp_path, monkeypatch, ["030826 HERV -clodro n1 A.ome.tiff"])
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    edited = discover.read_table(prj / "layout.csv")
    edited.loc[0, ["include", "note"]] = ["no", "out of focus"]
    discover.write_layout(edited, prj / "layout.csv")

    cli.main(["init", "tumour_exp"])
    assert discover.read_table(prj / "layout.csv").loc[0, "note"] == "out of focus"
    assert discover.read_table(prj / "layout.new.csv").loc[0, "include"] == "yes"
    cli.main(["init", "tumour_exp", "--force-layout"])
    assert discover.read_table(prj / "layout.csv").loc[0, "include"] == "yes"
    assert list((prj / "_work" / "backups").glob("layout*.csv"))  # the edited version is backed up


def test_excluded_tumours_and_their_reasons_are_written_with_the_run(tmp_path, monkeypatch):
    from tumour_synth import SynthTumour, write_ome

    from kerncount import cli, discover

    _, prj = _workspace(tmp_path, monkeypatch, ["030826 HERV -clodro n1 B.ome.tiff"])
    t = SynthTumour(semi_axes_um=(15, 35, 30), voxel_um=(2.36, 1.25, 1.25), shape=(24, 128, 128), plateau=90)
    write_ome(tmp_path / "data" / "tumour_exp" / KO_SESSION / "030826 HERV -clodro n1 A.ome.tiff", t.stack(),
              t.voxel_um)
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    layout = discover.read_table(prj / "layout.csv")
    layout.loc[layout["file"].str.endswith("n1 B.ome.tiff"), ["include", "note"]] = ["no", "stack corrupted"]
    discover.write_layout(layout, prj / "layout.csv")

    cli.main(["run", "tumour_exp", "--no-review"])
    assert list(pd.read_csv(prj / "tumours.csv")["tumour_id"]) == ["TLR7 KO / HERV / -clodro / n1A"]
    ex = pd.read_csv(prj / "excluded_tumours.csv")
    assert list(ex["tumour_id"]) == ["TLR7 KO / HERV / -clodro / n1B"]
    assert list(ex["exclusion_reason"]) == ["stack corrupted"]
    assert list(ex["tumour_line"]) == ["HERV"]


# ---------------------------------------------------------------- real folder (slow)

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7"


@pytest.mark.slow
@pytest.mark.skipif(not (REAL / KO_SESSION).is_dir(), reason="real tumour data not present")
def test_real_folder_layout_is_valid_with_48_tumours(tmp_path):
    from kerncount.tumour import layout as lay

    s = _settings(tmp_path)
    layout = lay.scan(REAL, s)
    layout = layout[layout["imaging_session"] == KO_SESSION].reset_index(drop=True)  # other sessions may be added
    _assert_real_layout(layout)
    assert lay.unmatched_files(REAL / KO_SESSION, layout.assign(
        file=layout["file"].str.split("/").str[-1], legacy_mask=layout["legacy_mask"].str.split("/").str[-1],
        legacy_stats=layout["legacy_stats"].str.split("/").str[-1])) == []
    assert sorted(os.listdir(REAL / KO_SESSION)) == sorted(REAL_NAMES)  # the hard-coded list is the real folder
    v = lay.check(layout, None, s)
    assert len(v.included) == 48
