"""Structure-agnostic ingestion (ticket 13): the folder structure below a session folder is irrelevant.

Seams: `layout.scan` of a nested tree (genotype from the top-level session folder only; stacks at any depth; .czi
twins as provenance, .czi-only tumours reported; `_s1`, `1n`, optional date prefix; culture batch from the file-name
date; IDs unique across culture batches), `layout.validate` with an unavailable (dangling-link) session folder, and the
statistics with culture batch as a fixed block effect (ADR 0003 amendment 2026-11-20)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import settings as settings_mod


def _settings(tmp_path: Path, extra: str = "") -> dict:
    (tmp_path / "settings.toml").write_text('[experiment]\ntype = "tumour_volume"\n' + extra)
    return settings_mod.load(tmp_path)


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


# A miniature of the real `5d WT` folder (2026-11-20): <DDMMYY> OBS/<DDMMYY> <line> <±clodro> [n<k>]/<files>
WT_TREE = [
    "190626 OBS/190626 cntrl -clodro n1/190626 cntrl -clodro n1 B.czi",
    "190626 OBS/190626 cntrl -clodro n1/190626 cntrl -clodro n1 B.ome.tiff",
    "190626 OBS/190626 cntrl +clodro n1/190626 cntrl +clodro 1n A.ome.tiff",
    "190626 OBS/190626 cntrl +clodro n1/cntrl +clodro n1 A.czi",  # czi without date prefix, section written n1
    "190626 OBS/190626 HERV -clodro n1/190626 HERV -clodro n1 A.czi",
    "190626 OBS/190626 HERV -clodro n1/190626 HERV -clodro n1 A.ome.tiff",
    "190626 OBS/190626 N4 +clodro n3/190626 N4 +clodro n3 B.czi",  # czi only: reported, not measured
    "100726 OBS/100726 cntrl -clodro/100726 cntrl -clodro n1 A_s1.ome.tiff",
    "100726 OBS/100726 cntrl -clodro/100726 cntrl -clodro n1 B.ome.tiff",  # same ID as 190626 n1 B -> batch suffix
    "100726 OBS/100726 HERV -clodro/100726 HERV -clodro n1 A_s1.ome.tiff",
]


def _wt(tmp_path: Path, extra_sessions: dict[str, list[str]] | None = None) -> tuple[Path, dict, pd.DataFrame]:
    from kerncount.tumour import layout as lay

    raw = tmp_path / "raw"
    for rel in WT_TREE:
        _touch(raw / "5d WT" / rel)
    for session, names in (extra_sessions or {}).items():
        for n in names:
            _touch(raw / session / n)
    s = _settings(tmp_path)
    return raw, s, lay.scan(raw, s)


def test_a_nested_session_folder_gives_one_row_per_stack_with_genotype_from_the_top_folder(tmp_path):
    raw, s, layout = _wt(tmp_path, {"5d K07": ["030826 HERV -clodro n1 A.ome.tiff"]})
    wt = layout[layout["imaging_session"] == "5d WT"]
    assert len(wt) == 6 and set(wt["genotype"]) == {"WT"}  # 6 .ome.tiff; .czi twins are no extra rows
    assert all(f.startswith("5d WT/") and f.endswith(".ome.tiff") for f in wt["file"])
    ids = dict(zip(wt["file"], wt["tumour_id"]))
    assert ids["5d WT/190626 OBS/190626 cntrl +clodro n1/190626 cntrl +clodro 1n A.ome.tiff"] == \
        "WT / cntrl / +clodro / n1A"
    assert ids["5d WT/100726 OBS/100726 HERV -clodro/100726 HERV -clodro n1 A_s1.ome.tiff"] == \
        "WT / HERV / -clodro / n1A / 2026-07-10"
    assert ids["5d WT/190626 OBS/190626 HERV -clodro n1/190626 HERV -clodro n1 A.ome.tiff"] == \
        "WT / HERV / -clodro / n1A / 2026-06-19"
    ko = layout[layout["imaging_session"] == "5d K07"].iloc[0]
    assert ko["tumour_id"] == "TLR7 KO / HERV / -clodro / n1A" and ko["culture_batch"] == "2026-08-03"


def test_culture_batch_is_the_file_name_date_and_ids_stay_unique_across_batches(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _wt(tmp_path)
    by = dict(zip(layout["file"].str.rsplit("/", n=1).str[-1], layout["culture_batch"]))
    assert by["190626 cntrl -clodro n1 B.ome.tiff"] == "2026-06-19"
    assert by["100726 cntrl -clodro n1 B.ome.tiff"] == "2026-07-10"
    v = lay.check(layout, raw, s)
    assert v.included["tumour_id"].is_unique and len(v.included) == 6
    assert "WT / cntrl / -clodro / n1B / 2026-07-10" in set(v.included["tumour_id"])
    assert "WT / cntrl / -clodro / n1A" in set(v.included["tumour_id"])  # no clash: plain ID


def test_name_variants_are_normalised_and_noted(tmp_path):
    _, _, layout = _wt(tmp_path)
    note = dict(zip(layout["file"].str.rsplit("/", n=1).str[-1], layout["note"]))
    assert "scene suffix '_s1' ignored" in note["100726 cntrl -clodro n1 A_s1.ome.tiff"]
    assert "spacing" not in note["100726 cntrl -clodro n1 A_s1.ome.tiff"]
    assert "section label '1n' read as 'n1'" in note["190626 cntrl +clodro 1n A.ome.tiff"]
    assert note["190626 cntrl -clodro n1 B.ome.tiff"] == ""


def test_czi_twins_are_provenance_and_czi_only_tumours_are_reported(tmp_path):
    from kerncount.tumour import layout as lay

    raw, s, layout = _wt(tmp_path)
    czi = dict(zip(layout["file"].str.rsplit("/", n=1).str[-1], layout["czi"]))
    assert czi["190626 cntrl +clodro 1n A.ome.tiff"] == "5d WT/190626 OBS/190626 cntrl +clodro n1/cntrl +clodro n1 A.czi"
    assert czi["190626 HERV -clodro n1 A.ome.tiff"].endswith("190626 HERV -clodro n1 A.czi")
    assert czi["100726 cntrl -clodro n1 B.ome.tiff"] == ""  # 190626 n1 B's czi is not its twin (other batch)
    only = lay.czi_only(raw, layout)
    assert only == ["5d WT/190626 OBS/190626 N4 +clodro n3/190626 N4 +clodro n3 B.czi"]
    assert lay.unmatched_files(raw, layout) == []  # czi files are either twins or reported as czi-only
    v = lay.validate(layout, raw, s)
    assert any("190626 N4 +clodro n3 B.czi" in w and "CZI" in w for w in v.warnings)


def test_deeper_folders_are_never_read_for_factors(tmp_path):
    from kerncount.tumour import layout as lay

    raw = tmp_path / "raw"
    _touch(raw / "5d WT" / "KO stuff" / "HERV +clodro" / "190626 cntrl -clodro n1 A.ome.tiff")
    s = _settings(tmp_path)
    row = lay.scan(raw, s).iloc[0]
    assert (row["genotype"], row["tumour_line"], row["clodronate"]) == ("WT", "cntrl", "-clodro")


def test_a_symlinked_session_folder_is_scanned_and_an_unmounted_one_blocks_with_a_clear_message(tmp_path):
    from kerncount.tumour import layout as lay

    ext = tmp_path / "drive" / "5d WT"
    _touch(ext / "190626 OBS" / "x" / "190626 HERV -clodro n1 A.ome.tiff")
    raw = tmp_path / "raw"
    _touch(raw / "5d K07" / "030826 HERV -clodro n1 A.ome.tiff")
    os.symlink(ext, raw / "5d WT")
    s = _settings(tmp_path)
    layout = lay.scan(raw, s)
    assert set(layout["genotype"]) == {"WT", "TLR7 KO"}
    ext.rename(tmp_path / "drive" / "gone")  # drive unmounted: dangling link
    v = lay.validate(layout, raw, s)
    (p,) = [p for p in v.problems if "5d WT" in p]
    assert "not available" in p and str(ext) in p and "1 tumour" in p
    assert not any("not found in the raw folder" in p for p in v.problems)
    assert len(lay.scan(raw, s)) == 1  # scanning does not crash


@pytest.mark.slow
@pytest.mark.skipif(not (Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7" / "5d WT").is_dir(),
                    reason="WT data (external drive) not mounted")
def test_real_5d_wt_folder_gives_46_tumours_in_two_culture_batches(tmp_path):
    from kerncount.tumour import layout as lay

    raw = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7"
    s = _settings(tmp_path)
    layout = lay.scan(raw, s)
    wt = layout[layout["imaging_session"] == "5d WT"]
    assert len(wt) == 46 and set(wt["genotype"]) == {"WT"}
    assert wt["culture_batch"].value_counts().to_dict() == {"2026-06-19": 34, "2026-07-10": 12}
    v = lay.check(wt, None, s)
    assert len(v.included) == 46 and v.included["tumour_id"].is_unique
    groups = v.included.groupby(["tumour_line", "clodronate"]).size().to_dict()
    assert groups == {("cntrl", "-clodro"): 11, ("cntrl", "+clodro"): 6, ("HERV", "-clodro"): 12,
                      ("HERV", "+clodro"): 6, ("N4", "-clodro"): 6, ("N4", "+clodro"): 5}
    assert (wt["czi"] != "").sum() == 34 and lay.czi_only(raw, layout) == []


# ---------------------------------------------------------------- culture batch as a fixed block effect


def _batch_table(shift: float = 0.0, seed: int = 3) -> tuple[pd.DataFrame, dict]:
    """WT in two culture batches (190626: all six groups; 100726: cntrl and HERV -clodro only, every value x
    exp(shift)), KO in one (030826). True effects: N4 = 1.5 x cntrl, HERV = cntrl, +clodro = 0.8 x -clodro."""
    from kerncount.tumour import table

    rng = np.random.default_rng(seed)
    rows, meas = [], []
    plan = [("WT", "2026-06-19", [(ln, c) for ln in ("cntrl", "HERV", "N4") for c in ("-clodro", "+clodro")]),
            ("WT", "2026-07-10", [("cntrl", "-clodro"), ("HERV", "-clodro")]),
            ("TLR7 KO", "2026-08-03", [(ln, c) for ln in ("cntrl", "HERV", "N4") for c in ("-clodro", "+clodro")])]
    for g, batch, cells in plan:
        for line, c in cells:
            for k in range(1, 6):
                tid = f"{g} / {line} / {c} / n{k}A / {batch}"
                mu = np.log(1.5) * (line == "N4") + np.log(0.8) * (c == "+clodro") + shift * (batch == "2026-07-10")
                rows.append(dict(tumour_id=tid, imaging_session=f"5d {g}", session_label="5d", file_date=batch,
                                 culture_batch=batch, acquisition_day="", genotype=g, tumour_line=line, clodronate=c,
                                 section=f"n{k}", hemisphere="A", tumour=f"n{k}A", file=f"{tid}.ome.tiff",
                                 legacy_mask="", legacy_stats="", czi="", include="yes", note=""))
                meas.append(dict(tumour_id=tid, main_mass_volume_um3=float(np.exp(20 + mu + rng.normal(0, 0.2)))))
    s = settings_mod._merge(settings_mod.DEFAULTS, {"experiment": {"type": "tumour_volume"}})
    return table.build(pd.DataFrame(rows), pd.DataFrame(meas), s), s


def _results(t: pd.DataFrame, s: dict) -> pd.DataFrame:
    from kerncount.tumour import config, stats

    ro = config.Readout("v", "main_mass_volume_mm3", "v", "confirmatory")
    return stats.readout_results(t, config.factors(s), config.families(s), ro).set_index(["family", "comparison"])


def test_culture_batch_is_the_default_batch():
    from kerncount.tumour import config

    t, s = _batch_table()
    assert config.batch(s) == "culture_batch"
    assert list(t["batch"]) == list(t["culture_batch"])


def test_a_batch_shift_inside_one_genotype_leaves_within_genotype_contrasts_unchanged(tmp_path):
    base = _results(*_batch_table(0.0))
    shifted = _results(*_batch_table(np.log(2.0)))
    for fam in ("F1", "F2", "F3", "F4", "F5", "F7"):
        a, b = base.loc[fam], shifted.loc[fam]
        tested = a[np.isfinite(a["p"])]
        assert len(tested), fam
        for col in ("estimate", "ci_low", "ci_high", "p"):
            np.testing.assert_allclose(b.loc[tested.index, col], tested[col], rtol=1e-9, err_msg=f"{fam} {col}")
    wt_n4 = shifted.loc[("F1", "N4 vs cntrl | WT, -clodro")]
    assert 1.2 < wt_n4["estimate"] < 1.9 and np.isfinite(wt_n4["p"]) and not wt_n4["confounded"]
    f6 = shifted.loc["F6"]
    assert f6["confounded"].all() and f6["note"].str.contains("confounded with culture batch").all()
    assert not shifted.loc["F4", "confounded"].any()


def test_without_the_block_effect_the_planted_shift_would_bias_the_wt_line_effect():
    """Sanity check of the simulation: the raw ratio of group means is pulled by the shifted 100726 tumours."""
    t, _ = _batch_table(np.log(2.0))
    wt = t[(t["genotype"] == "WT") & (t["clodronate"] == "-clodro")]
    gm = wt.groupby("tumour_line")["main_mass_volume_mm3"].apply(lambda v: np.exp(np.log(v).mean()))
    assert gm["N4"] / gm["cntrl"] < 1.2  # true 1.5: half of the cntrl tumours are inflated x2


def test_sensitivity_models_and_normalised_values_use_batch_adjusted_values():
    from kerncount.tumour import config, sensitivity, stats
    from kerncount.tumour.families import contrasts

    t, s = _batch_table(np.log(2.0))
    facs, fams = config.factors(s), [f for f in config.families(s) if f.id == "F1"]
    ro = config.Readout("v", "main_mass_volume_mm3", "v", "confirmatory")
    model = stats.fit(t, facs, ro.column)
    adj = stats.batch_adjusted(t, model)
    res = stats.readout_results(t, facs, fams, ro, model)
    sens = sensitivity.readout_sensitivity(adj, facs, fams, ro, res, 199, 1).set_index("comparison")
    n4 = next(c for c in contrasts(facs, fams[0]) if c.comparison == "N4 vs cntrl | WT, -clodro")
    w = sensitivity.welch(sensitivity.groups(adj, facs, ro.column), n4)
    np.testing.assert_allclose(w.ratio, res.set_index("comparison").loc["N4 vs cntrl | WT, -clodro", "estimate"])
    assert np.isfinite(sens.loc["N4 vs cntrl | WT, -clodro", "welch_p"])
    pct = stats.normalised(t, model, facs[1])
    wt_cntrl = (t["genotype"] == "WT") & (t["tumour_line"] == "cntrl") & (t["clodronate"] == "-clodro")
    by_batch = pct[wt_cntrl].groupby(t.loc[wt_cntrl, "culture_batch"]).apply(lambda v: np.exp(np.log(v).mean()))
    assert all(70 < x < 140 for x in by_batch)  # both batches scatter around 100 % (the shift is removed)


# ---------------------------------------------------------------- acquisition report: differing settings, batches


def test_acquisition_report_names_the_groups_of_a_differing_setting_and_tabulates_culture_batches():
    from test_tumour_acquisition import _xml

    from kerncount.tumour import acquisition
    from kerncount.tumour.config import Factor

    rows, tum = [], []
    for i, (line, batch, volt) in enumerate([("cntrl", "2026-06-19", 597.0), ("cntrl", "2026-07-10", 554.0),
                                             ("HERV", "2026-07-10", 554.0), ("HERV", "2026-06-19", 597.0)]):
        rows.append({"tumour_id": f"t{i}", "file": f"f{i}", **acquisition.parse(_xml(volt=volt))})
        tum.append({"tumour_id": f"t{i}", "tumour_line": line, "culture_batch": batch, "voxel_x_um": 0.62})
        # (the readout table also carries measured columns named like acquisition settings, e.g. voxel_x_um)
    rows[1]["voxel_x_um"] = "1.2478516"  # one binned stack
    line = Factor("tumour_line", ("cntrl", "HERV"), "cntrl", "x")
    md = "\n".join(acquisition.report_lines(pd.DataFrame(rows), pd.DataFrame(tum), [line], "culture_batch"))
    assert "detector_voltage = 554.0: 2 tumours (cntrl 1, HERV 1; culture batch 2026-07-10)" in md
    assert "detector_voltage = 597.0: 2 tumours (cntrl 1, HERV 1; culture batch 2026-06-19)" in md
    assert "voxel_x_um = 1.2478516: 1 tumour (cntrl 1; culture batch 2026-07-10)" in md
    assert "Pixel size differs" in md
    assert "Tumours per group and culture batch" in md and "| cntrl | 1 | 1 |" in md
    assert "imaging day is recorded, never a batch" in md


def test_web_layout_view_and_results_stay_usable_when_the_drive_is_not_mounted(tmp_path):
    from kerncount import discover, gui_tumour
    from kerncount.tumour import layout as lay

    ext = tmp_path / "drive" / "5d WT"
    _touch(ext / "190626 OBS" / "x" / "190626 HERV -clodro n1 A.ome.tiff")
    raw = tmp_path / "raw"
    _touch(raw / "5d K07" / "030826 HERV -clodro n1 A.ome.tiff")
    os.symlink(ext, raw / "5d WT")
    prj = tmp_path / "prj"
    prj.mkdir()
    s = _settings(prj, f'[project]\nraw_dir = "{raw}"\n')
    discover.write_layout(lay.scan(raw, s), prj / "layout.csv", ",")
    (prj / "stats_long.csv").write_text("x\n1\n")
    ext.rename(tmp_path / "drive" / "gone")
    info = gui_tumour.layout_info(prj)
    assert any("not available" in p for p in info["problems"]) and len(info["rows"]) == 2
    assert gui_tumour.results_info(prj) is not None
