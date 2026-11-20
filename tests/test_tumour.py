"""tumour_volume experiment type (ADR 0002-0004): layout from flat files, main-mass volume, readout table, CLI."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import experiment_types, settings as settings_mod

KO_SESSION = "ready 5d KO7 030826"


def _tumour_settings(tmp_path: Path) -> dict:
    (tmp_path / "settings.toml").write_text('[experiment]\ntype = "tumour_volume"\n')
    return settings_mod.load(tmp_path)


def _touch(folder: Path, *names: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for n in names:
        (folder / n).write_bytes(b"")


# ---------------------------------------------------------------- seam 1: raw folder -> layout rows


def test_scan_layout_reads_factors_from_filenames_and_session_folder(tmp_path):
    raw = tmp_path / "raw"
    _touch(raw / KO_SESSION,
           "030826 HERV -clodro n3 A.ome.tiff", "030826 HERV -clodro n3 Af.tif",
           "Statistics for 030826 HERV -clodro n3 Af.csv",
           "030826 HERV -clodro n3 B.ome.tiff", "030826 HERV -clodro n3 Bf.tif",
           "Statistics for 030826 HERV -clodro n3 Bf.txt",
           "030826 cntrl +clodro n3 A.ome.tiff", "030826 N4 -clodro n1 B.ome.tiff")
    _touch(raw / "ready 5d WT 120926", "120926 HERV -clodro n3 A.ome.tiff")
    s = _tumour_settings(tmp_path)

    layout = experiment_types.get("tumour_volume").scan_layout(raw, s)

    assert len(layout) == 5  # legacy Fiji masks / statistics are not tumours
    row = layout.set_index("tumour_id").loc["TLR7 KO / HERV / -clodro / n3A"]
    assert row["genotype"] == "TLR7 KO" and row["tumour_line"] == "HERV" and row["clodronate"] == "-clodro"
    assert row["section"] == "n3" and row["hemisphere"] == "A" and row["tumour"] == "n3A"
    assert row["imaging_session"] == KO_SESSION and row["file_date"] == "2026-08-03"
    assert row["file"] == f"{KO_SESSION}/030826 HERV -clodro n3 A.ome.tiff"
    assert row["legacy_mask"] == f"{KO_SESSION}/030826 HERV -clodro n3 Af.tif"
    assert row["legacy_stats"] == f"{KO_SESSION}/Statistics for 030826 HERV -clodro n3 Af.csv"
    assert row["include"] == "yes"
    b = layout.set_index("tumour_id").loc["TLR7 KO / HERV / -clodro / n3B"]
    assert b["legacy_stats"].endswith("n3 Bf.txt")
    wt = layout[layout["genotype"] == "WT"]
    assert list(wt["tumour_id"]) == ["WT / HERV / -clodro / n3A"]  # same n3A, other genotype: still unique
    assert layout["tumour_id"].is_unique
    n4 = layout[layout["tumour_line"] == "N4"].iloc[0]  # the tumour line N4 is not a section label
    assert n4["section"] == "n1" and n4["tumour_id"] == "TLR7 KO / N4 / -clodro / n1B"


def test_scan_layout_normalises_the_clodron_typo_with_a_note(tmp_path):
    raw = tmp_path / "raw"
    _touch(raw / KO_SESSION, "030826 HERV -clodron n6 B.ome.tiff")
    layout = experiment_types.get("tumour_volume").scan_layout(raw, _tumour_settings(tmp_path))
    row = layout.iloc[0]
    assert row["clodronate"] == "-clodro" and row["tumour_id"] == "TLR7 KO / HERV / -clodro / n6B"
    assert "-clodron" in row["note"]


def test_factors_are_declared_in_settings_not_hard_coded(tmp_path):
    """A fourth factor (here: timepoint from a filename token) needs only a settings entry."""
    raw = tmp_path / "raw"
    _touch(raw / KO_SESSION, "030826 HERV -clodro 10d n1 A.ome.tiff")
    (tmp_path / "settings.toml").write_text(
        '[experiment]\ntype = "tumour_volume"\n'
        '[[tumour.factors]]\nname = "genotype"\nlevels = ["WT", "TLR7 KO"]\nreference = "WT"\nrole = "row"\n'
        'source = "session"\naliases = { KO7 = "TLR7 KO" }\n'
        '[[tumour.factors]]\nname = "tumour_line"\nlevels = ["cntrl", "HERV", "N4"]\nreference = "cntrl"\nrole = "x"\n'
        '[[tumour.factors]]\nname = "clodronate"\nlevels = ["-clodro", "+clodro"]\nreference = "-clodro"\n'
        'role = "panel"\n'
        '[[tumour.factors]]\nname = "timepoint"\nlevels = ["5d", "10d"]\nreference = "5d"\nrole = "colour"\n')
    layout = experiment_types.get("tumour_volume").scan_layout(raw, settings_mod.load(tmp_path))
    row = layout.iloc[0]
    assert row["timepoint"] == "10d"
    assert row["tumour_id"] == "TLR7 KO / HERV / -clodro / 10d / n1A"


# ---------------------------------------------------------------- seam 2: z-stack -> main-mass volume

from tumour_synth import Satellite, SynthTumour, write_ome  # noqa: E402

VOX = (2.36, 1.25, 1.25)  # anisotropic (z, y, x) µm; coarser xy than the microscope keeps the tests fast
SHAPE = (40, 256, 256)
VOLUME_TOL = 0.05  # measured main mass within 5 % of the analytic ellipsoid volume
BRIGHTNESS_TOL = 0.03  # uniform x0.4 / x2.5 brightness changes the volume by less than 3 %


def _measure(tmp_path, tumour: SynthTumour, name="t.ome.tiff", voxel_um="from tumour"):
    from kerncount.tumour import config, segment

    path = write_ome(tmp_path / name, tumour.stack(), tumour.voxel_um if voxel_um == "from tumour" else voxel_um)
    return segment.measure_file(path, config.tumour_settings({}))


@pytest.mark.parametrize("semi_axes", [(15, 30, 25), (30, 70, 60), (40, 140, 120)])
def test_main_mass_volume_matches_the_analytic_volume_across_sizes(tmp_path, semi_axes):
    """Small tumours in a large field included: the threshold must not depend on how much of the stack is tumour."""
    t = SynthTumour(semi_axes_um=semi_axes, voxel_um=VOX, shape=SHAPE, plateau=90)
    res = _measure(tmp_path, t)
    assert res.readouts["main_mass_volume_um3"] == pytest.approx(t.volume_um3, rel=VOLUME_TOL)
    assert res.readouts["voxel_z_um"] == pytest.approx(2.36) and res.readouts["n_layers"] == SHAPE[0]


def test_uniform_brightness_scaling_does_not_change_the_volume(tmp_path):
    """Same tumour recorded x0.4 dimmer and x2.5 brighter (clipped at 255, so it saturates): same main mass."""
    vols, sats = {}, {}
    for b in (0.4, 1.0, 2.5):
        t = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=VOX, shape=SHAPE, plateau=100, brightness=b)
        r = _measure(tmp_path, t, name=f"b{b}.ome.tiff").readouts
        vols[b], sats[b] = r["main_mass_volume_um3"], r["saturated_frac_main_mass"]
    assert sats[2.5] > 0.05 and sats[1.0] == 0  # the bright version really clips
    for b in (0.4, 2.5):
        assert vols[b] == pytest.approx(vols[1.0], rel=BRIGHTNESS_TOL), (b, vols, sats)


def test_satellites_do_not_merge_into_the_main_mass(tmp_path):
    cz, cy, cx = (s * v / 2 for s, v in zip(SHAPE, VOX))
    sats = [Satellite((cz, cy, cx + 60 + 8 + 8), 8.0),  # 8 µm gap to the main-mass surface
            Satellite((cz, cy - 70 - 10 - 10, cx), 10.0)]
    t = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=VOX, shape=SHAPE, plateau=90, satellites=sats)
    res = _measure(tmp_path, t)
    assert res.readouts["main_mass_volume_um3"] == pytest.approx(t.volume_um3, rel=VOLUME_TOL)
    assert res.readouts["n_satellites"] == 2
    assert int(res.labels.max()) == 3  # 1 = main mass, 2.. = satellites by size


def test_a_dark_core_is_filled_in_3d(tmp_path):
    """A necrotic / dim core enclosed by tumour belongs to the main mass (ADR 0002: 3D hole filling)."""
    t = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=VOX, shape=SHAPE, plateau=90, hollow_frac=0.6)
    assert _measure(tmp_path, t).readouts["main_mass_volume_um3"] == pytest.approx(t.volume_um3, rel=VOLUME_TOL)


def test_voxel_size_is_read_per_file(tmp_path):
    """The z-step varies between stacks (1.75 / 2.0 / 2.36 µm): the same tumour gives the same volume."""
    a = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=(1.75, 1.25, 1.25), shape=(54, 256, 256), plateau=90)
    b = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=VOX, shape=SHAPE, plateau=90)
    va = _measure(tmp_path, a, "a.ome.tiff").readouts
    vb = _measure(tmp_path, b, "b.ome.tiff").readouts
    assert va["voxel_z_um"] == pytest.approx(1.75)
    assert va["main_mass_volume_um3"] == pytest.approx(vb["main_mass_volume_um3"], rel=VOLUME_TOL)


def test_missing_voxel_size_is_an_error(tmp_path):
    from kerncount.tumour import segment

    t = SynthTumour(semi_axes_um=(15, 30, 25), voxel_um=VOX, shape=(20, 64, 64), plateau=90)
    with pytest.raises(segment.MissingCalibration, match="voxel size"):
        _measure(tmp_path, t, voxel_um=None)


# ---------------------------------------------------------------- seam 3: layout + measurements -> readout table


def test_readout_table_has_one_row_per_included_tumour_with_factors_and_volume(tmp_path):
    s = _tumour_settings(tmp_path)
    layout = pd.DataFrame({
        "tumour_id": ["TLR7 KO / HERV / -clodro / n3A", "TLR7 KO / HERV / -clodro / n3B",
                      "TLR7 KO / cntrl / +clodro / n1A"],
        "imaging_session": KO_SESSION, "file_date": "2026-08-03", "genotype": "TLR7 KO",
        "tumour_line": ["HERV", "HERV", "cntrl"], "clodronate": ["-clodro", "-clodro", "+clodro"],
        "section": ["n3", "n3", "n1"], "hemisphere": ["A", "B", "A"], "tumour": ["n3A", "n3B", "n1A"],
        "file": ["a.ome.tiff", "b.ome.tiff", "c.ome.tiff"], "legacy_mask": "", "legacy_stats": "",
        "include": ["yes", "yes", "no"], "note": ["", "", "out of focus"]})
    measurements = pd.DataFrame({
        "tumour_id": layout["tumour_id"][:2], "voxel_z_um": 2.36, "voxel_y_um": 0.6239, "voxel_x_um": 0.6239,
        "n_layers": [50, 44], "main_mass_volume_um3": [2.5e8, 1.0e7], "threshold_method": "intermedians",
        "threshold": [40.0, 25.5], "saturated_frac": [0.01, 0.0], "saturated_frac_main_mass": [0.05, 0.0]})

    table = experiment_types.get("tumour_volume").readout_table(tmp_path, layout, measurements, s)

    assert list(table["tumour_id"]) == list(layout["tumour_id"][:2])  # excluded tumour is not a datapoint
    assert list(table["main_mass_volume_mm3"]) == pytest.approx([0.25, 0.01])
    for col in ("genotype", "tumour_line", "clodronate", "imaging_session", "voxel_z_um", "n_layers", "threshold",
                "saturated_frac", "main_mass_volume_um3"):
        assert col in table.columns, col
    on_disk = pd.read_csv(tmp_path / "tumours.csv")
    assert list(on_disk["tumour_id"]) == list(table["tumour_id"])


# ---------------------------------------------------------------- seam 4: kerncount init --type / run end-to-end

WT_SESSION = "ready 5d WT 120926"
SMALL = dict(voxel_um=VOX, shape=(30, 160, 160), plateau=90)


def _experiment(tmp_path, monkeypatch) -> dict[str, SynthTumour]:
    """data/tumour_exp/<KO session>/ and <WT session>/ with flat synthetic stacks + one Fiji legacy file."""
    tumours = {
        f"{KO_SESSION}/030826 HERV -clodro n1 A.ome.tiff": SynthTumour(semi_axes_um=(20, 50, 40), seed=1, **SMALL),
        f"{KO_SESSION}/030826 HERV -clodro n1 B.ome.tiff": SynthTumour(semi_axes_um=(15, 35, 30), seed=2, **SMALL),
        f"{KO_SESSION}/030826 cntrl +clodro n1 A.ome.tiff": SynthTumour(semi_axes_um=(20, 60, 50), seed=3, **SMALL),
        f"{WT_SESSION}/120926 N4 -clodro n2 B.ome.tiff": SynthTumour(semi_axes_um=(18, 45, 45), seed=4, **SMALL),
    }
    raw = tmp_path / "data" / "tumour_exp"
    for rel, t in tumours.items():
        write_ome(raw / rel, t.stack(), t.voxel_um)
    (raw / KO_SESSION / "030826 HERV -clodro n1 Af.tif").write_bytes(b"legacy")
    monkeypatch.chdir(tmp_path)
    return tumours


def _snapshot(folder: Path) -> dict[str, tuple[int, int]]:
    return {str(p.relative_to(folder)): (p.stat().st_size, p.stat().st_mtime_ns) for p in folder.rglob("*")}


def test_init_type_and_run_measure_and_plot_a_tumour_experiment(tmp_path, monkeypatch):
    from kerncount import cli

    tumours = _experiment(tmp_path, monkeypatch)
    raw, prj = tmp_path / "data" / "tumour_exp", tmp_path / "analysis" / "tumour_exp"
    before = _snapshot(raw)

    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    assert settings_mod.load(prj)["experiment"]["type"] == "tumour_volume"
    layout = pd.read_csv(prj / "layout.csv")
    assert len(layout) == 4 and set(layout["genotype"]) == {"WT", "TLR7 KO"}

    cli.main(["run", "tumour_exp"])
    table = pd.read_csv(prj / "tumours.csv")
    assert len(table) == 4
    for _, row in table.iterrows():
        truth = tumours[row["file"]].volume_um3
        assert row["main_mass_volume_um3"] == pytest.approx(truth, rel=VOLUME_TOL), row["tumour_id"]
    assert (prj / "plots" / "main_mass_volume_absolute.png").exists()
    assert _snapshot(raw) == before  # data/ untouched
    work = sum(p.stat().st_size for p in (prj / "_work").rglob("*") if p.is_file())
    assert work < 2_000_000  # cached masks are compressed and small

    cached = _snapshot(prj / "_work" / "tumour")
    cli.main(["run", "tumour_exp"])  # re-run uses the cache: nothing re-segmented
    assert _snapshot(prj / "_work" / "tumour") == cached
    pd.testing.assert_frame_equal(pd.read_csv(prj / "tumours.csv"), table)


def test_group_plot_works_with_one_genotype(tmp_path, monkeypatch):
    from kerncount import cli

    _experiment(tmp_path, monkeypatch)
    import shutil

    shutil.rmtree(tmp_path / "data" / "tumour_exp" / WT_SESSION)
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    cli.main(["run", "tumour_exp", "--no-review"])
    prj = tmp_path / "analysis" / "tumour_exp"
    assert set(pd.read_csv(prj / "tumours.csv")["genotype"]) == {"TLR7 KO"}
    assert (prj / "plots" / "main_mass_volume_absolute.png").exists()


def test_group_plot_keeps_long_titles_and_axis_labels_inside_the_figure(monkeypatch):
    """Real run (ticket 11): the confirmatory suptitle and the normalised y label were cut off at the figure edge."""
    from kerncount.tumour import plots
    from kerncount.tumour.config import Factor

    figs = []
    monkeypatch.setattr(plots, "_save", lambda fig, out: figs.append(fig) or out)
    facs = [Factor("tumour_line", ("cntrl", "HERV", "N4"), "cntrl", "x"),
            Factor("clodronate", ("-clodro", "+clodro"), "-clodro", "panel")]
    t = pd.DataFrame({"tumour_line": ["cntrl", "HERV", "N4"] * 4, "clodronate": ["-clodro"] * 6 + ["+clodro"] * 6,
                      "v": np.linspace(0.1, 0.5, 12), "flagged": [True] + [False] * 11})
    long = ("% of model-estimated cntrl in the same cell of the other factors; bars = fold-change 95 % CI; "
            "stars: F1 (Holm)")
    plots.group_plot(t, facs, "v", "Cross-sectional area, mid-section (mm²), % of cntrl", Path("x"), suptitle=long,
                     hollow="flagged")
    fig = figs[0]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    for text in [fig._suptitle, *[ax.yaxis.label for ax in fig.axes if ax.get_ylabel()]]:
        e = text.get_window_extent(r)
        assert e.x0 >= -1 and e.x1 <= W + 1 and e.y0 >= -1 and e.y1 <= H + 1, (text.get_text(), e, W, H)


# ---------------------------------------------------------------- real data smoke test (slow)

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7" / "5d K07"  # renamed from KO_SESSION 2026-11-20
REAL_STEM = "030826 cntrl +clodro n1 B"


@pytest.mark.slow
@pytest.mark.skipif(not (REAL / f"{REAL_STEM}.ome.tiff").exists(), reason="real tumour data not present")
def test_real_stack_gives_a_plausible_main_mass_volume():
    """One real stack: calibrated from its OME metadata, a plausible volume, roughly comparable with Fiji's largest
    object (methods differ by design: no rolling ball, intermedians vs Triangle, 3D vs 2D fill), so only printed."""
    from kerncount.tumour import config, segment

    res = segment.measure_file(REAL / f"{REAL_STEM}.ome.tiff", config.tumour_settings({}))
    r = res.readouts
    assert r["voxel_z_um"] == pytest.approx(2.0) and r["voxel_x_um"] == pytest.approx(0.6239258)
    assert 1e6 < r["main_mass_volume_um3"] < 1e9  # 0.001-1 mm³: a 5-day tumour in a 640 x 640 x 86 µm stack
    fiji = pd.read_csv(REAL / f"Statistics for {REAL_STEM}f.csv")["Volume (micron^3)"].max()
    ratio = r["main_mass_volume_um3"] / fiji
    print(f"\nmain mass {r['main_mass_volume_um3']:.4g} µm³, Fiji largest object {fiji:.4g} µm³, ratio {ratio:.3f}, "
          f"threshold {r['threshold']:.1f} (background {r['background_level']:.1f}, tumour {r['tumour_level']:.1f}), "
          f"satellites {r['n_satellites']}, saturated in main mass {r['saturated_frac_main_mass']:.3%}")
    assert 0.25 < ratio < 4
