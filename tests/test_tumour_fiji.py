"""Fiji legacy readout and method concordance (ticket 07, ADR 0002): parsing Fiji's 3D Objects Counter tables,
the macro re-implementation on synthetic stacks, the threshold scale of the primary segmentation, concordance and
sweep tables. Real-data checks (2 stacks) are slow and skip without the data."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tumour_synth import SynthTumour, write_ome

HEADER = ("Volume (micron^3),Surface (micron^2),Nb of obj. voxels,Nb of surf. voxels,IntDen,Mean,StdDev,Median,Min,"
          "Max,X,Y,Z,Mean dist. to surf. (micron),SD dist. to surf. (micron),Median dist. to surf. (micron),XM,YM,ZM,"
          "BX,BY,BZ,B-width,B-height,B-depth")
ROW_SMALL = ("3847005.250,611277.125,4187404,569675,1071909568.000,255.984,0.998,255,255,255,281.580,575.427,19.949,"
             "104.041,28.938,101.178,280.854,573.146,19.778,20,219,1,489,642,44")
ROW_BIG = ("18981228.000,1783628.250,20660766,1962261,4294967296.000,207.880,49.705,255,255,255,1175.361,865.184,"
           "22.298,216.552,73.006,218.248,1439.211,1064.566,27.112,605,380,1,1027,1176,48")


# ---------------------------------------------------------------- seam 1: Fiji statistics -> objects -> readouts


def test_statistics_with_header_give_one_object_per_row():
    from kerncount.tumour import fiji

    st = fiji.parse_stats(f"{HEADER}\n{ROW_SMALL}\n{ROW_BIG}\n")
    assert st.problems == []
    assert [o["voxels"] for o in st.objects] == [4187404, 20660766]
    assert st.objects[1]["surface_um2"] == 1783628.25 and st.objects[1]["b_depth"] == 48


def test_spreadsheet_damage_is_skipped_and_reported_never_guessed():
    """The real `N4 -clodro n3 B` file: `;;` endings, CRLF, a blank `;;` line and decimal-comma sum rows."""
    from kerncount.tumour import fiji

    text = (f"{HEADER};;\r\n{ROW_SMALL};;\r\n{ROW_BIG};;\r\n;;\r\n3847005,25;;611.277\r\n18981228;;1783628,25\r\n"
            "22828233,25;;2.394.905\r\n")
    st = fiji.parse_stats(text, "s.csv")
    assert len(st.objects) == 2
    assert len(st.problems) == 3 and all("s.csv line" in p for p in st.problems)
    assert "3847005,25;;611.277" in st.problems[0]


def test_header_less_tab_separated_txt_is_read_in_fijis_column_order():
    from kerncount.tumour import fiji

    st = fiji.parse_stats(ROW_BIG.replace(",", "\t") + "\n" + ROW_SMALL.replace(",", "\t") + "\n")
    assert st.problems == [] and st.objects[0]["volume_um3"] == 18981228.0 and st.objects[0]["bz"] == 1


def test_truncated_header_and_a_leading_row_number_column():
    from kerncount.tumour import fiji

    short = ",".join(HEADER.split(",")[:16])  # the header Excel cut off; rows still carry all 25 values
    assert fiji.parse_stats(f"{short}\n{ROW_SMALL}\n").objects[0]["b_depth"] == 44
    st = fiji.parse_stats(f" ,{HEADER}\n1,{ROW_SMALL}\n2,{ROW_BIG}\n")
    assert [o["voxels"] for o in st.objects] == [4187404, 20660766] and st.problems == []
    assert fiji.parse_stats("").problems  # empty file: reported


def test_readouts_take_the_largest_object_not_the_first_row():
    from kerncount.tumour import fiji

    objs = fiji.parse_stats(f"{HEADER}\n{ROW_SMALL}\n{ROW_BIG}\n").objects
    r = fiji.readouts(objs, n_layers=48)
    assert r["fiji_main_volume_mm3"] == pytest.approx(18981228e-9)
    assert r["fiji_main_surface_um2"] == 1783628.25 and r["fiji_main_surface_mm2"] == pytest.approx(1.78362825)
    assert r["fiji_n_objects"] == 2 and r["fiji_total_volume_mm3"] == pytest.approx((18981228 + 3847005.25) / 1e9)
    assert r["fiji_main_z_start"] == 1 and r["fiji_main_z_layers"] == 48 and r["fiji_main_full_depth"] == "yes"
    assert fiji.readouts(objs, min_voxels=5_000_000)["fiji_n_objects"] == 1
    empty = fiji.readouts([])
    assert math.isnan(empty["fiji_main_volume_mm3"]) and empty["fiji_n_objects"] == 0


# ---------------------------------------------------------------- seam 2: macro re-implementation (synthetic)


def test_triangle_matches_imagej_on_a_hand_worked_histogram():
    from kerncount.tumour import fiji

    # peak at 2, long tail to the right: ImageJ reverses the histogram, the line runs from the peak (2, 100) to
    # the empty bin after the last value (8, 0); the furthest bin below the line is 4 -> level 4 + 1 (reversal)
    h = np.zeros(10, int)
    h[[1, 2, 3, 4, 5, 6, 7]] = [10, 100, 60, 10, 8, 6, 4]
    assert fiji.triangle(h) == 5
    assert fiji.triangle(h[::-1].copy()) == 4  # mirrored: same line, no reversal, split 5 - 1
    assert fiji.triangle(np.zeros(256, int)) == 0


def test_rolling_ball_removes_flat_background_and_keeps_small_spots():
    from kerncount.tumour import fiji

    ball, shrink = fiji._ball(50)
    assert shrink == 4 and ball.shape == (19, 19)  # ImageJ: radius 12.5 on the 4x shrunken layer, trimmed 32 %
    y, x = np.mgrid[:256, :256]
    layer = 30 + 0.05 * x + 120 * np.exp(-((y - 128) ** 2 + (x - 100) ** 2) / (2 * 4.0 ** 2))  # ramp + small spot
    out = fiji.subtract_background(np.clip(layer, 0, 255).astype(np.uint8)[None])[0]
    assert out[128, 100] > 100  # the spot survives
    assert np.percentile(out[np.hypot(y - 128, x - 100) > 30], 99) <= 2  # the ramp is gone


def test_macro_mask_finds_a_bright_object_and_objects_3d_gives_fiji_style_numbers():
    from kerncount.tumour import fiji

    t = SynthTumour(semi_axes_um=(12, 20, 20), voxel_um=(2.0, 1.0, 1.0), shape=(24, 96, 96), plateau=150, seed=3)
    truth = t.occupancy() > 0.5
    mask, _ = fiji.macro_mask(t.stack())
    assert np.count_nonzero(mask & truth) / np.count_nonzero(truth) > 0.95  # the object is found ...
    # ... and, as in the real data, Triangle sits near the background and keeps the blurred tails: larger, deeper
    assert mask.sum() > truth.sum() and mask.any(axis=(1, 2)).sum() > truth.any(axis=(1, 2)).sum()
    labels, objs = fiji.objects_3d(mask, (2.0, 1.0, 1.0), min_voxels=1000)
    assert len(objs) == 1 and objs[0]["voxels"] == int(np.count_nonzero(labels == 1))


def test_objects_3d_counts_26_connected_objects_and_fijis_border_voxel_surface():
    from kerncount.tumour import fiji

    m = np.zeros((6, 10, 10), bool)
    m[1:4, 2:6, 2:7] = True  # 3 x 4 x 5 box
    m[4, 6, 7] = True  # touches the box only diagonally (26-connected) -> same object
    m[0:2, 8:10, 0:2] = True  # separate 2 x 2 x 2 cube at the stack border, 8 voxels
    labels, objs = fiji.objects_3d(m, (2.0, 0.5, 0.5), min_voxels=1)
    assert [o["voxels"] for o in objs] == [61, 8]
    only_box = np.zeros_like(m)
    only_box[1:4, 2:6, 2:7] = True
    o = fiji.objects_3d(only_box, (2.0, 0.5, 0.5), 1)[1][0]
    vz, vy, vx = 2.0, 0.5, 0.5
    assert o["surface_um2"] == pytest.approx(2 * (4 * 5 * vy * vx + 3 * 5 * vz * vx + 3 * 4 * vz * vy))
    assert o["surface_voxels"] == 60 - 1 * 2 * 3 and o["bz"] == 2 and o["b_depth"] == 3
    assert o["volume_um3"] == pytest.approx(60 * vz * vy * vx)
    assert fiji.objects_3d(m, (2.0, 0.5, 0.5), min_voxels=10)[1][0]["voxels"] == 61  # small cube dropped


def test_legacy_source_order_statistics_then_mask_then_re_implementation(tmp_path):
    import tifffile

    from kerncount.tumour import fiji, legacy

    t = SynthTumour(semi_axes_um=(12, 20, 20), voxel_um=(2.0, 1.0, 1.0), shape=(24, 96, 96), plateau=150, seed=3)
    orig = write_ome(tmp_path / "t.ome.tiff", t.stack(), t.voxel_um)
    stats = tmp_path / "Statistics for tf.csv"
    stats.write_text(f"{HEADER};;\n{ROW_SMALL};;\n1,5;;2\n")
    r = legacy.tumour_readouts(orig, stats, None, min_voxels=1000)
    assert r["fiji_source"] == "Fiji statistics" and r["fiji_main_volume_um3"] == 3847005.25
    assert "not a row of numbers" in r["fiji_note"]

    mask, _ = fiji.macro_mask(t.stack())
    mpath = tmp_path / "tf.tif"
    tifffile.imwrite(mpath, (mask * 255).astype(np.uint8))
    r2 = legacy.tumour_readouts(orig, tmp_path / "missing.csv", mpath, min_voxels=1000)
    r3 = legacy.tumour_readouts(orig, None, None, min_voxels=1000)
    assert r2["fiji_source"] == "Fiji mask" and r3["fiji_source"] == "re-implementation"
    assert r2["fiji_main_volume_um3"] == r3["fiji_main_volume_um3"] > 0
    (tmp_path / "bad.tif").write_bytes(b"legacy")
    r4 = legacy.tumour_readouts(orig, None, tmp_path / "bad.tif", min_voxels=1000)
    assert r4["fiji_source"] == "re-implementation" and "unreadable" in r4["fiji_note"]
    off = legacy.tumour_readouts(orig, None, None, min_voxels=1000, fallback=False)
    assert off["fiji_source"] == "" and math.isnan(off["fiji_main_volume_um3"])


# ---------------------------------------------------------------- seam 3: threshold scale -> threshold, sweep


def test_threshold_scale_multiplies_the_automatic_threshold_and_sweep_equals_single_runs():
    from kerncount.tumour import segment

    t = SynthTumour(semi_axes_um=(15, 30, 25), voxel_um=(2.36, 1.25, 1.25), shape=(24, 96, 96), plateau=90, seed=5)
    st = t.stack()
    base = segment.segment(st, t.voxel_um).readouts
    hi = segment.segment(st, t.voxel_um, threshold_scale=1.2).readouts
    assert hi["threshold"] == pytest.approx(1.2 * base["threshold"], abs=2e-3)
    assert hi["threshold_auto"] == pytest.approx(base["threshold"]) and hi["threshold_scale"] == 1.2
    assert "threshold_scale" not in base  # the primary readouts are unchanged
    sw = segment.sweep(st, t.voxel_um, [0.8, 1.2])
    assert sw[1.2]["main_mass_volume_um3"] == hi["main_mass_volume_um3"]
    assert sw[0.8]["main_mass_volume_um3"] > base["main_mass_volume_um3"] > sw[1.2]["main_mass_volume_um3"]


# ---------------------------------------------------------------- seam 4: concordance and sweep tables


def _res(est: dict[str, float], p_adj: dict[str, float], family="F1", role="confirmatory") -> pd.DataFrame:
    return pd.DataFrame([dict(readout="main_mass_volume", readout_role="confirmatory", family=family,
                              family_role=role, family_label="", comparison=k, estimate=v, ci_low=v / 2,
                              ci_high=v * 2, p=p_adj[k], p_adj=p_adj[k], significance="", N=12, df=40,
                              confounded=False, note="" if np.isfinite(p_adj[k]) else "descriptive only")
                         for k, v in est.items()])


def test_concordance_flags_direction_and_significance_agreement_per_contrast():
    from kerncount.tumour import concordance

    new = _res({"a": 2.0, "b": 1.5, "c": 0.8, "d": 1.2, "e": 1.1},
               {"a": 0.001, "b": 0.01, "c": 0.4, "d": 0.03, "e": float("nan")})
    old = _res({"a": 3.0, "b": 1.3, "c": 1.1, "d": 0.7, "e": 1.0},
               {"a": 0.0001, "b": 0.2, "c": 0.6, "d": 0.01, "e": float("nan")})
    c = concordance.compare(new, old, "fiji_main_volume_mm3").set_index("comparison")
    assert list(c["direction_agree"]) == ["yes", "yes", "no", "no", "no"]  # e: 1.1 vs exactly 1 -> no direction
    assert list(c["significance_agree"]) == ["yes", "no", "yes", "no", ""]  # d: both significant, opposite signs
    s = concordance.summary(c.reset_index()).iloc[0]
    assert s["direction_agree"] == "2/5" and s["significance_agree"] == "2/4" and s["conclusion_differs"] == "yes"


def test_sweep_results_rerun_the_statistics_on_each_scale_for_the_given_readouts():
    from kerncount.tumour import concordance
    from kerncount.tumour.config import Factor, Family, Readout

    line = Factor("line", ("cntrl", "HERV"), "cntrl", "x")
    rng = np.random.default_rng(0)
    base = pd.DataFrame({"tumour_id": [f"t{i}" for i in range(16)], "line": ["cntrl"] * 8 + ["HERV"] * 8,
                         "v": np.exp(rng.normal(0, 0.2, 16)) * np.repeat([1.0, 2.0], 8)})
    tabs = {s: base.assign(v=base["v"] * np.where(base["line"] == "HERV", s, 1.0)) for s in (0.8, 1.0, 1.2)}
    sw = concordance.sweep_results(tabs, [line], [Family("F1", "simple", ("line",), "confirmatory")],
                                   [Readout("vol", "v", "Volume", "confirmatory")])
    assert list(sw["scale"]) == [0.8, 1.0, 1.2]
    e = sw.set_index("scale")["estimate"]
    assert e[1.2] / e[1.0] == pytest.approx(1.2) and e[0.8] / e[1.0] == pytest.approx(0.8)
    assert set(sw["significant"]) <= {"yes", "no"}


def test_concordance_outputs_follow_the_declared_confirmatory_readouts(tmp_path):
    """Readout-generic: a confirmatory area readout with its own `legacy` column is compared, volume is not."""
    from kerncount.tumour import concordance
    from kerncount.tumour.config import factors, families
    from test_tumour_stats import _settings, synth_table

    rng = np.random.default_rng(2)
    t = synth_table(rng, genotypes=("TLR7 KO",)).rename(columns={"volume": "area_mm2"})
    t["fiji_area"] = t["area_mm2"] * 1.8 * np.exp(rng.normal(0, 0.1, len(t)))
    t["main_mass_volume_mm3"] = t["area_mm2"]
    t["fiji_main_volume_mm3"] = t["area_mm2"]
    s = _settings(readouts=[{"name": "area", "column": "area_mm2", "label": "Area", "role": "confirmatory",
                             "legacy": "fiji_area"},
                            {"name": "main_mass_volume", "column": "main_mass_volume_mm3", "role": "exploratory"}])
    report = tmp_path / "stats_report.md"
    report.write_text("# r\n")
    written = concordance.analyse(tmp_path, t, s, factors(s), families(s), report)
    assert (tmp_path / "plots" / "area_vs_fiji.png").exists()
    assert not (tmp_path / "plots" / "main_mass_volume_vs_fiji.png").exists()
    cmp = pd.read_csv(tmp_path / "concordance.csv")
    assert set(cmp["readout"]) == {"area"} and set(cmp["fiji_column"]) == {"fiji_area"}
    md = report.read_text()
    assert "Method concordance" in md and "Threshold sweep" in md and "threshold_sweep = false" in md
    assert str(tmp_path / "concordance.csv") in written


# ---------------------------------------------------------------- end to end: run with legacy + sweep


def test_run_adds_fiji_columns_and_the_threshold_sweep(tmp_path, monkeypatch):
    from kerncount import cli, settings as settings_mod
    from test_tumour import _experiment

    _experiment(tmp_path, monkeypatch)
    prj = tmp_path / "analysis" / "tumour_exp"
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    p = prj / "settings.toml"
    p.write_text(p.read_text().replace("threshold_sweep = false", "threshold_sweep = true")
                 .replace("fiji_min_voxels = 200000", "fiji_min_voxels = 1000"))
    assert settings_mod.load(prj)["tumour"]["threshold_sweep"] is True
    cli.main(["run", "tumour_exp"])
    table = pd.read_csv(prj / "tumours.csv")
    assert (table["fiji_main_volume_mm3"] > 0).all()
    assert set(table["fiji_source"]) == {"re-implementation"}  # the only 'legacy' file is not a TIFF
    assert table["fiji_note"].fillna("").str.contains("unreadable").sum() == 1
    sw = pd.read_csv(prj / "threshold_sweep.csv")
    assert sorted(set(sw["scale"])) == [0.8, 0.9, 1.0, 1.1, 1.2]
    # confirmatory readouts since the 2026-11-20 amendment: mid-section area and perimeter (ticket 05)
    for name in ("cross_sectional_area", "perimeter"):
        assert (prj / "plots" / f"{name}_threshold_sweep.png").exists()
        assert (prj / "plots" / f"{name}_vs_fiji.png").exists()
    assert (table["fiji_cross_sectional_area_mm2"] > 0).all()
    assert set(sw["readout"]) >= {"cross_sectional_area", "perimeter"}
    cached = {q.name: q.stat().st_mtime_ns for q in (prj / "_work" / "tumour").glob("*.json")}
    assert any(n.endswith(".sweep.json") for n in cached) and any(n.endswith(".fiji.json") for n in cached)
    cli.main(["run", "tumour_exp"])  # cached: nothing recomputed
    assert {q.name: q.stat().st_mtime_ns for q in (prj / "_work" / "tumour").glob("*.json")} == cached


# ---------------------------------------------------------------- real data (slow, 2 stacks)

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7" / "ready 5d KO7 030826"


@pytest.mark.slow
@pytest.mark.parametrize("stem", ["030826 cntrl +clodro n1 B", "030826 N4 -clodro n3 B"])
def test_re_implementation_reproduces_fijis_mask_and_statistics_on_real_stacks(stem):
    """Fiji's own outputs vs (a) objects counted on Fiji's mask and (b) the macro re-run from the original stack.
    Measured 2026-11-19: both stacks voxel-identical (Dice 1.0, volume ratio 1.0); statistics exact."""
    from kerncount.tumour import fiji, segment

    orig = REAL / f"{stem}.ome.tiff"
    if not orig.exists():
        pytest.skip("real tumour data not present")
    vox = segment.voxel_size_um(orig)
    st = fiji.read_stats(REAL / f"Statistics for {stem}f.csv")
    ref = max(st.objects, key=lambda o: o["voxels"])
    fm = fiji.read_mask(REAL / f"{stem}f.tif")
    _, objs = fiji.objects_3d(fm, vox)
    assert objs[0]["voxels"] == ref["voxels"] and objs[0]["surface_um2"] == pytest.approx(ref["surface_um2"], rel=1e-6)
    assert objs[0]["surface_voxels"] == ref["surface_voxels"] and objs[0]["b_depth"] == ref["b_depth"]
    mask, _ = fiji.macro_mask(segment.read_stack(orig))
    d = fiji.dice(mask, fm)
    _, mine = fiji.objects_3d(mask, vox)
    ratio = mine[0]["volume_um3"] / ref["volume_um3"]
    print(f"\n{stem}: Dice {d:.5f}, main-object volume ratio {ratio:.5f}, statistics problems {st.problems}")
    assert d > 0.99 and ratio == pytest.approx(1.0, abs=0.01)
