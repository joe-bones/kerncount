"""Ticket 05 (ADR 0002/0003 amendments 2026-11-20): mid-section, cross-sectional area, perimeter, shape and burden
readouts. Synthetic shapes at the real pixel size (0.6239 µm); real-data checks (2 stacks) are slow."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

PX = 0.6239  # real xy pixel size (µm)
PERIMETER_TOL = 0.005  # outline length within 0.5 % of the analytic perimeter (discs r >= 20 µm, ellipses)
AREA_TOL = 0.005


def _ellipse(a_um: float, b_um: float, px: float = PX, sub: int = 4, angle: float = 0.0) -> np.ndarray:
    """Digitised filled ellipse (semi-axes a along x, b along y, rotated by `angle`), pixel = inside at its centre
    after `sub`-fold supersampling majority (i.e. area-correct digitisation)."""
    n = int(2 * max(a_um, b_um) / px) + 20
    c = n * px / 2
    off = (np.arange(sub) + 0.5) / sub
    y = ((np.arange(n)[:, None] + off[None]).ravel() * px - c)
    x = y.copy()
    yy, xx = np.meshgrid(y, x, indexing="ij")
    ca, sa = math.cos(angle), math.sin(angle)
    u, v = xx * ca + yy * sa, -xx * sa + yy * ca
    inside = (u / a_um) ** 2 + (v / b_um) ** 2 <= 1
    return inside.reshape(n, sub, n, sub).mean(axis=(1, 3)) >= 0.5


def _ramanujan(a: float, b: float) -> float:
    h = ((a - b) / (a + b)) ** 2
    return math.pi * (a + b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))


# ---------------------------------------------------------------- seam 1: outline length of a digitised section


@pytest.mark.parametrize("r", [20.0, 80.0, 300.0, 500.0])
def test_outline_length_of_digitised_discs_matches_2_pi_r(r):
    from kerncount.tumour import shape

    m = _ellipse(r, r)
    p = shape.outline_length(m, (PX, PX), smooth_um=2.0)
    assert p == pytest.approx(2 * math.pi * r, rel=PERIMETER_TOL)
    assert shape.section_area(m, (PX, PX)) == pytest.approx(math.pi * r * r, rel=AREA_TOL)


@pytest.mark.parametrize("r", [80.0, 250.0, 500.0])
def test_at_the_default_8_um_outline_scale_discs_from_80_um_radius_stay_within_tolerance(r):
    """Real tumour sections are r ~ 250-400 µm; the default scale's shrinkage (~sigma²/2r) is negligible there."""
    from kerncount.tumour import shape

    assert shape.ShapeOptions().outline_smooth_um == 8.0
    p = shape.outline_length(_ellipse(r, r), (PX, PX), smooth_um=8.0)
    assert p == pytest.approx(2 * math.pi * r, rel=PERIMETER_TOL)


@pytest.mark.parametrize("a,b,angle", [(400.0, 150.0, 0.0), (300.0, 120.0, 0.6), (250.0, 200.0, 1.1)])
def test_outline_length_of_digitised_ellipses_matches_ramanujan(a, b, angle):
    from kerncount.tumour import shape

    p = shape.outline_length(_ellipse(a, b, angle=angle), (PX, PX), smooth_um=2.0)
    assert p == pytest.approx(_ramanujan(a, b), rel=PERIMETER_TOL)


def test_pixel_counting_overestimates_but_the_estimator_does_not():
    """Why not count boundary pixel edges: a digitised disc's staircase is 4/pi = 27 % too long at every scale."""
    from kerncount.tumour import shape

    m = _ellipse(300.0, 300.0)
    edges = (np.count_nonzero(m[:, 1:] != m[:, :-1]) + np.count_nonzero(m[1:] != m[:-1])) * PX
    assert edges / (2 * math.pi * 300) == pytest.approx(4 / math.pi, rel=0.01)
    assert shape.outline_length(m, (PX, PX), 2.0) / (2 * math.pi * 300) == pytest.approx(1.0, rel=PERIMETER_TOL)


# ---------------------------------------------------------------- seam 2: mid-section selection


def test_mid_section_is_centred_on_the_largest_cross_section_with_fixed_physical_thickness():
    from kerncount.tumour import shape

    areas = np.array([0, 1, 3, 6, 9, 10, 9.5, 7, 4, 2, 1, 0, 0, 0, 0, 0], float)
    ms = shape.mid_section(areas, z_um=2.36, thickness_um=20.0)
    assert (ms.centre, ms.first, ms.last, ms.fits) == (5, 1, 9, True)  # +-4 layers = +-9.4 µm
    assert shape.mid_section(np.array([0, 5, 5, 5, 0.0]), 2.0, 2.0).centre == 2  # plateau: its middle
    ms = shape.mid_section(areas, z_um=1.75, thickness_um=20.0)
    assert (ms.first, ms.last, ms.fits) == (0, 10, True)  # +-5 layers = +-8.75 µm
    ms = shape.mid_section(areas, z_um=1.0, thickness_um=20.0)
    assert (ms.first, ms.last, ms.fits) == (0, 15, False)  # needs layers -5..15: clipped and flagged


# ---------------------------------------------------------------- seam 3: segment -> confirmatory readouts

from tumour_synth import SynthTumour, write_ome  # noqa: E402

VOX = (2.36, PX, PX)
SEMI = (40.0, 160.0, 110.0)  # (z, y, x) µm
SHAPE = (40, 600, 460)
READOUT_TOL = 0.02  # vs the analytic mid-section mean of the blurred, thresholded synthetic ellipsoid


def _truth(t: SynthTumour, centre_layer: int, half: int) -> tuple[float, float]:
    """Analytic mean cross-sectional area and perimeter of the ellipsoid over layers centre +- half."""
    cz = t._centre()[0]
    a_z, b, c = t.semi_axes_um
    areas, perims = [], []
    for k in range(centre_layer - half, centre_layer + half + 1):
        z = (k + 0.5) * t.voxel_um[0]
        s = math.sqrt(max(1 - ((z - cz) / a_z) ** 2, 0))
        areas.append(math.pi * b * c * s * s)
        perims.append(_ramanujan(b * s, c * s))
    return float(np.mean(areas)), float(np.mean(perims))


def _seg(t: SynthTumour, **kw):
    from kerncount.tumour import segment

    return segment.segment(t.stack(), t.voxel_um, **kw)


def test_cross_sectional_area_and_perimeter_match_the_analytic_mid_section_of_an_off_centre_ellipsoid():
    t = SynthTumour(semi_axes_um=SEMI, centre_um=(31.0, 190.0, 140.0), voxel_um=VOX, shape=SHAPE, plateau=90)
    r = _seg(t).readouts
    centre = int(31.0 // VOX[0])  # layer holding the ellipsoid's centre
    assert abs(r["mid_section_centre_layer"] - (centre + 1)) <= 1  # 1-based in the table
    assert r["mid_section_layers"] == 9 and r["mid_section_fits"] == "yes"
    area, perim = _truth(t, r["mid_section_centre_layer"] - 1, 4)
    assert r["cross_sectional_area_um2"] == pytest.approx(area, rel=READOUT_TOL)
    assert r["perimeter_um"] == pytest.approx(perim, rel=READOUT_TOL)


BRIGHTNESS_TOL = 0.03


def test_confirmatory_readouts_do_not_move_with_uniform_brightness():
    """Same tumour recorded x0.4 dimmer and x2.5 brighter (clipped at 255): area and perimeter within 3 %."""
    r = {b: _seg(SynthTumour(semi_axes_um=SEMI, voxel_um=VOX, shape=SHAPE, plateau=100, brightness=b)).readouts
         for b in (0.4, 1.0, 2.5)}
    assert r[2.5]["saturated_frac_main_mass"] > 0.05  # the bright version really clips
    for b in (0.4, 2.5):
        for k in ("cross_sectional_area_um2", "perimeter_um"):
            assert r[b][k] == pytest.approx(r[1.0][k], rel=BRIGHTNESS_TOL), (b, k)


def test_confirmatory_readouts_do_not_depend_on_the_z_step():
    """1.75 vs 2.36 µm z-steps (both in the real data): same mid-section thickness in µm, same readouts."""
    out = {}
    for vz, nz in ((1.75, 54), (2.36, 40)):
        t = SynthTumour(semi_axes_um=SEMI, voxel_um=(vz, PX, PX), shape=(nz, 600, 460), plateau=90)
        out[vz] = _seg(t).readouts
    assert out[1.75]["mid_section_layers"] == 11 and out[2.36]["mid_section_layers"] == 9
    for k in ("cross_sectional_area_um2", "perimeter_um"):
        assert out[1.75][k] == pytest.approx(out[2.36][k], rel=0.01), k


def test_the_threshold_sweep_recomputes_the_confirmatory_readouts():
    from kerncount.tumour import segment

    t = SynthTumour(semi_axes_um=(30, 70, 60), voxel_um=(2.36, 1.25, 1.25), shape=(30, 128, 128), plateau=90, seed=3)
    st = t.stack()
    sw = segment.sweep(st, t.voxel_um, [0.8, 1.2])
    hi = segment.segment(st, t.voxel_um, threshold_scale=1.2).readouts
    for k in ("cross_sectional_area_um2", "perimeter_um", "mid_section_centre_layer"):
        assert sw[1.2][k] == hi[k], k
    assert sw[0.8]["cross_sectional_area_um2"] > hi["cross_sectional_area_um2"]


def test_mid_section_that_does_not_fit_in_the_stack_is_clipped_and_reported():
    t = SynthTumour(semi_axes_um=SEMI, centre_um=(4.0, 190.0, 140.0), voxel_um=VOX, shape=SHAPE, plateau=90)
    r = _seg(t).readouts
    assert r["mid_section_fits"] == "no" and r["mid_section_first_layer"] == 1
    assert r["mid_section_layers"] < 9


# ---------------------------------------------------------------- seam 4: exploratory shape readouts (2D, mid-section)


def _rough_disc(r: float, amp: float, waves: int) -> np.ndarray:
    n = int(2 * (r + amp) / PX) + 20
    c = n * PX / 2
    y, x = np.mgrid[:n, :n] * PX - c + PX / 2
    rho, th = np.hypot(y, x), np.arctan2(y, x)
    return rho <= r + amp * np.sin(waves * th)


def _lobed(r: float) -> np.ndarray:
    """Three discs of radius r touching a smaller centre: lobed (concave), not elongated."""
    n = int(6 * r / PX) + 20
    c = n * PX / 2
    y, x = np.mgrid[:n, :n] * PX - c
    m = np.hypot(y, x) <= 0.8 * r
    for a in (0, 2 * math.pi / 3, 4 * math.pi / 3):
        m |= np.hypot(y - 1.6 * r * math.sin(a), x - 1.6 * r * math.cos(a)) <= r
    return m


def test_shape_descriptors_separate_elongation_lobes_and_a_ragged_edge():
    from kerncount.tumour import shape

    o = shape.ShapeOptions()
    disc = shape.section_shape(_ellipse(250, 250), (PX, PX), o)
    ellipse = shape.section_shape(_ellipse(400, 150), (PX, PX), o)
    rough = shape.section_shape(_rough_disc(250, 6, 120), (PX, PX), o)  # 6 µm ripples every 13 µm
    lobed = shape.section_shape(_lobed(120), (PX, PX), o)
    assert disc["elongation"] == pytest.approx(1.0, abs=0.01) and ellipse["elongation"] == pytest.approx(400 / 150, rel=0.01)
    assert rough["elongation"] == pytest.approx(1.0, abs=0.02)
    # edge roughness: only the ragged outline is longer than its heavily smoothed self
    assert disc["edge_roughness"] == pytest.approx(1.0, abs=0.01)
    assert ellipse["edge_roughness"] == pytest.approx(1.0, abs=0.01)
    assert rough["edge_roughness"] > 1.3
    # solidity: convex shapes ~1 (also when elongated), lobes clearly less
    assert disc["solidity"] > 0.995 and ellipse["solidity"] > 0.995
    assert lobed["solidity"] < 0.9
    # circularity (2D "sphericity", reference only) confounds elongation and lobes
    assert disc["circularity"] == pytest.approx(1.0, abs=0.01)
    assert max(ellipse["circularity"], lobed["circularity"]) < 0.9


def test_footprint_slab_volume_burden_and_invasion_from_labels():
    from kerncount.tumour import shape

    vox = (2.0, 1.0, 1.0)
    lab = np.zeros((30, 80, 80), np.uint8)
    lab[5:25, 10:50, 10:40] = 1  # 20 layers x 40 x 30 box main mass
    lab[10, 20:30, 45:55] = 0  # (outside the box: no effect)
    lab[12, 60:70, 60:70] = 2  # satellite 10 x 10 x 1
    r = shape.readouts(lab, vox, shape.ShapeOptions(mid_section_um=10, slab_um=20))
    assert r["footprint_area_um2"] == pytest.approx(40 * 30)
    assert r["cross_sectional_area_um2"] == pytest.approx(40 * 30)
    assert r["slab_volume_um3"] == pytest.approx(11 * 2.0 * 40 * 30)  # +-5 layers of 2 µm around the centre layer
    main, sat = 20 * 40 * 30 * 2.0, 100 * 2.0
    assert r["tumour_burden_um3"] == pytest.approx(main + sat)
    assert r["invasion_fraction"] == pytest.approx(sat / (main + sat))


def test_mesh_surface_area_of_a_sphere_on_anisotropic_voxels():
    """Exploratory 3D surface (demoted): marching cubes on an isotropic 2 µm grid, within 3 % of 4 pi r^2."""
    from kerncount.tumour import shape

    vox = (2.36, PX, PX)
    r = 40.0
    z = (np.arange(60) + 0.5) * vox[0] - 70
    y = (np.arange(160) + 0.5) * vox[1] - 50
    m = (z[:, None, None] ** 2 + y[None, :, None] ** 2 + y[None, None, :] ** 2) <= r * r
    assert shape.surface_area(m, vox) == pytest.approx(4 * math.pi * r * r, rel=0.03)


# ---------------------------------------------------------------- seam 5: Fiji legacy counterparts from the lab's mask


def test_fiji_mid_section_area_and_perimeter_come_from_the_labs_mask(tmp_path):
    """A Fiji-like 0/255 mask spanning every layer (as the lab's per-layer Triangle masks do) with its widest layer
    at 14: the same mid-section rule, applied to Fiji's largest object, gives the analytic discs' mean."""
    import tifffile

    from kerncount.tumour import legacy

    vox = (2.36, PX, PX)
    nz, n = 30, 560
    radii = 150 - 0.15 * ((np.arange(nz) - 14) * vox[0]) ** 2  # widest (150 µm) in layer 14, >= ~70 µm everywhere
    c = n * PX / 2
    y, x = (np.mgrid[:n, :n] + 0.5) * PX - c
    rho = np.hypot(y, x)
    mask = np.stack([rho <= r for r in radii]).astype(np.uint8) * 255
    mask[:, 5:12, 5:12] = 255  # a small separate object (not the main one)
    stack = write_ome(tmp_path / "t.ome.tiff", np.zeros((nz, n, n), np.uint8), vox)
    tifffile.imwrite(tmp_path / "tf.tif", mask)

    r = legacy.tumour_readouts(stack, None, tmp_path / "tf.tif", min_voxels=1000)

    used = radii[10:19]  # 9 layers = 20 µm at 2.36 µm
    assert r["fiji_mid_section_centre_layer"] == 15 and r["fiji_mid_section_fits"] == "yes"
    assert r["fiji_cross_sectional_area_um2"] == pytest.approx(np.mean(np.pi * used ** 2), rel=AREA_TOL)
    assert r["fiji_perimeter_um"] == pytest.approx(np.mean(2 * np.pi * used), rel=PERIMETER_TOL)
    assert r["fiji_section_source"] == "Fiji mask"


# ---------------------------------------------------------------- seam 6: per-tumour cache follows the readouts


def test_cached_measurements_from_before_the_readouts_or_with_other_scales_are_re_measured(tmp_path, monkeypatch):
    import pandas as pd

    from kerncount.tumour import measure

    t = SynthTumour(semi_axes_um=(20, 50, 40), voxel_um=(2.36, 1.25, 1.25), shape=(30, 128, 128), plateau=90)
    write_ome(tmp_path / "raw" / "t.ome.tiff", t.stack(), t.voxel_um)
    layout = pd.DataFrame({"tumour_id": ["a"], "file": ["t.ome.tiff"], "include": ["yes"]})
    s = {"project": {"raw_dir": str(tmp_path / "raw")}, "tumour": {}}
    monkeypatch.setattr(measure, "CACHE_VERSION", 1)  # a cache written before ticket 05 ...
    old = measure.measure(tmp_path, layout, s)
    meta = next((tmp_path / measure.WORK).glob("*.json"))
    import json

    j = json.loads(meta.read_text())
    j["readouts"] = {k: v for k, v in j["readouts"].items() if k not in ("cross_sectional_area_um2", "perimeter_um")}
    meta.write_text(json.dumps(j))  # ... lacked the new columns
    monkeypatch.undo()
    new = measure.measure(tmp_path, layout, s)
    assert "cross_sectional_area_um2" in new.columns and new["cross_sectional_area_um2"].iloc[0] > 0
    assert measure.measure(tmp_path, layout, s).equals(new)  # cached now
    thin = measure.measure(tmp_path, layout, {**s, "tumour": {"mid_section_um": 5.0}})
    assert thin["mid_section_layers"].iloc[0] == 3 and new["mid_section_layers"].iloc[0] == 9
    assert old["main_mass_volume_um3"].iloc[0] == new["main_mass_volume_um3"].iloc[0]


# ---------------------------------------------------------------- seam 7: analyse -> shape panel


def test_analyse_draws_the_shape_panel_and_analyses_area_and_perimeter_as_confirmatory(tmp_path):
    from kerncount.experiment_types import RunOptions, get
    from test_tumour_stats import _settings, synth_table

    rng = np.random.default_rng(4)
    t = synth_table(rng, genotypes=("TLR7 KO",)).rename(columns={"volume": "cross_sectional_area_mm2"})
    t["perimeter_mm"] = 2 * np.sqrt(np.pi * t["cross_sectional_area_mm2"]) * np.exp(rng.normal(0.1, 0.05, len(t)))
    for col, lo, hi in (("elongation", 1, 2), ("solidity", 0.8, 1), ("edge_roughness", 1, 1.3),
                        ("circularity", 0.5, 1), ("invasion_fraction", 0, 0.1)):
        t[col] = rng.uniform(lo, hi, len(t))
    t["n_satellites"] = rng.integers(0, 4, len(t))
    s = _settings()
    get("tumour_volume").analyse(tmp_path, t, s, RunOptions())
    assert (tmp_path / "plots" / "shape_panel.png").exists()
    report = (tmp_path / "stats_report.md").read_text()
    assert "cross_sectional_area (confirmatory)" in report and "perimeter (confirmatory)" in report
    assert "elongation (exploratory)" in report


# ---------------------------------------------------------------- real data (slow, 2 stacks)

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7" / "ready 5d KO7 030826"


@pytest.mark.slow
@pytest.mark.parametrize("stem", ["030826 cntrl +clodro n1 B", "030826 HERV -clodro n2 A"])
def test_real_mid_section_area_and_perimeter_beside_fijis(stem):
    """Real stacks: the mid-section fits, area and perimeter are plausible, and the Fiji counterparts from the lab's
    mask (same mid-section rule) are of the same order (methods differ: per-layer Triangle keeps more haze)."""
    import time

    from kerncount.tumour import config, legacy, segment

    orig = REAL / f"{stem}.ome.tiff"
    if not orig.exists():
        pytest.skip("real tumour data not present")
    t0 = time.time()
    r = segment.measure_file(orig, config.tumour_settings({})).readouts
    f = legacy.tumour_readouts(orig, None, REAL / f"{stem}f.tif", fallback=False)
    a, p = r["cross_sectional_area_um2"], r["perimeter_um"]
    print(f"\n{stem}: centre layer {r['mid_section_centre_layer']} ({r['mid_section_layers']} layers, fits "
          f"{r['mid_section_fits']}), area {a / 1e6:.4f} mm², perimeter {p / 1e3:.3f} mm, circularity "
          f"{r['circularity']:.3f}, roughness {r['edge_roughness']:.3f}, solidity {r['solidity']:.3f}, elongation "
          f"{r['elongation']:.2f}, footprint {r['footprint_area_um2'] / 1e6:.4f} mm², slab {r['slab_volume_um3'] / 1e9:.5f} "
          f"mm³, surface {r['main_mass_surface_um2'] / 1e6:.4f} mm², invasion {r['invasion_fraction']:.4f} | Fiji: "
          f"centre {f['fiji_mid_section_centre_layer']}, area {f['fiji_cross_sectional_area_um2'] / 1e6:.4f} mm² "
          f"(ratio {a / f['fiji_cross_sectional_area_um2']:.3f}), perimeter {f['fiji_perimeter_um'] / 1e3:.3f} mm "
          f"(ratio {p / f['fiji_perimeter_um']:.3f}); {time.time() - t0:.0f} s")
    assert r["mid_section_fits"] == "yes" and f["fiji_section_source"] == "Fiji mask"
    assert 0.02e6 < a < 2e6 and 2 * math.sqrt(math.pi * a) <= p  # a perimeter is never shorter than a circle's
    assert 0.5 < a / f["fiji_cross_sectional_area_um2"] < 2 and 0.5 < p / f["fiji_perimeter_um"] < 2
