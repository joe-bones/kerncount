"""Fluorescence intensity (Ki-67): interleaved channel files, per-nucleus measurement, flat field, positivity,
tuning API, metadata, plate-map control wells, and single-file projects unchanged (docs/intensity.md)."""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile

from kerncount import design, discover, flatfield, imageio, intensity, platemap, positivity, settings
from kerncount import settings_form

RNG = np.random.default_rng(0)


def _touch_tifs(folder: Path, numbers, value=1):
    folder.mkdir(parents=True, exist_ok=True)
    for k in numbers:
        tifffile.imwrite(folder / f"{k}.tif", np.full((8, 8), value, np.uint16))


def _inter_settings(order=("Ki67", "DAPI")):
    s = settings._merge(settings.DEFAULTS, {"channels": {"files": "interleaved", "file_order": list(order),
                                                         "segment": "DAPI"}})
    return s


# ---------------------------------------------------------------- interleaved files

def test_group_files_ignores_gaps_and_reports_leftover():
    nums = [1, 2, 4, 5, 6, 7, 9]  # 3 missing, 8 missing
    fields, left = discover.group_files([(k, Path(f"{k}.tif")) for k in reversed(nums)], 2)
    assert [f for f, _, _ in fields] == [1, 2, 3]
    assert [n for _, _, n in fields] == [(1, 2), (4, 5), (6, 7)]
    assert [p.name for p in left] == ["9.tif"]


def test_records_and_scan_interleaved(tmp_path):
    raw = tmp_path / "raw"
    n = raw / "05.06.26 exp"
    _touch_tifs(n / "cntrl", range(2, 26))  # 2..25: 12 fields, numbering does not start at 1
    _touch_tifs(n / "30e", [1, 2] + list(range(4, 26)))  # skips 3
    _touch_tifs(n / "lyovec", range(1, 8))  # odd: leftover 7
    s = _inter_settings()
    lay = discover.scan(raw, files="interleaved", n_channels=2)
    assert dict(zip(lay["treatment_folder"], lay["images"])) == {"30e": 12, "cntrl": 12, "lyovec": 3}
    assert "leftover 7.tif skipped" in lay.loc[lay["treatment_folder"] == "lyovec", "note"].iloc[0]
    recs, warns = discover.records(lay, raw, s)
    c = [r for r in recs if r.treatment == "Control"]
    assert len(c) == 12 and c[0].numbers == (2, 3) and c[-1].numbers == (24, 25) and c[0].image == 1
    e = [r for r in recs if r.treatment == "miR-30e"]
    assert e[0].numbers == (1, 2) and e[1].numbers == (4, 5) and e[1].image == 2
    assert e[1].paths[0].name == "4.tif" and e[1].path == e[1].paths[0] and e[1].sources == e[1].paths
    assert any("leftover" in w and "7.tif" in w for w in warns)
    # single-file mode is untouched: one record per file, numbered by the file
    recs1, _ = discover.records(discover.scan(raw), raw)
    assert len([r for r in recs1 if r.treatment == "Control"]) == 24 and recs1[0].paths == ()


def test_date_from_day_month_and_file_year(tmp_path):
    raw = tmp_path / "raw"
    _touch_tifs(raw / "Ki67 0611 dapicount" / "cntrl", [1, 2])
    lay = discover.scan(raw, files="interleaved", n_channels=2)
    yr = discover._file_year(raw)
    assert lay["date"].iloc[0] == f"{yr}-11-06" and lay["n"].iloc[0] == f"n1_{yr}-11-06"
    assert "check it" in lay["note"].iloc[0]


def test_rgb_collapse_and_field_stack(tmp_path):
    a = RNG.integers(0, 4000, (16, 20)).astype(np.uint16)
    b = RNG.integers(0, 4000, (16, 20)).astype(np.uint16)
    tifffile.imwrite(tmp_path / "1.tif", np.stack([a, a, a], -1), photometric="rgb")
    tifffile.imwrite(tmp_path / "2.tif", np.stack([b, b, b], -1), photometric="rgb")
    assert imageio.read(tmp_path / "1.tif").data.shape == (1, 16, 20)
    im = imageio.read_field([tmp_path / "1.tif", tmp_path / "2.tif"])
    assert im.data.shape == (2, 16, 20) and np.array_equal(im.data[1], b) and im.max_value == 65535
    tifffile.imwrite(tmp_path / "3.tif", np.stack([a, b, a], -1), photometric="rgb")  # real colour image
    with pytest.raises(ValueError, match="different channels"):
        imageio.read_field([tmp_path / "1.tif", tmp_path / "3.tif"])


def test_detect_interleaved(tmp_path):
    folder = tmp_path / "raw" / "n" / "cntrl"
    folder.mkdir(parents=True)
    for k in range(1, 9):
        v = 300 if k % 2 else 3000  # odd = dim marker, even = bright DAPI
        img = RNG.normal(v, 20, (32, 32)).clip(0).astype(np.uint16)
        tifffile.imwrite(folder / f"{k}.tif", np.stack([img] * 3, -1), photometric="rgb")
    d = discover.detect_interleaved(tmp_path / "raw")
    assert d["interleaved"] and d["bright_position"] == 1


# ---------------------------------------------------------------- synthetic nuclei

def _disks(shape, centres, r):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    lab = np.zeros(shape, np.int32)
    for i, (y, x) in enumerate(centres, 1):
        lab[(yy - y) ** 2 + (xx - x) ** 2 <= r * r] = i
    return lab


def _grid(shape, step, r, jitter=0):
    cs = [(y + RNG.integers(-jitter, jitter + 1), x + RNG.integers(-jitter, jitter + 1))
          for y in range(step, shape[0] - step + 1, step) for x in range(step, shape[1] - step + 1, step)]
    return _disks(shape, cs, r), cs


def test_alignment_check_finds_shift_and_mispairing():
    from scipy import ndimage as ndi

    cs = [(RNG.uniform(10, 246), RNG.uniform(10, 310)) for _ in range(60)]
    lab = _disks((256, 320), cs, 7)
    dapi = np.where(lab > 0, 3000.0, 200.0) + RNG.normal(0, 20, lab.shape)
    ki = np.where(lab > 0, 600.0, 150.0) * (1 - 0.3 * np.linspace(0, 1, 320))[None] + RNG.normal(0, 20, lab.shape)
    a0 = intensity.alignment(dapi, ki)
    assert a0["shift_px"] < 0.5 and a0["corr"] > 0.5
    moved = ndi.shift(ki, (3, -2), order=1, mode="nearest")
    a = intensity.alignment(dapi, moved)
    assert abs(a["dy"] + 3) < 0.4 and abs(a["dx"] - 2) < 0.4 and abs(a["shift_px"] - 3.6) < 0.4
    other = _disks((256, 320), [(RNG.uniform(10, 246), RNG.uniform(10, 310)) for _ in range(60)], 7)
    wrong = np.where(other > 0, 600.0, 150.0) + RNG.normal(0, 20, lab.shape)
    assert intensity.alignment(dapi, wrong)["corr"] < 0.2
    b = intensity.brightness_check(np.stack([ki, dapi]), lab, seg_ch=1)
    assert b["nuclear_brightest"] and not intensity.brightness_check(np.stack([ki, dapi]), lab, seg_ch=0)["nuclear_brightest"]


def test_erode_rule():
    assert list(intensity.erode_amount(np.array([50, 150, 202, 2000]))) == [1, 1, 2, 2]  # 16 px across = 201.1 px²
    assert list(intensity.erode_amount(np.array([50, 2000]), 3)) == [3, 3]


def test_eroded_mask_excludes_bright_rim():
    lab, _ = _grid((120, 120), 40, 10)
    img = np.full(lab.shape, 100.0, np.float32)
    img[lab > 0] = 500
    rim =(lab > 0) & (np.pad(lab, 1)[2:, 1:-1] == 0) | (lab > 0) & (np.pad(lab, 1)[:-2, 1:-1] == 0) | \
        (lab > 0) & (np.pad(lab, 1)[1:-1, 2:] == 0) | (lab > 0) & (np.pad(lab, 1)[1:-1, :-2] == 0)
    img[rim] = 5000  # a bright halo on the outline only
    df, _ = intensity.measure_nuclei(img, lab, 65535, 20)
    assert (df["erode_px"] == 2).all() and (df["measured_px"] < df["area"] * 0.85).all()
    assert np.allclose(df["mean"], 500) and np.allclose(df["p90"], 500)
    tiny = _disks((40, 40), [(20, 20)], 2)  # 13 px: shrinking would leave too little, so it is lowered
    d2, _ = intensity.measure_nuclei(np.ones((40, 40), np.float32), tiny, 65535, 5)
    assert d2["erode_px"].iloc[0] < 1 or d2["measured_px"].iloc[0] >= 9


def test_local_background_follows_gradient_and_ignores_neighbours():
    shape = (160, 200)
    lab, cs = _grid(shape, 40, 9)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    bg = 200 + 2.0 * xx  # strong left-right gradient
    img = bg.astype(np.float32).copy()
    sig = RNG.uniform(50, 400, len(cs))
    for i, s_ in enumerate(sig, 1):
        img[lab == i] += s_
    df, info = intensity.measure_nuclei(img, lab, 65535, 18)
    true_bg = np.array([200 + 2.0 * x for _, x in cs])
    assert np.abs(df["bg"].to_numpy() - true_bg).max() < 6  # local, not the image median (~400)
    assert np.abs((df["mean"] - df["bg"]).to_numpy() - sig).max() < 6
    assert (df["bg_source"] == "ring").all()
    # crowded: neighbours almost touching -> ring too small -> wider ring or image background, never nuclei
    lab2, cs2 = _grid((100, 100), 13, 6)
    img2 = np.where(lab2 > 0, 5000.0, 100.0).astype(np.float32)
    d2, _ = intensity.measure_nuclei(img2, lab2, 65535, 12, bg_min=40)
    assert (d2["bg"] < 200).all() and set(d2["bg_source"]) <= {"ring", "wide ring", "image"}


def test_saturated_nuclei_flagged_and_excluded_but_counted():
    lab, cs = _grid((100, 100), 30, 8)
    img = np.full(lab.shape, 100.0, np.float32)
    img[lab == 1] = 1000
    img[tuple(np.array(cs[1]))] = 65535
    df, _ = intensity.measure_nuclei(img, lab, 65535, 16)
    assert df.loc[df["label"] == 2, "saturated"].item() and not df.loc[df["label"] == 1, "saturated"].item()
    nuc = df.assign(n="n1", treatment="A", image=1, channel="Ki67", category="nucleus")
    out = intensity.derive(nuc, {}, 3.0)
    assert out["counted"].all() and out["measurable"].sum() == len(out) - 1


def test_noise_estimate_handles_correlated_noise():
    from scipy import ndimage as ndi

    white = RNG.normal(0, 10, (400, 400))
    corr = ndi.gaussian_filter(RNG.normal(0, 10, (400, 400)), 1.0)
    m = np.ones(white.shape, bool)
    s1, a1 = intensity.background_noise(white.astype(np.float32), m, 20)
    s2, a2 = intensity.background_noise(corr.astype(np.float32), m, 20)
    assert abs(s1 - 10) < 1 and a1 < 1.5
    assert a2 > 4  # neighbouring pixels share noise: fewer independent pixels per nucleus


# ---------------------------------------------------------------- flat field

def _vignetted_stack(n_img=80, shape=(256, 320), seed=1):
    rng = np.random.default_rng(seed)
    H, W = shape
    yy, xx = np.mgrid[0:H, 0:W]
    V = 1 - 0.3 * (((yy - 0.5 * H) / (0.5 * H)) ** 2 + ((xx - 0.6 * W) / (0.5 * W)) ** 2) / 2
    V = V / V.mean()
    imgs = []
    for _ in range(n_img):
        s = np.zeros(shape)
        for y, x, a in zip(rng.uniform(-8, H + 8, 120), rng.uniform(-8, W + 8, 120), rng.lognormal(np.log(800), 0.4, 120)):
            y0, y1, x0, x1 = int(max(y - 8, 0)), int(min(y + 8, H)), int(max(x - 8, 0)), int(min(x + 8, W))
            s[y0:y1, x0:x1][(yy[y0:y1, x0:x1] - y) ** 2 + (xx[y0:y1, x0:x1] - x) ** 2 < 49] += a
        imgs.append(s * V + 400 + 80 * V + rng.normal(0, 20, shape))
    return imgs, V


def test_flatfield_recovers_synthetic_vignetting():
    imgs, V = _vignetted_stack()
    ff = flatfield.estimate(imgs)
    F = ff.flat()
    assert np.abs(F / V - 1).mean() < 0.02 and np.abs(F / V - 1).max() < 0.08
    assert abs(F.mean() - 1) < 0.01
    mag = flatfield.magnitude(ff)
    assert mag["darkest_side"] == "left" and 5 < mag["edge_dim_pct"] < 20
    assert "dimmer than the centre" in flatfield.describe(mag)
    d = flatfield.FlatField.from_dict(json.loads(json.dumps(ff.to_dict())))
    assert np.allclose(d.flat(np.array([3.0, 100.0]), np.array([5.0, 200.0])), ff.flat(np.array([3.0, 100.0]),
                                                                                          np.array([5.0, 200.0])))
    flat_imgs = [np.full((64, 64), 100.0) + RNG.normal(0, 1, (64, 64)) for _ in range(5)]
    assert flatfield.estimate(flat_imgs[:2]).source == "off"  # too few images: no correction


def test_flatfield_from_reference_images(tmp_path):
    _, V = _vignetted_stack(n_img=1)
    for i in range(3):
        tifffile.imwrite(tmp_path / f"Ki67_ref{i}.tif", (1000 * V + 100 + RNG.normal(0, 5, V.shape)).astype(np.uint16))
    tifffile.imwrite(tmp_path / "dark.tif", np.full(V.shape, 100, np.uint16))
    refs, darks = flatfield.reference_files(tmp_path, "Ki67")
    assert len(refs) == 3 and len(darks) == 1
    ff = flatfield.from_reference(refs, V.shape, darks)
    assert ff.source == "reference" and np.abs(ff.flat() / V - 1).max() < 0.03


# ---------------------------------------------------------------- positivity

def _background_only_nuclei(n_img=6, shape=(300, 300), corr=1.0, signal=0.0):
    from scipy import ndimage as ndi

    out = []
    for i in range(n_img):
        lab, cs = _grid(shape, 30, 8, jitter=3)
        noise = ndi.gaussian_filter(RNG.normal(0, 30, shape), corr) if corr else RNG.normal(0, 30, shape)
        yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
        img = (500 + 0.5 * xx + noise).astype(np.float32)
        if signal:
            img[lab > 0] += signal
        df, _ = intensity.measure_nuclei(img, lab, 65535, 16)
        out.append(df.assign(n="n1", treatment="A", image=i + 1, channel="Ki67", category="nucleus"))
    return intensity.derive(pd.concat(out, ignore_index=True), {}, 3.0)


@pytest.mark.parametrize("corr", [0, 1.5])
def test_lod_rule_false_positive_rate_on_pure_background(corr):
    nuc = _background_only_nuclei(corr=corr)
    s = settings.DEFAULTS
    dec = positivity.decide(nuc, s, "Ki67", "n1", [])
    assert dec["mean"]["rule"] == "lod" and dec["p90"]["rule"] == "lod"
    fp = (nuc["mean_bgsub_ff"].to_numpy() > dec["mean"]["thr"]).mean()
    fp90 = (nuc["p90_bgsub_ff"].to_numpy() > dec["p90"]["thr"]).mean()
    assert len(nuc) > 400 and fp < 0.01 and fp90 < 0.02
    stained = _background_only_nuclei(corr=corr, signal=60)  # 2 x the per-pixel noise SD
    d2 = positivity.decide(stained, s, "Ki67", "n1", [])
    assert (stained["mean_bgsub_ff"].to_numpy() > d2["mean"]["thr"]).mean() > 0.95


def _two_groups():
    a = pd.DataFrame(dict(treatment="secondary only", mean_bgsub_ff=RNG.normal(0, 10, 400), p90_bgsub_ff=RNG.normal(20, 10, 400)))
    b = pd.DataFrame(dict(treatment="Control", mean_bgsub_ff=RNG.normal(100, 30, 400), p90_bgsub_ff=RNG.normal(150, 30, 400)))
    d = pd.concat([a, b], ignore_index=True)
    return d.assign(thr_lod=5.0, thr90_lod=40.0)


def test_negative_control_percentile_cutoff():
    d = _two_groups()
    s = settings._merge(settings.DEFAULTS, {"experiment": {"negative_control": "Secondary_Only"}})
    dec = positivity.decide(d, s, "Ki67", "n1", [])
    ctrl = d["treatment"] == "secondary only"
    assert dec["mean"]["rule"] == "negative_control"
    assert np.isclose(dec["mean"]["value"], np.percentile(d.loc[ctrl, "mean_bgsub_ff"], 99))
    assert abs((d.loc[ctrl, "mean_bgsub_ff"] > dec["mean"]["value"]).mean() - 0.01) < 0.006
    assert positivity.is_negative_control("secondary-only", s) and not positivity.is_negative_control("Control", s)


def test_locked_manual_cutoff_precedence(tmp_path):
    d = _two_groups()
    d = d[d["treatment"] == "Control"]  # no control wells in this n
    s = settings.DEFAULTS
    assert positivity.decide(d, s, "Ki67", "n1", positivity.load_locks(tmp_path))["mean"]["rule"] == "lod"
    positivity.save_lock(tmp_path, "Ki67", 80.0, ["fa", "fb"])
    positivity.save_lock(tmp_path, "Ki67", 90.0, ["fc"], n="n2")
    locks = positivity.load_locks(tmp_path)
    dec = positivity.decide(d, s, "ki67", "n1", locks)
    assert dec["mean"]["rule"] == "manual" and dec["mean"]["value"] == 80.0 and "2 fields" in dec["mean"]["detail"]
    assert dec["p90"]["rule"] == "lod"  # a mean lock does not set the 90th-percentile cut-off
    assert positivity.decide(d, s, "Ki67", "n2", locks)["mean"]["value"] == 90.0  # n-specific lock wins
    s2 = settings._merge(settings.DEFAULTS, {"positivity": {"manual_cutoff": {"Ki67": 70}}})
    assert positivity.decide(d, s2, "Ki67", "n1", locks)["mean"]["value"] == 70.0  # settings beat the lock
    both = _two_groups()
    s3 = settings._merge(settings.DEFAULTS, {"experiment": {"negative_control": "secondary only"},
                                             "positivity": {"manual_cutoff": 70}})
    assert positivity.decide(both, s3, "Ki67", "n1", locks)["mean"]["rule"] == "negative_control"  # control wins
    data = json.loads((tmp_path / "cutoffs.json").read_text())
    assert data["locks"][0]["locked_at"] and data["locks"][0]["fields"] == ["fa", "fb"]


def test_sweep_is_monotone_and_includes_the_cutoff():
    d = _two_groups()
    sw = positivity.sweep(d, np.full(len(d), 50.0))
    for t, g in sw.groupby("treatment"):
        g = g.sort_values("cutoff")
        assert (np.diff(g["frac_positive"].to_numpy()) <= 1e-12).all()
    one = sw[np.isclose(sw["factor"], 1.0) & (sw["treatment"] == "all")]
    assert len(one) == 1 and np.isclose(one["cutoff"].item(), 50.0)
    assert np.isclose(one["frac_positive"].item(), (d["mean_bgsub_ff"] > 50).mean())


def test_registry_labels_channel_readouts():
    from kerncount import metrics

    m = metrics.describe_column("readout_frac_positive_Ki67")
    assert m.kind == "fraction" and "Ki67-positive" in m.label
    m90 = metrics.describe_column("readout_frac_positive_p90_Ki67")
    assert m90.kind == "fraction" and "90th" in m90.label
    assert metrics.describe_column("intensity_Ki67_nucleus_mean_bgsub").kind == "positive"
    from kerncount.classify import READOUTS, channel_metric_info

    assert "readout_frac_positive_Ki67" not in list(READOUTS)  # readouts.include lists only the static ones
    assert channel_metric_info("intensity_Ki67_positive_mean_bgsub")["kind"] == "positive"


# ---------------------------------------------------------------- metadata, plate maps, settings

def test_metadata_columns_with_layout_overrides():
    s = settings._merge(settings.DEFAULTS, {"experiment": {"cell_line": "LN229", "magnification": 20,
                                                           "seeding_density": 30000}})
    t = pd.DataFrame(dict(n=["n1", "n1", "n2"], treatment=["A", "B", "A"], image=[1, 1, 1]))
    lay = pd.DataFrame(dict(n=["n1", "n1", "n2"], treatment=["A", "B", "A"], n_folder=["x", "x", "y"],
                            cell_line=["", "", "U251"], seeding_density=["", "15000", ""]))
    out = design.assign(t, s, layout=lay)
    assert list(out.columns[-8:]) == design.metadata_columns()
    assert out["cell_line"].tolist() == ["LN229", "LN229", "U251"]
    assert out["meta_seeding_density"].tolist() == ["30000", "15000", "30000"]
    assert out["meta_magnification"].tolist() == ["20"] * 3 and out["meta_seeding_unit"].iloc[0] == "cells/well"
    assert design.output_columns(out) == design.COLUMNS + design.metadata_columns()
    empty = design.assign(t.iloc[0:0], s)
    assert "cell_line" in empty and "meta_stain" in empty
    assert design.assign(t, settings.DEFAULTS)["cell_line"].tolist() == ["", "", ""]


def test_pixel_size_from_camera_and_magnification():
    s = settings._merge(settings.DEFAULTS, {"experiment": {"camera_pixel_um": 6.5, "magnification": "20x"}})
    assert np.isclose(settings.pixel_size_um(s), 0.325)
    assert settings.pixel_size_um(settings._merge(s, {"segmentation": {"pixel_size_um": 0.5}})) == 0.5
    assert settings.pixel_size_um(settings.DEFAULTS) is None


def test_platemap_reserves_negative_control_wells(tmp_path):
    conds = ["Control", "LyoVec", "miR-a"]
    p = platemap.generate(conds, 24, 2, n_plates=2, seed=3, negative_control="secondary only", negative_control_wells=1)
    for m in p:
        assert (m["condition"] == "secondary only").sum() == 1 and (m["condition"] == "Control").sum() == 2
    assert p[0].attrs["conditions"][-1] == "secondary only"
    with pytest.raises(ValueError):
        platemap.generate(conds, 24, 2, negative_control="x", negative_control_wells=3)
    settings.write_template(tmp_path, tmp_path)
    s_path = tmp_path / "settings.toml"
    s_path.write_text(settings_form.update_text(s_path.read_text(), settings_form.coerce(
        {"experiment": {"negative_control": "secondary only"}})))
    out = platemap.plan_project(tmp_path, plates=1, conditions=conds, seed=4, reference_wells=2)
    m = pd.read_csv(out / "plate_1.csv")
    assert (m["condition"] == "secondary only").sum() == 2
    assert "negative control 'secondary only'" in (out / "README.txt").read_text()


def test_settings_sections_and_form():
    keys = list(settings.DEFAULTS)
    assert keys.index("experiment") == keys.index("project") + 1
    for sec in ("experiment", "channels", "intensity", "positivity"):
        assert set(settings_form.FIELDS[sec]) == set(settings.DEFAULTS[sec])
    assert [x["section"] for x in settings_form.schema()][:3] == ["project", "experiment", "channels"]
    typed = settings_form.coerce({"channels": {"files": "interleaved", "file_order": "Ki67, DAPI", "segment": "DAPI"},
                                  "experiment": {"cell_line": "LN229", "seeding_density": "30000"}})
    out = settings_form.update_text(settings.TEMPLATE.format(raw_dir="x"), typed)
    import tomllib

    d = tomllib.loads(out)
    assert d["channels"]["file_order"] == ["Ki67", "DAPI"] and d["experiment"]["seeding_density"] == 30000
    s = settings._merge(settings.DEFAULTS, d)
    assert settings.channel_index(s, "dapi", 2) == 1 and settings.channel_index(s, 1, 2) == 0
    assert settings.channel_names(s, 2) == ["Ki67", "DAPI"]


def test_single_file_segmentation_hash_unchanged_by_new_settings():
    import hashlib

    legacy = {"segmentation": {**settings.DEFAULTS["segmentation"], "pixel_size_um": ""},
              "channels": {"names": [], "segment": 1},
              "intensity": {"enabled": True, "measure": "all", "regions": ["nucleus", "cell"], "ring_px": "auto",
                            "subtract_background": True}}
    blob = json.dumps(legacy, sort_keys=True, default=str)
    old = hashlib.sha1(blob.encode()).hexdigest()[:12]  # the formula used before this change
    new = settings.section_hash(settings.seg_hash_input(settings.DEFAULTS), "segmentation", "channels", "intensity")
    assert new == old
    changed = settings._merge(settings.DEFAULTS, {
        "intensity": {"flatfield": "off", "erode_px": 2, "bg_ring_px": 9, "regions": ["nucleus", "cell", "fiji_cell"]},
        "positivity": {"lod_k": 4}, "experiment": {"cell_line": "X", "magnification": 20},
        "channels": {"align": "correct"}})
    assert settings.section_hash(settings.seg_hash_input(changed), "segmentation", "channels", "intensity") == old
    inter = settings._merge(settings.DEFAULTS, {"channels": {"files": "interleaved", "file_order": ["K", "D"]}})
    assert settings.section_hash(settings.seg_hash_input(inter), "segmentation", "channels", "intensity") != old


def test_single_channel_project_writes_no_intensity_outputs(tmp_path):
    class P:
        settings = settings.DEFAULTS
        path = tmp_path
        work = tmp_path / "_work"

        def seg_channel(self, n):
            return 0

    fields = [(None, None, pd.DataFrame(), {"n_channels": 1}, {})]
    assert intensity.process(P(), fields) is None
    assert not (tmp_path / "intensity").exists()


# ---------------------------------------------------------------- end to end (classical detector) + tuning API

def _synthetic_interleaved(raw: Path):
    """2 conditions x 3 fields; marker first, DAPI second. 'hi' nuclei are 80% stained, 'lo' 20%."""
    rng = np.random.default_rng(5)
    shape = (256, 256)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    for cond, frac in (("cntrl", 0.8), ("lyovec", 0.2)):
        folder = raw / "05.06.26 test" / cond
        folder.mkdir(parents=True)
        k = 1
        for f in range(3):
            dapi = np.full(shape, 300.0)
            ki = np.full(shape, 200.0) * (1 - 0.15 * (xx / shape[1]))  # dimmer to the right
            for y in range(24, 240, 36):
                for x in range(24, 240, 36):
                    cy, cx = y + rng.integers(-4, 5), x + rng.integers(-4, 5)
                    m = (yy - cy) ** 2 + (xx - cx) ** 2 <= 81
                    dapi[m] += 3000
                    if rng.uniform() < frac:
                        ki[m] += 400
            if cond == "lyovec" and f == 2:  # a marker file from another field: must be caught as mis-paired
                ki = np.full(shape, 200.0)
                for _ in range(40):
                    cy, cx = rng.uniform(10, 246, 2)
                    ki[(yy - cy) ** 2 + (xx - cx) ** 2 <= 81] += 400
            for img in (ki, dapi):
                img = (img + rng.normal(0, 15, shape)).clip(0).astype(np.uint16)
                tifffile.imwrite(folder / f"{k}.tif", np.stack([img] * 3, -1), photometric="rgb")
                k += 1


@pytest.fixture(scope="module")
def ki_project(tmp_path_factory):
    from kerncount import pipeline

    root = tmp_path_factory.mktemp("ki")
    raw, prj = root / "raw", root / "prj"
    _synthetic_interleaved(raw)
    prj.mkdir()
    settings.write_template(prj, raw)
    sp = prj / "settings.toml"
    sp.write_text(settings_form.update_text(sp.read_text(), settings_form.coerce({
        "channels": {"files": "interleaved", "file_order": "Ki67, DAPI", "segment": "DAPI"},
        "segmentation": {"backend": "classical"}, "experiment": {"cell_line": "LN229"},
        "intensity": {"regions": ["nucleus", "cell", "fiji_cell"]}})))
    s = settings.load(prj)
    discover.write_layout(discover.scan(raw, files="interleaved", n_channels=2), prj / "layout.csv")
    pipeline.run(prj, no_review=True)
    return prj


def test_end_to_end_interleaved_project(ki_project):
    prj = ki_project
    qc = pd.read_csv(prj / "qc_details.csv")
    assert len(qc) == 6 and (qc["auto_count"] >= 30).all()
    for c in ("readout_frac_positive_Ki67", "readout_frac_positive_p90_Ki67", "intensity_Ki67_nucleus_mean_bgsub",
              "intensity_Ki67_positive_mean_bgsub", "intensity_Ki67_fiji_cell", "n_measured_Ki67", "channel_check"):
        assert c in qc, c
    by = qc.groupby("treatment")["readout_frac_positive_Ki67"].mean()
    assert abs(by["Control"] - 0.8) < 0.15 and abs(by["LyoVec"] - 0.2) < 0.15
    cut = pd.read_csv(prj / "intensity" / "cutoffs.csv")
    assert set(cut["rule"]) == {"lod"} and set(cut["measure"]) == {"mean", "p90"}
    man = json.loads((prj / "intensity" / "manifest.json").read_text())
    paths = {f["path"] for f in man["files"]}
    assert all((prj / "intensity" / p).exists() for p in paths)
    for kind in ("cutoffs", "sweep", "distribution", "density", "alignment", "flatfield_table", "summary"):
        assert any(f["kind"] == kind for f in man["files"]), kind
    al = pd.read_csv(prj / "intensity" / "alignment.csv")
    bad = (al["treatment"] == "LyoVec") & (al["image"] == 3)
    assert (al.loc[~bad, "shift_px"] < 1).all() and (al.loc[~bad, "corr"] > 0.3).all()
    assert al.loc[bad, "corr"].item() < 0.2 and "mis-paired" in al.loc[bad, "flags"].item()
    row = qc[(qc["treatment"] == "LyoVec") & (qc["image"] == 3)].iloc[0]
    assert np.isnan(row["readout_frac_positive_Ki67"]) and row["n_measured_Ki67"] == 0 and row["auto_count"] > 30
    ff = pd.read_csv(prj / "intensity" / "flatfield.csv")
    assert ff["source"].tolist() == ["off"]  # 6 fields are too few to estimate an illumination pattern
    nuc = pd.read_csv(prj / "_work" / "nuclei_intensity.csv")
    assert {"mean", "total", "p90", "area", "bg", "saturated", "mean_bgsub_ff", "positive"} <= set(nuc.columns)
    used = sorted((prj / "generations").glob("g*/counts_used.csv"))
    assert used


def test_tuning_api_blinded_and_fast(ki_project):
    from kerncount import tuning

    prj = ki_project
    assert tuning.channels(prj) == ["Ki67"]
    smp = tuning.sample(prj, k=4, seed=1)
    ids = [f["field_id"] for f in smp["fields"]]
    assert len(ids) == 4 and len(set(ids)) == 4
    assert tuning.sample(prj, k=4, seed=1) == smp
    blob = json.dumps(smp) + json.dumps(tuning.nuclei(prj, ids[0])) + json.dumps(tuning.preview(prj, 50.0)) \
        + json.dumps(tuning.suggest_range(prj))
    for leak in ("Control", "LyoVec", "cntrl", "lyovec", ".tif"):
        assert leak not in blob
    nuc = tuning.nuclei(prj, ids[0])
    assert nuc["nuclei"] and len(nuc["nuclei"][0]["contour"]) >= 4
    png = tuning.render(prj, ids[0], "Ki67", cutoff=100.0)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    t0 = time.perf_counter()
    for c in (10.0, 50.0, 100.0, 200.0):
        tuning.render(prj, ids[0], "Ki67", cutoff=c)
        tuning.preview(prj, c, "Ki67")
    assert (time.perf_counter() - t0) / 4 < 0.2
    lo, hi = tuning.preview(prj, 1e9)["frac_positive"], tuning.preview(prj, -1e9)["frac_positive"]
    assert lo == 0 and hi == 1
    with pytest.raises(tuning.TuningError):
        tuning.render(prj, "f_nope", "Ki67", 1.0)
    lock = tuning.lock(prj, "Ki67", 150.0, ids)
    assert lock["cutoff"] == 150.0 and lock["fields"] == ids and lock["locked_at"]
    assert positivity.find_lock(positivity.load_locks(prj), "Ki67", "n1_2026-06-05")["cutoff"] == 150.0
    (prj / "cutoffs.json").unlink()
