"""Per-image readouts from classified objects (classify.compute_readouts): synthetic object tables, known answers."""
import math

import numpy as np
import pandas as pd
import pytest

from kerncount import classify, settings, settings_form


def obj(category, area, x, y, solidity=0.97, edge=False):
    return dict(category=category, area=float(area), x=float(x), y=float(y), solidity=solidity,
                edge_top=edge, edge_left=False, edge_bottom=False, edge_right=False)


def grid(n, spacing, area, category="nucleus", **kw):
    """n nuclei on a square grid with the given centre spacing (all well inside the image)."""
    side = math.ceil(math.sqrt(n))
    return [obj(category, area, spacing * (1 + i % side), spacing * (1 + i // side), **kw) for i in range(n)]


def run(rows, shape=(1000, 1000), **kw):
    return classify.compute_readouts(pd.DataFrame(rows), shape, **kw)


def test_registry_contract():
    assert classify.READOUTS
    for name, meta in classify.READOUTS.items():
        assert name.startswith("readout_")
        assert set(meta) == {"label", "kind", "unit", "description", "min_nuclei"}
        assert meta["kind"] in ("fraction", "positive", "count")
        assert isinstance(meta["min_nuclei"], int) and meta["min_nuclei"] >= 0
        assert (meta["kind"] == "fraction") == name.startswith("readout_frac_"), name
    # the stdlib-only settings form keeps its own copy of the names: it must match
    assert ["readout_" + n for n in settings_form.READOUT_NAMES] == list(classify.READOUTS)


def test_known_values():
    rows = grid(16, 100, 400)  # 16 whole nuclei, area 400, nearest neighbour 100 px
    rows += [obj("out_of_focus", 900, 800, 800)] * 4
    rows += [obj("debris", 30, 100 + 15, 100)]  # small, touching a nucleus -> micronucleus candidate
    rows += [obj("debris", 30, 900, 100)]  # small but far from any nucleus -> dust
    rows += [obj("debris", 70, 200, 215)]  # near, but too big (> 1/9 of a nucleus) -> fragment, not micronucleus
    ro, df = run(rows)
    assert ro["readout_frac_out_of_focus"] == pytest.approx(4 / 20)
    assert ro["readout_micronuclei_per_nucleus"] == pytest.approx(1 / 16)
    assert ro["readout_nuclear_area_median"] == pytest.approx(400)
    assert ro["readout_nuclear_area_cv"] == pytest.approx(0)
    assert ro["readout_frac_irregular"] == 0
    assert ro["readout_coverage"] == pytest.approx(16 * 400 / 1e6)
    diam = math.sqrt(4 * 400 / math.pi)
    assert ro["readout_nn_distance_median"] == pytest.approx(100 / diam, rel=1e-4)
    assert "readout_nuclear_area_median_um2" not in ro  # no pixel size known
    assert df["micronucleus_candidate"].sum() == 1
    assert df.loc[df["category"] == "nucleus", "nn_distance_px"].eq(100).all()
    assert df.loc[df["category"] != "nucleus", "nn_distance_px"].isna().all()


def test_size_shape_use_whole_nuclei_only():
    rows = grid(10, 100, 400)
    rows += [obj("nucleus", 4000, 50, 990, edge=True)]  # cut by the border: counted, but not a whole nucleus
    rows += [obj("merged", 5000, 700, 700, solidity=0.8)]  # suspected merge: not a single nucleus
    rows[0]["solidity"] = rows[1]["solidity"] = 0.85  # two irregular nuclei
    rows[2]["area"] = 800
    ro, _ = run(rows)
    areas = np.array([800] + [400] * 9, float)
    assert ro["readout_nuclear_area_median"] == 400
    assert ro["readout_nuclear_area_cv"] == pytest.approx(areas.std(ddof=1) / areas.mean(), rel=1e-4)
    assert ro["readout_frac_irregular"] == pytest.approx(0.2)
    # coverage counts every in-focus nucleus, including cut and merged ones
    assert ro["readout_coverage"] == pytest.approx((areas.sum() + 4000 + 5000) / 1e6)


def test_pixel_size_adds_um2():
    ro, _ = run(grid(9, 100, 400), pixel_size_um=0.5)
    assert ro["readout_nuclear_area_median_um2"] == pytest.approx(100)


@pytest.mark.parametrize("rows", [[], [obj("nucleus", 400, 100, 100)], [obj("debris", 20, 5, 5)]])
def test_empty_or_single_nucleus_gives_nan_not_crash(rows):
    ro, df = run(rows)
    assert set(ro) == set(classify.READOUTS) - {"readout_nuclear_area_median_um2"}
    for k, v in ro.items():
        if k == "readout_coverage":
            assert v == pytest.approx(sum(r["area"] for r in rows if r["category"] == "nucleus") / 1e6)
        else:
            assert np.isnan(v), k
    assert len(df) == len(rows)


def test_classify_empty_output_is_accepted():
    df, _ = classify.classify(pd.DataFrame(), {}, classify.FolderRef(np.nan, np.nan), settings.DEFAULTS)
    ro, _ = classify.compute_readouts(df, (100, 100))
    assert ro["readout_coverage"] == 0 and np.isnan(ro["readout_frac_out_of_focus"])


def test_min_nuclei_thresholds():
    ro, _ = run(grid(9, 100, 400))  # 9 nuclei: enough for median area / NN (5), not for CV / fractions (10)
    assert not np.isnan(ro["readout_nuclear_area_median"]) and not np.isnan(ro["readout_nn_distance_median"])
    assert np.isnan(ro["readout_nuclear_area_cv"]) and np.isnan(ro["readout_frac_irregular"])
    assert np.isnan(ro["readout_micronuclei_per_nucleus"])
    # out-of-focus fraction needs 10 nuclear objects in total, in focus or not: all-floating fields still score
    ro, _ = run([obj("out_of_focus", 900, 100 * i, 500) for i in range(1, 11)])
    assert ro["readout_frac_out_of_focus"] == 1.0


def test_magnification_invariance():
    """The same scene imaged at 2x magnification: relative readouts are unchanged, px² area scales by 4."""
    def scene(z):
        rows = grid(25, 90 * z, 500 * z * z)
        rows[3]["solidity"] = 0.8
        rows += [obj("out_of_focus", 1200 * z * z, 950 * z, 950 * z)] * 3
        rows += [obj("debris", 40 * z * z, 90 * z + 14 * z, 90 * z)]
        return run(rows, shape=(1000 * z, 1000 * z))[0]

    a, b = scene(1), scene(2)
    for k in a:
        if k == "readout_nuclear_area_median":
            assert b[k] == pytest.approx(4 * a[k])
        else:
            assert b[k] == pytest.approx(a[k], rel=1e-3), k


def test_include_subset_and_prefix():
    ro, _ = run(grid(12, 100, 400), include=["coverage", "readout_frac_irregular"])
    assert set(ro) == {"readout_coverage", "readout_frac_irregular"}
    assert classify.readout_names("all") == list(classify.READOUTS)


def test_settings_defaults_and_form():
    keys = list(settings.DEFAULTS)
    assert keys[keys.index("intensity") + 1] == "readouts"
    assert settings.DEFAULTS["readouts"] == {"enabled": True, "include": "all"}
    fkeys = list(settings_form.FORM_ITEMS)
    assert fkeys[fkeys.index("intensity") + 1] == "readouts"
    assert "readouts" in settings_form.SECTION_TITLES
    import tomllib
    assert tomllib.loads(settings.TEMPLATE.format(raw_dir="x"))["readouts"] == settings.DEFAULTS["readouts"]
    typed = settings_form.coerce({"readouts": {"enabled": False, "include": "coverage, readout_frac_out_of_focus"}})
    assert typed["readouts"] == {"enabled": False, "include": ["frac_out_of_focus", "coverage"]}
    assert settings_form.coerce({"readouts": {"include": " all "}})["readouts"]["include"] == "all"
    assert settings_form.coerce({"readouts": {"include": ""}})["readouts"]["include"] == "all"
    with pytest.raises(ValueError, match="Readouts to include"):
        settings_form.coerce({"readouts": {"include": "coverage, nonsense"}})
    text = settings.TEMPLATE.format(raw_dir="x")
    out = settings_form.update_text(text, typed)
    assert tomllib.loads(out)["readouts"]["include"] == ["frac_out_of_focus", "coverage"]


def test_readout_settings_do_not_change_segmentation_hash():
    a = settings.section_hash(settings.DEFAULTS, "segmentation", "channels", "intensity")
    s = {**settings.DEFAULTS, "readouts": {"enabled": False, "include": ["coverage"]}}
    assert settings.section_hash(s, "segmentation", "channels", "intensity") == a
