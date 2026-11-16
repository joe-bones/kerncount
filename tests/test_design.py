"""Plate design: well / position / timepoint assignment (kerncount.design) and its wiring."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import analysis, design, discover, settings_form
from kerncount import settings as st_mod


def _settings(**d):
    s = st_mod.load(Path("/nonexistent"))
    s["design"].update(d)
    return s


def _table(spec: dict) -> pd.DataFrame:
    """{(n, treatment): [image numbers]} -> rows."""
    return pd.DataFrame([dict(n=n, treatment=t, image=i) for (n, t), imgs in spec.items() for i in imgs])


def _wp(out, n="n1", t="A"):
    g = out[(out["n"] == n) & (out["treatment"] == t)]
    return list(zip(g["image"], g["well"], g["position"]))


def _row(desc, n="n1", t="A"):
    return desc[(desc["n"] == n) & (desc["treatment"] == t)].iloc[0]


def test_defaults_are_two_wells_auto_fields_fixed_positions():
    sp = design.spec(_settings())
    assert (sp.wells, sp.fields_per_well, sp.well_layout, sp.field_positions) == (2, "auto", "consecutive", "fixed")
    assert design.positions_meaningful(_settings())
    assert not design.positions_meaningful(_settings(field_positions="random"))


def test_standard_twelve_images_consecutive_and_alternating():
    t = _table({("n1", "A"): range(1, 13)})
    out = design.assign(t, _settings())
    assert list(out["well"]) == [1] * 6 + [2] * 6 and list(out["position"]) == list(range(1, 7)) * 2
    assert list(out["timepoint"]) == [""] * 12
    alt = design.assign(t, _settings(well_layout="alternating"))
    assert list(alt["well"]) == [1, 2] * 6 and list(alt["position"]) == [p for p in range(1, 7) for _ in (1, 2)]
    d = _row(design.describe(out))
    assert (d["images"], d["wells_present"], d["fields_per_well"], d["flags"]) == (12, 2, 6, "")
    assert not d["unbalanced"] and not d["single_well"]


def test_gap_in_image_numbers_keeps_positions_and_flags_unbalanced():
    # the real n1 / CL075: image 8 missing -> still 2 wells x 6 (auto), well 2 has 5 fields, position 2 empty
    out = design.assign(_table({("n1", "CL075"): [i for i in range(1, 13) if i != 8]}), _settings())
    assert _wp(out, t="CL075")[6:8] == [(7, 2, 1), (9, 2, 3)]
    d = _row(design.describe(out), t="CL075")
    assert d["fields_per_well"] == 6 and d["missing_images"] == "8" and d["unbalanced"]
    assert (d["fields_min"], d["fields_max"]) == (5, 6) and "unbalanced" in d["flags"]
    line = design.summary_line(design.describe(out), _settings())
    assert "2 wells x 6 fields (auto)" in line and "1 condition unbalanced (n1 / CL075: 11 images)" in line


def test_eleven_images_last_missing_is_still_six_per_well():
    d = _row(design.describe(design.assign(_table({("n1", "A"): range(1, 12)}), _settings())))
    assert d["fields_per_well"] == 6 and d["missing_images"] == "12" and d["unbalanced"]


@pytest.mark.parametrize("n_img, wells, k", [(18, 2, 9), (24, 2, 12), (18, 3, 6), (6, 2, 3)])
def test_auto_fields_adapt_to_image_count(n_img, wells, k):
    out = design.assign(_table({("n1", "A"): range(1, n_img + 1)}), _settings(wells=wells))
    assert out.groupby("well").size().tolist() == [k] * wells
    assert out["position"].max() == k
    d = _row(design.describe(out))
    assert d["fields_per_well"] == k and d["flags"] == ""


def test_fixed_fields_images_beyond_design_stay_in_last_well_and_are_flagged():
    out = design.assign(_table({("n1", "A"): range(1, 15)}), _settings(fields_per_well=6))
    assert _wp(out)[-3:] == [(12, 2, 6), (13, 2, 7), (14, 2, 8)]
    d = _row(design.describe(out))
    assert d["beyond_design"] == "13; 14" and "beyond the design" in d["flags"] and d["fields_source"] == "fixed"
    alt = design.assign(_table({("n1", "A"): range(1, 15)}), _settings(fields_per_well=6, well_layout="alternating"))
    assert _wp(alt)[-2:] == [(13, 1, 7), (14, 2, 7)]
    assert "positions above 6" in _row(design.describe(alt))["flags"]


def test_single_image_and_single_well():
    d = design.describe(design.assign(_table({("n1", "A"): [1], ("n1", "B"): range(1, 7)}), _settings()))
    a, b = _row(d, t="A"), _row(d, t="B")
    assert a["images"] == 1 and "single image" in a["flags"] and a["single_well"]
    assert b["wells_present"] == 2 and b["fields_per_well"] == 3  # 6 images, 2 wells -> 3 + 3
    one = design.assign(_table({("n1", "A"): range(1, 13)}), _settings(wells=1))
    assert set(one["well"]) == {1} and list(one["position"]) == list(range(1, 13))
    r = _row(design.describe(one))
    assert r["single_well"] and not r["unbalanced"] and r["flags"] == "single well"


def test_conditions_in_one_n_with_different_image_counts():
    out = design.assign(_table({("n1", "A"): range(1, 13), ("n1", "B"): range(1, 19), ("n2", "A"): range(1, 11)}),
                        _settings())
    d = design.describe(out)
    assert list(d["fields_per_well"]) == [6, 9, 5]
    assert _wp(out, t="B")[9] == (10, 2, 1)
    assert "5-9 fields" in design.summary_line(d, _settings())


def test_rows_are_not_dropped_or_reordered_and_index_kept():
    t = _table({("n1", "A"): range(1, 13), ("n1", "B"): range(1, 13)}).sample(frac=1, random_state=0)
    t.index = t.index * 10 + 3
    t["value"] = np.arange(len(t))
    out = design.assign(t, _settings())
    assert out.index.equals(t.index) and (out["value"] == t["value"]).all() and (out["image"] == t["image"]).all()
    assert all(((i - 1) // 6 + 1, (i - 1) % 6 + 1) == (w, p) for i, w, p in zip(out["image"], out["well"], out["position"]))
    assert "well" not in t  # input untouched


def test_old_names_still_work():
    s = _settings()
    s["design"].update(plates=3, fields_per_plate=4, plate_layout="alternating")
    sp = design.spec(s)
    assert (sp.wells, sp.fields_per_well, sp.well_layout) == (3, 4, "alternating") and len(sp.notes) == 3
    out = design.assign(_table({("n1", "A"): range(1, 13)}), s)
    assert list(out["well"][:4]) == [1, 2, 3, 1]


@pytest.mark.parametrize("bad", [dict(wells=0), dict(wells="two"), dict(well_layout="random"),
                                 dict(field_positions="maybe"), dict(fields_per_well="x")])
def test_invalid_settings_raise(bad):
    with pytest.raises(ValueError):
        design.spec(_settings(**bad))


def test_wells_file_most_specific_rule_wins(tmp_path):
    (tmp_path / "wells.csv").write_text(
        "n,treatment,image,well,position\n"
        + "".join(f",,{i},{1 if i % 2 else 2},{(i + 1) // 2}\n" for i in range(1, 13))  # generic: odd/even
        + "n1,,1,2,9\n"            # n prefix matches n1_2026-01-01
        + ",B,1,3,1\n"             # treatment only
        + "n1,B,1,1,4\n"           # n + treatment beats both
        + "n2,A,2,3,\n")           # blank position: numbered in image order within the well
    s = _settings(well_layout="file")
    t = _table({("n1_2026-01-01", "A"): range(1, 14), ("n1_2026-01-01", "B"): [1], ("n2", "B"): [1], ("n2", "A"): [1, 2]})
    out = design.assign(t, s, tmp_path)
    assert _wp(out, "n1_2026-01-01", "A")[:3] == [(1, 2, 9), (2, 2, 1), (3, 1, 2)]
    assert _wp(out, "n1_2026-01-01", "A")[-1] == (13, 2, 6)  # not in the file (13 images -> 7 per well): consecutive fallback, flagged
    assert _wp(out, "n1_2026-01-01", "B") == [(1, 1, 4)]
    assert _wp(out, "n2", "B") == [(1, 3, 1)]
    assert _wp(out, "n2", "A") == [(1, 1, 1), (2, 3, 1)]
    d = design.describe(out)
    assert "13 not in wells.csv" in _row(d, "n1_2026-01-01", "A")["flags"]
    assert _row(d, "n2", "A")["fields_source"] == "file"


def test_wells_file_missing_or_bad(tmp_path):
    t = _table({("n1", "A"): range(1, 13)})
    with pytest.raises(ValueError, match="wells.csv"):
        design.assign(t, _settings(well_layout="file"), tmp_path)
    (tmp_path / "wells.csv").write_text("image,position\n1,1\n")
    with pytest.raises(ValueError, match="well"):
        design.assign(t, _settings(well_layout="file"), tmp_path)


def test_wells_auto_uses_file_if_present_else_two(tmp_path):
    t = _table({("n1", "A"): range(1, 13)})
    assert design.assign(t, _settings(wells="auto"), tmp_path)["well"].max() == 2
    (tmp_path / "wells.csv").write_text("image,well,position\n" + "".join(f"{i},{(i - 1) // 4 + 1},{(i - 1) % 4 + 1}\n" for i in range(1, 13)))
    out = design.assign(t, _settings(wells="auto"), tmp_path)
    assert out.groupby("well").size().tolist() == [4, 4, 4]
    assert _row(design.describe(out))["wells_expected"] == 3


def test_template_round_trip(tmp_path):
    t = _table({("n1", "A"): range(1, 13), ("n1", "B"): range(1, 19), ("n2", "A"): [1, 2, 3, 5]})
    s = _settings()
    path = design.write_wells_template(tmp_path, t, s)
    assert path.name == "wells.csv"
    tmpl = pd.read_csv(path, keep_default_na=False)
    assert (tmpl["n"] == "").sum() == 18  # one generic row per image number ...
    assert len(tmpl) < 34  # ... plus rows only where a condition differs
    before = design.assign(t, s)
    after = design.assign(t, _settings(well_layout="file"), tmp_path)
    assert (before[["well", "position"]] == after[["well", "position"]]).all().all()
    assert design.write_wells_template(tmp_path, t, s).name == "wells.new.csv"  # never overwrites by default


@pytest.mark.parametrize("name, tp", [
    ("06.10.26 U251 dapi count, MH, 72h, output", "72h"), ("05.06.26,U251,DAPI count, MH, output", ""),
    ("exp 48 h", "48h"), ("24hr run", "24h"), ("U251 72 hrs", "72h"), ("day 3", "d3"), ("U251_d5", "d5"), ("x_72h_y", "72h"),
    ("2 days", "d2"), ("U251 dapi", ""), ("29.09.26", "")])
def test_parse_timepoint(name, tp):
    assert discover.parse_timepoint(name) == tp


def test_timepoint_from_layout_old_and_new(tmp_path):
    t = _table({("n1", "A"): [1, 2], ("n2", "A"): [1]})
    new = pd.DataFrame(dict(n=["n1", "n2"], n_folder=["a 72h", "b"], timepoint=["48h", ""]))
    assert list(design.assign(t, _settings(), layout=new)["timepoint"]) == ["48h", "48h", ""]  # column wins
    old = new.drop(columns="timepoint")  # layout.csv written before the column existed: parse the folder
    assert list(design.assign(t, _settings(), layout=old)["timepoint"]) == ["72h", "72h", ""]
    old.assign(treatment_folder="x").to_csv(tmp_path / "layout.csv", index=False)
    assert list(design.assign(t, _settings(), tmp_path)["timepoint"]) == ["72h", "72h", ""]
    kept = design.assign(t.assign(timepoint="d2"), _settings())  # no layout at all: an existing column is kept
    assert set(kept["timepoint"]) == {"d2"}


def test_scan_writes_timepoint(tmp_path):
    for folder in ("06.10.26 U251, 72h, output", "05.06.26 U251"):
        d = tmp_path / folder / "cntrl"
        d.mkdir(parents=True)
        (d / "1.tif").write_bytes(b"")
    lay = discover.scan(tmp_path)
    assert "timepoint" in lay and dict(zip(lay["n"], lay["timepoint"])) == {"n1_2026-06-05": "", "n2_2026-10-06": "72h"}


def test_per_n_uses_well_column(tmp_path):
    rng = np.random.default_rng(0)
    df = pd.DataFrame([dict(n="n1", treatment=t, image=i, value=rng.normal(50, 5))
                       for t in ("Control", "LyoVec", "X") for i in range(1, 13)])
    analysis.per_n(df.assign(well=1), "count", "nuclei", "n1", tmp_path, _settings())
    assert not (tmp_path / "count_twoway_anova.csv").exists()  # one well: no treatment x well ANOVA
    analysis.per_n(df.assign(well=[1 + i % 3 for i in df["image"]]), "count", "nuclei", "n1", tmp_path, _settings())
    aov = pd.read_csv(tmp_path / "count_twoway_anova.csv")
    assert aov.loc[aov["source"] == "well", "df"].iloc[0] == 2  # three wells from the column, not the settings


def test_generation_counts_used_has_design_columns(tmp_path):
    from kerncount import generations
    rng = np.random.default_rng(0)
    rows = [dict(n=n, treatment=t, image=i, count=int(rng.integers(40, 60)), source="auto", excluded=False)
            for n in ("n1", "n2") for t in ("Control", "LyoVec", "X") for i in range(1, 13) if not (t == "X" and i == 8)]
    tab = pd.DataFrame(rows)
    pd.DataFrame(dict(n=["n1", "n2"], n_folder=["a", "b 72h"])).to_csv(tmp_path / "layout.csv", index=False)
    gen = dict(name="g1", table=tab, source="auto", changes=0, fingerprint="x")
    out = generations.build(tmp_path, gen, tab[["n", "treatment", "image"]].assign(qc_level="high"), {}, _settings())
    used = pd.read_csv(out / "counts_used.csv", keep_default_na=False)
    assert list(used.columns) == ["n", "treatment", "image", "count", "source", "excluded", "qc_level",
                                  "well", "position", "timepoint"]
    assert used.loc[(used["n"] == "n2") & (used["image"] == 7), "well"].tolist() == [2, 2, 2]
    assert set(used.loc[used["n"] == "n2", "timepoint"]) == {"72h"}
    assert list(pd.read_excel(out / "results.xlsx", sheet_name="Counts used").columns)[-3:] == ["well", "position", "timepoint"]


def test_settings_form_design_fields():
    typed = settings_form.coerce({"design": dict(wells="auto", fields_per_well="auto", well_layout="file",
                                                 field_positions="random")})
    assert typed["design"] == dict(wells="auto", fields_per_well="auto", well_layout="file", field_positions="random")
    assert settings_form.coerce({"design": dict(fields_per_well="9", wells="3")})["design"] == dict(fields_per_well=9, wells=3)
    with pytest.raises(ValueError):
        settings_form.coerce({"design": dict(field_positions="sometimes")})
    assert set(settings_form.FIELDS["design"]) == set(st_mod.DEFAULTS["design"])
