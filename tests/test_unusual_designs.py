"""Unusual designs (docs/design.md, "Unusual designs"): every case runs end to end without a traceback, and what
drops out of the statistics is said in the outputs and the browser, never silently.

Fast tests check each rule on small tables; the end-to-end runs (marked slow) push synthetic projects through the
whole chain: discover -> classify (case a: real images) -> design check -> models / stats -> readout ladder ->
duplicates -> sampling + planner -> plate maps -> the GUI endpoints (tests/unusual.py)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import design, discover, generations, gui, intensity, metrics, models, settings
from unusual import chain, image_project, project, rows, run_chain

S = settings._merge(settings.DEFAULTS, {"analysis": {"reference": "LyoVec", "normalise_to": ["LyoVec", "Control"]}})


def _assigned(data):
    return design.assign(pd.DataFrame(data)[["n", "treatment", "image"]], S)


# ---------------------------------------------------------------- a. interleaved channel files

def test_a_interleaved_gaps_odd_count_and_one_channel_folder(tmp_path):
    raw = tmp_path / "raw" / "05.06.26 x"
    for folder, numbers in (("cntrl", [1, 2, 5, 6, 9, 10, 11, 12]),  # gaps: still 4 pairs
                            ("lyovec", range(1, 10)),                 # odd: 9.tif unpaired
                            ("21", [])):                              # empty
        (raw / folder).mkdir(parents=True)
        for k in numbers:
            (raw / folder / f"{k}.tif").write_bytes(b"x")
    lay = discover.scan(tmp_path / "raw", files="interleaved", n_channels=2)
    got = dict(zip(lay["treatment_folder"], lay["images"]))
    assert got == {"cntrl": 4, "lyovec": 4, "21": 0}
    note = dict(zip(lay["treatment_folder"], lay["note"]))
    assert "leftover 9.tif skipped" in note["lyovec"] and "no images found" in note["21"]
    s = settings._merge(settings.DEFAULTS, {"channels": {"files": "interleaved", "file_order": ["Ki67", "DAPI"]}})
    recs, warns = discover.records(lay, tmp_path / "raw", s)
    assert [r.numbers for r in recs if r.treatment == "Control"] == [(1, 2), (5, 6), (9, 10), (11, 12)]
    assert any("leftover" in w and "9.tif" in w for w in warns)


def _fi(corr=0.9, brightest=True, shift=0.1):
    return {"alignment": {"Ki67": {"corr": corr, "shift_px": shift, "dy": 0, "dx": 0}},
            "brightness": {"nuclear_brightest": brightest}}


class _R:
    def __init__(self, n, t, i):
        self.n, self.treatment, self.image = n, t, i


def test_a_folder_that_fails_the_channel_check_is_suspect():
    """A folder with only one channel (DAPI paired with DAPI of the next field) or shifted pairs: most fields fail."""
    fin = [(_R("n1", "A", i), {}, _fi(), {}) for i in range(1, 5)]
    fin += [(_R("n1", "B", i), {}, _fi(corr=0.05), {}) for i in range(1, 4)] + [(_R("n1", "B", 4), {}, _fi(), {})]
    fin += [(_R("n1", "C", i), {}, _fi(brightest=i % 2 == 0), {}) for i in range(1, 5)]
    fin += [(_R("n1", "D", 1), {}, _fi(corr=0.05), {})] + [(_R("n1", "D", i), {}, _fi(), {}) for i in range(2, 6)]
    fin += [(_R("n1", "E", i), {}, _fi(shift=20 if i > 1 else 0.3), {}) for i in range(1, 5)]  # pairs off by a field
    assert intensity.suspect_folders(fin) == {("n1", "B"): (3, 4), ("n1", "C"): (2, 4), ("n1", "E"): (3, 4)}


def test_a_attention_lists_unused_files_and_suspect_folders(tmp_path):
    prj = tmp_path
    (prj / "layout.csv").write_text("n_folder,treatment_folder,n,date,timepoint,treatment,include,images,note\n"
                                    "a,lyovec,n1,,,LyoVec,yes,4,9 files are not a multiple of 2 channels: leftover 9.tif skipped\n"
                                    "a,cntrl,n1,,,Control,yes,4,\n")
    (prj / "intensity").mkdir()
    (prj / "intensity" / "alignment.csv").write_text(
        "n,treatment,image,channel,files,shift_px,corr,flags\n"
        "n1,miR-21,1,Ki67,1.tif; 2.tif,0.1,0.05,folder fails the channel check (2 of 2 fields): only one channel, or files mis-paired?\n"
        "n1,miR-21,2,Ki67,3.tif; 4.tif,0.1,0.04,folder fails the channel check (2 of 2 fields): only one channel, or files mis-paired?\n")
    items = {a["key"]: a for a in gui.attention_items(prj, S, {}, None, [], [])}
    assert "leftover 9.tif" in items["files"]["text"] and items["files"]["level"] == "warn"
    assert "Whole folders fail the channel check (n1 / miR-21)" in items["mispaired"]["text"]


# ---------------------------------------------------------------- b. a missing reference / treatment

def test_b_missing_reference_and_treatment_are_flagged():
    data = [r for r in rows(treatments=("Control", "LyoVec", "miR-1", "miR-2"))
            if not (r["n"] == "n2" and r["treatment"] == "LyoVec") and not (r["n"] == "n3" and r["treatment"] == "miR-2")]
    a = _assigned(data)
    gone = design.absent(a, S).set_index(["n", "treatment"])
    assert set(gone.index) == {("n2", "LyoVec"), ("n3", "miR-2")}
    assert gone.loc[("n2", "LyoVec"), "flags"].startswith("reference LyoVec absent from this n (present in n1, n3)")
    assert "condition absent from this n (present in n1, n2)" in gone.loc[("n3", "miR-2"), "flags"]
    assert (gone["images"] == 0).all()
    line = design.summary_line(pd.concat([design.describe(a), gone.reset_index()]), S)
    assert "2 conditions missing from an n (n2 / LyoVec, n3 / miR-2)" in line and "4 conditions" in line


def test_b_cell_lines_are_never_compared_for_missing_conditions():
    data = rows(ns=("n1", "n2"), treatments=("Control", "LyoVec", "miR-1"), cell_line="U251") + \
        rows(ns=("n3",), treatments=("Control", "LyoVec"), cell_line="LN229")
    a = _assigned(data).assign(cell_line=lambda d: np.where(d["n"] == "n3", "LN229", "U251"))
    assert design.absent(a, S).empty


# ---------------------------------------------------------------- c. unbalanced fields

def test_c_one_well_or_trailing_images_missing_is_flagged_against_the_rest_of_the_n():
    data = [r for r in rows() if not (r["n"] == "n1" and r["treatment"] == "miR-1" and r["image"] > 6)]
    data = [r for r in data if not (r["n"] == "n2" and r["treatment"] == "miR-1" and r["image"] not in (1, 7, 8, 9))]
    d = design.describe(_assigned(data)).set_index(["n", "treatment"])
    f1 = d.loc[("n1", "miR-1"), "flags"]
    assert "fields per well read as 3" in f1 and "6 in the other conditions" in f1 and "wells.csv" in f1
    assert "fields per well read as 5" in d.loc[("n2", "miR-1"), "flags"]
    assert d.loc[("n1", "Control"), "flags"] == ""


def test_c_a_well_with_no_usable_field_and_a_readout_without_values_are_noted():
    t = pd.DataFrame(rows(ns=("n1",), treatments=("LyoVec", "miR-1")))
    t["well"] = np.where(t["image"] <= 6, 1, 2)
    t["excluded"] = (t["treatment"] == "LyoVec") & (t["well"] == 2)
    t["source"] = np.where(t["excluded"], "auto (excluded: QC low)", "auto")
    t["readout_frac_positive_Ki67"] = np.where(t["treatment"] == "miR-1", np.nan, 0.3)
    reg = metrics.registry(["count", "readout_frac_positive_Ki67"], S)
    notes = generations.coverage_notes(t, reg)
    assert any(x.startswith("n1 / LyoVec: well 2 has no usable field (all 6 excluded: QC low)") for x in notes)
    assert any("miR-1 has no value in any field of any n (e.g. no measurable nuclei); not testable" in x for x in notes)
    t["excluded"] = t["treatment"] == "LyoVec"
    assert any("n1 / LyoVec: no usable field (all 12 excluded" in x for x in generations.coverage_notes(t, {}))


# ---------------------------------------------------------------- d. empty and all-zero conditions

def test_d_empty_folder_is_a_missing_condition():
    lay = pd.DataFrame([dict(n_folder="a", treatment_folder="cl075", n="n1", treatment="CL075", include="yes"),
                        dict(n_folder="a", treatment_folder="cntrl", n="n1", treatment="Control", include="yes")])
    gone = design.absent(_assigned(rows(ns=("n1",), treatments=("Control",))), S, lay)
    assert gone["treatment"].tolist() == ["CL075"] and "folder a/cl075 is empty" in gone["flags"].iloc[0]


def test_d_all_zero_counts_and_all_or_nothing_fractions_are_explained():
    d = pd.DataFrame(dict(treatment=["A"] * 3 + ["B"] * 3 + ["C"] * 3, value=[0, 0, 0, 5, 6, 7, 0, 0, 0]))
    cnt = models.boundary_notes(d, metrics.COUNT)
    assert len(cnt) == 1 and cnt[0].startswith("A, C: 0 nuclei in every field") and "offset" in cnt[0]
    frac = metrics.registry(["readout_frac_positive_Ki67"], S)["readout_frac_positive_Ki67"]
    f = pd.DataFrame(dict(treatment=["A"] * 2 + ["B"] * 2, value=[1.0, 1.0, 0.2, 0.3]))
    assert models.boundary_notes(f, frac) == [
        "A: 100% in every field; the empirical logit keeps it finite, so its odds ratios are arbitrary: read the "
        "percentage-point difference"]
    # back-transformed means never go below 0 because of the count offset
    assert float(models.Scale("count", 0.5).back(np.log(0.3))) == 0.0


# ---------------------------------------------------------------- end to end (slow): every case through the chain

def _no_nan_headline(out):
    """Every all-n model row either has an estimate or says why not."""
    sl = out["stats_long"]
    m = sl[(sl["level"] == "all_n") & (sl["family"].str.startswith("vs ") | (sl["family"] == "estimated means"))]
    bad = m[m["estimate"].isna() & m["note"].fillna("").eq("")]
    assert bad.empty, bad[["metric", "family", "comparison"]].to_string()


def _ok(out):
    assert out["manifest"]["files"] and out["stats"]["rows"]
    assert all("error" not in p or p.get("ok") is False for p in out["planner"])
    assert out["plans"]["plans"]
    _no_nan_headline(out)


@pytest.mark.slow
def test_e2e_a_interleaved_files_with_gaps_an_odd_file_and_a_one_channel_folder(tmp_path):
    def f(k=4, c=20, fr=0.3):
        return [(c, fr)] * k

    spec = {"01.06.26 72h": {"cntrl": [(20, .3), None, (21, .3), None, None, (19, .3), (22, .3)],  # gaps
                             "lyovec": f() + ["Ki67"],                                          # unpaired last file
                             "21": f(c=12, fr=.6)},
            "02.06.26 72h": {"cntrl": f(c=22), "lyovec": f(c=21), "21": ["DAPI"] * 8}}           # DAPI files only
    out = run_chain(image_project(tmp_path, spec))
    _ok(out)
    qc = pd.read_csv(tmp_path / "analysis" / "A" / "qc_details.csv")
    qc["n"] = qc["n"].str[:2]
    ctl = qc[(qc["n"] == "n1") & (qc["treatment"] == "Control")]
    assert ctl["image"].tolist() == [1, 2, 3, 4] and ctl["source_file_Ki67"].str.endswith(("1.tif", "4.tif", "8.tif",
                                                                                          "10.tif")).all()
    assert len(qc[(qc["n"] == "n1") & (qc["treatment"] == "LyoVec")]) == 4  # the unpaired 9.tif is skipped ...
    one = qc[(qc["n"] == "n2") & (qc["treatment"] == "miR-21")]
    assert one["readout_frac_positive_Ki67"].isna().all()                    # ... the DAPI-only folder gives no marker
    assert one["channel_check"].str.contains("folder fails the channel check").all()
    assert "WARNING: n2_2026-06-02 / miR-21" in out["log"] and "only one channel" in out["log"]
    att = {a["key"]: a for a in out["overview"]["attention"]}
    assert "leftover" in att["files"]["text"] and "Whole folders fail the channel check" in att["mispaired"]["text"]
    assert any(x.startswith("% Ki67-positive nuclei: miR-21 has no value in any field in n2_2026-06-02")
               for x in out["notes"])


@pytest.mark.slow
def test_e2e_b_missing_reference_and_treatment(tmp_path):
    data = [r for r in rows(treatments=("Control", "LyoVec", "miR-1", "miR-2"))
            if not (r["n"] == "n2" and r["treatment"] == "LyoVec") and not (r["n"] == "n3" and r["treatment"] == "miR-2")]
    out = chain(project(tmp_path, data))
    _ok(out)
    d = out["design"].set_index(["n", "treatment"])
    assert d.loc[("n2", "LyoVec"), "images"] == 0 and "reference LyoVec absent" in d.loc[("n2", "LyoVec"), "flags"]
    assert "WARNING: n2 / LyoVec" in out["log"] and "2 conditions missing from an n" in out["log"]
    att = {a["key"]: a for a in out["overview"]["attention"]}
    assert "n2 / LyoVec" in att["absent"]["text"] and "n3 / miR-2" in att["absent"]["text"]
    c = pd.read_csv(out["gen"] / "all_n" / "count_contrasts.csv").set_index(["family", "comparison"])
    assert c.loc[("vs LyoVec", "miR-1 vs LyoVec"), "N"] == 2 and c.loc[("vs Control", "miR-1 vs Control"), "N"] == 3


@pytest.mark.slow
def test_e2e_c_unbalanced_fields_and_an_empty_well(tmp_path):
    d = rows()
    d = [r for r in d if not (r["n"] == "n1" and r["treatment"] == "miR-1" and r["image"] > 6)]
    d = [r for r in d if not (r["n"] == "n2" and r["treatment"] == "miR-1" and r["image"] not in (1, 7, 8, 9))]
    d = [r for r in d if not (r["n"] == "n3" and r["treatment"] == "Control")]
    d += rows(ns=("n3",), treatments=("Control",), images=range(1, 19), seed=4)
    for r in d:
        r["qc"] = "low" if (r["n"], r["treatment"]) == ("n3", "LyoVec") and r["image"] > 6 else "high"
    out = chain(project(tmp_path, d, {"analysis": {"exclude_low_qc": True}}))
    _ok(out)
    des = out["design"].set_index(["n", "treatment"])
    assert "fields per well read as 3" in des.loc[("n1", "miR-1"), "flags"]
    assert des.loc[("n3", "Control"), "fields_per_well"] == 9 and des.loc[("n3", "Control"), "flags"] == ""
    assert any("n3 / LyoVec: well 2 has no usable field" in x for x in out["notes"])


@pytest.mark.slow
def test_e2e_d_empty_and_all_zero_conditions(tmp_path):
    d = rows(treatments=("Control", "LyoVec", "miR-1", "miR-2"))
    for r in d:
        if r["treatment"] == "miR-1":
            r["count"], r["frac"] = 0, np.nan
        elif r["treatment"] == "miR-2":
            r["frac"] = 0.0
        elif r["treatment"] == "Control":
            r["frac"] = 1.0
    out = chain(project(tmp_path, d, empty=(("n1", "CL075"),)))
    _ok(out)
    des = out["design"].set_index(["n", "treatment"])
    assert des.loc[("n1", "CL075"), "images"] == 0 and "is empty" in des.loc[("n1", "CL075"), "flags"]
    notes = " | ".join(out["notes"])
    assert "miR-1 has no value in any field of any n" in notes and "not testable" in notes
    assert "miR-1: 0 nuclei in every field" in notes and "Control: 100% in every field" in notes
    assert "miR-2: 0% in every field" in notes
    means = pd.read_csv(out["gen"] / "all_n" / "count_model_means.csv").set_index("treatment")
    assert means.loc["miR-1", "estimate"] >= 0 and means.loc["miR-1", "ci_low"] >= 0
    c = pd.read_csv(out["gen"] / "all_n" / "readout_frac_positive_Ki67_contrasts.csv")
    assert np.isfinite(c["pp_difference"]).all() and "miR-1" not in set(c["group"])


@pytest.mark.slow
def test_e2e_e_single_condition_single_n_and_three_cell_lines(tmp_path):
    out = chain(project(tmp_path / "one", rows(ns=("n1",), treatments=("LyoVec", "miR-1"))))
    _ok(out)
    assert any("needs ≥ 2 n" in x for x in out["notes"])
    assert out["planner"][0]["ok"] is False and "n1" in {a["key"] for a in out["overview"]["attention"]}
    d = rows(ns=("n1", "n2"), cell_line="U251") + rows(ns=("n3", "n4"), cell_line="LN229", seed=2) \
        + rows(ns=("n5",), cell_line="T98G", seed=3)
    out = chain(project(tmp_path / "lines", d))
    _ok(out)
    folders = {f["cell_line"]: f["n_count"] for f in out["manifest"]["all_n_folders"]}
    assert folders == {"LN229": 2, "T98G": 1, "U251": 2}
    assert {p["cell_line"]: p["ok"] for p in out["planner"]} == {"LN229": True, "T98G": False, "U251": True}
    assert out["design"]["flags"].eq("").all()  # cell lines are never compared: nothing is "missing"


@pytest.mark.slow
def test_e2e_f_copied_condition_plus_missing_reference(tmp_path):
    data = [r for r in rows() if not (r["n"] == "n2" and r["treatment"] == "LyoVec")]
    out = chain(project(tmp_path, data, copies={("n3", "Control"): ("n3", "LyoVec")}))
    _ok(out)
    dup = pd.read_csv(tmp_path / "analysis" / "E" / "duplicates.csv")
    assert len(dup) == 24 and dup["cross_condition"].all() and set(dup["n"]) == {"n3"}
    att = [a["key"] for a in out["overview"]["attention"]]
    assert att[0] == "duplicates" and "absent" in att and "acceptance" in att
    acc = pd.read_csv(out["gen"] / "acceptance.csv")
    assert acc.loc[(acc["n"] == "n3") & (acc["criterion"] == "overall"), "result"].item() == "fail"
    assert any("DUPLICATE RAW IMAGES" in x for x in out["notes"])
    assert "comparisons vs LyoVec: only 1 n (n1) has a LyoVec that is neither missing nor a duplicate - needs ≥ 2 n " \
           "for a trustworthy test" in out["notes"]
    assert any(x.startswith("comparisons vs Control: only 2 n (n1, n2)") for x in out["notes"])
