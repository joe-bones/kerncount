"""Ticket 10: tumour experiments in the web interface (lock / unlock of the analysis plan, headline statistics without
sensitivity rows, blinded review and corrections, layout, settings, results). Tiny synthetic projects only: nothing
is segmented here."""
import json
import tomllib
from pathlib import Path

import pytest

from kerncount import gui, gui_tumour, settings, settings_form
from kerncount.tumour import config, lock_log
from kerncount.tumour.results import lock as plan_lock

LEVELS = ["WT", "TLR7 KO", "cntrl", "HERV", "N4", "-clodro", "+clodro"]


def _write(p: Path, text: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _stats_long(prj: Path, h: str):
    _write(prj / "stats_long.csv",
           "level,n,metric,family,test,comparison,estimate,estimate_type,ci_low,ci_high,p,p_adj,significance,N,model,"
           "note,estimand,cell_line\n"
           f"all_n,,cross_sectional_area,estimated means,m,TLR7 KO / cntrl / -clodro,0.2,mean,0.1,0.3,,,,6,M,"
           f"confirmatory; lock {h},e,\n"
           f"all_n,,cross_sectional_area,F1,t,HERV vs cntrl | TLR7 KO / -clodro,1.5,ratio,1.1,2.0,0.01,0.02,*,12,M,"
           f"confirmatory family; lock {h},e,\n"
           f"all_n,,cross_sectional_area,F1 (Welch-type),w,HERV vs cntrl | TLR7 KO / -clodro,1.5,ratio,1.0,2.1,0.03,"
           f"0.06,,12,M,lock {h},e,\n"
           f"all_n,,cross_sectional_area,F1 (permutation),p,HERV vs cntrl | TLR7 KO / -clodro,1.5,ratio,,,0.02,0.04,,12,"
           f"M,lock {h},e,\n"
           f"all_n,,cross_sectional_area,F1 (without flagged),t,HERV vs cntrl | TLR7 KO / -clodro,1.4,ratio,1.0,2.0,"
           f"0.04,0.08,,10,M,lock {h},e,\n"
           f"all_n,,footprint_area,F3,t,x,1.1,ratio,0.9,1.3,0.4,0.8,,12,M,exploratory family; lock {h},e,\n")


def _tumour_project(ws: Path, exp="T", results=False, extra=()) -> Path:
    """A tumour project set up like `kerncount init --type tumour_volume`, with dummy stack files (only their names
    are read here) and, optionally, results of a run (stats_long.csv with the plan's lock hash, plots, review)."""
    raw = ws / "data" / exp
    sess = raw / "030826 KO7"
    for name in ("030826 cntrl -clodro n1 A.ome.tiff", "030826 cntrl -clodro n1 B.ome.tiff",
                 "030826 HERV +clodro n2 A.ome.tiff", "030826 N4 -clodron n1 A.ome.tiff", *extra):
        _write(sess / name, "x")
    prj = ws / "analysis" / exp
    prj.mkdir(parents=True, exist_ok=True)
    settings.write_template(prj, raw, template=config.TEMPLATE)
    from kerncount.tumour import layout

    s = settings.load(prj)
    lay = layout.scan(raw, s)
    from kerncount import discover

    discover.write_layout(lay, prj / "layout.csv", ",")
    if results:
        h, _ = plan_lock(s)
        _stats_long(prj, h)
        _write(prj / "tumours.csv", "tumour_id,genotype,flagged\nTLR7 KO / cntrl / -clodro / n1A,TLR7 KO,no\n")
        for name in ("cross_sectional_area_absolute", "cross_sectional_area_normalised", "cross_sectional_area_forest",
                     "perimeter_absolute", "shape_panel"):
            _write(prj / "plots" / f"{name}.png", "png")
            _write(prj / "plots" / f"{name}.pdf", "pdf")
        _write(prj / "results.xlsx", "xlsx")
        _write(prj / "stats_report.md", f"# Tumour statistics report\n- lock {h}\n")
        # blinded review: the key (code -> tumour ID, i.e. every level) sits OUTSIDE review/
        _write(prj / "review_key.csv", "code,tumour_id\nBCDFG,TLR7 KO / cntrl / -clodro / n1A\n"
               "HJKLM,TLR7 KO / HERV / +clodro / n2A\n")
        _write(prj / "review" / "corrections.csv", "code,objects,qc_flags,manual_1,manual_2,manual_3,note\n"
               "BCDFG,3,,,,,\nHJKLM,1,edge contact,,,,\n")
        _write(prj / "review" / "BCDFG.jpg", "jpg")
        _write(prj / "review" / "HJKLM.jpg", "jpg")
        _write(prj / "review" / "overview_1.jpg", "jpg")
        _write(prj / "corrections_used.csv", "code,tumour_id,main_mass_object,source,generation\n"
               "BCDFG,TLR7 KO / cntrl / -clodro / n1A,1,auto,g1\n")
    return prj


# ---------------------------------------------------------------- lock / unlock of the analysis plan (ADR 0003)

def test_plan_is_not_locked_before_results_exist(tmp_path):
    st = lock_log.status(tmp_path, "aaaaaaaaaaaa")
    assert not st["results"] and not st["locked"] and st["results_lock"] == ""
    lock_log.check_change(tmp_path, "aaaaaaaaaaaa", "bbbbbbbbbbbb")  # no results: any change is fine


def test_results_lock_the_plan_until_an_explicit_logged_unlock(tmp_path):
    _stats_long(tmp_path, "aaaaaaaaaaaa")
    assert lock_log.results_lock(tmp_path) == "aaaaaaaaaaaa"
    st = lock_log.status(tmp_path, "aaaaaaaaaaaa")
    assert st["results"] and st["locked"]
    lock_log.check_change(tmp_path, "aaaaaaaaaaaa", "aaaaaaaaaaaa")  # no change to the plan: always allowed
    with pytest.raises(lock_log.PlanLocked, match="locked"):
        lock_log.check_change(tmp_path, "aaaaaaaaaaaa", "bbbbbbbbbbbb")
    with pytest.raises(ValueError, match="reason"):
        lock_log.unlock(tmp_path, "  ", "aaaaaaaaaaaa")
    lock_log.unlock(tmp_path, "lab decided perimeter is exploratory", "aaaaaaaaaaaa")
    assert not lock_log.status(tmp_path, "aaaaaaaaaaaa")["locked"]
    lock_log.check_change(tmp_path, "aaaaaaaaaaaa", "bbbbbbbbbbbb")
    lock_log.record_change(tmp_path, "aaaaaaaaaaaa", "bbbbbbbbbbbb")
    log = (tmp_path / lock_log.LOG).read_text(encoding="utf-8")
    assert "unlock" in log and "lab decided perimeter is exploratory" in log and "bbbbbbbbbbbb" in log
    # a run with the new plan closes the unlock: the new results are locked again
    lock_log.record_run(tmp_path, "bbbbbbbbbbbb")
    _stats_long(tmp_path, "bbbbbbbbbbbb")
    st = lock_log.status(tmp_path, "bbbbbbbbbbbb")
    assert st["locked"] and st["results_lock"] == "bbbbbbbbbbbb"
    lines = lock_log.report_lines(tmp_path)
    text = "\n".join(lines)
    assert "aaaaaaaaaaaa" in text and "bbbbbbbbbbbb" in text and "lab decided perimeter is exploratory" in text
    assert [e["action"] for e in lock_log.entries(tmp_path)] == ["unlock", "plan changed", "run"]


def test_relock_and_runs_without_an_unlock_are_recorded(tmp_path):
    _stats_long(tmp_path, "aaaaaaaaaaaa")
    lock_log.unlock(tmp_path, "typo in a level", "aaaaaaaaaaaa")
    lock_log.relock(tmp_path, "aaaaaaaaaaaa")
    assert lock_log.status(tmp_path, "aaaaaaaaaaaa")["locked"]
    with pytest.raises(lock_log.PlanLocked):
        lock_log.check_change(tmp_path, "aaaaaaaaaaaa", "cccccccccccc")
    # settings.toml edited outside the web interface: the run says so in the log (and so in the report)
    lock_log.record_run(tmp_path, "cccccccccccc")
    last = lock_log.entries(tmp_path)[-1]
    assert last["action"] == "run" and "without a logged unlock" in last["reason"]
    assert any("without a logged unlock" in x for x in lock_log.report_lines(tmp_path))
    # a run with the same plan and no unlock pending writes nothing
    n = len(lock_log.entries(tmp_path))
    _stats_long(tmp_path, "cccccccccccc")
    lock_log.record_run(tmp_path, "cccccccccccc")
    assert len(lock_log.entries(tmp_path)) == n


# ---------------------------------------------------------------- headline statistics without sensitivity rows

def test_headline_statistics_leave_out_sensitivity_rows(tmp_path):
    _stats_long(tmp_path, "aaaaaaaaaaaa")
    cols, rows, _ = gui.read_csv(tmp_path / "stats_long.csv")
    head = gui_tumour.headline(cols, rows)
    fams = [r[cols.index("family")] for r in head]
    assert fams == ["estimated means", "F1", "F3"]
    assert gui_tumour.is_sensitivity("F2 (Welch-type)") and gui_tumour.is_sensitivity("F7 (permutation)")
    assert gui_tumour.is_sensitivity("F1 (without flagged)") and not gui_tumour.is_sensitivity("F1")


def test_a_rebuild_with_a_changed_plan_shows_the_lock_change_in_the_report(tmp_path):
    """results.analyse records the run in lock_log.csv and prints the lock history in the report / About sheet."""
    import numpy as np
    import pandas as pd

    from kerncount.tumour import results
    from test_tumour_stats import _settings, _vol_table

    t = _vol_table(np.random.default_rng(10))
    vol = [{"name": "main_mass_volume", "column": "main_mass_volume_mm3", "role": "confirmatory"}]
    s1 = _settings(readouts=vol, sensitivity_permutations=99)
    results.analyse(tmp_path, t, s1)
    h1 = plan_lock(s1)[0]
    assert lock_log.results_lock(tmp_path) == h1 and not (tmp_path / lock_log.LOG).exists()
    lock_log.unlock(tmp_path, "volume is exploratory after all", h1)
    s2 = _settings(readouts=[{**vol[0], "role": "exploratory"}], sensitivity_permutations=99)
    h2 = plan_lock(s2)[0]
    results.analyse(tmp_path, t, s2)
    assert lock_log.results_lock(tmp_path) == h2 != h1
    report = (tmp_path / "stats_report.md").read_text(encoding="utf-8")
    assert "Analysis plan lock history" in report and f"{h1} -> {h2}" in report
    assert "volume is exploratory after all" in report
    about = pd.read_excel(tmp_path / "results.xlsx", "About", header=None)[0].astype(str).tolist()
    assert any("lock history" in x for x in about)
    assert lock_log.status(tmp_path, h2)["locked"]


# ---------------------------------------------------------------- the web interface on tumour projects

from test_gui import Server  # noqa: E402


@pytest.fixture()
def tsrv(tmp_path):
    ws = tmp_path / "ws"
    _tumour_project(ws, "T", results=True)
    _tumour_project(ws, "U", results=False, extra=["030826 XYZ -clodro n3 A.ome.tiff"])  # an unreadable name
    s = Server(ws)
    yield s
    s.close()


def test_state_types_and_set_up_with_a_type(tsrv):
    code, st = tsrv.call("/api/state")
    assert code == 200 and {"cell_count", "tumour_volume"} <= {t["name"] for t in st["types"]}
    e = {x["name"]: x for x in st["experiments"]}
    assert e["T"]["type"] == "tumour_volume" and e["T"]["results"] and not e["U"]["results"]
    assert gui.build_args("init", "E", {"type": "tumour_volume"}) == ["init", "--type", "tumour_volume", "--", "E"]
    assert gui.build_args("init", "E", {}) == ["init", "--", "E"]
    with pytest.raises(ValueError, match="unknown experiment type"):
        gui.build_args("init", "E", {"type": "--force"})


def test_tumour_settings_form_and_rerun_advice(tsrv):
    code, p = tsrv.call("/api/project", exp="T")
    assert code == 200 and p["type"] == "tumour_volume"
    secs = {s["section"]: s for s in p["schema"]}
    assert set(secs) == {"project", "experiment", "tumour"}
    assert [i["key"] for i in secs["experiment"]["items"]] == ["type", "notes"]
    items = {i["key"]: i for i in secs["tumour"]["items"]}
    for k in ("smooth_sigma_um", "threshold", "mid_section_um", "outline_smooth_um", "threshold_sweep"):
        assert items[k]["help"] or k == "threshold_sweep"
    assert items["mid_section_um"]["default"] == config.DEFAULTS["mid_section_um"] and items["mid_section_um"]["reseg"]
    # not part of the plan: allowed although results exist (the plan is locked)
    code, r = tsrv.call("/api/settings", body={"experiment": "T", "values": {
        "tumour": {"threshold_sweep": True, "sweep_scales": "0.85, 1.15"}}})
    assert code == 200, r
    assert [s["action"] for s in r["steps"]] == ["run"] and not r["steps"][0]["slow"]
    s = settings.load(tsrv.state.workspace / "analysis" / "T")
    assert s["tumour"]["threshold_sweep"] is True and s["tumour"]["sweep_scales"] == [0.85, 1.15]
    code, r = tsrv.call("/api/settings", body={"experiment": "T", "values": {"tumour": {"mid_section_um": "30"}}})
    assert code == 200 and r["steps"][0]["slow"]
    # cell-count projects keep their form
    assert "tumour" not in {s["section"] for s in settings_form.schema()}


def test_locked_plan_needs_a_logged_unlock(tsrv):
    prj = tsrv.state.workspace / "analysis" / "T"
    code, plan = tsrv.call("/api/tumour/plan", exp="T")
    assert code == 200 and plan["lock"]["locked"] and plan["lock"]["results_lock"] == plan["lock"]["settings_lock"]
    assert [f["name"] for f in plan["factors"]] == ["genotype", "tumour_line", "clodronate"]
    assert {r["name"]: r["role"] for r in plan["readouts"]}["perimeter"] == "confirmatory"
    assert {f["id"]: f["role"] for f in plan["families"]}["F6"] == "descriptive"
    change = {"experiment": "T", "readouts": [{"name": "perimeter", "role": "exploratory"}]}
    code, r = tsrv.call("/api/tumour/plan", body=change)
    assert code == 400 and "locked" in r["error"]
    # the same change through the raw settings.toml editor is blocked too
    text = (prj / "settings.toml").read_text(encoding="utf-8")
    edited = text.replace('column = "perimeter_mm"\nlabel = "Perimeter, mid-section (mm)"\nrole = "confirmatory"',
                          'column = "perimeter_mm"\nlabel = "Perimeter, mid-section (mm)"\nrole = "exploratory"')
    assert edited != text
    code, r = tsrv.call("/api/settings", body={"experiment": "T", "text": edited})
    assert code == 400 and "locked" in r["error"]
    assert (prj / "settings.toml").read_text(encoding="utf-8") == text
    assert tsrv.call("/api/settings_reset", body={"experiment": "U"})[0] == 200  # no results: free
    # unlock needs a reason; then the change is saved and logged
    assert tsrv.call("/api/tumour/unlock", body={"experiment": "T", "reason": ""})[0] == 400
    code, r = tsrv.call("/api/tumour/unlock", body={"experiment": "T", "reason": "lab meeting 2026-11-20"})
    assert code == 200 and not r["lock"]["locked"]
    code, r = tsrv.call("/api/tumour/plan", body={**change, "factors": [
        {"name": "tumour_line", "levels": "cntrl, HERV, N4", "reference": "cntrl", "role": "panel"}]})
    assert code == 200 and r["changed"], r
    s = settings.load(prj)
    assert {x["name"]: x["role"] for x in s["tumour"]["readouts"]}["perimeter"] == "exploratory"
    assert {x["name"]: x["role"] for x in s["tumour"]["factors"]}["tumour_line"] == "panel"
    assert "# Factors: one column per factor" in (prj / "settings.toml").read_text(encoding="utf-8")  # comments kept
    st = r["lock"]
    assert st["changed"] and st["settings_lock"] == plan_lock(s)[0] != st["results_lock"]
    assert [e["action"] for e in lock_log.entries(prj)] == ["unlock", "plan changed"]
    # an invalid plan is refused (reference not a level)
    bad = {"experiment": "T", "factors": [{"name": "clodronate", "levels": "-clodro, +clodro", "reference": "none"}]}
    assert tsrv.call("/api/tumour/plan", body=bad)[0] == 400
    code, r = tsrv.call("/api/tumour/relock", body={"experiment": "T"})
    assert code == 200 and r["lock"]["locked"]


def test_review_and_corrections_are_blind(tsrv):
    prj = tsrv.state.workspace / "analysis" / "T"
    code, rv = tsrv.call("/api/tumour/review", exp="T")
    assert code == 200 and [r["code"] for r in rv["rows"]] == ["BCDFG", "HJKLM"]
    assert rv["rows"][0]["image"] == "review/BCDFG.jpg" and rv["overviews"] == ["review/overview_1.jpg"]
    assert rv["manual"] == ["manual_1", "manual_2", "manual_3"]
    body = {"experiment": "T", "mtime": rv["mtime"], "changes": [
        {"code": "BCDFG", "col": "manual_1", "value": "2"},
        {"code": "HJKLM", "col": "manual_1", "value": "X:  cut  in half"},
        {"code": "HJKLM", "col": "note", "value": "edge"}]}
    code, r = tsrv.call("/api/tumour/corrections", body=body)
    assert code == 200, r
    cols, rows, _ = gui.read_csv(prj / "review" / "corrections.csv")
    assert rows[0][cols.index("manual_1")] == "2" and rows[1][cols.index("manual_1")] == "x: cut in half"
    assert rows[1][cols.index("note")] == "edge"
    for bad in ({"code": "BCDFG", "col": "manual_1", "value": "7"},      # only 3 objects on the image
                {"code": "BCDFG", "col": "manual_1", "value": "0"},
                {"code": "BCDFG", "col": "manual_1", "value": "maybe"},
                {"code": "ZZZZZ", "col": "manual_1", "value": "1"},
                {"code": "BCDFG", "col": "objects", "value": "9"},
                {"code": "BCDFG", "col": "manual_9", "value": "1"}):
        assert tsrv.call("/api/tumour/corrections", body={"experiment": "T", "changes": [bad]})[0] == 400, bad
    # filling the last round adds the next one
    tsrv.call("/api/tumour/corrections", body={"experiment": "T", "changes": [
        {"code": "BCDFG", "col": "manual_3", "value": "1"}]})
    assert "manual_4" in gui.read_csv(prj / "review" / "corrections.csv")[0]
    # nothing the Review / Corrections pages receive names a group, a session, a date or a tumour ID
    code, rv = tsrv.call("/api/tumour/review", exp="T")
    blob = json.dumps(rv)
    for word in LEVELS + ["KO7", "030826", "n1A", "n2A", "TLR7"]:
        assert word not in blob, word


def test_layout_page_shows_validation_and_saves_edits(tsrv):
    code, lay = tsrv.call("/api/tumour/layout", exp="U")
    assert code == 200
    cols = lay["columns"]
    for c in ("tumour_id", "imaging_session", "genotype", "tumour_line", "clodronate", "section", "hemisphere",
              "include", "note"):
        assert c in cols, c
    assert [f["name"] for f in lay["factors"]] == ["genotype", "tumour_line", "clodronate"]
    assert any("XYZ" in p for p in lay["problems"]), lay["problems"]
    assert {g["group"]: g["tumours"] for g in lay["groups"]}["TLR7 KO / cntrl / -clodro"] == 2
    rows = [dict(zip(cols, r)) for r in lay["rows"]]
    for r in rows:
        if "XYZ" in r["file"]:
            r["include"], r["note"] = "no", "unreadable name, checked with the lab"
    code, r = tsrv.call("/api/tumour/layout", body={"experiment": "U", "rows": rows})
    assert code == 200 and not r["layout"]["problems"], r
    rows[0]["file"] = "other.ome.tiff"
    assert tsrv.call("/api/tumour/layout", body={"experiment": "U", "rows": rows})[0] == 400


def test_results_overview_and_headline_statistics(tsrv):
    code, res = tsrv.call("/api/tumour/results", exp="T")
    assert code == 200 and res["xlsx"] == "results.xlsx" and res["report"] == "stats_report.md"
    ro = {r["name"]: r for r in res["readouts"]}
    assert res["readouts"][0]["role"] == "confirmatory"
    assert [p["kind"] for p in ro["cross_sectional_area"]["plots"]] == ["absolute", "normalised", "forest"]
    assert ro["cross_sectional_area"]["plots"][0]["pdf"] == "plots/cross_sectional_area_absolute.pdf"
    assert res["extra"][0]["path"] == "plots/shape_panel.png" and "tumours.csv" in res["tables"]
    code, st = tsrv.call("/api/tumour/stats", exp="T")
    fam = st["columns"].index("family")
    assert code == 200 and st["hidden"] == 3 and not any("(" in r[fam] for r in st["rows"])
    code, st = tsrv.call("/api/tumour/stats", exp="T", sensitivity="1")
    assert st["hidden"] == 0 and len(st["rows"]) == 6
    code, ov = tsrv.call("/api/tumour/overview", exp="T")
    assert code == 200 and ov["tumours"] == 1 and ov["generation"] == "g1" and ov["lock"]["locked"]
    code, ov = tsrv.call("/api/tumour/overview", exp="U")
    assert code == 200 and ov["attention"][0]["key"] == "layout" and ov["attention"][0]["level"] == "error"


def test_page_scripts_load_the_tumour_views_and_the_blind_ones_call_only_blind_endpoints():
    import re

    web = Path(gui.__file__).with_name("web")
    assert '<script src="/static/tumour.js"></script>' in (web / "index.html").read_text()
    js = (web / "tumour.js").read_text()
    blind = js.split("// ---------------------------------------------------------------- blinded review")[1].split(
        "// ---------------------------------------------------------------- results and statistics")[0]
    assert set(re.findall(r'"(/api/[^"?]+)', blind)) == {"/api/tumour/review", "/api/tumour/corrections"}
    assert "review_key" not in blind.replace("<code>review_key.csv</code>", "")
    for sec in ("home", "layout", "review", "corrections", "results", "stats"):
        assert f"T.{sec} = async function" in js, sec
