import tomllib

import pytest

from kerncount import gui, settings, settings_form


def test_settings_form_keeps_comments_and_types(tmp_path):
    settings.write_template(tmp_path, tmp_path / "raw")
    text = (tmp_path / "settings.toml").read_text(encoding="utf-8")
    typed = settings_form.coerce({
        "counting": {"edge_rule": "exclude", "min_area_frac": "0,25", "exclude_out_of_focus": False},
        "intensity": {"measure": "2, 3", "regions": ["ring", "nucleus"]},
        "segmentation": {"nucleus_diameter_px": "auto", "pixel_size_um": "0.65"},
        "analysis": {"normalise_to": "LyoVec, Control"},
    })
    out = settings_form.update_text(text, typed)
    d = tomllib.loads(out)
    assert d["counting"]["edge_rule"] == "exclude" and d["counting"]["min_area_frac"] == 0.25
    assert d["counting"]["exclude_out_of_focus"] is False
    assert d["intensity"]["measure"] == [2, 3] and d["intensity"]["regions"] == ["nucleus", "ring"]
    assert d["segmentation"]["pixel_size_um"] == 0.65 and d["segmentation"]["nucleus_diameter_px"] == "auto"
    assert d["analysis"]["normalise_to"] == ["LyoVec", "Control"]
    assert "# kerncount settings" in out and "unbiased" in out  # comments survive
    assert out.count("\n") == text.count("\n")


def test_settings_form_appends_missing_keys_and_rejects_bad_values():
    out = settings_form.update_text("[counting]\nedge_rule = \"half\"\n", settings_form.coerce(
        {"counting": {"focus_frac": "0.4"}, "review": {"jpeg_quality": "70"}}))
    d = tomllib.loads(out)
    assert d["counting"] == {"edge_rule": "half", "focus_frac": 0.4} and d["review"]["jpeg_quality"] == 70
    for bad in ({"counting": {"edge_rule": "sometimes"}}, {"review": {"jpeg_quality": "101"}},
                {"counting": {"min_area_frac": "abc"}}, {"intensity": {"regions": []}}):
        with pytest.raises(ValueError):
            settings_form.coerce(bad)


def _counts(tmp_path, delim=","):
    p = tmp_path / "counts.csv"
    p.write_text(delim.join("n treatment image auto_count qc_ok qc_level manual_1 manual_2".split()) + "\n"
                 + delim.join(["n1", "Control", "1", "50", "yes", "high", "", ""]) + "\n"
                 + delim.join(["n1", "Control", "2", "60", "no", "low", "", ""]) + "\n", encoding="utf-8")
    return p


@pytest.mark.parametrize("delim", [",", ";"])
def test_save_counts_writes_manual_values_and_pads(tmp_path, delim):
    p = _counts(tmp_path, delim)
    ch = lambda img, col, v: dict(n="n1", treatment="Control", image=str(img), col=col, value=v)  # noqa: E731
    gui.save_counts(tmp_path, [ch(1, "manual_1", "55"), ch(2, "manual_1", "X")], None)
    cols, rows, d = gui.read_csv(p)
    assert d == delim and cols[-2:] == ["manual_1", "manual_2"]
    assert rows[0][6] == "55" and rows[1][6] == "x"
    gui.save_counts(tmp_path, [ch(1, "manual_2", "40")], None)  # filling the last column adds a fresh one
    cols, rows, _ = gui.read_csv(p)
    assert cols[-1] == "manual_3" and rows[0][7] == "40"
    gui.save_counts(tmp_path, [ch(1, "manual_1", "")], None)  # clearing
    assert gui.read_csv(p)[1][0][6] == ""


def test_save_counts_rejects_bad_input(tmp_path):
    p = _counts(tmp_path)
    for change in (dict(n="n1", treatment="Control", image="1", col="manual_1", value="4.5"),
                   dict(n="n1", treatment="Control", image="1", col="auto_count", value="4"),
                   dict(n="n1", treatment="Control", image="9", col="manual_1", value="4"),
                   dict(n="n1", treatment="Control", image="1", col="manual_7", value="4")):
        with pytest.raises(ValueError):
            gui.save_counts(tmp_path, [change], None)
    with pytest.raises(ValueError, match="changed by another program"):
        gui.save_counts(tmp_path, [], p.stat().st_mtime - 10)


def test_paths_cannot_leave_project(tmp_path):
    with pytest.raises(ValueError):
        gui.inside(tmp_path, "../x")
    assert gui.inside(tmp_path, "a/b") == (tmp_path / "a" / "b").resolve()


def test_build_args():
    assert gui.build_args("regen", "e", {"rebuild": "g2, g3; rm"}) == ["regen", "e", "--rebuild", "g2", "g3"]
    assert gui.build_args("sampling", "e", {"generation": "g2", "iterations": "200"}) == \
        ["sampling", "e", "--generation", "g2", "--iterations", "200"]


# ================================================================ HTTP server: endpoints, security, blinding

import json  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import threading  # noqa: E402
import urllib.error  # noqa: E402
import urllib.request  # noqa: E402
from http.server import ThreadingHTTPServer  # noqa: E402
from pathlib import Path  # noqa: E402
from urllib.parse import urlencode  # noqa: E402

WEB = Path(gui.__file__).with_name("web")


class Server:
    """The real request handler on a free port, talking HTTP (token, Host and Origin checks included)."""

    def __init__(self, ws: Path):
        self.state = gui.State(ws)
        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), gui.make_handler(self.state))
        self.state.httpd = self.httpd
        self.port = self.httpd.server_address[1]
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def call(self, endpoint, body=None, token=True, host=None, origin=None, **params):
        if token and body is None:
            params["t"] = gui.TOKEN
        url = f"http://127.0.0.1:{self.port}{endpoint}" + ("?" + urlencode(params) if params else "")
        req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                     method="GET" if body is None else "POST")
        if body is not None:
            req.add_header("Content-Type", "application/json")
            if token:
                req.add_header("X-Token", gui.TOKEN)
        if host:
            req.add_header("Host", host)
        if origin:
            req.add_header("Origin", origin)
        try:
            r = urllib.request.urlopen(req, timeout=60)
            code, raw, ct = r.status, r.read(), r.headers.get("Content-Type", "")
        except urllib.error.HTTPError as e:
            code, raw, ct = e.code, e.read(), e.headers.get("Content-Type", "")
        return code, (json.loads(raw) if ct.startswith("application/json") else raw)

    def close(self):
        self.httpd.shutdown()


def _write(p: Path, text: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _results_project(ws: Path, exp="E") -> Path:
    """A small project with one generation (manifest, stats_long, acceptance), duplicates, a design flag and a
    planner built from simulated data - enough for every results endpoint without running the pipeline."""
    import numpy as np

    from kerncount import planner
    from test_planner import _settings, simulate

    (ws / "data" / exp / "05.06.26 a" / "cntrl").mkdir(parents=True)
    prj = ws / "analysis" / exp
    prj.mkdir(parents=True)
    settings.write_template(prj, ws / "data" / exp)
    _write(prj / "layout.csv", "n_folder,treatment_folder,n,date,timepoint,treatment,include,images,note\n"
           "05.06.26 a,cntrl,n1_2026-06-05,2026-06-05,72h,Control,yes,2,\n05.06.26 a,lyo,n1_2026-06-05,2026-06-05,72h,LyoVec,yes,2,\n")
    _write(prj / "qc_details.csv", "n,treatment,image,auto_count,count_low,count_high,qc_level,qc_flags,source_file\n"
           "n1_2026-06-05,Control,1,50,49,51,high,,/raw/cntrl/1.tif\nn1_2026-06-05,Control,2,40,30,60,low,blur,/raw/cntrl/2.tif\n"
           "n1_2026-06-05,LyoVec,1,52,50,53,high,,/raw/lyo/1.tif\nn1_2026-06-05,LyoVec,2,41,40,42,low,,/raw/lyo/2.tif\n")
    _write(prj / "design_check.csv", "n,treatment,images,wells_expected,fields_per_well,fields_source,unbalanced,flags\n"
           "n1_2026-06-05,Control,2,2,1,auto,False,\nn1_2026-06-05,LyoVec,2,2,1,auto,True,missing image(s) 3\n")
    _write(prj / "duplicates.csv", "group,n,treatment,image,source_file,hash,copies,cross_condition,shared_with\n"
           "1,n1_2026-06-05,Control,1,/raw/cntrl/1.tif,ab,2,True,n1_2026-06-05 / LyoVec / 1\n")
    g = prj / "generations" / "g1"
    _write(g / "generation.json", json.dumps({"name": "g1", "source": "auto", "changes": 0, "built": "2026-11-16T10:00:00"}))
    _write(g / "acceptance.csv", "n,criterion,value,limit,result,detail\nn1_2026-06-05,no_duplicate_images,1,0,fail,1 image shared\n"
           "n1_2026-06-05,overall,,,fail,failed\n")
    _write(g / "all_n" / "n_colours.csv", "n,colour\nn1_2026-06-05,#1f77b4\n")
    _write(g / "all_n" / "count.png", "png")
    _write(g / "stats_long.csv", "level,n,metric,family,test,comparison,estimate,estimate_type,ci_low,ci_high,p,p_adj,"
           "significance,N,model,note,estimand,cell_line\nall_n,,count,vs Control,mixed,LyoVec vs Control,1.04,ratio,,,,,,1,"
           "descriptive,needs ≥ 2 n,ratio of mean nuclei per field,\nn,n1_2026-06-05,count,omnibus,anova,all,,,,,0.5,0.5,ns,,classical,nan,,\n")
    man = {"schema": "kerncount.generation.manifest/1", "generation": "g1", "n_count": 1,
           "metrics": [{"id": "count", "label": "Nuclei per field", "unit": "nuclei per field", "kind": "count"}],
           "all_n_folders": [{"cell_line": None, "path": "all_n", "n_count": 1}], "model_status": {},
           "duplicates": {"summary": ["n1 / Control 1 = n1 / LyoVec 1"]}, "acceptance": {"enforce": False},
           "notes": ["only one n"], "layers": {},
           "files": [{"path": "all_n/count.png", "level": "all_n", "n": None, "treatment": None, "metric": "count",
                      "kind": "absolute", "title": "All n: count", "description": "d", "cell_line": None}]}
    _write(g / "manifest.json", json.dumps(man))
    plan = planner.build(simulate(np.random.default_rng(3), n=4), _settings())
    _write(prj / "sampling" / "g1" / "planner.json", json.dumps(planner._jsonable(plan)))
    _write(prj / "sampling" / "g1" / "REPORT.txt", "Recommendations for the next experiment\n=====\nplan 5 n's\n\n\nHow many images\n")
    return prj


@pytest.fixture(scope="module")
def srv(tmp_path_factory):
    ws = tmp_path_factory.mktemp("ws")
    _results_project(ws)
    s = Server(ws)
    yield s
    s.close()


def test_page_token_host_and_origin_checks(srv):
    code, page = srv.call("/", token=False)
    assert code == 200 and gui.TOKEN.encode() in page
    assert srv.call("/api/state", token=False)[0] == 403
    assert srv.call("/api/state", t="wrong", token=False)[0] == 403
    assert srv.call("/api/state")[0] == 200
    assert srv.call("/api/state", host="evil.example:80")[0] == 403          # DNS rebinding
    assert srv.call("/api/state", origin="http://evil.example")[0] == 403    # other websites
    assert srv.call("/api/tuning/render", token=False, exp="E", field="f1")[0] == 403
    assert srv.call("/api/field_overlay", token=False, exp="E", n="n1", treatment="Control", image="1")[0] == 403
    for path in ("/api/platemap", "/api/planner/evaluate", "/api/tuning/lock", "/api/tuning/sample", "/api/settings"):
        assert srv.call(path, body={"experiment": "E"}, token=False)[0] == 403, path


def test_paths_are_confined_to_the_project(srv):
    assert srv.call("/static/../gui.py", token=False)[0] == 404
    for p in ("../../x", "/etc/passwd", "../E/settings.toml/../../../.."):
        assert srv.call("/file", exp="E", path=p)[0] in (400, 404), p
        assert srv.call("/api/table", exp="E", path=p)[0] == 400, p
    assert srv.call("/file", exp="E", path="settings.toml")[0] == 200
    assert srv.call("/api/manifest", exp="E", gen="../g1")[0] == 400
    assert srv.call("/api/stats", exp="E", gen="g1/../../x")[0] == 400
    assert srv.call("/api/overview", exp="../E")[0] == 400   # unknown experiment


def test_static_files_and_javascript_syntax(srv):
    for f in ("app.js", "results.js", "intensity.js", "plan.js", "app.css"):
        assert srv.call(f"/static/{f}", token=False)[0] == 200, f
    page = (WEB / "index.html").read_text()
    for f in ("app.js", "results.js", "intensity.js", "plan.js"):
        assert f"/static/{f}" in page
    node = shutil.which("node")
    if node:
        for f in WEB.glob("*.js"):
            r = subprocess.run([node, "--check", str(f)], capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
    for f in WEB.glob("*.js"):  # in-page confirmations only
        text = f.read_text()
        assert "window.confirm" not in text and "window.alert" not in text and "alert(" not in text, f.name


def test_manifest_stats_and_overview(srv):
    code, m = srv.call("/api/manifest", exp="E", gen="g1")
    assert code == 200 and m["manifest"]["files"][0]["path"] == "all_n/count.png"
    assert m["colours"] == {"n1_2026-06-05": "#1f77b4"}
    code, st = srv.call("/api/stats", exp="E", gen="g1")
    assert code == 200 and st["columns"][0] == "level" and len(st["rows"]) == 2
    assert st["rows"][1][st["columns"].index("note")] == ""  # "nan" is shown as empty
    code, ov = srv.call("/api/overview", exp="E")
    assert code == 200 and ov["latest"] == "g1" and ov["meta"]["reference"]
    keys = {a["key"] for a in ov["attention"]}
    assert {"duplicates", "acceptance", "design", "n1", "qc"} <= keys
    assert [a["level"] for a in ov["attention"]][0] == "error"  # most serious first
    assert ov["ns"][0]["timepoints"] == ["72h"] and ov["ns"][0]["acceptance"] == "fail"
    code, rv = srv.call("/api/review", exp="E")
    img = rv[0]["treatments"][0]["images"][0]
    assert img["sources"] == {"all": "/raw/cntrl/1.tif"} and img["duplicate"] == "n1_2026-06-05 / LyoVec / 1"


def test_planner_endpoint(srv):
    code, p = srv.call("/api/planner", exp="E")
    assert code == 200 and p["ok"] and p["gen"] == "g1" and "draws" not in p and "_np" not in p
    assert p["recommendations"].startswith("Recommendations") and "How many images" not in p["recommendations"]
    body = {"experiment": "E", "gen": "g1", "n": 5, "wells": 2, "fields": 6, "effect": "0.7", "vary": "fields"}
    code, r = srv.call("/api/planner/evaluate", body=body)
    assert code == 200
    x = r["result"]
    assert x["n"] == 5 and x["power"]["low"] <= x["power"]["mid"] <= x["power"]["high"]
    assert r["vary"] == "fields" and [c["x"] for c in r["curve"]][:3] == [1, 2, 3]
    assert {"power", "ci_halfwidth_pct", "mde_ratio", "cost"} <= set(r["curve"][0])
    code, r = srv.call("/api/planner/evaluate", body={**body, "effect": "miR-0", "reference_wells": "auto", "vary": "n"})
    assert code == 200 and r["result"]["effect_label"] == "miR-0" and r["result"]["reference_wells"] > 2
    for bad in ({"n": 1}, {"effect": "-3"}, {"vary": "plates"}, {"alpha": 0.9}, {"gen": "../g1"}):
        assert srv.call("/api/planner/evaluate", body={**body, **bad})[0] == 400, bad


def test_platemap_endpoint(srv):
    body = {"experiment": "E", "plates": 2, "plate_format": "12", "wells": 2, "reference_wells": "auto", "seed": "7",
            "negative_control": "secondary only", "negative_control_wells": 1, "avoid_edge": True}
    code, r = srv.call("/api/platemap", body=body)
    assert code == 200 and r["name"].startswith("platemap_")
    code, pl = srv.call("/api/plans", exp="E")
    plan = pl["plans"][0]
    assert plan["name"] == r["name"] and len(plan["plates"]) == 2 and "seed 7" in plan["readme"]
    assert plan["edge"]["columns"][0] == "condition" and pl["conditions"] == ["Control", "LyoVec"]
    assert srv.call("/file", exp="E", path=plan["plates"][0])[0] == 200
    for bad in ({"plate_format": "10"}, {"plates": 0}, {"seed": "x"}, {"negative_control_wells": 3}):
        assert srv.call("/api/platemap", body={**body, **bad})[0] == 400, bad
    assert gui.build_args("platemap", "E", {"plates": 2, "plate_format": "24", "reference_wells": "auto", "seed": "3",
                                            "avoid_edge": True, "negative_control": "neg"}) == \
        ["platemap", "E", "--plates", "2", "--wells", "2", "--format", "24", "--reference-wells", "auto",
         "--avoid-edge", "--seed", "3", "--negative-control", "neg", "--control-wells", "2"]
    assert gui.build_args("wells-template", "E", {"force": True}) == ["wells-template", "E", "--force"]


# one form item per kind: (section, key, value sent by the page, value expected in settings.toml)
KIND_CASES = {
    "bool": ("counting", "exclude_out_of_focus", False, False),
    "int": ("review", "jpeg_quality", "70", 70),
    "float": ("counting", "min_area_frac", "0,3", 0.3),
    "choice": ("counting", "edge_rule", "exclude", "exclude"),
    "multi": ("intensity", "regions", ["ring", "nucleus"], ["nucleus", "ring"]),
    "list": ("analysis", "normalise_to", ["LyoVec", "Control"], ["LyoVec", "Control"]),
    "num_or": ("experiment", "magnification", "20", 20),
    "measure": ("intensity", "measure", "2, 3", [2, 3]),
    "names": ("readouts", "include", ["coverage", "frac_out_of_focus"], ["frac_out_of_focus", "coverage"]),
    "table": ("analysis", "classes", '{ mirna = ["miR-*"] }', {"mirna": ["miR-*"]}),
    "text": ("experiment", "cell_line", "LN229", "LN229"),
}


def test_settings_round_trip_for_every_form_item_kind(srv):
    kinds = {meta["kind"] for keys in settings_form.FORM_ITEMS.values() for meta in keys.values()}
    assert kinds == set(KIND_CASES), "a new form-item kind needs a round-trip case here and a widget in app.js"
    values: dict = {}
    for sec, key, sent, _ in KIND_CASES.values():
        values.setdefault(sec, {})[key] = sent
    code, r = srv.call("/api/settings", body={"experiment": "E", "values": values})
    assert code == 200, r
    code, p = srv.call("/api/project", exp="E")
    for kind, (sec, key, _, want) in KIND_CASES.items():
        assert settings_form.FORM_ITEMS[sec][key]["kind"] == kind
        assert p["settings"][sec][key] == want, kind
    assert all(sec["items"] for sec in p["schema"])  # the schema uses "items"
    # what to re-run: regions re-detect -> a slow run (which also covers the statistics settings)
    assert r["steps"][0]["action"] == "run" and r["steps"][0].get("slow")
    assert {c["key"] for c in r["changed"]} >= {"regions", "classes", "cell_line"}
    # special values of num_or / names / table
    code, r = srv.call("/api/settings", body={"experiment": "E", "values": {
        "segmentation": {"nucleus_diameter_px": "auto"}, "readouts": {"include": "all"}, "analysis": {"classes": "auto"}}})
    p = srv.call("/api/project", exp="E")[1]["settings"]
    assert p["segmentation"]["nucleus_diameter_px"] == "auto" and p["readouts"]["include"] == "all"
    assert p["analysis"]["classes"] == "auto"
    assert srv.call("/api/settings", body={"experiment": "E", "values": {"review": {"jpeg_quality": "500"}}})[0] == 400


def test_rerun_advice():
    base = settings.load(Path("/nonexistent"))
    after = json.loads(json.dumps(base))
    after["analysis"]["primary_test"] = "rm_anova"
    r = gui.rerun_advice(base, after)
    assert [s["action"] for s in r["steps"]] == ["regen"] and r["steps"][0]["options"] == {"rebuild": "all"}
    after["planner"]["max_n"] = 20
    assert [s["action"] for s in gui.rerun_advice(base, after)["steps"]] == ["regen", "sampling"]
    after["positivity"]["lod_k"] = 4.0
    assert gui.rerun_advice(base, after)["steps"][0]["action"] == "run"
    assert gui.rerun_advice(base, base) == {"changed": [], "steps": []}


# ---------------------------------------------------------------- blinded tuning on a real (small) run

@pytest.fixture(scope="module")
def ki_srv(tmp_path_factory):
    from kerncount import discover, pipeline
    from test_intensity import _synthetic_interleaved

    ws = tmp_path_factory.mktemp("kiws")
    raw = ws / "data" / "Ki"
    _synthetic_interleaved(raw)
    prj = ws / "analysis" / "Ki"
    prj.mkdir(parents=True)
    settings.write_template(prj, raw)
    sp = prj / "settings.toml"
    sp.write_text(settings_form.update_text(sp.read_text(), settings_form.coerce({
        "channels": {"files": "interleaved", "file_order": "Ki67, DAPI", "segment": "DAPI"},
        "segmentation": {"backend": "classical"}})))
    discover.write_layout(discover.scan(raw, files="interleaved", n_channels=2), prj / "layout.csv")
    pipeline.run(prj, no_review=True)
    s = Server(ws)
    yield s, prj
    s.close()


def test_tuning_endpoints_are_blinded(ki_srv):
    srv, prj = ki_srv
    code, info = srv.call("/api/tuning/info", exp="Ki")
    assert code == 200 and info["channels"] == ["Ki67"] and info["locks"] == []
    ns = info["scopes"]
    code, rng = srv.call("/api/tuning/range", exp="Ki", channel="Ki67", measure="mean")
    assert code == 200 and rng["min"] < rng["max"]
    code, smp = srv.call("/api/tuning/sample", body={"experiment": "Ki", "channel": "Ki67", "k": 4, "seed": 3})
    assert code == 200 and len(smp["fields"]) == 4 and set(smp) == {"seed", "channel", "fields"}
    ids = [f["field_id"] for f in smp["fields"]]
    code, png = srv.call("/api/tuning/render", exp="Ki", field=ids[0], channel="Ki67", cutoff="100", measure="mean")
    assert code == 200 and png[:8] == b"\x89PNG\r\n\x1a\n"
    code, pv = srv.call("/api/tuning/preview", exp="Ki", cutoff="100", channel="Ki67", measure="mean", fields=",".join(ids))
    assert code == 200 and 0 <= pv["frac_positive"] <= 1 and pv["sample"]["nuclei"] > 0
    code, lk = srv.call("/api/tuning/lock", body={"experiment": "Ki", "channel": "Ki67", "cutoff": 150.0,
                                                  "measure": "mean", "field_ids": ids})
    assert code == 200 and lk["lock"]["cutoff"] == 150.0 and lk["lock"]["fields"] == 4
    code, info2 = srv.call("/api/tuning/info", exp="Ki")
    assert info2["locks"][0]["cutoff"] == 150.0 and info2["history"][-1]["action"] == "lock"
    # nothing a tuning response returns may name a condition, an n or a file
    blob = json.dumps([smp, pv, lk, rng, {k: v for k, v in info2.items() if k != "scopes"}])
    for secret in ("Control", "LyoVec", "cntrl", "lyovec", ".tif", *ns):
        assert secret not in blob, secret
    stored = json.loads((prj / "cutoffs.json").read_text())
    assert stored["locks"][0]["field_files"]  # the record keeps the files; the page never sees them
    assert srv.call("/api/tuning/render", exp="Ki", field="f_nope", channel="Ki67", cutoff="1")[0] == 400
    assert srv.call("/api/tuning/lock", body={"experiment": "Ki", "channel": "Ki67", "cutoff": 1, "field_ids": []})[0] == 400
    code, un = srv.call("/api/tuning/unlock", body={"experiment": "Ki", "channel": "Ki67", "measure": "mean"})
    assert code == 200 and srv.call("/api/tuning/info", exp="Ki")[1]["locks"] == []
    assert srv.call("/api/tuning/unlock", body={"experiment": "Ki", "channel": "Ki67", "measure": "mean"})[0] == 400


def test_review_overlay_and_intensity_page(ki_srv):
    srv, prj = ki_srv
    code, rv = srv.call("/api/review", exp="Ki")
    n = rv[0]["n"]
    t = next(x for x in rv[0]["treatments"] if x["treatment"] == "Control")
    img = t["images"][0]
    assert set(img["sources"]) >= {"all", "Ki67", "DAPI"} and "Ki67" in img["marker"]
    code, png = srv.call("/api/field_overlay", exp="Ki", n=n, treatment="Control", image=str(img["image"]), channel="Ki67")
    assert code == 200 and png[:4] == b"\x89PNG"
    code, it = srv.call("/api/intensity", exp="Ki")
    assert code == 200 and it["manifest"]["channels"] == ["Ki67"] and "locks" in it
    code, ov = srv.call("/api/overview", exp="Ki")
    assert ov["intensity"] and "mispaired" in {a["key"] for a in ov["attention"]}
