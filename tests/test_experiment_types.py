"""Experiment types (ADR 0004): `[experiment] type` in settings.toml selects how an experiment is discovered (layout),
measured, tabulated (readout table) and analysed (statistical design). Omitted = the cell-count / intensity path."""
import pandas as pd
import pytest

from kerncount import experiment_types, settings


def test_type_omitted_selects_cell_count(tmp_path):
    (tmp_path / "settings.toml").write_text('[project]\nraw_dir = "raw"\n')
    t = experiment_types.of(settings.load(tmp_path))
    assert t.name == "cell_count"


def _workspace(tmp_path, monkeypatch, exp="expA"):
    import numpy as np
    import tifffile

    d = tmp_path / "data" / exp / "01.02.26 run" / "cntrl"
    d.mkdir(parents=True)
    tifffile.imwrite(d / "1.tif", np.zeros((16, 16), np.uint16))
    monkeypatch.chdir(tmp_path)
    return tmp_path / "analysis" / exp


def _set_type(prj, name):
    sp = prj / "settings.toml"
    text = sp.read_text()
    old = 'type = "cell_count"'
    sp.write_text(text.replace(old, f'type = "{name}"') if old in text
                  else text.replace("[experiment]", f'[experiment]\ntype = "{name}"', 1))


def test_unknown_type_names_the_valid_types(tmp_path, monkeypatch):
    from kerncount import cli

    prj = _workspace(tmp_path, monkeypatch)
    cli.main(["init", "expA"])
    _set_type(prj, "tumor_volumes")
    for cmd in (["init", "expA", "--force-layout"], ["run", "expA"]):
        with pytest.raises(SystemExit) as e:
            cli.main(cmd)
        msg = str(e.value)
        assert "tumor_volumes" in msg and "valid types: cell_count" in msg, msg


class _Spheres(experiment_types.ExperimentType):
    """A toy type: one row per raw file; its 'measurement' is the file size; analysis writes the mean."""
    name, label = "spheres", "Toy spheres"

    def scan_layout(self, raw, settings):
        files = sorted(p for p in raw.rglob("*.tif"))
        return pd.DataFrame({"unit": [p.stem for p in files], "file": [p.name for p in files]})

    def measure(self, project_dir, layout, settings, options):
        return pd.DataFrame({"unit": layout["unit"], "size": 1.0})

    def readout_table(self, project_dir, layout, measurements, settings):
        table = layout.merge(measurements, on="unit")
        table.to_csv(project_dir / "units.csv", index=False)
        return table

    def analyse(self, project_dir, table, settings, options):
        (project_dir / "mean.txt").write_text(f"{table['size'].mean()} rebuild={options.rebuild}")


def test_selected_type_drives_init_and_run(tmp_path, monkeypatch):
    from kerncount import cli

    monkeypatch.setitem(experiment_types.TYPES, "spheres", _Spheres())
    prj = _workspace(tmp_path, monkeypatch)
    cli.main(["init", "expA"])
    assert "treatment" in pd.read_csv(prj / "layout.csv").columns  # default type until selected
    _set_type(prj, "spheres")
    cli.main(["init", "expA", "--force-layout"])
    assert list(pd.read_csv(prj / "layout.csv").columns) == ["unit", "file"]
    cli.main(["run", "expA", "--rebuild", "g2"])
    assert pd.read_csv(prj / "units.csv").to_dict("records") == [{"unit": 1, "file": "1.tif", "size": 1.0}]
    assert (prj / "mean.txt").read_text() == "1.0 rebuild=['g2']"
    assert not (prj / "counts.csv").exists() and not (prj / "generations").exists()


def test_type_setting_in_defaults_template_and_form():
    import tomllib

    from kerncount import settings_form

    assert settings.DEFAULTS["experiment"]["type"] == "cell_count"
    assert tomllib.loads(settings.TEMPLATE.format(raw_dir="x"))["experiment"]["type"] == "cell_count"
    item = settings_form.FORM_ITEMS["experiment"]["type"]
    assert item["kind"] == "choice" and "cell_count" in item["options"]
    assert settings_form.coerce({"experiment": {"type": "cell_count"}})["experiment"]["type"] == "cell_count"


def _synthetic_plates(raw):
    """Single-channel nuclei, 2 n x 2 treatments x 2 fields: the smallest project that builds every generation output."""
    import numpy as np
    import tifffile

    rng = np.random.default_rng(1)
    yy, xx = np.mgrid[0:128, 0:128]
    for n, day in enumerate(("05.06.26", "12.06.26")):
        for cond, count in (("cntrl", 12), ("lyovec", 6)):
            folder = raw / f"{day} test" / cond
            folder.mkdir(parents=True)
            for k in (1, 2):
                img = np.full((128, 128), 100.0)
                for _ in range(count + n + k):
                    cy, cx = rng.uniform(12, 116, 2)
                    img[(yy - cy) ** 2 + (xx - cx) ** 2 <= 36] += 800
                tifffile.imwrite(folder / f"{k}.tif", (img + rng.normal(0, 10, img.shape)).clip(0).astype(np.uint16))


def test_run_without_type_is_the_original_pipeline_call(tmp_path, monkeypatch):
    """Fast guard (the end-to-end comparison below is slow): with no type, `kerncount run` makes exactly the call
    it made before experiment types existed."""
    from kerncount import cli, pipeline

    calls = []
    monkeypatch.setattr(pipeline, "run", lambda prj, **kw: calls.append((prj, kw)))
    prj = _workspace(tmp_path, monkeypatch)
    cli.main(["init", "expA"])
    _set_type(prj, "")  # an empty type is the default too
    cli.main(["run", "expA", "--rebuild", "g2", "g3", "--resegment", "--no-review"])
    cli.main(["run", "expA"])
    assert calls == [(prj.resolve(), dict(rebuild=["g2", "g3"], resegment=True, no_review=True)),
                     (prj.resolve(), dict(rebuild=[], resegment=False, no_review=False))]


@pytest.mark.slow
def test_existing_project_runs_identically(tmp_path):
    """A project written before experiment types existed (no `type` line) gives the same outputs through
    `kerncount run` as the original pipeline.run. p-values from Monte-Carlo integration (Dunnett) vary between any
    two runs, so numbers are compared to 1e-2 relative."""
    from kerncount import cli, discover, pipeline, settings_form

    raw = tmp_path / "raw"
    _synthetic_plates(raw)
    prjs = {}
    for name in ("via_cli", "direct"):
        prj = prjs[name] = tmp_path / name
        prj.mkdir()
        settings.write_template(prj, raw)
        sp = prj / "settings.toml"
        text = settings_form.update_text(sp.read_text(), settings_form.coerce({"segmentation": {"backend": "classical"}}))
        sp.write_text("\n".join(l for l in text.splitlines() if not l.startswith("type =")) + "\n")
        assert "type" not in settings.tomllib.loads(sp.read_text())["experiment"]
        discover.write_layout(discover.scan(raw), prj / "layout.csv")
    cli.main(["run", str(prjs["via_cli"]), "--no-review"])
    pipeline.run(prjs["direct"], no_review=True)

    def outputs(prj):
        return sorted(p.relative_to(prj).as_posix() for p in prj.rglob("*") if p.is_file() and "_work" not in p.parts)

    files = outputs(prjs["direct"])
    assert outputs(prjs["via_cli"]) == files and "counts.csv" in files
    assert any(f.startswith("generations/g1/all_n/") for f in files)
    for f in files:
        if f.endswith(".csv"):
            pd.testing.assert_frame_equal(pd.read_csv(prjs["via_cli"] / f, sep=None, engine="python"),
                                          pd.read_csv(prjs["direct"] / f, sep=None, engine="python"),
                                          check_exact=False, rtol=1e-2, obj=f)
