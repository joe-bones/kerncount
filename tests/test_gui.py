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
