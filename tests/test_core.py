from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import analysis, discover, master, segment


@pytest.mark.parametrize("raw, expected", [
    ("hsamiRNA 21-5p", "miR-21-5p"), ("hsa mi RNA 155-5p", "miR-155-5p"), ("21 5p", "miR-21-5p"),
    ("hsa 30a-3p", "miR-30a-3p"), ("30a 3p", "miR-30a-3p"), ("LyoVec veh. cntrl", "LyoVec"), ("lyo vec", "LyoVec"),
    ("010125 cntrl", "Control"), ("-PDL cntrl", "Control"), ("CL075", "CL075"),
])
def test_canonical_treatment(raw, expected):
    assert discover.canonical_treatment(raw)[0] == expected


def test_dates_and_image_numbers():
    assert str(discover.parse_date("05.06.26,cells,DAPI count")) == "2026-06-05"
    assert str(discover.parse_date("06.10.26 cells dapi count, 72h")) == "2026-10-06"
    assert discover.image_number(Path("`10.tif")) == 10
    assert discover.image_number(Path("03.tif")) == 3
    assert discover.image_number(Path("5 zoom.tif")) is None


def _auto(n=3):
    return pd.DataFrame(dict(n="n1", treatment="A", image=range(1, n + 1), auto_count=[10, 20, 30][:n],
                             qc_ok="yes", qc_level="high"))


def test_master_keeps_manual_and_pads(tmp_path):
    p = tmp_path / "counts.csv"
    df, _ = master.merge(p, _auto(), [], ",", tmp_path / "bk")
    assert master.manual_columns(df) == ["manual_1", "manual_2", "manual_3"]
    df.loc[1, "manual_1"] = "25"
    df.loc[2, "manual_2"] = "x"
    df.loc[0, "manual_3"] = "11"
    df.to_csv(p, index=False, sep=";")  # as German Excel would save it
    auto2 = _auto()
    auto2.loc[0, "auto_count"] = 12  # re-segmentation changed auto counts
    df2, _ = master.merge(p, auto2, [], ",", tmp_path / "bk")
    assert df2.loc[1, "manual_1"] == "25" and df2.loc[2, "manual_2"] == "x"
    assert master.manual_columns(df2)[-1] == "manual_4"  # always an empty column at the end
    assert list((tmp_path / "bk").iterdir())


def test_generation_chain():
    df = _auto().assign(manual_1=["", "25", ""], manual_2=["", "", "x"], manual_3="", manual_4=["99", "", ""])
    gens = master.chain(df)
    assert [g["name"] for g in gens] == ["g1", "g2", "g3"]  # stops at first empty column (manual_3)
    g1, g2, g3 = (g["table"] for g in gens)
    assert g1["count"].tolist() == [10, 20, 30]
    assert g2["count"].tolist() == [10, 25, 30] and not g2["excluded"].any()
    assert g3["count"].tolist() == [10, 25, 30] and g3["excluded"].tolist() == [False, False, True]
    assert len({g["fingerprint"] for g in gens}) == 3


def test_bad_manual_value():
    df = _auto().assign(manual_1=["abc", "", ""])
    with pytest.raises(ValueError, match="manual value"):
        master.chain(df)


def test_merge_oversplit_joins_lobes_not_neighbours():
    img = np.zeros((120, 200), np.float32)
    yy, xx = np.mgrid[:120, :200]
    lobe1 = (yy - 60) ** 2 + (xx - 40) ** 2 < 15 ** 2
    lobe2 = (yy - 60) ** 2 + (xx - 62) ** 2 < 15 ** 2
    img[lobe1 | lobe2] = 1000  # one bilobed nucleus: no dark seam
    a = (yy - 60) ** 2 + (xx - 130) ** 2 < 15 ** 2
    b = (yy - 60) ** 2 + (xx - 161) ** 2 < 15 ** 2
    img[a | b] = 1000
    img[:, 145:147] = 100  # two touching nuclei with a dark seam
    labels = np.zeros(img.shape, np.int32)
    labels[lobe1 & (xx <= 51)] = 1
    labels[lobe2 & (xx > 51)] = 2
    labels[a & (xx <= 145)] = 3
    labels[b & (xx >= 147)] = 4
    out, _, merged = segment.merge_oversplit(labels, img, expected_area=np.pi * 22 ** 2, probs=None)
    assert merged == 1 and out.max() == 3


def test_holm_and_stars():
    assert analysis.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])
    assert analysis.stars(0.00001) == "****" and analysis.stars(0.2) == "ns"


def test_significance_marks_only_where_significant(tmp_path):
    from kerncount import settings as st_mod
    rng = np.random.default_rng(0)
    rows = []
    for n in ["n1", "n2", "n3", "n4"]:
        for t, mu in [("LyoVec", 100), ("Strong", 30), ("Null", 100)]:
            for i in range(1, 13):
                rows.append(dict(n=n, treatment=t, image=i, value=rng.normal(mu, 8)))
    df = pd.DataFrame(rows)
    s = st_mod.DEFAULTS
    st = analysis.per_n(df, "count", "nuclei", "n1", tmp_path, s)
    marks = analysis.significance_marks(st, "LyoVec", "Dunnett")
    assert "Strong" in marks and "Null" not in marks
    _, st_all = analysis.all_n(df, "count", "nuclei", tmp_path, s)
    marks_all = analysis.significance_marks(st_all, "LyoVec", "vs LyoVec (")
    assert "Strong" in marks_all and "Null" not in marks_all
    assert (tmp_path / "count.png").exists() and (tmp_path / "count_pct_of_LyoVec.png").exists()


def test_workspace_init_by_name(tmp_path, monkeypatch):
    import tifffile
    from kerncount import cli, settings as st_mod

    for exp in ("expA", "expB"):
        d = tmp_path / "data" / exp / "01.02.26 run" / "cntrl"
        d.mkdir(parents=True)
        tifffile.imwrite(d / "1.tif", np.zeros((16, 16), np.uint16))
    monkeypatch.chdir(tmp_path)
    cli.main(["init", "expA"])
    prj = tmp_path / "analysis" / "expA"
    assert (prj / "layout.csv").exists()
    s = st_mod.load(prj)
    assert s["project"]["raw_dir"] == "../../data/expA"  # relative: the workspace can be moved
    assert cli.resolve_project("expA") == prj.resolve()
    cli.main(["init", "--all"])
    assert (tmp_path / "analysis" / "expB" / "settings.toml").exists()


def _sampling_settings(**design):
    from kerncount import settings
    s = settings.load(Path("/nonexistent"))
    s["design"].update(design)
    return s


def test_sampling_design_accepts_old_plate_names():
    from kerncount import sampling
    s = _sampling_settings()
    assert sampling.design(s)[:3] == (2, 6, "consecutive") and not sampling.design(s)[3]
    s["design"].update(plates=3, fields_per_plate=4, plate_layout="alternating")
    w, k, lay, notes = sampling.design(s)
    assert (w, k, lay) == (3, 4, "alternating") and len(notes) == 3
    assert [sampling.well_of(i, 6, 2, "consecutive") for i in (1, 6, 7, 12)] == [1, 1, 2, 2]
    assert [sampling.well_of(i, 6, 2, "alternating") for i in (1, 2, 3)] == [1, 2, 1]


def test_sampling_batch_checks_and_missing_reference(tmp_path):
    from kerncount import sampling
    rng = np.random.default_rng(1)
    rows = []
    for n, strength in [("n1", 0.5), ("n2", 0.5), ("n3", 0.15), ("n4", 0.5)]:
        for t, base in [("Control", 80), ("LyoVec", 80), ("miR-a", 80 * strength), ("miR-b", 80 * strength)]:
            if n == "n2" and t == "LyoVec":
                continue  # reference missing in one n
            for img in range(1, 13):
                f = 4.0 if (n == "n4" and t == "miR-a" and img > 6) else 1.0  # one odd well
                rows.append(dict(n=n, treatment=t, image=img, value=max(base * f * rng.lognormal(0, 0.15), 1)))
    df = pd.DataFrame(rows)
    s = _sampling_settings()
    summary = sampling.run(df, s, tmp_path, iterations=20)
    assert set(summary["reference"]) == {"LyoVec", "Control"}
    eff = pd.read_csv(tmp_path / "effects.csv")
    assert eff.loc[eff["reference"] == "LyoVec", "n_pairs"].max() == 3  # n2 left out vs LyoVec ...
    assert eff.loc[eff["reference"] == "Control", "n_pairs"].max() == 4  # ... but kept vs Control
    by_n = pd.read_csv(tmp_path / "batch_by_n.csv").set_index("n")
    assert by_n.loc["n2", "missing_references"] == "LyoVec"
    assert by_n.loc["n3", "shift_vs_Control_flag"] == "check" and by_n.loc["n3", "shift_vs_Control"] < 1
    cond = pd.read_csv(tmp_path / "batch_by_condition.csv").fillna("")
    flagged = cond[cond["well_flag"] == "check"][["n", "treatment"]].values.tolist()
    assert ["n4", "miR-a"] in flagged
    rep = (tmp_path / "REPORT.txt").read_text()
    assert "n2 has no LyoVec" in rep and "plate" not in rep.replace("one plate per n", "").replace(
        "within a plate", "").replace("plate, day", "").replace("Plate-wide", "")


def test_default_reference_is_control_and_both_controls_get_tests(tmp_path):
    from kerncount import settings as st_mod
    s = st_mod.DEFAULTS
    assert analysis.control_groups(s["analysis"])[:2] == ["Control", "LyoVec"]
    rng = np.random.default_rng(1)
    rows = [dict(n=n, treatment=t, image=i, value=rng.normal(mu, 5))
            for n in ["n1", "n2", "n3"] for t, mu in [("Control", 100), ("LyoVec", 80), ("X", 40)] for i in range(1, 13)]
    df = pd.DataFrame(rows)
    st = analysis.per_n(df, "count", "nuclei", "n1", tmp_path, s)
    assert {"Dunnett vs Control (fields, descriptive)", "Dunnett vs LyoVec (fields, descriptive)"} <= set(st["test"])
    assert (tmp_path / "count.png").exists() and (tmp_path / "count_vs_LyoVec.png").exists()
    assert (tmp_path / "count_twoway.png").exists()  # treatment x well
    _, st_all = analysis.all_n(df, "count", "nuclei", tmp_path, s)
    assert st_all["test"].str.contains("vs Control (", regex=False).any() and st_all["test"].str.contains("vs LyoVec (", regex=False).any()
    aov = pd.read_csv(tmp_path / "count_twoway_anova.csv")
    assert list(aov["source"]) == ["treatment", "n", "treatment x n (interaction)", "residual"]


def test_two_way_one_block_matches_scipy_tukey():
    from scipy import stats
    rng = np.random.default_rng(2)
    groups = {"A": rng.normal(10, 2, 6), "B": rng.normal(12, 2, 8), "C": rng.normal(15, 2, 7)}
    df = pd.DataFrame([dict(treatment=t, blk="b1", value=v) for t, vs in groups.items() for v in vs])
    anova, pairs, _ = analysis.two_way(df, "blk", ["A", "B", "C"])
    f, p = stats.f_oneway(*groups.values())
    assert anova.loc[anova["source"] == "treatment", "F"].iloc[0] == pytest.approx(f)
    ref = stats.tukey_hsd(*groups.values())
    got = pairs.set_index("comparison")["p_adj"]
    assert got["A vs B"] == pytest.approx(ref.pvalue[0, 1], abs=1e-4)
    assert got["A vs C"] == pytest.approx(ref.pvalue[0, 2], abs=1e-4)
    assert got["B vs C"] == pytest.approx(ref.pvalue[1, 2], abs=1e-4)


def test_two_way_handles_missing_group_and_interaction():
    rng = np.random.default_rng(3)
    rows = []
    for n, shift in [("n1", 0), ("n2", 20), ("n3", -10)]:
        for t, mu in [("Control", 100), ("LyoVec", 100), ("X", 60)]:
            if n == "n2" and t == "LyoVec":
                continue  # this n has no LyoVec
            rows += [dict(treatment=t, n=n, value=v) for v in rng.normal(mu + shift, 5, 12)]
    anova, pairs, cm = analysis.two_way(pd.DataFrame(rows), "n", ["Control", "LyoVec", "X"])
    p = pairs.set_index("comparison")
    assert p.at["Control vs LyoVec", "matched_levels"] == 2 and p.at["Control vs X", "matched_levels"] == 3
    assert p.at["Control vs X", "p_adj"] < 1e-4 and p.at["Control vs LyoVec", "p_adj"] > 0.05
    assert anova.loc[anova["source"] == "treatment", "p"].iloc[0] < 1e-4
    assert cm.loc["LyoVec"].isna().sum() == 1


def test_dunnett_matches_scipy():
    from scipy import stats
    rng = np.random.default_rng(4)
    k, n = 5, 6
    ctrl = rng.normal(10, 2, n)
    others = [rng.normal(10 + d, 2, n) for d in (0, 1, 2.5, 4)]
    ref = stats.dunnett(*others, control=ctrl)
    mse = np.mean([np.var(g, ddof=1) for g in [ctrl] + others])
    t = np.array([(g.mean() - ctrl.mean()) / np.sqrt(2 * mse / n) for g in others])
    R = np.full((4, 4), 0.5) + 0.5 * np.eye(4)
    padj, _ = analysis._dunnett_padj(t, R, k * n - k)
    assert padj == pytest.approx(ref.pvalue, abs=2e-3)


def test_all_n_prism_standard_drives_stars(tmp_path):
    from kerncount import settings as st_mod
    rng = np.random.default_rng(6)
    rows = []
    for n, day in [("n1", 1.0), ("n2", 1.5), ("n3", 0.7), ("n4", 1.2)]:
        for t, mu in [("Control", 100), ("LyoVec", 100), ("Drop", 50), ("Null", 100)]:
            if n == "n2" and t == "LyoVec":
                continue
            rows += [dict(n=n, treatment=t, image=i, value=v) for i, v in enumerate(rng.normal(mu * day, 6, 12), 1)]
    _, st = analysis.all_n(pd.DataFrame(rows), "count", "nuclei", tmp_path, st_mod.DEFAULTS)
    d = st[st["test"].str.startswith("Dunnett vs Control")].set_index("comparison")
    assert d.at["Drop vs Control", "p_adj"] < 0.05 and d.at["Null vs Control", "p_adj"] > 0.05
    assert d.at["LyoVec vs Control", "N"] == 3  # n2 has no LyoVec but still counts for the others
    assert st["test"].str.startswith("Dunnett vs LyoVec").any()
    assert st["test"].str.contains("ratio paired t-test", regex=False).any()
    assert "Drop" in analysis.significance_marks(st, "Control", "Dunnett vs Control (")


# ---------------------------------------------------------------- I5: generations swapped in by rename

def test_swap_in_keeps_user_files_and_drops_known_outputs(tmp_path):
    from kerncount import fileio

    final, tmp = tmp_path / "g1", tmp_path / "g1.tmp"
    for rel, txt in (("stats_long.csv", "old"), ("all_n/old_plot.png", "x"), ("my_notes.txt", "mine"),
                     ("_kept/earlier.txt", "kept before")):
        (final / rel).parent.mkdir(parents=True, exist_ok=True)
        (final / rel).write_text(txt)
    tmp.mkdir()
    (tmp / "stats_long.csv").write_text("new")
    kept = fileio.swap_in(tmp, final, known={"stats_long.csv", "all_n/old_plot.png"})
    assert kept == ["my_notes.txt"] and not tmp.exists()
    assert (final / "stats_long.csv").read_text() == "new"
    assert not (final / "all_n" / "old_plot.png").exists()  # a known output that is no longer written
    assert (final / "_kept" / "my_notes.txt").read_text() == "mine"
    assert (final / "_kept" / "earlier.txt").read_text() == "kept before"
    assert not any(p.name.startswith(".g1") for p in tmp_path.iterdir())


def test_swap_in_reports_a_file_open_in_excel(tmp_path, monkeypatch):
    import os

    from kerncount import fileio

    final, tmp = tmp_path / "g2", tmp_path / "g2.tmp"
    final.mkdir()
    (final / "results.xlsx").write_text("old")
    tmp.mkdir()
    (tmp / "results.xlsx").write_text("new")
    real = os.replace

    def locked(src, dst):  # Windows: renaming a folder with an open file fails
        if Path(src) == final:
            raise PermissionError(13, "in use")
        return real(src, dst)

    monkeypatch.setattr(fileio.os, "replace", locked)
    with pytest.raises(fileio.FolderInUse, match="open in Excel"):
        fileio.swap_in(tmp, final, known={"results.xlsx"})
    assert (final / "results.xlsx").read_text() == "old" and (tmp / "results.xlsx").read_text() == "new"


def test_failed_generation_build_leaves_the_old_one(tmp_path, monkeypatch):
    from kerncount import generations

    final = tmp_path / "generations" / "g1"
    final.mkdir(parents=True)
    (final / "stats_long.csv").write_text("old")

    def boom(out, *a, **k):
        (out / "half.csv").write_text("half")
        raise RuntimeError("model crashed")

    monkeypatch.setattr(generations, "_build_into", boom)
    with pytest.raises(RuntimeError):
        generations.build(tmp_path, dict(name="g1"), None, None, {})
    assert (final / "stats_long.csv").read_text() == "old"
    assert not (tmp_path / "generations" / "g1.tmp").exists()


def test_atomic_write_keeps_permissions(tmp_path):
    from kerncount import fileio

    p = tmp_path / "settings.toml"
    p.write_text("a = 1\n")
    p.chmod(0o640)
    fileio.atomic_write(p, "a = 2\n")
    assert p.read_text() == "a = 2\n" and (p.stat().st_mode & 0o777) == 0o640
    assert [x.name for x in tmp_path.iterdir()] == ["settings.toml"]  # no temp file left


# ---------------------------------------------------------------- workspace paths, backups, model cache

def test_resolve_raw_prefers_data_and_keeps_symlinks(tmp_path, monkeypatch):
    from kerncount import cli, settings

    ws = tmp_path / "ws"
    real = tmp_path / "elsewhere" / "exp"
    real.mkdir(parents=True)
    (ws / "data").mkdir(parents=True)
    (ws / "data" / "exp").symlink_to(real)
    (ws / "exp").mkdir()  # a same-named folder in the working directory must not win over data/exp
    monkeypatch.chdir(ws)
    monkeypatch.delenv("KERNCOUNT_WORKSPACE", raising=False)
    raw = cli.resolve_raw("exp")
    assert raw == ws.resolve() / "data" / "exp" and raw.is_symlink()  # the link, not its target
    assert cli.resolve_raw("./exp") == ws.resolve() / "exp"  # a path is used as given
    prj = ws / "analysis" / "exp"
    prj.mkdir(parents=True)
    settings.write_template(prj, raw)
    assert settings.load(prj)["project"]["raw_dir"] == "../../data/exp"


def test_overwrites_keep_timestamped_backups(tmp_path):
    from kerncount import design, settings

    prj = tmp_path / "prj"
    prj.mkdir()
    assert settings.write_template(prj, tmp_path / "raw") is None  # nothing to back up yet
    (prj / "settings.toml").write_text((prj / "settings.toml").read_text() + "\n# my edit\n")
    bk = settings.write_template(prj, tmp_path / "raw")
    assert bk.parent == prj / "_work" / "backups" and "# my edit" in bk.read_text()
    assert "# my edit" not in (prj / "settings.toml").read_text()
    bk2 = settings.backup(prj, "settings.toml")
    assert bk2 != bk and bk2.exists()  # same second: a second copy, never overwritten
    table = pd.DataFrame(dict(n="n1", treatment="A", image=range(1, 13)))
    s = settings.load(prj)
    design.write_wells_template(prj, table, s)
    (prj / "wells.csv").write_text((prj / "wells.csv").read_text() + "# hand edit\n")
    assert design.write_wells_template(prj, table, s).name == "wells.new.csv"  # default: never overwrite
    design.write_wells_template(prj, table, s, overwrite=True)
    kept = list((prj / "_work" / "backups").glob("wells_*.csv"))
    assert len(kept) == 1 and "# hand edit" in kept[0].read_text()


def test_stardist_weights_fall_back_to_a_user_cache(tmp_path, monkeypatch):
    import sys
    import types

    from kerncount import segment

    pkg, user = tmp_path / "pkg_models", tmp_path / "home_models"
    pkg.mkdir()
    pkg.chmod(0o500)  # installed read-only
    monkeypatch.setattr(segment, "MODEL_DIR", pkg)
    monkeypatch.setattr(segment, "USER_MODEL_DIR", user)
    try:
        assert segment.model_cache_dir() == user
        dl = tmp_path / "keras_cache" / "2D_versatile_fluo"
        dl.mkdir(parents=True)
        (dl / "config.json").write_text("{}")
        loaded = []

        class Fake:
            def __init__(self, config=None, name=None, basedir=None):
                self.logdir, self.basedir = dl, basedir
                loaded.append(basedir)

            @classmethod
            def from_pretrained(cls, name):
                return cls()

        mod = types.ModuleType("stardist.models")
        mod.StarDist2D = Fake
        monkeypatch.setitem(sys.modules, "stardist", types.ModuleType("stardist"))
        monkeypatch.setitem(sys.modules, "stardist.models", mod)
        monkeypatch.setattr(segment, "_MODEL", None)
        segment._stardist_model()  # downloads (fake) and saves a copy in the user cache
        assert (user / "2D_versatile_fluo" / "config.json").exists() and loaded == [None]
        monkeypatch.setattr(segment, "_MODEL", None)
        segment._stardist_model()  # next time: loaded from the user cache, no download
        assert loaded[-1] == str(user)
    finally:
        pkg.chmod(0o700)


def test_concurrent_lock_writes_do_not_lose_updates(tmp_path):
    """cutoffs.json is read-modify-written under a per-file lock with an atomic replace: 16 threads locking
    different channels at once keep all 16 locks, and no temporary file is left behind."""
    import json
    import threading

    from kerncount import positivity

    ts = [threading.Thread(target=positivity.save_lock, args=(tmp_path, f"C{i}", float(i), ["a", "b"]))
          for i in range(16)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    data = json.loads((tmp_path / "cutoffs.json").read_text())
    assert sorted(lk["channel"] for lk in data["locks"]) == sorted(f"C{i}" for i in range(16))
    assert all(lk["fields"] == 2 for lk in data["locks"])
    assert [p.name for p in tmp_path.iterdir()] == ["cutoffs.json"]
