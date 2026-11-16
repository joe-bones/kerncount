"""Randomised plate maps: fit checks, reproducibility, balance across plates, edge avoidance, rendering, CLI."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kerncount import platemap

CONDS = ["Control", "LyoVec", "CL075", "miR-a", "miR-b", "miR-c", "miR-d", "miR-e"]


def test_fits_and_columns():
    plates = platemap.generate(CONDS, 24, 2, n_plates=2, seed=1)
    assert len(plates) == 2
    for p in plates:
        assert list(p.columns) == platemap.COLUMNS and len(p) == 24
        used = p[p["condition"] != platemap.EMPTY]
        assert used.groupby("condition").size().to_dict() == {c: 2 for c in CONDS}
        assert (p.loc[p["condition"] == platemap.EMPTY, "replicate"] == 0).all()
        assert sorted(used[used["condition"] == "CL075"]["replicate"]) == [1, 2]
        assert p["plate_well"].is_unique and set(p["row"]) == set("ABCD") and set(p["col"]) == set(range(1, 7))


def test_dict_wells_and_overfull_errors():
    counts = {c: (5 if c in ("Control", "LyoVec") else 2) for c in CONDS}
    p = platemap.generate(CONDS, 24, counts, seed=2)[0]
    assert (p["condition"] == "Control").sum() == 5 and (p["condition"] == platemap.EMPTY).sum() == 2
    with pytest.raises(ValueError, match="only 24"):
        platemap.generate(CONDS, 24, 4)
    with pytest.raises(ValueError, match="no entry"):
        platemap.generate(CONDS, 24, {"Control": 2})
    with pytest.raises(ValueError, match="plate format"):
        platemap.generate(CONDS, 30, 1)


def test_seed_reproducible_and_randomised():
    a = platemap.generate(CONDS, 24, 2, n_plates=3, seed=7)
    b = platemap.generate(CONDS, 24, 2, n_plates=3, seed=7)
    c = platemap.generate(CONDS, 24, 2, n_plates=3, seed=8)
    assert all(x.equals(y) for x, y in zip(a, b))
    assert not all(x.equals(y) for x, y in zip(a, c))
    assert not a[0]["condition"].equals(a[1]["condition"])  # every plate randomised separately
    assert isinstance(platemap.generate(CONDS, 24, 2)[0].attrs["seed"], int)  # no seed: a recorded random one


def test_avoid_edge():
    p = platemap.generate(["A", "B", "C", "D"], 24, 2, n_plates=3, avoid_edge=True, seed=1)
    for x in p:  # 8 wells needed, 8 inner wells: no edge well used
        assert not x.loc[x["condition"] != platemap.EMPTY, "edge"].any()
    p = platemap.generate(CONDS, 24, 2, n_plates=2, avoid_edge=True, seed=1)
    used = p[0][p[0]["condition"] != platemap.EMPTY]
    assert used.loc[~used["edge"]].shape[0] == len(platemap.inner_wells(24))  # all inner wells used first
    assert "edge wells are used" in p[0].attrs["notes"][0]


def test_balance_spreads_edge_wells_across_plates():
    def spread(balance):
        sds = []
        for seed in range(5):
            eu = platemap.edge_usage(platemap.generate(CONDS, 24, 2, n_plates=4, seed=seed, balance=balance))
            sds.append(eu["edge_wells"].std())
        return np.mean(sds)

    assert spread(True) < spread(False)
    eu = platemap.edge_usage(platemap.generate(CONDS, 24, 2, n_plates=4, seed=0))
    assert eu["edge_wells"].max() - eu["edge_wells"].min() <= 2
    assert (eu["rows_used"] >= 3).all()


@pytest.mark.parametrize("fmt", [6, 12, 24, 48, 96])
def test_every_format_generates_and_renders(fmt, tmp_path):
    conds = CONDS[: min(len(CONDS), fmt // 2)]
    plates = platemap.generate(conds, fmt, 2 if fmt > 6 else 1, n_plates=2, avoid_edge=fmt >= 48, seed=3)
    assert len(plates[1]) == fmt
    paths = platemap.render(plates[1], tmp_path / f"p{fmt}.png")
    assert all(p.exists() and p.stat().st_size > 1000 for p in paths)
    assert {p.suffix for p in paths} == {".png", ".pdf"}


def test_plan_project_writes_folder(tmp_path):
    from kerncount import settings as st_mod

    prj = tmp_path / "exp"
    prj.mkdir()
    st_mod.write_template(prj, tmp_path / "raw")
    pd.DataFrame(dict(n=["n1"] * 4, treatment=["Control", "LyoVec", "miR-x", "CL075"], image=[1] * 4,
                      auto_count=[1] * 4)).to_csv(prj / "counts.csv", index=False)
    out = platemap.plan_project(prj, plates=2, plate_format=12, wells=2, reference_wells="auto", seed=5)
    assert out.parent == prj / "plans" and out.name.startswith("platemap_")
    for f in ("plate_1.csv", "plate_2.png", "plate_2.pdf", "all_plates.csv", "edge_usage.csv", "README.txt"):
        assert (out / f).exists(), f
    p1 = pd.read_csv(out / "plate_1.csv")
    # auto: 2 x sqrt(3) = 3.5 -> 3 per reference; 2 x 2 + 2 x 3 = 10 of 12 wells
    assert (p1["condition"] == "Control").sum() == 3 and (p1["condition"] == "miR-x").sum() == 2
    assert "--seed 5" in (out / "README.txt").read_text()
    assert platemap.project_conditions(prj, st_mod.load(prj))[:2] == ["Control", "LyoVec"]


def test_cli_platemap(tmp_path, monkeypatch):
    from kerncount import cli
    from kerncount import settings as st_mod

    prj = tmp_path / "analysis" / "exp"
    prj.mkdir(parents=True)
    (tmp_path / "data" / "exp").mkdir(parents=True)
    st_mod.write_template(prj, tmp_path / "data" / "exp")
    pd.DataFrame(dict(n=["n1"] * 3, treatment=["Control", "LyoVec", "miR-x"], image=[1] * 3)).to_csv(
        prj / "counts.csv", index=False)
    monkeypatch.setenv("KERNCOUNT_WORKSPACE", str(tmp_path))
    cli.main(["platemap", "exp", "--plates", "2", "--format", "12", "--avoid-edge", "--seed", "1"])
    assert len(list((prj / "plans").glob("platemap_*/plate_2.png"))) == 1
    with pytest.raises(SystemExit):
        cli.main(["platemap", "exp", "--format", "6", "--wells", "3", "--reference-wells", "3"])
