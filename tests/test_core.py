from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellcount import analysis, discover, master, segment


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
