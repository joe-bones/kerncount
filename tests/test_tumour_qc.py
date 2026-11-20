"""Ticket 06: QC flags (truncation = xy-edge contact or mid-section outside the stack; saturation), flow table,
hollow markers, sensitivity without flagged tumours, blinded review sheets, corrections, A-B hemisphere check."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

VOX = (2.0, 0.6239, 0.6239)  # (z, y, x) µm


def _ellipsoid(shape, centre_vox, semi_vox) -> np.ndarray:
    z, y, x = np.ogrid[: shape[0], : shape[1], : shape[2]]
    (cz, cy, cx), (a, b, c) = centre_vox, semi_vox
    return ((z - cz) / a) ** 2 + ((y - cy) / b) ** 2 + ((x - cx) / c) ** 2 <= 1


def _labels(*objects: np.ndarray) -> np.ndarray:
    lab = np.zeros(objects[0].shape, np.uint8)
    for i, o in enumerate(objects, 1):
        lab[o & (lab == 0)] = i
    return lab


def _qc(labels, stack=None, **tumour):
    from kerncount.tumour import qc

    stack = np.full(labels.shape, 100, np.uint8) if stack is None else stack
    return qc.tumour_qc(labels, stack, VOX, {"tumour": tumour})


# ---------------------------------------------------------------- seam 1: cached labels (+ stack) -> per-tumour QC

def test_a_tumour_inside_the_field_with_a_centred_mid_section_is_not_flagged():
    shape = (40, 200, 200)
    main = _ellipsoid(shape, (20, 100, 100), (30, 60, 50))  # extends beyond the stack in z: never a flag by itself
    q = _qc(_labels(main), mid_section_um=20.0)
    assert q["edge_contact"] == "no" and q["edge_outline_frac"] == 0.0
    assert q["mid_section_fits"] == "yes" and q["truncated"] == "no" and q["flagged"] == "no"
    assert q["qc_flags"] == ""


def test_a_tumour_cut_by_the_image_edge_is_flagged_with_the_fraction_of_its_outline_on_the_edge():
    shape = (30, 200, 200)
    r = 60
    main = _ellipsoid(shape, (15, 100, 0), (100, r, r))  # disc cut in half by the left edge in every layer
    q = _qc(_labels(main), mid_section_um=20.0)
    assert q["edge_contact"] == "yes" and q["truncated"] == "yes" and q["flagged"] == "yes"
    half_disc = 2 * r / (math.pi * r + 2 * r)  # straight cut / (half circle + cut) = 0.389
    assert q["edge_outline_frac"] == pytest.approx(half_disc, abs=0.02)
    assert q["edge_outline_frac_mid"] == pytest.approx(half_disc, abs=0.02)
    assert "image edge" in q["qc_flags"]


def test_a_tumour_grazing_the_image_edge_is_flagged_without_claiming_0_percent_of_its_outline():
    """Real data (cntrl -clodro n1A): contact over a few pixels rounds to 0 % — the flag must not read '0%'."""
    shape = (30, 200, 200)
    main = _ellipsoid(shape, (15, 100, 139), (100, 60, 60))  # reaches the last column at one point
    assert main[:, :, -1].any()
    q = _qc(_labels(main), mid_section_um=20.0)
    assert q["edge_contact"] == "yes" and q["edge_outline_frac"] < 0.005
    assert "(<1% of its outline)" in q["qc_flags"] and "(0%" not in q["qc_flags"]


def test_a_mid_section_that_does_not_fit_in_the_imaged_layers_is_flagged():
    shape = (30, 200, 200)
    main = _ellipsoid(shape, (2, 100, 100), (25, 60, 60))  # largest cross-section in layer 3
    q = _qc(_labels(main), mid_section_um=20.0)
    assert q["edge_contact"] == "no"
    assert q["mid_section_fits"] == "no" and q["truncated"] == "yes" and q["flagged"] == "yes"
    assert "mid-section" in q["qc_flags"]


def test_satellites_do_not_count_for_edge_contact_and_a_satellite_on_the_edge_is_ignored():
    shape = (30, 200, 200)
    main = _ellipsoid(shape, (15, 100, 100), (40, 40, 40))
    sat = _ellipsoid(shape, (15, 5, 5), (3, 6, 6))
    q = _qc(_labels(main, sat), mid_section_um=20.0)
    assert q["edge_contact"] == "no" and q["flagged"] == "no"


def test_saturation_of_the_main_mass_in_the_mid_section_is_reported_with_a_warning_above_the_threshold():
    shape = (40, 120, 120)
    main = _ellipsoid(shape, (20, 60, 60), (40, 40, 40))
    stack = np.full(shape, 120, np.uint8)
    stack[:, :, 60:] = 255  # right half of the main mass saturated (and everything outside the mid-section too)
    q = _qc(_labels(main), stack, mid_section_um=20.0)
    assert q["saturated_frac_mid_section"] == pytest.approx(0.5, abs=0.03)
    assert q["saturation_warning"] == "yes" and "saturat" in q["qc_flags"]
    assert q["flagged"] == "no"  # a warning, not a truncation
    q = _qc(_labels(main), stack, mid_section_um=20.0, qc_saturation_warn=0.6)
    assert q["saturation_warning"] == "no"


@pytest.mark.parametrize("vz, outer", [(2.36, (10 - 3.5 * 2.36) / 2.36), (1.75, (10 - 5.5 * 1.75) / 1.75)])
def test_mid_section_saturation_weights_the_outermost_layers_like_the_readouts(vz, outer):
    """Review fix: QC on the mid-section uses the same exact-thickness layer weights as the readouts."""
    from kerncount.tumour import qc

    shape = (40, 80, 80)
    main = np.zeros(shape, bool)
    main[:, 20:60, 20:60] = True  # z-invariant: all layers tie, centre = middle layer 20
    stack = np.full(shape, 100, np.uint8)
    half = int(10 / vz + 0.5)  # outermost layer with a weight > 0
    stack[20 + half] = 255
    q = qc.tumour_qc(_labels(main), stack, (vz, 0.6239, 0.6239), {"tumour": {"mid_section_um": 20.0}})
    assert q["saturated_frac_mid_section"] == pytest.approx(outer / (20 / vz), abs=1e-4)


# ---------------------------------------------------------------- seam 2: group plots draw flagged tumours hollow

def _facs():
    from kerncount.tumour.config import factors

    return factors({"tumour": {}})


def _small_table(n=4, flagged_every=3):
    import pandas as pd

    rows = []
    for g in ("WT", "TLR7 KO"):
        for line in ("cntrl", "HERV", "N4"):
            for c in ("-clodro", "+clodro"):
                for i in range(n):
                    rows.append({"tumour_id": f"{g} / {line} / {c} / n{i // 2 + 1}{'AB'[i % 2]}",
                                 "imaging_session": f"s {g}", "genotype": g, "tumour_line": line, "clodronate": c,
                                 "section": f"n{i // 2 + 1}", "hemisphere": "AB"[i % 2],
                                 "value": float(np.exp(np.random.default_rng(len(rows)).normal(0, 0.3))),
                                 "flagged": len(rows) % flagged_every == 0})
    return pd.DataFrame(rows)


def test_group_plot_draws_flagged_tumours_as_hollow_markers(tmp_path, monkeypatch):
    from matplotlib.colors import to_rgba

    from kerncount.tumour import plots

    figs = []
    real_save = plots._save
    monkeypatch.setattr(plots, "_save", lambda fig, out: (figs.append(fig), real_save(fig, out))[1])
    t = _small_table()
    plots.group_plot(t, _facs(), "value", "value", tmp_path / "v", hollow="flagged")
    filled = hollow = 0
    for ax in figs[0].axes:
        for col in ax.collections:
            if type(col).__name__ != "PathCollection":  # scatter dots only (not the mean / CI lines)
                continue
            fc = col.get_facecolors()
            n = len(col.get_offsets())
            if n == 0:
                continue
            if len(fc) == 0 or fc[0][3] == 0 or tuple(fc[0]) == to_rgba("white"):
                hollow += n
            else:
                filled += n
    assert hollow == int(t["flagged"].sum()) and filled == int((~t["flagged"]).sum())
    assert "hollow" in figs[0]._suptitle.get_text()


# ---------------------------------------------------------------- seam 3: flow table per group (all factors)

def test_flow_table_counts_files_exclusions_with_reasons_flagged_and_analysed_per_group():
    import pandas as pd

    from kerncount.tumour import qc

    t = _small_table(n=4, flagged_every=5)
    excluded = pd.DataFrame([
        {"tumour_id": "WT / HERV / +clodro / n3A", "genotype": "WT", "tumour_line": "HERV", "clodronate": "+clodro",
         "exclusion_reason": "section torn"},
        {"tumour_id": "WT / HERV / +clodro / n3B", "genotype": "WT", "tumour_line": "HERV", "clodronate": "+clodro",
         "exclusion_reason": "section torn"},
        {"tumour_id": "TLR7 KO / N4 / -clodro / n3A", "genotype": "TLR7 KO", "tumour_line": "N4",
         "clodronate": "-clodro", "exclusion_reason": "review (manual_1): not a tumour"},
    ])
    flow = qc.flow_table(t, excluded, _facs())
    assert list(flow.columns[:3]) == ["genotype", "tumour_line", "clodronate"]
    assert list(flow.columns[3:]) == ["files", "tumours", "excluded", "exclusion_reasons", "analysed", "flagged",
                                      "analysed_without_flagged"]
    assert len(flow) == 12 + 1  # every group + total
    g = flow[(flow.genotype == "WT") & (flow.tumour_line == "HERV") & (flow.clodronate == "+clodro")].iloc[0]
    sel = (t.genotype == "WT") & (t.tumour_line == "HERV") & (t.clodronate == "+clodro")
    assert (g.files, g.tumours, g.excluded, g.analysed) == (6, 6, 2, 4)
    assert g.exclusion_reasons == "section torn (2)"
    assert g.flagged == int(t.loc[sel, "flagged"].sum())
    assert g.analysed_without_flagged == 4 - g.flagged
    total = flow.iloc[-1]
    assert total.genotype == "all" and total.files == len(t) + 3 and total.excluded == 3
    assert total.flagged == int(t["flagged"].sum()) and total.analysed == len(t)


# ---------------------------------------------------------------- seam 4: results with QC (analyse on a readout table)

VOLUME_CONFIRMATORY = [{"name": "main_mass_volume", "column": "main_mass_volume_mm3", "label": "Main-mass volume (mm³)",
                        "role": "confirmatory"}]


def _settings(**tumour) -> dict:
    from kerncount import settings as settings_mod

    tumour.setdefault("readouts", VOLUME_CONFIRMATORY)
    return settings_mod._merge(settings_mod.DEFAULTS, {"experiment": {"type": "tumour_volume"}, "tumour": tumour})


def _vol_table(seed=0, n=6, sigma=0.05):
    """True HERV/cntrl ratio 2 everywhere; in WT -clodro HERV half the tumours are flagged and read 3x too large."""
    import pandas as pd

    rng = np.random.default_rng(seed)
    rows = []
    for g in ("WT", "TLR7 KO"):
        for line in ("cntrl", "HERV", "N4"):
            for c in ("-clodro", "+clodro"):
                for i in range(n):
                    v = (2.0 if line == "HERV" else 1.0) * float(np.exp(rng.normal(0, sigma)))
                    bad = (g, line, c) == ("WT", "HERV", "-clodro") and i % 2 == 0
                    rows.append({"tumour_id": f"{g} / {line} / {c} / n{i // 2 + 1}{'AB'[i % 2]}",
                                 "imaging_session": f"s {g}", "genotype": g, "tumour_line": line, "clodronate": c,
                                 "section": f"n{i // 2 + 1}", "hemisphere": "AB"[i % 2],
                                 "main_mass_volume_mm3": v * (3.0 if bad else 1.0), "flagged": bad,
                                 "truncated": bad, "qc_flags": "main mass touches the image edge" if bad else ""})
    return pd.DataFrame(rows)


def test_sensitivity_without_flagged_tumours_is_reported_beside_the_primary_with_flow_and_hemisphere_plot(tmp_path):
    import pandas as pd

    from kerncount.experiment_types import RunOptions, get

    t = _vol_table()
    pd.DataFrame([{"tumour_id": "WT / N4 / +clodro / n9A", "genotype": "WT", "tumour_line": "N4",
                   "clodronate": "+clodro", "exclusion_reason": "out of focus"}]).to_csv(
        tmp_path / "excluded_tumours.csv", index=False)
    get("tumour_volume").analyse(tmp_path, t, _settings(), RunOptions())
    xl = pd.ExcelFile(tmp_path / "results.xlsx")
    assert {"QC flow", "Without flagged"} <= set(xl.sheet_names)
    side = pd.read_excel(xl, "Without flagged")
    row = side[(side.family == "F1") & side.comparison.str.contains("HERV vs cntrl")
               & side.comparison.str.contains("WT") & side.comparison.str.contains("-clodro")].iloc[0]
    assert row["estimate"] > 2.5  # primary: inflated by the flagged tumours (kept in the primary analysis)
    assert row["without_flagged_estimate"] == pytest.approx(2.0, rel=0.08)
    assert row["N"] - row["without_flagged_N"] == 3
    flow = pd.read_csv(tmp_path / "flow.csv")
    total = flow.iloc[-1]
    assert (total.files, total.excluded, total.flagged, total.analysed_without_flagged) == (len(t) + 1, 1, 3,
                                                                                            len(t) - 3)
    st = pd.read_csv(tmp_path / "stats_long.csv")
    wf = st[st.family == "F1 (without flagged)"]
    assert len(wf) and (wf["significance"].fillna("") == "").all()  # sensitivity rows never carry stars
    report = (tmp_path / "stats_report.md").read_text()
    assert "## QC and problem tumours" in report and "out of focus" in report and "without flagged" in report
    assert (tmp_path / "plots" / "main_mass_volume_hemispheres_AB.png").exists()


def test_without_flagged_tumours_nothing_is_dropped_and_the_section_says_so(tmp_path):
    from kerncount.experiment_types import RunOptions, get

    t = _vol_table().assign(flagged=False)
    get("tumour_volume").analyse(tmp_path, t, _settings(), RunOptions())
    assert "no tumour is flagged" in (tmp_path / "stats_report.md").read_text()


# ---------------------------------------------------------------- seam 5: run -> blinded review sheets + corrections

KO_SESSION = "ready 5d KO7 030826"
WT_SESSION = "ready 5d WT 120926"
SVOX = (2.36, 1.25, 1.25)
SMALL = dict(voxel_um=SVOX, shape=(30, 160, 160), plateau=90)
GROUP_WORDS = ["HERV", "cntrl", "clodro", "TLR7", "KO7", "N4 ", "030826", "120926", "ready"]


def _run_experiment(tmp_path, monkeypatch):
    from tumour_synth import Satellite, SynthTumour, write_ome

    from kerncount import cli

    edge = 160 * 1.25
    tumours = {
        f"{KO_SESSION}/030826 HERV -clodro n1 A.ome.tiff": SynthTumour(
            semi_axes_um=(20, 45, 40), centre_um=(35, 70, 90), seed=1,
            satellites=[Satellite((35, 160, 160), 18.0)], **SMALL),
        f"{KO_SESSION}/030826 HERV -clodro n1 B.ome.tiff": SynthTumour(semi_axes_um=(15, 35, 30), seed=2, **SMALL),
        f"{KO_SESSION}/030826 cntrl +clodro n1 A.ome.tiff": SynthTumour(  # cut by the right image edge
            semi_axes_um=(20, 60, 50), centre_um=(35, 100, edge - 10), seed=3, **SMALL),
        f"{WT_SESSION}/120926 N4 -clodro n2 B.ome.tiff": SynthTumour(semi_axes_um=(18, 45, 45), seed=4, **SMALL),
    }
    raw = tmp_path / "data" / "tumour_exp"
    for rel, t in tumours.items():
        write_ome(raw / rel, t.stack(), t.voxel_um)
    monkeypatch.chdir(tmp_path)
    cli.main(["init", "tumour_exp", "--type", "tumour_volume"])
    cli.main(["run", "tumour_exp"])
    return tumours, tmp_path / "analysis" / "tumour_exp"


def _read(path):
    import pandas as pd

    return pd.read_csv(path, dtype=str, keep_default_na=False)


def test_run_flags_a_synthetic_tumour_cut_by_the_image_edge_and_writes_blinded_review_sheets(tmp_path, monkeypatch):
    import pandas as pd

    tumours, prj = _run_experiment(tmp_path, monkeypatch)
    table = pd.read_csv(prj / "tumours.csv")
    flagged = dict(zip(table["file"], table["flagged"]))
    assert flagged[f"{KO_SESSION}/030826 cntrl +clodro n1 A.ome.tiff"] == "yes"
    assert list(flagged.values()).count("yes") == 1
    assert {"edge_outline_frac", "mid_section_fits", "saturated_frac_mid_section"} <= set(table.columns)

    key = _read(prj / "review_key.csv")
    assert set(key["tumour_id"]) == set(table["tumour_id"]) and key["code"].is_unique
    review = prj / "review"
    files = sorted(p for p in review.rglob("*") if p.is_file())
    assert {p.name for p in files} >= {f"{c}.jpg" for c in key["code"]} | {"corrections.csv", "overview_1.jpg"}
    for p in files:  # blind: no factor level, session or file name in any name or content of the review folder
        blob = p.name.encode() + p.read_bytes()
        assert not [w for w in GROUP_WORDS if w.encode() in blob], p.name
    assert not (review / "review_key.csv").exists()
    corr = _read(review / "corrections.csv")
    assert list(corr["code"]) == sorted(key["code"])  # random order, not grouped
    assert {"manual_1", "manual_2", "manual_3"} <= set(corr.columns)
    size = sum(p.stat().st_size for p in files)
    assert size < 2_000_000


def test_corrections_re_pick_the_main_mass_or_exclude_and_rebuild_without_re_segmenting(tmp_path, monkeypatch):
    import pandas as pd

    from kerncount import cli
    from kerncount.tumour import segment

    tumours, prj = _run_experiment(tmp_path, monkeypatch)
    before = pd.read_csv(prj / "tumours.csv").set_index("file")
    key = _read(prj / "review_key.csv").merge(pd.read_csv(prj / "tumours.csv")[["tumour_id", "file"]])
    code = dict(zip(key["file"], key["code"]))
    sat_file = f"{KO_SESSION}/030826 HERV -clodro n1 A.ome.tiff"
    bad_file = f"{WT_SESSION}/120926 N4 -clodro n2 B.ome.tiff"
    corr = _read(prj / "review" / "corrections.csv")
    corr.loc[corr.code == code[sat_file], "manual_1"] = "2"
    corr.loc[corr.code == code[bad_file], "manual_1"] = "x: air bubble over the tumour"
    corr.to_csv(prj / "review" / "corrections.csv", index=False)

    def no_resegment(*a, **k):
        raise AssertionError("re-segmented")
    monkeypatch.setattr(segment, "measure_file", no_resegment)
    cli.main(["run", "tumour_exp"])

    after = pd.read_csv(prj / "tumours.csv").set_index("file")
    assert bad_file not in after.index and len(after) == 3
    sat_vol = tumours[sat_file].satellites[0].volume_um3
    assert after.loc[sat_file, "main_mass_volume_um3"] == pytest.approx(sat_vol, rel=0.1)
    assert before.loc[sat_file, "main_mass_volume_um3"] > 3 * sat_vol
    other = f"{KO_SESSION}/030826 HERV -clodro n1 B.ome.tiff"
    assert after.loc[other, "main_mass_volume_um3"] == before.loc[other, "main_mass_volume_um3"]
    ex = _read(prj / "excluded_tumours.csv")
    assert "air bubble over the tumour" in ex["exclusion_reason"].iloc[0] and "manual_1" in ex["exclusion_reason"].iloc[0]
    used = _read(prj / "corrections_used.csv").set_index("code")
    assert used.loc[code[sat_file], "source"].startswith("manual_1") and used.loc[code[other], "source"] == "auto"
    assert "g2" in (prj / "stats_report.md").read_text()
    assert pd.read_csv(prj / "flow.csv").iloc[-1]["excluded"] == 1

    corr.loc[corr.code == code[sat_file], "manual_2"] = "99"  # no such object: the run stops and says why
    corr.to_csv(prj / "review" / "corrections.csv", index=False)
    with pytest.raises(SystemExit, match=code[sat_file]):
        cli.main(["run", "tumour_exp"])


def test_the_threshold_sweep_follows_review_corrections(tmp_path, monkeypatch):
    """Review fix: a re-picked main mass is followed at every sweep scale (the object overlapping the corrected x1
    main mass most), and a tumour excluded at review stays out of the sweep."""
    import pandas as pd

    from kerncount import cli, settings as settings_mod
    from kerncount.tumour import sweep

    tumours, prj = _run_experiment(tmp_path, monkeypatch)
    p = prj / "settings.toml"
    p.write_text(p.read_text().replace("threshold_sweep = false", "threshold_sweep = true"))
    key = _read(prj / "review_key.csv").merge(pd.read_csv(prj / "tumours.csv")[["tumour_id", "file"]])
    code = dict(zip(key["file"], key["code"]))
    sat_file = f"{KO_SESSION}/030826 HERV -clodro n1 A.ome.tiff"
    bad_file = f"{WT_SESSION}/120926 N4 -clodro n2 B.ome.tiff"
    corr = _read(prj / "review" / "corrections.csv")
    corr.loc[corr.code == code[sat_file], "manual_1"] = "2"
    corr.loc[corr.code == code[bad_file], "manual_1"] = "x: air bubble over the tumour"
    corr.to_csv(prj / "review" / "corrections.csv", index=False)
    cli.main(["run", "tumour_exp"])

    table = pd.read_csv(prj / "tumours.csv")
    tabs = sweep.tables(prj, table, settings_mod.load(prj))  # cached by the run: read back, nothing recomputed
    sat_vol = tumours[sat_file].satellites[0].volume_um3
    for scale, t in tabs.items():
        t = t.set_index("file")
        assert bad_file not in t.index, scale
        assert t.loc[sat_file, "main_mass_volume_um3"] == pytest.approx(sat_vol, rel=0.35), scale


# ---------------------------------------------------------------- real data (slow): QC + review image of real stacks

REAL = Path(__file__).resolve().parents[2] / "data" / "tumour_TLR7" / "5d K07"  # renamed from KO_SESSION 2026-11-20
REAL_STEMS = ["030826 cntrl +clodro n1 B", "030826 HERV -clodro n2 A"]


@pytest.mark.slow
@pytest.mark.skipif(not (REAL / f"{REAL_STEMS[0]}.ome.tiff").exists(), reason="real tumour data not present")
@pytest.mark.parametrize("stem", REAL_STEMS)
def test_real_stack_qc_and_review_image_are_small_and_blind(tmp_path, stem):
    import time

    from kerncount.tumour import qc, review
    from kerncount.tumour.config import tumour_settings
    from kerncount.tumour.segment import measure_file, read_stack, voxel_size_um

    path = REAL / f"{stem}.ome.tiff"
    settings = {"tumour": {}}
    seg = measure_file(path, tumour_settings(settings))
    vox = voxel_size_um(path)
    stack = read_stack(path)
    t0 = time.time()
    q, centre = qc.assess(seg.labels, stack, vox, settings)
    out = review.render(stack, seg.labels, centre, "TEST1", tmp_path / "TEST1.jpg", (vox[1], vox[2]))
    took = time.time() - t0
    print(stem, q, f"{took:.1f} s", out.stat().st_size)
    assert set(q) == set(qc.QC_COLUMNS)
    assert q["edge_contact"] in ("yes", "no") and 0 <= q["saturated_frac_mid_section"] <= 1
    assert out.stat().st_size < 400_000  # review images stay small (48 tumours: a few MB)
    assert not [w for w in GROUP_WORDS if w.encode() in out.read_bytes()]
    assert took < 60


def test_hemisphere_plot_keeps_title_and_axis_labels_inside_the_figure(monkeypatch):
    """Real run (ticket 11): the long confirmatory label pushed the title and y label off the figure."""
    from pathlib import Path

    import pandas as pd

    from kerncount.tumour import plots, qc

    figs = []
    monkeypatch.setattr(plots, "_save", lambda fig, out: figs.append(fig) or out)
    pairs = pd.DataFrame({"group": ["g1", "g1", "g2"], "section": ["n1", "n2", "n1"], "A": [0.2, 0.3, 0.4],
                          "B": [0.25, 0.2, 0.5], "flagged": [False, True, False]})
    qc.hemisphere_plot(pairs, "Cross-sectional area, mid-section (mm²)", Path("x"))
    fig = figs[0]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    ax = fig.axes[0]
    for text in (ax.title, ax.yaxis.label, ax.xaxis.label):
        e = text.get_window_extent(r)
        assert e.x0 >= -1 and e.x1 <= W + 1 and e.y0 >= -1 and e.y1 <= H + 1, (text.get_text(), e, W, H)
