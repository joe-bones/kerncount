"""Readout ladder (layers.py, intensity.field_layers), the metric N column, and the generation input state."""
import json

import numpy as np
import pandas as pd
import pytest

from kerncount import classify, design, generations, intensity, layers, metrics, models, pipeline
from kerncount import settings as st_mod


# ---------------------------------------------------------------- field layers on synthetic images

def _field(rng, signal=400.0, bg=200.0, flat=None, spots=30, shape=(256, 256)):
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    sig = np.zeros(shape)
    for _ in range(spots):
        cy, cx = rng.uniform(15, shape[0] - 15, 2)
        sig[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] = signal
    F = np.ones(shape) if flat is None else flat
    return sig * F + bg + rng.normal(0, 10, shape), sig


def test_field_layers_raw_corrected_stained_and_saturation():
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[:256, :256]
    F = 1.2 - 0.4 * xx / 255  # vignetting, mean ~1
    F = F / F.mean()
    img, sig = _field(rng, flat=F)
    out = intensity.field_layers(img, 65535, F, np.full(F.shape, 200.0))  # BaSiC: observed = signal x F + D
    assert out["field_raw"] == pytest.approx(img.mean(), rel=1e-6)  # raw = the image as acquired
    # corrected: (img - D) / F minus the background mode -> mean of the true signal
    assert out["field_corrected"] == pytest.approx(sig.mean(), rel=0.15)
    assert out["stained_area_mean"] == pytest.approx(400, rel=0.1)
    assert out["stained_frac"] == pytest.approx((sig > 0).mean(), rel=0.15)
    assert out["field_total"] == pytest.approx(out["stained_area_mean"] * out["stained_frac"] * img.size, rel=1e-6)
    # no stain: total 0, no stained-area mean
    blank = 200 + rng.normal(0, 10, (256, 256))
    o2 = intensity.field_layers(blank, 65535)
    assert o2["field_total"] == 0.0 and np.isnan(o2["stained_area_mean"]) and abs(o2["field_corrected"]) < 2
    # saturated field: every layer left empty
    sat = img.copy()
    sat[:40] = 65535
    o3 = intensity.field_layers(sat, 65535)
    assert np.isnan(o3["field_raw"]) and np.isnan(o3["field_corrected"]) and o3["saturated_frac"] > 0.1


def test_background_mode_ignores_bright_objects():
    rng = np.random.default_rng(1)
    v = np.concatenate([rng.normal(100, 5, 60000), rng.normal(600, 50, 40000)])
    m, sd = intensity.background_mode(v)
    assert m == pytest.approx(100, abs=2) and sd == pytest.approx(5, rel=0.25)


# ---------------------------------------------------------------- registry: labels, kinds, N column

def test_registry_labels_and_n_column():
    cols = ["count", "intensity_Ki67_field_raw", "intensity_Ki67_field_corrected", "intensity_Ki67_nucleus_mean_raw",
            "intensity_Ki67_nucleus_mean_bgsub", "intensity_Ki67_stained_area_mean", "intensity_Ki67_field_total",
            "readout_frac_positive_Ki67", "readout_frac_positive_p90_Ki67", "n_measured_Ki67", "intensity_C2_nucleus"]
    reg = metrics.registry(cols)
    assert "n_measured_Ki67" not in reg
    for c in cols[1:7]:
        info = classify.READOUTS.get(c)
        assert reg[c].label == info["label"] and reg[c].unit == info["unit"] and reg[c].kind == "positive"
    assert reg["intensity_Ki67_field_raw"].label == "Ki67 field MFI, raw"
    assert reg["intensity_C2_nucleus"].label == "C2 nucleus mean intensity"  # legacy columns keep the fallback
    assert reg["readout_frac_positive_Ki67"].n_column == "n_measured_Ki67"
    assert reg["readout_frac_positive_p90_Ki67"].n_column == "n_measured_Ki67"
    assert metrics.n_columns(reg) == ["n_measured_Ki67"]
    assert metrics.normalise({"readout_frac_positive_Ki67": "x"})["readout_frac_positive_Ki67"].n_column == "n_measured_Ki67"


def test_fraction_model_uses_n_measured_not_count():
    met = metrics.describe_column("readout_frac_positive_Ki67")
    f = pd.DataFrame(dict(count=[100, 100], n_measured_Ki67=[50, 20], value=[0.5, 0.5]))
    assert models.nuclei(f, met).tolist() == [50, 20]
    assert models.nuclei(f.drop(columns="n_measured_Ki67"), met).tolist() == [100, 100]
    rows = []
    for n in ("n1", "n2"):
        for t, p in (("LyoVec", 0.5), ("X", 0.3)):
            for i in range(1, 13):
                rows.append(dict(n=n, treatment=t, image=i, well=1 + (i > 6), position=(i - 1) % 6 + 1, timepoint="",
                                 count=200, n_measured_Ki67=40, value=p))
    w, notes, has_n = models.aggregate_wells(pd.DataFrame(rows), met, True)
    assert has_n and (w["nuclei"] == 6 * 40).all()
    assert any("n_measured_Ki67" in x for x in notes)


# ---------------------------------------------------------------- concordance

def _eff(calls, n=None, level="all_n"):
    rows = []
    for (k, lab, _), (est, call) in zip(layers.LADDER, calls):
        d = "higher" if call.startswith("higher") else "lower" if call.startswith("lower") else "no change"
        rows.append(dict(layer=lab, layer_key=k, layer_index=layers.LADDER.index((k, lab, _)) + 1, extra=False,
                         metric=layers.column("Ki67", k), reference="LyoVec", cell_line=None, level=level, n=n,
                         treatment="miR-1", estimate=est, ci_low=np.nan, ci_high=np.nan, direction=d, call=call,
                         significance="", effect_type="", interval="", p_adj=np.nan, N=2))
    return pd.DataFrame(rows)


def test_concordance_names_the_step_that_changes_the_conclusion():
    c = layers.concordance(_eff([(95, "lower (ns)"), (80, "lower *"), (85, "lower *"), (70, "lower *"),
                                 (-12, "lower *")]))
    r = c.iloc[0]
    assert r["direction_agrees"] and not r["conclusion_agrees"] and r["changes_at"] == "1->2"
    assert "flat-field and background correction" in r["summary"] and "Headline (5 % positive): lower *" in r["summary"]
    c = layers.concordance(_eff([(90, "lower"), (80, "lower"), (88, "lower"), (75, "lower"), (-10, "lower")]),
                           "significance needs ≥ 2 n")
    assert c.iloc[0]["conclusion_agrees"] and "needs ≥ 2 n" in c.iloc[0]["summary"]
    c = layers.concordance(_eff([(90, "lower *"), (85, "lower *"), (95, "lower *"), (90, "lower *"),
                                 (4, "higher (ns)")]))
    assert not c.iloc[0]["direction_agrees"] and c.iloc[0]["changes_at"] == "4->5"
    assert "positivity cut-off" in c.iloc[0]["summary"]


def test_layer_channels_need_two_layers():
    assert layers.channels(["count", "readout_frac_positive_Ki67", "intensity_Ki67_nucleus_mean_bgsub"]) == ["Ki67"]
    assert layers.channels(["count", "readout_frac_positive_Ki67"]) == []
    assert layers.channels(["readout_frac_positive_p90_Ki67", "intensity_Ki67_field_raw"]) == []


# ---------------------------------------------------------------- generation: layers, metadata, inputs

def _settings(**over):
    s = st_mod.load(__import__("pathlib").Path("/nonexistent"))
    for sec, kv in over.items():
        s[sec].update(kv)
    return s


def _ladder_table(n_list=("n1", "n2"), shift=1.0):
    rng = np.random.default_rng(4)
    rows = []
    for n in n_list:
        for t, eff in (("Control", 1.0), ("LyoVec", 1.0), ("miR-1", 0.6)):
            for i in range(1, 13):
                sig = 500 * eff * shift * np.exp(rng.normal(0, 0.1))
                pos = 0.6 * eff + rng.normal(0, 0.03)
                rows.append(dict(n=n, treatment=t, image=i, count=int(rng.integers(80, 120)), source="auto",
                                 excluded=False, qc_level="high", n_measured_Ki67=60,
                                 source_file=f"/raw/{n}/{t}/{2 * i - 1}.tif; /raw/{n}/{t}/{2 * i}.tif",
                                 source_file_Ki67=f"/raw/{n}/{t}/{2 * i - 1}.tif",
                                 intensity_Ki67_field_raw=300 + 0.2 * sig, intensity_Ki67_field_corrected=0.2 * sig,
                                 intensity_Ki67_nucleus_mean_raw=300 + sig, intensity_Ki67_nucleus_mean_bgsub=sig,
                                 intensity_Ki67_stained_area_mean=0.9 * sig, intensity_Ki67_field_total=900 * sig,
                                 readout_frac_positive_Ki67=float(np.clip(pos, 0, 1))))
    return pd.DataFrame(rows)


@pytest.mark.slow
def test_generation_writes_layers_metadata_and_source_files(tmp_path):
    tab = _ladder_table()
    s = _settings(experiment={"cell_line": "LN229", "magnification": 20})
    s["analysis"].update(reference="LyoVec", normalise_to=["Control"])
    img = tab.drop(columns=["count", "source", "excluded"])
    reg = metrics.registry(img.columns, s)
    reg = {k: m for k, m in reg.items() if k != "count"}
    gen = dict(name="g1", table=tab[["n", "treatment", "image", "count", "source", "excluded", "qc_level"]],
               source="auto", changes=0, fingerprint="x")
    out = generations.build(tmp_path, gen, img, reg, s)
    used = pd.read_csv(out / "counts_used.csv", keep_default_na=False)
    assert list(used.columns)[:11] == ["n", "treatment", "image", "count", "source", "excluded", "qc_level", "well",
                                       "position", "timepoint", "cell_line"]
    assert (used["cell_line"] == "LN229").all() and (used["meta_magnification"].astype(str) == "20").all()
    assert {"source_file", "source_file_Ki67"} <= set(used.columns)
    for f in ("layers_Ki67_effects.csv", "layers_Ki67_concordance.csv", "layers_Ki67_absolute.csv",
              "layers_Ki67_vs_LyoVec.png", "layers_Ki67_vs_Control.png", "layers_Ki67_absolute.png"):
        assert (out / "all_n" / f).exists(), f
        assert (out / "n1" / f).exists(), f
    eff = pd.read_csv(out / "all_n" / "layers_Ki67_effects.csv")
    m = eff[(eff["treatment"] == "miR-1") & (eff["reference"] == "LyoVec")].set_index("layer_key")
    assert m.loc["nucleus_mean_bgsub", "estimate"] == pytest.approx(60, abs=8)  # % of reference
    assert 60 < m.loc["field_raw", "estimate"] < 100  # raw layers compress toward 100%
    assert m.loc["frac_positive", "effect_type"] == "pp_difference"
    assert m.loc["frac_positive", "estimate"] == pytest.approx(-24, abs=4)
    conc = pd.read_csv(out / "all_n" / "layers_Ki67_concordance.csv")
    assert {"summary", "changes_at", "conclusion_agrees"} <= set(conc.columns)
    man = json.loads((out / "manifest.json").read_text())
    kinds = {f["path"]: f["kind"] for f in man["files"]}
    assert kinds["all_n/layers_Ki67_vs_LyoVec.png"] == "layers" and kinds["n1/layers_Ki67_concordance.csv"] == "layers"
    assert [x["key"] for x in man["layers"]["Ki67"]["layers"]] == [k for k, _, _ in layers.LADDER]
    assert man["layers"]["Ki67"]["references"] == ["LyoVec", "Control"]


def test_input_state_describes_cutoff_and_metadata_changes(tmp_path):
    class P:
        path = tmp_path

    tab = _ladder_table(("n1",))
    keys = tab[["n", "treatment", "image"]]
    s = _settings(experiment={"cell_line": "LN229"})
    dt = design.assign(keys, s, None)
    (tmp_path / "intensity").mkdir()
    pd.DataFrame([dict(n="n1", channel="Ki67", measure="mean", rule="lod", value=120.0)]).to_csv(
        tmp_path / "intensity" / "cutoffs.csv", index=False)
    a = pipeline.input_state(P, tab, dt)
    assert pipeline.describe_input_changes(a, pipeline.input_state(P, tab, dt)) == []
    pd.DataFrame([dict(n="n1", channel="Ki67", measure="mean", rule="manual", value=150.0)]).to_csv(
        tmp_path / "intensity" / "cutoffs.csv", index=False)
    b = pipeline.input_state(P, tab.assign(readout_frac_positive_Ki67=tab["readout_frac_positive_Ki67"] / 2), dt)
    ch = pipeline.describe_input_changes(a, b)
    assert any("cut-off (mean) in n1: lod 120 -> manual 150" in x for x in ch)
    dt2 = design.assign(keys, _settings(experiment={"cell_line": "U251"}), None)
    ch = pipeline.describe_input_changes(a, pipeline.input_state(P, tab, dt2))
    assert ch and "cell line n1: LN229 -> U251" in ch[0]
    # a DAPI-only table has no intensity state: unaffected by intensity rules
    assert pipeline.input_state(P, tab[["n", "treatment", "image"]], dt)["intensity_fingerprint"] == ""
