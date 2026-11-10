"""Write one generation folder (g1, g2, ...) with outputs per treatment, per n and for all n."""
from __future__ import annotations

import json
import re
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from . import analysis

README = """Generation {name}
==============
Source: {source_text}
Built:  {built}

Layout
  counts_used.csv            the value used for every image in this generation (and where it came from)
  results.xlsx               everything below in one readable workbook
  <n>/<treatment>/           per-field values and plot for one treatment in one n
  <n>/                       one experiment: summary, Prism table, stats, plot per metric
  all_n/                     all experiments combined: summary, Prism tables (raw and % of control),
                             stats, plots (dot colour = n, see n_colours.csv)

Prism import: *_prism*.csv are Column tables (one column per treatment). In Prism use
File > Import, or open the CSV and paste into a Column table.

Statistics: see the header of each *_stats.csv. Within one n, fields are technical replicates
(descriptive only). Across n, each n is one biological replicate.
"""


def _slug(s: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", s).strip("_")


def _autosize(writer: pd.ExcelWriter):
    for ws in writer.book.worksheets:
        ws.freeze_panes = "A2"
        for col in ws.columns:
            width = max((len(str(c.value)) if c.value is not None else 0) for c in col)
            ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)


def _fields_plot(sub: pd.DataFrame, out: Path, title: str):
    fig, ax = plt.subplots(figsize=(5, 3))
    ok = sub[~sub["excluded"]]
    ax.bar(ok["image"].astype(str), ok["count"], color="#9ecae1", edgecolor="black")
    bad = sub[sub["excluded"]]
    if len(bad):
        ax.bar(bad["image"].astype(str), [0] * len(bad), color="none")
        for x in bad["image"].astype(str):
            ax.text(x, 0, "x", ha="center", va="bottom", color="red")
    ax.set_xlabel("image")
    ax.set_ylabel("nuclei")
    ax.set_title(title, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def build(project: Path, gen: dict, image_table: pd.DataFrame, intensity_metrics: dict[str, str], settings: dict) -> Path:
    out = project / "generations" / gen["name"]
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    a = settings["analysis"]

    t = gen["table"].copy()
    t = t.merge(image_table, on=["n", "treatment", "image"], how="left", suffixes=("", "_img"))
    if a.get("exclude_low_qc"):
        drop = (t["qc_level"] == "low") & (t["source"] == "auto")
        t.loc[drop, "excluded"] = True
        t.loc[drop, "source"] = "auto (excluded: QC low)"
    used_cols = ["n", "treatment", "image", "count", "source", "excluded", "qc_level"]
    t[used_cols].to_csv(out / "counts_used.csv", index=False)

    metrics = {"count": "Nuclei per field"} | intensity_metrics
    live = t[~t["excluded"]]
    per_n_summ, per_n_stats, all_summ, all_stats = {}, {}, {}, {}

    for n in sorted(t["n"].unique()):
        nd = out / _slug(n)
        nd.mkdir()
        for tr in analysis.treatment_order(t.loc[t["n"] == n, "treatment"], a["treatment_order"]):
            td = nd / _slug(tr)
            td.mkdir()
            sub = t[(t["n"] == n) & (t["treatment"] == tr)].sort_values("image")
            keep = ["image", "count", "source", "excluded", "qc_level"] + [m for m in intensity_metrics if m in sub]
            sub[keep].to_csv(td / "fields.csv", index=False)
            _fields_plot(sub, td / "count_fields.png", f"{n} / {tr}")
        for m, label in metrics.items():
            if m == "count":
                vals = live.assign(value=live["count"])
            elif m in live:
                vals = live.assign(value=live[m])
            else:
                continue
            if vals["value"].notna().sum() == 0:
                continue
            st = analysis.per_n(vals, m, label, n, nd, settings)
            per_n_stats.setdefault(m, []).append(st)
            per_n_summ.setdefault(m, []).append(pd.read_csv(nd / f"{m}_summary.csv"))

    ad = out / "all_n"
    ad.mkdir()
    for m, label in metrics.items():
        col = "count" if m == "count" else m
        if col not in live or live[col].notna().sum() == 0:
            continue
        summ, st = analysis.all_n(live.assign(value=live[col]), m, label, ad, settings)
        all_summ[m], all_stats[m] = summ, st

    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl") as xw:
        about = pd.DataFrame({"": [f"Generation {gen['name']}", f"Source: {gen['source']}",
                                   f"Manual changes in this step: {gen['changes']}",
                                   f"Built {datetime.now():%Y-%m-%d %H:%M}",
                                   f"Reference group: {a['reference']}",
                                   "Within-n stats treat fields as technical replicates (descriptive).",
                                   "All-n stats: each n = one biological replicate; ratio paired t-tests vs reference, Holm-adjusted."]})
        about.to_excel(xw, sheet_name="About", index=False, header=False)
        t[used_cols].to_excel(xw, sheet_name="Counts used", index=False)
        if "count" in all_summ:
            all_summ["count"].to_excel(xw, sheet_name="Count - all n", index=False)
            all_stats["count"].to_excel(xw, sheet_name="Count - stats all n", index=False)
        if "count" in per_n_summ:
            pd.concat(per_n_summ["count"]).to_excel(xw, sheet_name="Count - per n", index=False)
            pd.concat(per_n_stats["count"]).to_excel(xw, sheet_name="Count - stats per n", index=False)
        inten = [m for m in intensity_metrics if m in all_summ]
        if inten:
            pd.concat([all_summ[m].assign(metric=m) for m in inten]).to_excel(xw, sheet_name="Intensity - all n", index=False)
            pd.concat([all_stats[m].assign(metric=m) for m in inten]).to_excel(xw, sheet_name="Intensity - stats all n", index=False)
            pd.concat([pd.concat(per_n_summ[m]).assign(metric=m) for m in inten if m in per_n_summ]).to_excel(
                xw, sheet_name="Intensity - per n", index=False)
        _autosize(xw)

    src = "fully automatic (auto_count)" if gen["source"] == "auto" else \
        f"previous generation + manual values from column {gen['source']} of counts.csv ({gen['changes']} cells)"
    (out / "README.txt").write_text(README.format(name=gen["name"], source_text=src, built=f"{datetime.now():%Y-%m-%d %H:%M}"),
                                    encoding="utf-8")
    (out / "generation.json").write_text(json.dumps(dict(name=gen["name"], source=gen["source"], fingerprint=gen["fingerprint"],
                                                         changes=gen["changes"], built=datetime.now().isoformat(timespec="seconds")),
                                                    indent=2), encoding="utf-8")
    return out
