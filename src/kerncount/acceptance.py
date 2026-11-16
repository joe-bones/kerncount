"""Per-n acceptance criteria ([acceptance] in settings.toml), evaluated on the counts of one generation.

Criteria (each may be switched off with ""):
    vehicle_within             the vehicle / control ratio of mean nuclei per field in that n lies within 1 +- value
                               (the pair comes from [analysis] vehicle_check; "auto" = the second control vs the first)
    positive_control_max_ratio positive control / primary reference must be <= value (if value < 1) or >= value
                               (if value > 1): the positive control must show its expected effect
    min_images_per_condition   every condition in that n has at least this many usable images
    max_low_qc_fraction        at most this fraction of that n's images has QC level low

enforce = false: report only. enforce = true: n's failing any applicable criterion are left out of every all-n
statistic (they keep their per-n outputs).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

COLUMNS = ["n", "criterion", "value", "limit", "result", "detail"]


def _off(v) -> bool:
    return v is None or (isinstance(v, str) and v.strip() == "")


def evaluate(used: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """used: counts_used rows (n, treatment, image, count, excluded, qc_level). One row per n x criterion plus an
    'overall' row per n. result is pass | fail | not applicable."""
    from .models import config

    acc = settings.get("acceptance", {})
    cfg = config(settings, used["treatment"].unique())
    rows = []
    live = used[~used["excluded"].astype(bool)]
    for n in sorted(used["n"].astype(str).unique()):
        u = used[used["n"].astype(str) == n]
        lv = live[live["n"].astype(str) == n]
        means = lv.groupby("treatment")["count"].mean()
        res = []

        v = acc.get("vehicle_within", "")
        if not _off(v):
            pair = cfg.vehicle_pairs[0] if cfg.vehicle_pairs else None
            if pair and pair[0] in means and pair[1] in means and means[pair[1]] > 0:
                ratio = float(means[pair[0]] / means[pair[1]])
                ok = abs(ratio - 1) <= float(v)
                res.append(("vehicle_within", ratio, f"{1 - float(v):.3g}-{1 + float(v):.3g}", "pass" if ok else "fail",
                            f"{pair[0]} / {pair[1]} mean nuclei per field"))
            else:
                res.append(("vehicle_within", np.nan, f"1 ± {float(v):g}", "not applicable",
                            "vehicle check pair not present in this n" if pair else "no vehicle check pair configured"))

        v = acc.get("positive_control_max_ratio", "")
        if not _off(v):
            pc, ref = cfg.positive_control, cfg.primary_ref
            if ref not in means:
                ref = next((c for c in cfg.controls if c in means), None)
            if pc and ref and pc in means and means[ref] > 0:
                ratio = float(means[pc] / means[ref])
                lim = float(v)
                ok = ratio <= lim if lim < 1 else ratio >= lim
                res.append(("positive_control_max_ratio", ratio, f"{'<=' if lim < 1 else '>='} {lim:g}",
                            "pass" if ok else "fail", f"{pc} / {ref} mean nuclei per field"))
            else:
                res.append(("positive_control_max_ratio", np.nan, str(v), "not applicable",
                            "no positive control configured" if not pc else "positive control or reference missing in this n"))

        v = acc.get("min_images_per_condition", "")
        if not _off(v):
            per = lv.groupby("treatment").size().reindex(sorted(u["treatment"].unique()), fill_value=0)
            worst = per.idxmin() if len(per) else ""
            res.append(("min_images_per_condition", float(per.min()) if len(per) else np.nan, f">= {int(v)}",
                        "pass" if len(per) and per.min() >= int(v) else "fail", f"fewest: {worst}"))

        v = acc.get("max_low_qc_fraction", "")
        if not _off(v):
            frac = float((u["qc_level"].astype(str) == "low").mean()) if len(u) else np.nan
            res.append(("max_low_qc_fraction", frac, f"<= {float(v):g}", "pass" if frac <= float(v) else "fail",
                        f"{int((u['qc_level'].astype(str) == 'low').sum())} of {len(u)} images QC low"))

        for crit, value, limit, result, detail in res:
            rows.append(dict(n=n, criterion=crit, value=value, limit=limit, result=result, detail=detail))
        fails = [x[0] for x in res if x[3] == "fail"]
        rows.append(dict(n=n, criterion="overall", value=np.nan, limit="", result="fail" if fails else "pass",
                         detail=("failed: " + ", ".join(fails)) if fails else "all applicable criteria passed"))
    return pd.DataFrame(rows, columns=COLUMNS)


def failing(table: pd.DataFrame) -> list[str]:
    o = table[table["criterion"] == "overall"]
    return o.loc[o["result"] == "fail", "n"].tolist()


def summary(table: pd.DataFrame, settings: dict) -> dict:
    enforce = bool(settings.get("acceptance", {}).get("enforce", False))
    fail = failing(table)
    return dict(enforce=enforce, n_checked=int(table["n"].nunique()), failing=fail,
                excluded=fail if enforce else [],
                criteria={c: settings.get("acceptance", {}).get(c) for c in
                          ("vehicle_within", "positive_control_max_ratio", "min_images_per_condition", "max_low_qc_fraction")})
