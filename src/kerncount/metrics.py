"""Metric registry: which per-field columns are analysed, and how (a shared contract).

Every analysed column becomes a `Metric` with a `kind` that decides the statistical scale:

    count      nuclei per field                      -> log(value + COUNT_OFFSET)
    positive   intensities, densities, areas, ...    -> log(value)       (results are ratios)
    fraction   a proportion 0..1 per field           -> logit            (results are odds ratios)

Recognised columns:
    count                 nuclei per field (always first)
    intensity_*           mean intensity per cell (positive)
    nuclei_per_mm2        density (positive; only present when the pixel size is known)
    readout_*             per-image readouts; label / kind / unit come from `kerncount.classify.READOUTS` when it
                          exists, otherwise the kind is inferred from the name (readout_frac_* = fraction, else positive)

`registry(columns, settings)` is what the pipeline uses; any registered column flows through every statistic.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

KINDS = ("count", "positive", "fraction")
DENOMINATORS = ("count", "count_plus_fraction")


@dataclass(frozen=True)
class Metric:
    id: str  # used in file names (<id>_stats.csv, ...)
    column: str  # column in the per-field table
    label: str  # human label for plots and tables
    unit: str
    kind: str  # count | positive | fraction
    description: str = ""
    # fractions only: the nuclei the fraction is taken of, from the field's counted nuclei c and the fraction f:
    # "count" (N = c, e.g. % Ki-67-positive of counted nuclei) or "count_plus_fraction" (N = c / (1 - f), when the
    # fraction's own nuclei are not counted, e.g. out of focus / (counted + out of focus))
    denominator: str = "count"

    def as_dict(self) -> dict:
        return asdict(self)


COUNT = Metric("count", "count", "Nuclei per field", "nuclei per field", "count", "Number of counted nuclei in one field.")


def _readout_info() -> dict:
    try:  # added by the readouts work; absent in older versions
        from .classify import READOUTS  # type: ignore[attr-defined]
    except Exception:
        return {}
    return READOUTS if isinstance(READOUTS, dict) else {}


def _pretty(s: str) -> str:
    return s.replace("_", " ").strip()


def describe_column(col: str) -> Metric | None:
    """The Metric for one column, or None if the column is not an analysable metric."""
    if col == "count":
        return COUNT
    if col.startswith("intensity_"):
        return Metric(col, col, f"{_pretty(col[len('intensity_'):])} mean intensity", "intensity (a.u.)", "positive",
                      "Mean fluorescence intensity per counted cell in one field (background subtracted if set).")
    if col == "nuclei_per_mm2":
        return Metric(col, col, "Nuclei per mm²", "nuclei/mm²", "positive", "Counted nuclei per mm² of imaged area.")
    if col.startswith("readout_"):
        info = _readout_info().get(col, {})
        kind = str(info.get("kind", "")).strip().lower()
        if kind not in KINDS:
            kind = "fraction" if col.startswith("readout_frac_") else "positive"
        unit = info.get("unit") or ("fraction" if kind == "fraction" else "")
        den = str(info.get("denominator") or ("count_plus_fraction" if col == "readout_frac_out_of_focus" else "count"))
        return Metric(col, col, str(info.get("label") or _pretty(col[len("readout_"):])), str(unit), kind,
                      str(info.get("description", "")), den if den in DENOMINATORS else "count")
    return None


def registry(columns, settings: dict | None = None) -> dict[str, Metric]:
    """{metric id: Metric} for the analysable columns among `columns` (count always first, then in column order).
    `settings` is accepted for future per-project overrides and is currently unused."""
    out = {"count": COUNT}
    for c in columns:
        m = describe_column(str(c))
        if m is not None and m.id not in out:
            out[m.id] = m
    return out


def normalise(metrics) -> dict[str, Metric]:
    """Accept a registry, a legacy {column: label} dict, or None; always include count."""
    out = {"count": COUNT}
    for k, v in (metrics or {}).items():
        if isinstance(v, Metric):
            out[k] = v
        else:
            m = describe_column(k) or Metric(k, k, str(v), "", "positive")
            out[k] = Metric(m.id, m.column, str(v) if v else m.label, m.unit, m.kind, m.description, m.denominator)
    return out


def estimand(metric: Metric, unit: str = "well", what: str = "effect") -> str:
    """Plain-language estimand of the all-n model for a metric (`what` = effect | mean).
    unit "well": count/positive -> ratio of mean <unit>; fraction -> difference in percentage points.
    unit "field": geometric means of fields (the older field-level model)."""
    u = metric.unit or metric.label
    if metric.kind == "fraction":
        if what == "mean":
            return "% of nuclei (pooled per well; typical well, logit model)"
        return "difference in percentage points (pooled per well; logit model)"
    if unit == "field":
        return f"geometric mean {u}" if what == "mean" else f"ratio of geometric means of fields ({u})"
    return f"mean {u}" if what == "mean" else f"ratio of mean {u}"
