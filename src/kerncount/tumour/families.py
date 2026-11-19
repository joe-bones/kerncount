"""Factor table + family spec -> contrasts (pure; ADR 0003 amendment).

A **cell** is one combination of factor levels (a group), keyed by a tuple of levels in factor order. A contrast is a
set of weights on cell means of log(readout); exp(estimate) is a fold-change (simple effect) or a ratio of ratios
(interaction). Families are declared generically in settings and never refer to particular factor names here:

- simple effect of factor X: every non-reference level of X vs the reference, within every cell of the other factors;
- interaction of X1 × … × Xk: product of (non-reference − reference) over the listed factors, within every cell of the
  rest. For two factors this is (b / ref in a) ÷ (b / ref in ref), the last factor being the inner ratio.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from .config import Factor, Family

Cell = tuple[str, ...]


@dataclass(frozen=True)
class Contrast:
    family: Family
    level: str  # the non-reference level(s) of the family's factors this contrast is about ("HERV", "TLR7 KO, HERV")
    within: tuple[str, ...]  # levels of the other factors (factor order) the contrast is taken within
    comparison: str  # readable label
    weights: dict[Cell, float]  # non-zero weights on cell means (log scale)


def cells(facs: list[Factor]) -> list[Cell]:
    return list(itertools.product(*(f.levels for f in facs)))


def unknown_factors(facs: list[Factor], fam: Family) -> list[str]:
    names = {f.name for f in facs}
    return [x for x in fam.factors if x not in names]


def _label(fam_facs: list[Factor], levels: tuple[str, ...]) -> str:
    if len(fam_facs) == 1:
        return f"{levels[0]} vs {fam_facs[0].reference}"
    if len(fam_facs) == 2:
        outer, inner = fam_facs
        a, b = levels
        return f"({b} / {inner.reference} in {a}) ÷ ({b} / {inner.reference} in {outer.reference})"
    return "interaction " + " × ".join(f"{lv}/{f.reference}" for f, lv in zip(fam_facs, levels))


def contrasts(facs: list[Factor], fam: Family) -> list[Contrast]:
    """Every contrast of one family over the declared factor table (all declared levels, present or not)."""
    if unknown_factors(facs, fam):
        raise ValueError(f"family {fam.id}: unknown factors {unknown_factors(facs, fam)}")
    idx = {f.name: i for i, f in enumerate(facs)}
    fam_facs = [facs[idx[n]] for n in fam.factors]
    rest = [f for f in facs if f.name not in fam.factors]
    out = []
    for within in itertools.product(*(f.levels for f in rest)):
        for levels in itertools.product(*([lv for lv in f.levels if lv != f.reference] for f in fam_facs)):
            weights: dict[Cell, float] = {}
            # each family factor at its level (+1) or reference (-1); weight = product of signs
            for pick in itertools.product((0, 1), repeat=len(fam_facs)):
                cell: list[str] = [""] * len(facs)
                for f, lv in zip(rest, within):
                    cell[idx[f.name]] = lv
                sign = 1.0
                for f, lv, p in zip(fam_facs, levels, pick):
                    cell[idx[f.name]] = lv if p == 0 else f.reference
                    sign *= 1.0 if p == 0 else -1.0
                weights[tuple(cell)] = sign
            comp = _label(fam_facs, levels) + (f" | {', '.join(within)}" if within else "")
            out.append(Contrast(fam, ", ".join(levels), tuple(within), comp, weights))
    return out
