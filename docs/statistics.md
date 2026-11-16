# Statistics in kerncount

This page explains every statistical result in a generation folder (`generations/gN/`), when to trust it, and
the two machine-readable contracts the browser app reads: `manifest.json` and `stats_long.csv`.

Words used here (see `CONTEXT.md`): an **n** is one independent experiment (one plate, one day) and is the unit
of replication. Each **treatment** (condition) sits in one or more **wells**; each well is imaged in several
**fields** (images). Nesting: field ⊂ well ⊂ treatment ⊂ n.

Absolute values are never replaced by percentages. Every table and plot set gives both: the absolute readout
(nuclei per field, intensity units, ...) and the relative one (ratio or % of each control).

## 1. What is measured: the metric registry

`kerncount.metrics.registry(columns, settings)` decides which per-field columns are analysed. Every registered
column goes through every statistic below.

| column | kind | analysis scale | results reported as |
|---|---|---|---|
| `count` (nuclei per field) | count | log(count + 0.5) | ratio; geometric mean in nuclei per field |
| `intensity_*` | positive | log(value) | ratio; geometric mean in intensity units |
| `nuclei_per_mm2` | positive | log(value) | ratio; geometric mean in nuclei/mm² |
| `readout_*` | from `classify.READOUTS`, else `readout_frac_*` = fraction, others positive | log or logit | ratio, or odds ratio for fractions; mean = proportion |

Why log: effects are multiplicative. Baselines differ between n's (62–90 nuclei per field in the current data),
and a treatment that removes 70% of cells removes 70% whatever the baseline. On the log scale such effects are
additive and the result is a ratio ("0.3 × LyoVec" = 70% fewer nuclei). The 0.5 added to counts only matters for
empty or near-empty fields. Positive metrics with exact zeros get half their smallest positive value added
(noted). Fractions use the empirical logit log((k + 0.5) / (N − k + 0.5)) with k = fraction × nuclei in the
field, N = nuclei in the field (clipped to [0.01, 0.99] when N is unknown).

## 2. The tests, from the most descriptive to the main result

### Within one n (folder `<n>/`)
Fields are **technical replicates**: they show how well one experiment was measured, not whether a treatment
works in general.
- One-way ANOVA and Dunnett vs each control (`<metric>_stats.csv`, stars on `<metric>.png`).
- Two-way ANOVA treatment × well (`<metric>_twoway*`).
- New: `<metric>_fields_by_well.png` (every field, colour = well, marker = position) and
  `<metric>_position_effect.png` (each field divided by its treatment mean, averaged by position, per well).

Trust: describes that experiment only. Never use these p-values as the result of the study.

### Across all n (folder `all_n/`)
Each n is one biological replicate. These outputs are always computed:

1. **Prism standard** (unchanged): repeated-measures one-way ANOVA on the per-n means (raw scale), then Dunnett
   vs each control with the pooled residual error. *Weakness*: one pooled error for every comparison. When the
   controls are very stable between n's but the treatment effects vary between n's (as with the miRNAs), the
   pooled error is too small for the treatments and the stars are overconfident.
2. **Ratio paired t-test** vs each control (Holm within each control): one test per pair, on the log ratio of
   per-n means, no pooling. Honest but weak with 3–4 n's (df = n's − 1).
3. **Two-way ANOVA** treatment × n with fields as replicates: optimistic (fields are not independent
   experiments); kept for comparison only.
4. **Hierarchical (mixed) model** (new, `models.py`) — the default primary test, see section 3.

### Which test draws the stars
`[analysis] primary_test = "mixed"` (default) or `"rm_anova"`. Every plot footer says which test drew its
stars. If the mixed model cannot be fitted for a metric, that metric falls back to the RM ANOVA automatically;
the footer, `model_status` in the manifest and the About sheet say so.

## 3. The hierarchical model

Fitted per metric to the field-level data (one row per field):

```
y = treatment                      (fixed; + timepoint when >= 2 timepoints are known)
  + n                              random intercept: day / plate / passage
  + n × treatment                  random, a SEPARATE variance for each treatment class
  + well within n × treatment      random (only if some condition has >= 2 wells)
  + field position within n        random n:position (only if [design] field_positions = "fixed")
  + field residual
```

- **Treatment classes** (`[analysis] classes`, e.g. `{ mirna = ["miR-*"] }`; `*` is a wildcard). Default: the
  control groups (`reference` + `normalise_to`) form the class `controls`, everything else `treatments`. With
  classes defined, unmatched controls stay in `controls` and other unmatched treatments go to `other`. A class
  with fewer than 3 n × treatment cells cannot have its own variance; one common variance is then used.
  Each class's n × treatment variance says how much that kind of treatment's effect changes from experiment to
  experiment. This is the fix for the overconfident pooled error.
- **Position**: a random n:position term, not fixed position effects. Every treatment is imaged at the same
  positions, so position hardly changes the treatment estimates; it removes a big chunk of field-to-field noise
  (24% of the variance in the current counts) and its variance is reported. It is nested in n because the
  pattern differs between days. With fixed positions the n-level shift and the position pattern share
  information, so read the `n` and `field position` rows of the variance table together.
- **Timepoint** is a fixed effect only when at least two timepoints are known. It changes only between n's, so it
  is confounded with n: its effect is a between-experiment comparison (noted in the model notes).
- **Estimation**: REML. `statsmodels.MixedLM` builds the design (`vc_formula`) and gives one starting point;
  `reml.py` maximises the REML likelihood over variances ≥ 0 from several starts (bounded L-BFGS-B, analytic
  gradient; the best likelihood wins). statsmodels alone fails when a variance is exactly zero, which is a
  legitimate and common estimate here (very stable controls).
- **Fallbacks** (recorded as `model`): class-specific variances → one common n × treatment variance → without
  position → RM ANOVA. Failure = an error, no convergence, or a singular fixed-effect covariance. A variance
  estimated at zero is *not* a failure; it is listed in the notes ("estimated as zero (boundary)").

### Degrees of freedom (why the p-values are honest with 3–6 n's)
The usual Wald z test treats the variance estimates as exact; with 4 n's that roughly triples the false-positive
rate. kerncount uses a t distribution with:

- **Satterthwaite df** for each comparison (computed from the expected REML information of the variances),
- **capped by a containment ceiling**: the df of the n × treatment variance(s) the comparison relies on. For a
  class with several treatments that is `cells − treatments − n's + 1` (its share of the RM ANOVA error df), for
  a one-treatment class `n's − 1`; for a comparison between two classes the two are added; never more than the RM
  ANOVA error df. This stops a variance that happens to be estimated as zero from buying the comparison many df
  (with stable controls, LyoVec vs Control gets df 2 in the current data),
- a comparison needs at least 2 n's that contain both groups; otherwise the estimate is shown with
  "not estimable" and no p-value.
- Estimated means use Satterthwaite capped at (n's with that treatment) − 1.

**Null simulation** (4 n's, 8 treatments, 2 wells × 6 positions, variance components like the real counts, no
true effects; 150 data sets per row; `vs LyoVec` family):

| scenario | per-comparison rate (p < 0.05) | family-wise (adjusted p < 0.05) |
|---|---|---|
| treatments n × treatment variance 0.1, controls 0 | 0.043 | 0.027 |
| no n × treatment variance anywhere | 0.035 | 0.033 |
| treatments 0.3, controls 0 | 0.051 | 0.060 |
| treatments 0.1, controls 0.03 | 0.041 | 0.040 |
| **misspecified**: one treatment 0.3, the rest of its class 0.02 | 0.038 overall, **0.16 for that treatment** | 0.113 |

So the test is close to nominal when the classes are right. It is **not** protected against one treatment that
varies much more than the rest of its class: the class variance is pooled. Check the per-n effect plot
(`<metric>_per_n_ratio_vs_<control>.png`) and leave-one-n-out; if one treatment is much more variable, give it
its own class (e.g. `classes = { mir132 = ["miR-132-5p"], mirna = ["miR-*"] }`). Giving every treatment its own
class is the most conservative setting (each comparison then rests on its own n-to-n variability, like the ratio
paired t-test but with field-level information). The quick version of this simulation runs in
`tests/test_models.py::test_null_type_one_error_is_roughly_nominal_at_four_n`.

### Multiplicity: comparison families
Within each family the p-values are adjusted by the **single-step multivariate t** method (Dunnett-type; the
correlations come from the model, each comparison uses its own df) and simultaneous 95% CIs are given
(`ci_low_adj`, `ci_high_adj`; dashed on the forest plots). Families:

| family | setting | comparisons |
|---|---|---|
| `vs <control>` | every group in `reference` + `normalise_to` | every other treatment vs that control |
| primary | `primary_family` (empty = `reference`) | the family whose stars go on the headline plots |
| `vehicle check (equivalence)` | `vehicle_check` ("auto" = every other control vs `reference`; "" = off; several with `;`) | TOST, Holm across pairs |
| `vs positive control X` | `positive_control` | every non-control treatment vs X |
| `within class C (all pairs)` | `among_classes` | all pairs of treatments in class C |

### Equivalence test for the vehicle (TOST)
"Not significant" for LyoVec vs Control does **not** mean "no vehicle effect"; with 3 n's almost nothing is
significant. The vehicle check asks the right question: is the ratio vehicle / control inside the margin
(`equivalence_margin = 0.2` → between 1/1.2 = 0.83 and 1.2)? Two one-sided t tests (TOST) with the model's SE and
df; the 90% CI is reported.

- **equivalent**: the whole 90% CI lies inside the margin (TOST p < 0.05);
- **not equivalent**: the whole 90% CI lies outside the margin;
- **inconclusive**: anything else — the data cannot tell. With 3–4 n's this is the expected answer.

### Outputs per metric (`all_n/`)
| file | content |
|---|---|
| `<m>_model_means.csv`, `<m>_model_absolute.png` | estimated marginal means per treatment **in absolute units** (geometric means, proportions for fractions) with 95% CI and % of each control; plot with per-n means (coloured) and fields (grey) |
| `<m>_contrasts.csv`, `<m>_forest_<family>.png` | every comparison: ratio (odds ratio), 95% CI, simultaneous CI, % of comparator, t, df, p, adjusted p, stars / equivalence verdict |
| `<m>_variance_components.csv/.png` | variance per level (n, n × treatment per class, well, position, field) on the log scale, SD, ≈ CV, % of total, approximate CI (Wald on log-variance from the expected information; flagged "too uncertain" or "at zero" when not meaningful; residual: chi-square) |
| `<m>_diagnostics.png` | residuals vs fitted, normal QQ, estimated random effects per n |
| `<m>_model_info.csv` | model used, formula, df rule, convergence, estimation, notes, failed attempts |
| `<m>_loo.csv`, `<m>_loo_summary.csv` | leave-one-n-out with the primary test (estimate, CI, p, adjusted p for every drop) and per comparison `stable = yes/no` (same direction **and** same significance whichever n is dropped) |
| `<m>_per_n_ratios.csv`, `<m>_per_n_ratio_vs_<control>.png` | the effect in every n: ratio of geometric means (model scale) and ratio of arithmetic means |
| `count_position_effect.png` | the position pattern per n |

Note on estimates: the model ratio is a ratio of **geometric** means of fields. When the treated fields are much
more uneven than the control fields (miRNA wells with a few dense and many sparse fields), it is smaller than the
ratio of arithmetic means used by the % tables and the ratio t-test. Both are in the tables.

## 4. Acceptance criteria per n (`[acceptance]`)
Checked on the counts of every n and written to `acceptance.csv` (n, criterion, value, limit, result = pass /
fail / not applicable, detail; plus one `overall` row per n):

| key | default | meaning |
|---|---|---|
| `vehicle_within` | 0.3 | vehicle / control (vehicle-check pair) mean nuclei per field within 1 ± 0.3 |
| `positive_control_max_ratio` | "" (off) | positive control / primary reference ≤ value (or ≥ value if value > 1) |
| `min_images_per_condition` | 3 | every condition has at least this many usable images |
| `max_low_qc_fraction` | 0.5 | at most this share of the n's images are QC low |
| `enforce` | false | true: failing n's are left out of every all-n statistic (listed in About and the manifest) |

## 5. `stats_long.csv`: every result in one table

One row per test result, for interactive filtering. Columns:

| column | meaning |
|---|---|
| `level` | `n` (within one experiment) or `all_n` |
| `n` | the n for level `n`, else empty |
| `metric` | metric id (see `manifest.json` `metrics`) |
| `family` | `vs <control>`, `vehicle check (equivalence)`, `vs positive control X`, `within class C (all pairs)`, `omnibus`, `all pairs (two-way)`, `estimated means`, `<family> (leave-one-n-out)`, `vs <control> (per-n ratio)` |
| `test` | plain name of the test (e.g. `Dunnett vs LyoVec (RM ANOVA, pooled error; Prism standard)`) |
| `comparison` | `A vs B`, a treatment (estimated means), a factor (omnibus) or `all` |
| `estimate` | the effect |
| `estimate_type` | `ratio`, `odds_ratio`, `difference` (raw-scale difference of means), `mean` (absolute estimated mean); empty for omnibus tests |
| `ci_low`, `ci_high` | 95% CI (90% for the equivalence test); empty if not available |
| `p`, `p_adj` | raw and family-adjusted p |
| `significance` | stars (`ns`, `*` … `****`) or the equivalence verdict |
| `N` | n's in the comparison (fields for per-n ratios) |
| `model` | `classical`, `classical (fields as replicates)`, `descriptive`, or the mixed model description |
| `note` | df, "without nX" for leave-one-out, "not estimable ...", etc. |

## 6. `manifest.json`: every plot and table of a generation

```jsonc
{
  "schema": "kerncount.generation.manifest/1",
  "generation": "g1", "source": "auto", "built": "2026-11-16T17:12:05",
  "metrics": [ {"id": "count", "label": "Nuclei per field", "unit": "nuclei per field", "kind": "count",
                "description": "..."} ],
  "primary": { "test": "mixed" | "rm_anova", "family": "vs LyoVec", "reference": "LyoVec",
               "fallback_metrics": [],              // metrics where the mixed model failed (RM ANOVA stars)
               "controls": ["Control", "LyoVec"], "classes": {"controls": [...], "mirna": [...]},
               "positive_control": "CL075" | null, "vehicle_check": ["LyoVec vs Control"],
               "equivalence_margin": 0.2 },
  "model_status": { "<metric id>": {
       "status": "ok" | "fallback_rm_anova" | "not_estimable",
       "model": "class-specific n×treatment variances + well + position",
       "spec": "full" | "common" | "no_position" | "", "formula": "...", "converged": true,
       "df_rule": "...", "response": "log(count + 0.5)", "primary_test_used": "...",
       "stars_from": "mixed" | "rm_anova", "notes": [...], "attempts": [...], "warnings": [...],
       "optimiser": "..." } },
  "acceptance": { "enforce": false, "n_checked": 4, "failing": [], "excluded": [],
                  "criteria": {"vehicle_within": 0.3, ...} },
  "notes": [ "..." ],
  "files": [ {
       "path": "all_n/count_forest_vs_LyoVec.png",   // relative to the generation folder, "/" separators
       "pdf": "all_n/count_forest_vs_LyoVec.pdf",    // only for plots that also have a PDF
       "level": "treatment" | "n" | "all_n" | "generation",
       "n": "n1_2026-06-05" | null, "treatment": "CL075" | null, "metric": "count" | null,
       "kind": "absolute" | "relative" | "model" | "contrast" | "robustness" | "diagnostic" | "design"
               | "acceptance" | "data",
       "title": "All n: Nuclei per field: each treatment compared with LyoVec",
       "description": "One plain sentence a biologist understands." } ]
}
```

Guarantees: every `.png`, `.csv` and `.xlsx` in the generation folder is listed exactly once (PDFs are attached
to their PNG), every listed path exists, `level`/`kind` take only the values above, and values that are not
numbers (NaN) are written as `null`. Files are described from their names, so a new output that follows the
`<metric>_<what>` naming is listed automatically (with a generic description if no rule matches).

## 7. Excel workbook (`results.xlsx`)
All previous sheets are kept. New sheets: `Model - means (absolute)`, `Model - ratios vs controls`,
`Comparison families`, `Variance components`, `Robustness - stability`, `Robustness - leave-one-out`,
`Robustness - effect per n`, `Model status`, `Acceptance`. Metrics other than intensities (readouts, densities)
get `Other metrics - ...` sheets like the intensity ones. The About sheet states the primary test, the classes,
the model per metric and the acceptance result.

## 8. When to trust what — short version
- **Headline stars (mixed model)**: trustworthy when the per-n effect plot shows every n pointing the same way and
  leave-one-n-out says `stable = yes`. With 3–4 n's, read the CI, not only the stars.
- **RM ANOVA stars**: comparable with Prism; overconfident for treatments whose effect varies between n's.
- **Ratio paired t-test**: the most cautious; with 3 n's it rarely reaches p < 0.05 even for real effects.
- **Within-n and two-way p-values**: describe the measurement, not the biology.
- **Vehicle check**: only "equivalent" supports "the vehicle has no effect"; "inconclusive" means more n's are
  needed.
- **Variance components**: the n-level variances rest on 3–4 levels; treat them as rough.
