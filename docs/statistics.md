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

| column | kind | analysis scale (model unit = well, default) | headline effect (estimand) |
|---|---|---|---|
| `count` (nuclei per field) | count | log(well mean + 0.5) | **ratio of mean nuclei per field**; means in nuclei per field |
| `intensity_*` | positive (label / unit from `classify.READOUTS.get(column)`, i.e. `classify.channel_metric_info` for per-channel columns) | log(well mean) | ratio of mean intensity; means in intensity units |
| `nuclei_per_mm2` | positive | log(well mean) | ratio of mean nuclei/mm² |
| `readout_*` | from `classify.READOUTS`, else `readout_frac_*` = fraction, others positive | log, or logit of the pooled well fraction | ratio; fractions: **difference in percentage points** (odds ratio and % of reference secondary); means in % of nuclei |

The estimand is written into every table (`estimand` column of contrasts, model means and `stats_long.csv`), the
manifest (`estimands`) and every model plot footer.

### Model unit and estimand (`[analysis] model_unit`)
Biologists read "0.4 × LyoVec" as "40% of the mean nuclei per field of LyoVec". The earlier model used one row per
field on log(count + 0.5); that compares **geometric** means of fields, which sparse, uneven miRNA fields drag down
(real data: ratios vs LyoVec 0.21–0.36 against 0.27–0.52 for ratios of mean nuclei per field).

- `model_unit = "well"` (default): one row per well. Count / positive metrics: the arithmetic mean of the well's
  fields, analysed as log(mean). Fractions: positives and nuclei pooled over the well's fields, K / N, analysed on
  the empirical logit. Every well is imaged at the same positions, so position effects cancel in well means.
- **Missing or excluded fields** (e.g. 11 images): a plain mean of the fields present would inherit the position
  pattern (bright / dense centre fields). Such a well gets a **position-adjusted mean** = Σ fields present /
  Σ r_p over their positions, where r_p is that n's position pattern (field / well mean, averaged over the n's
  complete wells, normalised to mean 1). This is a ratio estimator: it equals the plain mean for a complete well
  and is unbiased for a well missing fields when the pattern is multiplicative. Without a complete well in that n,
  or with `[design] field_positions = "random"`, the plain mean is used (noted). Fractions are pooled over the
  fields present without position adjustment (the pooled fraction already weights fields by their nuclei).
  Tested in `tests/test_estimands.py::test_position_adjusted_well_mean_with_missing_fields`.
- `model_unit = "field"`: the earlier field-level model (position and well random terms), kept as an option.
- Small-sample note: the log of a mean of 6 very uneven fields is itself slightly low (Jensen). In a simulation
  with treated fields of log-SD 1.1 and a true ratio 0.5, the well unit gave 0.40 on average, the ratio of the pooled
  means 0.44, the field unit 0.27.

### Fractions (`readout_frac_*`, % Ki-67-positive)
Model on the logit scale with one row per well: y = log((K + 0.5) / (N − K + 0.5)). This is a binomial-type model
with overdispersion: each well carries its **known binomial sampling variance** 1/(K + 0.5) + 1/(N − K + 0.5)
(added to V; the variance of the empirical logit), and the well-to-well, n and n × treatment variances are
estimated on top. Why not simply weight by nuclei: a nuclei weight lets a dense well dominate even when
well-to-well biology is the main noise; a pure binomial model ignores that overdispersion and is overconfident.
Here a well with few nuclei counts less, and a big well cannot buy precision the biology does not give.

**Headline: difference in percentage points** vs each reference, with its 95% (and simultaneous) CI, beside the
absolute % per condition (`percent`, `percent_ci_low/high` in `<m>_model_means.csv`). The pp difference is
expit(μ_a) − expit(μ_b) at the model means. Its CI is computed by back-transforming the model's uncertainty:
4000 draws of (μ_a, μ_b) from the bivariate t with the model covariance and the contrast's df, then the 2.5 / 97.5 %
quantiles of 100·(expit(μ_a) − expit(μ_b)) (simultaneous CI: the quantiles matching the family's critical value;
fixed seed, so reproducible). Why not the delta method: it was tried first and, with the real out-of-focus
fractions (0.2–3 %), its symmetric intervals ran below −p_b for every comparison and had to be clipped to the
impossible-value boundary; the back-transformed draws stay in the possible range, are asymmetric near 0 % and
include the uncertainty of both levels. Back-transforming only the two marginal CIs separately cannot give a CI
for their difference. The p-value and stars come from the logit-scale test, which tests the same null hypothesis
p_a = p_b and behaves better near 0 % than a test on the % scale (so a pp CI can just touch 0 when the logit test is
borderline). The odds ratio
(`estimate`, `ci_low/high`) and % of the reference (`pct_of_comparator` = 100 p_a / p_b) are in the tables, never
the headline: the out-of-focus fraction is 0 in most control images, and a ratio of near-zero fractions means
nothing. The model % is that of a typical well (logit-normal median); `observed_mean` gives the pooled fraction
averaged over n's for comparison. Nuclei N per field: the metric's own N column when it has one
(`Metric.n_column`, e.g. `n_measured_Ki67` for `readout_frac_positive_Ki67` and its p90 variant: the counted,
non-saturated nuclei the fraction is taken of), else counted nuclei (`count`), or count / (1 − f) for
`readout_frac_out_of_focus` (out of focus / (counted + out of focus); `Metric.denominator`). N is the binomial N of
the empirical logit, the weight when fields are pooled per well, and the `nuclei` reported. The registry's N columns
(`metrics.n_columns(registry)`) are carried into the generation next to the metric columns.

`Metric` fields: `id`, `column`, `label`, `unit`, `kind`, `description`, `denominator`, `n_column`.

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

Fitted per metric. Default unit = well (one row per well, section 1); with `model_unit = "field"` one row per field:

```
y = treatment                      (fixed; + timepoint when >= 2 timepoints are known)
  + n                              random intercept: day / plate / passage
  + n × treatment                  random, a SEPARATE variance for each treatment class
  + residual                       well unit: well-to-well (+ known binomial variance for fractions)
  field unit only:
  + well within n × treatment      random (only if some condition has >= 2 wells)
  + field position within n        random n:position (only if [design] field_positions = "fixed")
  + field residual
```

- **Treatment classes** (`[analysis] classes`). **Default `"auto"`: one class per treatment, each control
  (`reference` + `normalise_to`) included as its own class.** Each comparison then rests on the n-to-n
  variability of its own two groups, like the ratio paired t-test but using all wells. Pooling the controls (the
  previous "auto") let a vehicle that varies between n's borrow the stability of Control: "vs LyoVec" had a type
  I error of 0.13 at 4 n's when LyoVec varies and Control does not (row "vehicle varies" below; now nominal).
  Pooling several treatments in one class was shown to give a type I error of 0.16 for a treatment that varies
  more than the rest of its class (table below); `"auto"` removes that failure mode. Other settings: `{}` =
  `controls` vs everything else (`treatments`; the earlier default), or a table `{ mirna = ["miR-*"] }` (`*` is a
  wildcard; unmatched controls stay in `controls`, other unmatched treatments go to `other`). A class with
  several treatments needs >= 3 n × treatment cells, a one-treatment class >= 2 (2 n's); otherwise one common
  variance is used (noted). A treatment present in only one n has no n × treatment term (its single cell is
  absorbed by its fixed effect; its comparisons are "needs ≥ 2 n").
- **One well per condition** (well unit): the well-to-well residual is not identifiable next to n × treatment;
  it is held at ~0 and reported as "not identifiable".
- **Position** (field unit): a random n:position term, not fixed position effects. Every treatment is imaged at the same
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

### Background-subtracted intensities: some well means ≤ 0 (raw scale)
A positive-kind metric (intensities, areas, ...) is analysed on the log scale, so its effects are ratios. A
background-subtracted intensity of a dim stain can have well means at or below 0 (e.g. true mean 1, field SD 3).
The log of those wells is impossible; dropping them (as before) keeps only the brighter wells and biases the
effect toward 1. Now, when **any analysed well mean is ≤ 0** (field unit: any field < 0;
`models.needs_raw_scale`), the whole metric is analysed **on the raw scale as a difference** (identity link,
intensity units; `metrics` kind `difference`):

- effects are **differences** in intensity units (`estimate_type = "difference"`, 0 = no effect); means are raw
  means; variance components are on the raw scale (no CV); per-n effects are differences of the field means with
  their technical intervals; there is no % of reference;
- it is flagged everywhere: `model_status[...]` `raw_scale = true`, `effect_type = "difference"` and
  `raw_scale_warning`; the manifest's `estimands[m].headline = "difference"`; the model notes and
  `<m>_model_info.csv` start with "RAW SCALE (difference) ..."; every plot footer of the metric; `stats_long.csv`
  `note` of every model row; the layer ladder draws that layer on its own linear axis (orange panel, "DIFFERENCE,
  raw units", 0 = no effect) and says so in its footer; the classical ratio t-test notes the n's it had to leave
  out (a ratio of a mean ≤ 0 is impossible);
- the equivalence margin for such a metric is ± margin × the comparator mean (not testable when that mean ≤ 0).

Tested with a true mean of 1, field SD 3, 4 n's (`tests/test_estimands.py`): the +2 difference is recovered
with a CI covering it, and an all-positive copy of the same data stays on the log scale.

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

**Null simulation of the default** (`model_unit = "well"`, `classes = "auto"`): 8 treatments (Control, LyoVec,
CL075, 5 miRNAs), 2 wells × 6 fixed positions, Poisson-lognormal field counts (~70 per field; n variance 0.05,
well 0.03, position 0.1, field log-SD 0.45), no true effect on the mean nuclei per field; with 4 n's LyoVec is
missing in one n (as in the real data). `vs LyoVec` family, 7 comparisons per data set; 400 data sets per row,
re-run for the current default (every control its own class; 2026-11-17). Script: `tests/sim_null.py`
(`python tests/sim_null.py . 400 well:auto "" 3,4,6`; quick versions in
`tests/test_estimands.py::test_type_one_error_default_including_misspecification` and
`tests/test_models.py::test_null_type_one_error_is_roughly_nominal_at_four_n`).

| scenario (n × treatment variances) | n's | per-comparison (p < 0.05) | family-wise (adj. p < 0.05) | the odd treatment (miR-1) |
|---|---|---|---|---|
| treatments 0.1, controls 0 | 3 / 4 / 6 | 0.024 / 0.026 / 0.034 | 0.020 / 0.018 / 0.043 | 0.020 / 0.038 / 0.033 |
| none anywhere | 3 / 4 / 6 | 0.005 / 0.009 / 0.022 | 0.000 / 0.003 / 0.010 | 0.003 / 0.020 / 0.015 |
| treatments 0.3, controls 0 | 3 / 4 / 6 | 0.043 / 0.043 / 0.041 | 0.053 / 0.060 / 0.070 | 0.033 / 0.043 / 0.048 |
| **misspecified for pooled classes**: miR-1 0.3, other treatments 0.02 | 3 / 4 / 6 | 0.016 / 0.016 / 0.026 | 0.013 / 0.008 / 0.010 | **0.063 / 0.045 / 0.038** |
| controls unequal: LyoVec 0.1, Control 0, treatments 0.1 | 3 / 4 / 6 | 0.035 / 0.038 / 0.041 | 0.028 / 0.033 / 0.043 | 0.023 / 0.008 / 0.038 |
| uneven fields: miRNA field log-SD 1.0 vs 0.4, equal means | 3 / 4 / 6 | 0.014 / 0.014 / 0.040 | 0.015 / 0.015 / 0.035 | 0.025 / 0.018 / 0.048 |
| **vehicle varies for pooled controls**: LyoVec 0.3, Control 0, treatments 0 | 3 / 4 / 6 | **0.048 / 0.048 / 0.040** | 0.030 / 0.033 / 0.030 | 0.053 / 0.050 / 0.038 |

Comparison at 4 n's, 200 data sets (earlier settings):

| settings | misspecified: per-comparison / family-wise / miR-1 | uneven fields: per-comparison / family-wise | vehicle varies: per-comparison (4 n) |
|---|---|---|---|
| well unit, `classes = {}` (controls vs the rest) | 0.046 / 0.090 / **0.225** | 0.021 / 0.010 | |
| field unit, `classes = {}` (the previous default) | 0.043 / 0.105 / **0.230** | **0.172 / 0.240** | |
| well unit, `"auto"` with the controls pooled (the previous "auto"; independent review) | | | **0.133** |

So: (1) pooling the n × treatment variance of several treatments lets the most variable one through about one
time in four-five; one class per treatment ("auto") keeps it at nominal. (2) The field-level model tests geometric
means of fields: when treated fields are much more uneven, it reports "effects" on the mean nuclei per field that
do not exist (17% false positives); the well unit does not. (3) The default is somewhat conservative when there is
little n-to-n variation (rates 0.003–0.03), the price of letting every treatment carry its own n-to-n variance with
3–6 n's. The per-comparison rate stays at or below ~0.05 in every scenario. (4) When the vehicle varies between
n's and Control does not, pooling the two controls in one class gave "vs LyoVec" a type I error of 0.13 at 4 n's;
with every control in its own class it is 0.048 (and 0.038–0.041 instead of 0.050–0.057 in the milder "controls
unequal" scenario).

### Multiplicity: comparison families
Within each family the p-values are adjusted by the **single-step multivariate t** method (Dunnett-type; the
correlations come from the model, each comparison uses its own df; the numerical integration is kept within the
exact bounds raw p ≤ p_adj ≤ min(1, k × raw p), so a far-tail p_adj is never below the raw p nor exactly 0) and
simultaneous 95% CIs are given
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
significant. The vehicle check asks the right question: is the vehicle / control effect inside the margin? Two
one-sided t tests (TOST) with the model's SE and df, Holm-adjusted across the vehicle pairs (m pairs).

- **margin**: count and positive metrics: the ratio within `[1/(1 + equivalence_margin), 1 + equivalence_margin]`
  (0.2 → 0.83–1.20); **fractions** (% positive): the difference within **± `equivalence_margin_pp`
  percentage points** (default 5; odds ratios are not used: an odds-ratio margin means very different things at
  5% and at 50% positive); raw-scale metrics (section 3, background-subtracted intensities): ± margin × the
  comparator mean, in intensity units. Rows carry `margin_low` / `margin_high`.
- **equivalent**: the Holm-adjusted TOST p < 0.05, i.e. the whole (1 − 2·0.05/m) CI lies inside the margin (90%
  CI for one pair, 95% for two; the level is in `test` and `note`);
- **not equivalent**: that CI lies wholly outside the margin;
- **inconclusive**: anything else — the data cannot tell. With 3–4 n's this is the expected answer.

The verdict and the reported CI use the same adjusted level (earlier versions called a pair "equivalent" from
the unadjusted 90% CI even when its Holm-adjusted p was above 0.05). Fraction intervals and p come from the same
draws of the model means as the percentage-point CIs (4000 draws, fixed seed).

### Outputs per metric (`all_n/`)
| file | content |
|---|---|
| `<m>_model_means.csv`, `<m>_model_absolute.png` | estimated marginal means per treatment **in absolute units** (mean per field; proportions + `percent` for fractions) with 95% CI, `estimand`, `observed_mean`, `wells`, `fields`, `note` and % of each control; plot with per-n values (coloured) and fields (grey) |
| `<m>_contrasts.csv`, `<m>_forest_<family>.png` | every comparison: `estimand`, ratio (odds ratio for fractions), 95% CI, simultaneous CI, `pct_of_comparator`; fractions also `pp_difference` with `pp_ci_low/high(_adj)`, `pct_group`, `pct_comparator`; t, df, p, adjusted p, stars / equivalence verdict, note. Fraction forest plots show percentage points |
| `<m>_variance_components.csv/.png` | variance per level (n, n × treatment per class, well, position, field) on the log scale, SD, ≈ CV, % of total, approximate CI (Wald on log-variance from the expected information; flagged "too uncertain" or "at zero" when not meaningful; residual: chi-square) |
| `<m>_diagnostics.png` | residuals vs fitted, normal QQ, estimated random effects per n |
| `<m>_model_info.csv` | model used, formula, df rule, convergence, estimation, notes, failed attempts |
| `<m>_loo.csv`, `<m>_loo_summary.csv` | leave-one-n-out with the primary test (estimate, CI, p, adjusted p for every drop) and per comparison `stable = yes/no` (same direction **and** same significance whichever n is dropped); with 2 n's `stable = "not testable"` (dropping one leaves one n: no test). Each refit is a multi-start REML fit like the main one (every default start plus the full fit's variances; only the slow statsmodels start is skipped) |
| `<m>_per_n_ratios.csv`, `<m>_per_n_ratio_vs_<control>.png` | the effect in every n: `ratio` (model scale), `ratio_of_means`, `pp_difference` (fractions), `difference` (raw-scale metrics), technical interval `tech_ci_low/high` + `interval` (descriptive) |
| `<fraction>_density_adjusted_means.csv`, `_density_adjusted_contrasts.csv`, `_density.png` | density-adjusted model (section 3, density) |
| `count_position_effect.png` | the position pattern per n |

Note on estimates: with the default well unit the model ratio is a ratio of **mean nuclei per field** (averaged
over wells and n's on the log scale), the quantity of the % tables and the ratio t-test. With `model_unit =
"field"` it is a ratio of geometric means of fields, which is smaller when the treated fields are much more uneven.

### One n, and comparisons without 2 n's
The all-n level keeps **exactly the same files, columns and plot layout** whatever the number of n's.
- With one n, estimates are filled in from that n (N = 1): model means, ratios / pp differences, per-n effects,
  Dunnett differences, ratio t-test ratios. Every across-n CI, p-value and star is empty (NaN in tables) with the
  note **"needs ≥ 2 n"** (`note` columns; `model_status.across_n_inference = false`; plot footers). No stars are
  drawn. Field or well variation is never used to make one experiment look like a biological test: the two-way
  ANOVA treatment × n is written as empty tables and a plot without tests; the variance-components table shows
  only the well-to-well variance, labelled "descriptive: technical replicates".
- The same applies with 2+ n's to any comparison where fewer than 2 n's contain both groups ("needs ≥ 2 n
  containing both groups"), and to the ratio t-test rows with fewer than 2 pairs.
- The per-n level carries technical intervals, labelled **"descriptive: technical replicates"**: the 95% CI of
  the field mean (`ci_low`, `ci_high`, `interval` in `<n>/<m>_summary.csv`) and the technical 95% interval of each
  per-n effect (`tech_ci_low`, `tech_ci_high` in `all_n/<m>_per_n_ratios.csv`; from the field SEs, delta method on
  the log ratio, or the ratio-estimator SE of a pooled fraction).

### Cell lines are never pooled
If the per-image table carries a `cell_line` column (from the design layer; "" = unknown) with more than one
non-empty value, the all-n analysis runs once per cell line in **`all_n/<cell line>/`** (rows with an empty
cell line form `all_n/unknown/`). One line (or no column): `all_n/` as before. Per-n folders are unaffected.
One named line plus n's without a cell line stay one block (there is nothing to split them by), but it is
labelled **"U251 + unlabelled"** (`cell_line` in stats_long / manifest), never plain "U251", with a warning and a
note naming the unlabelled n's. The planner follows the same blocks (docs/planner.md).
Manifest: `cell_lines`, `all_n_folders` (cell_line, path, n_count), `cell_line` on every file entry,
`model_status` keyed `<cell line>/<metric>`; `stats_long.csv` has a `cell_line` column.

### Negative control
`[experiment] negative_control` (a treatment name; "" = none) only sets positivity cut-offs. Its images are kept
in `counts_used.csv` (excluded, source "negative control (cut-off only; not analysed)") and left out of every test,
plot and summary of the generation; listed in the About sheet, README and manifest (`negative_control`, `notes`).

### Density-adjusted estimates (fraction metrics)
`[analysis] density_covariate` (default `readout_coverage`; "" = off). For every fraction metric a second model
adds the well mean of the covariate (centred; position-adjusted like the response) as a fixed effect on the logit
scale. Adjusted means and pp differences are at the average density of all wells; the raw estimate stays the
headline. Outputs: `<m>_density_adjusted_means.csv` (with `raw_percent`), `<m>_density_adjusted_contrasts.csv`
(with `raw_pp_difference`, the slope, `density_overlap`, `density_flag`), `all_n/<m>_density.png` (wells, colour =
treatment, marker = n, lines = adjusted model per treatment) and `<n>/<m>_density.png` (fields with a descriptive
common-slope fit within that n). Stats_long families `... (density-adjusted)`.

**Overlap criterion**: for a comparison a vs b, overlap = min(share of a's wells whose density lies inside b's
min–max well range, share of b's wells inside a's range). If overlap < 0.5 the row is flagged **"conditions do not
overlap in density; adjusted estimate not trustworthy"**: the adjustment would then extrapolate the density slope
outside the data of one group, so it rests on the assumed straight line, not on comparable wells. Treatments that
kill or detach cells (the miRNAs) usually fail this check, which is the point: their fractions cannot be separated
from their crowding by this data.

### Duplicate raw images
`build_generations` hashes (SHA-256) every raw file listed in `qc_details.csv` (`source_file`), caching hashes by
path, size and mtime in `_work/file_hashes.json` (a re-run reads only new or changed files). Images sharing
content are written to `duplicates.csv` in the project folder (columns group, n, treatment, image, source_file,
hash, copies, cross_condition, shared_with; copied into the generation), a WARNING is logged, the affected n /
treatments are listed in the About sheet, README and manifest (`duplicates`, `notes`), and the acceptance
criterion `no_duplicate_images` fails for an n with an image shared across conditions.

## 4. Acceptance criteria per n (`[acceptance]`)
Checked on the counts of every n and written to `acceptance.csv` (n, criterion, value, limit, result = pass /
fail / not applicable, detail; plus one `overall` row per n):

| key | default | meaning |
|---|---|---|
| `vehicle_within` | 0.3 | vehicle / control (vehicle-check pair) mean nuclei per field within 1 ± 0.3 |
| `positive_control_max_ratio` | "" (off) | positive control / primary reference ≤ value (or ≥ value if value > 1) |
| `min_images_per_condition` | 3 | every condition has at least this many usable images |
| `max_low_qc_fraction` | 0.5 | at most this share of the n's images are QC low |
| `no_duplicate_images` | true | no raw image of that n is byte-identical to an image of another condition (not applicable when the files were not checked) |
| `enforce` | false | true: failing n's are left out of every all-n statistic (listed in About and the manifest) |

## 5. `stats_long.csv`: every result in one table

One row per test result, for interactive filtering. Columns:

| column | meaning |
|---|---|
| `level` | `n` (within one experiment) or `all_n` |
| `n` | the n for level `n`, else empty |
| `metric` | metric id (see `manifest.json` `metrics`) |
| `family` | `vs <control>`, `vehicle check (equivalence)`, `vs positive control X`, `within class C (all pairs)`, `omnibus`, `all pairs (two-way)`, `estimated means`, `<family> (leave-one-n-out)`, `<family> (density-adjusted)`, `vs <control> (per-n ratio)`, `vs <control> (per-n difference)` (fractions) |
| `test` | plain name of the test (e.g. `Dunnett vs LyoVec (RM ANOVA, pooled error; Prism standard)`) |
| `comparison` | `A vs B`, a treatment (estimated means), a factor (omnibus) or `all` |
| `estimate` | the effect |
| `estimate_type` | `ratio`, `odds_ratio`, `difference` (raw-scale difference of means: the RM ANOVA / Dunnett rows, and **every model row of a raw-scale metric**, section 3), `pp_difference` (difference in percentage points; the headline for fractions), `mean` (absolute estimated mean); empty for omnibus tests. A model comparison of a fraction has two rows with the same test: `pp_difference` (headline) and `odds_ratio` |
| `ci_low`, `ci_high` | 95% CI (90% for the equivalence test); empty if not available |
| `p`, `p_adj` | raw and family-adjusted p |
| `significance` | stars (`ns`, `*` … `****`) or the equivalence verdict |
| `N` | n's in the comparison (fields for per-n ratios) |
| `model` | `classical`, `classical (fields as replicates)`, `descriptive`, or the mixed model description |
| `note` | df, "without nX" for leave-one-out, "needs ≥ 2 n", density overlap / flag, technical-interval label, "RAW SCALE: difference, some well means <= 0" (raw-scale metrics), etc. |
| `estimand` | (new) what the estimate is, e.g. `ratio of mean nuclei per field`, `difference in percentage points (pooled per well; logit model)`, `difference in mean nuclei per field` |
| `cell_line` | (new) cell line of the all-n analysis / of the n; empty with one (unknown) cell line |

Per-n rows (`level = n`, family `(per-n ratio)` / `(per-n difference)`) are descriptive with a technical
interval in `ci_low/ci_high` (from the fields; label in `note`). With one n every all-n `p`, `ci_*` is empty and
`significance` is empty.

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
               "equivalence_margin": 0.2, "classes_mode": "auto" | "two" | "user" },
  "model_unit": "well" | "field",
  "estimands": { "count": {"effect": "ratio of mean nuclei per field", "mean": "mean nuclei per field",
                           "headline": "ratio" | "pp_difference"} },
  "density_covariate": "readout_coverage" | null,
  "n_count": 4,
  "cell_lines": [],                                 // non-empty cell lines; >= 2 = separate all-n folders
  "all_n_folders": [ {"cell_line": null, "path": "all_n", "n_count": 4} ],
  "negative_control": {"treatment": null, "present": false},
  "layers": { "Ki67": {                             // readout ladder per channel ({} for DAPI-only projects)
       "layers": [ {"index": 1, "key": "field_raw", "label": "1 Field, raw", "metric": "intensity_Ki67_field_raw",
                    "step": ""}, ... ],             // 5 entries in ladder order; "step" = what changes from the previous
       "extras": [ {"key": "stained_area_mean", "label": "Stained-area MFI", "metric": "..."}, ... ],
       "references": ["LyoVec", "Control"],
       "files": {"effects": "layers_Ki67_effects.csv", "absolute": "layers_Ki67_absolute.csv",
                 "concordance": "layers_Ki67_concordance.csv", "ladder_plot": "layers_Ki67_vs_[reference].png",
                 "absolute_plot": "layers_Ki67_absolute.png"} } },
  "duplicates": {"checked": true, "images": 24, "cross_condition": 24,
                 "affected": {"n3_2026-09-29": ["Control", "LyoVec"]}, "summary": ["..."]},
  "model_status": { "<metric id>" (or "<cell line>/<metric id>"): {
       "status": "ok" | "fallback_rm_anova" | "not_estimable",
       "model": "class-specific n×treatment variances + well + position",
       "spec": "full" | "common" | "no_position" | "", "formula": "...", "converged": true,
       "df_rule": "...", "response": "log(count + 0.5)", "primary_test_used": "...",
       "stars_from": "mixed" | "rm_anova" | "none", "notes": [...], "attempts": [...], "warnings": [...],
       "optimiser": "...", "unit": "well", "estimand": "...", "mean_estimand": "...", "n_count": 4,
       "across_n_inference": true, "density_model": "...", "cell_line": null, "metric": "count" } },
  "acceptance": { "enforce": false, "n_checked": 4, "failing": [], "excluded": [],
                  "criteria": {"vehicle_within": 0.3, ...} },
  "notes": [ "..." ],
  "files": [ {
       "path": "all_n/count_forest_vs_LyoVec.png",   // relative to the generation folder, "/" separators
       "pdf": "all_n/count_forest_vs_LyoVec.pdf",    // only for plots that also have a PDF
       "level": "treatment" | "n" | "all_n" | "generation",
       "n": "n1_2026-06-05" | null, "treatment": "CL075" | null, "metric": "count" | null,
       "cell_line": "U251" | null,
       "kind": "absolute" | "relative" | "model" | "contrast" | "robustness" | "diagnostic" | "design"
               | "acceptance" | "data" | "layers",
       "channel": "Ki67", "reference": "LyoVec" | null,   // only for kind "layers" (metric is null there)
       "title": "All n: Nuclei per field: each treatment compared with LyoVec",
       "description": "One plain sentence a biologist understands." } ]
}
```

Guarantees: every `.png`, `.csv` and `.xlsx` in the generation folder is listed exactly once (PDFs are attached
to their PNG), every listed path exists, `level`/`kind` take only the values above, and values that are not
numbers (NaN) are written as `null`. Files are described from their names, so a new output that follows the
`<metric>_<what>` naming is listed automatically (with a generic description if no rule matches).

### Layer comparison (readout ladder)

For every marker channel with at least two ladder layers, each generation writes `layers_<ch>_*` files in `all_n/`
(or `all_n/<cell line>/`) and in every `<n>/` (docs/intensity.md, section 6b): a ladder plot per reference
(rows = treatments, columns = layers, effect with its interval), an absolute version, the effects table and a
concordance table that names the analysis step that changes the conclusion. They are built from the per-layer
results above (the layers are ordinary registry metrics), so their rows in `stats_long.csv` are the usual ones of
each layer metric (no extra rows): all n = the model contrasts (one n: estimates, "needs ≥ 2 n"); per n = the
descriptive per-n effects with technical intervals. Manifest kind: `"layers"`.

### When a generation is rebuilt

`generation.json` records the counts fingerprint and, since this version, `inputs`: fingerprints of the per-image
intensity columns (`intensity_*`, `readout_frac_positive_*`, `n_measured_*`) and of the design / metadata, with a
readable summary (mean per column and n, cut-off rule and value per n and channel, cell line per n). On `run` /
`regen`:

- counts changed (automatic counts or a manual column): **not overwritten**, warning with the command to rebuild
  (`--rebuild gN`), as before; the warning also lists any intensity / design change;
- counts unchanged but intensity results changed (new cut-off lock, positivity, flat-field or measurement settings)
  or the design / metadata changed: **rebuilt automatically**, because nothing manual lives in a generation and a
  stale one would contradict `intensity/` (rewritten on every run). The log line and `generation.json`
  (`rebuilt_because`) say exactly what changed, e.g. `Ki67 cut-off (mean) in n1: lod 36.7 -> manual 120`;
- a generation built by an older kerncount (no `inputs`): left as it is, with a note and the rebuild command.

A generation is **built in `generations/gN.tmp/` and swapped in by renaming** (`fileio.swap_in`): a build that
fails leaves the old `gN` untouched (and no `.tmp`), and on Windows a file of `gN` open in Excel gives the error
"cannot replace gN: is a file in gN open in Excel?" with the new results waiting in `gN.tmp` (the old folder is
not touched). Files in the old `gN` that kerncount did not write (not in its manifest, e.g. your own notes or an
edited copy of a table) are never deleted: they are copied to **`gN/_kept/`** (path kept) and listed in the log;
earlier `_kept/` files stay.

## 7. Excel workbook (`results.xlsx`)
All previous sheets are kept (all-n sheets gain `metric` / `cell_line` columns). Added in this revision:
`Density-adjusted means`, `Density-adjusted effects`, `Duplicate images` (when any); the About sheet states the
model unit, estimands, classes, negative control, cell lines and duplicate images. Earlier sheets: `Model - means (absolute)`, `Model - ratios vs controls`,
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
- **Fractions**: read the percentage-point difference and the absolute % per condition; odds ratios and % of the
  reference are secondary (meaningless near 0 %). The density-adjusted estimate is only a check; when it is flagged
  ("conditions do not overlap in density") it extrapolates and should not be quoted.
- **One n**: estimates only; there is no across-n test (all-n CIs, p and stars say "needs ≥ 2 n").
- **Duplicate images** (`duplicates.csv`, acceptance `no_duplicate_images`): results for the affected conditions
  of that n are not interpretable until the real images are identified.

### Real data (4 n's, 8 conditions, 2 wells × 6 fields; g1, default settings, reference Control)
Nuclei per field vs LyoVec, mixed model (well unit, classes "auto"): miR-132-5p 0.29 [0.16, 0.55] *,
miR-30e-3p 0.33 [0.15, 0.70] ns (adj. p 0.058), miR-338-3p 0.41 [0.22, 0.77] ns (0.063), miR-340-3p 0.47
[0.31, 0.72] *, miR-375-3p 0.37 [0.21, 0.64] *, CL075 0.76 ns, Control 0.91 ns. The previous default (field unit,
two classes) gave 0.21–0.36 with ***/**** for every miRNA; ratio paired t-tests give 0.27–0.52, all ns after Holm
(LyoVec is missing in n2). n3's Control and LyoVec images are byte-identical (duplicates.csv), so n3 fails
`no_duplicate_images`.
