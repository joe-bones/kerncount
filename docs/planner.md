# Planning the next experiment

`kerncount sampling <experiment>` looks back ("could we have used fewer images?") and then plans
forward: how many **n's**, **wells** and **fields** the next experiment needs, how to split a plate
between treatments and references, and how likely more n's are to settle a comparison. It also
writes randomised **plate maps** (`kerncount platemap`).

Code: `src/kerncount/planner.py` (model, predictions, outputs) and `src/kerncount/platemap.py`.
Terms follow `CONTEXT.md`: one n = one plate on one day; field ⊂ well ⊂ condition ⊂ n.

## In plain language

1. Each image count varies for four reasons:
   - **field**: two fields of the same well differ;
   - **field position**: some spots in a well (e.g. near the rim) have more cells; the pattern is
     shared by all wells of one plate and changes from day to day;
   - **well**: two wells of the same condition differ;
   - **biological (n-to-n)**: the treatment effect itself differs between independent experiments
     (biology, transfection efficiency, passage). Only more n's reduce it.
2. The main result is the ratio treatment / reference per n, tested across n's with a paired t-test
   on log ratios (as in the analysis). Its noise per n is the biological part plus a measurement
   part that more wells and fields shrink.
3. With 3–4 n's every estimate of variation is uncertain, so every prediction comes as
   **pessimistic / expected / optimistic**. Plan with the pessimistic number when it matters.
4. The number of n's is usually the bottleneck: the t critical value is 12.7 at 2 n's, 4.30 at 3,
   3.18 at 4, 2.78 at 5, and extra images cannot remove the biological part.

## The model

All components are on the relative scale (value / condition mean). A CV of 30% is a log-scale SD of
about 0.3 (delta method), so they plug directly into the log ratio.

### Within an n: field, position, well

For every condition the values are divided by the condition mean, `r = y / mean`. Then

- within-well deviations `dev = r - well mean` give the **within-well mean square**;
- if `design.positions_meaningful(settings)` (`[design] field_positions = "fixed"`), the
  deviations are averaged per n and field position. Their sum of squares
  `SS_pos = Σ_n Σ_p m_np · mean(dev)²` (df `Σ_n (positions_n − 1)`) is the position part, and
  `E[MS_pos] = σ²_field + m0 σ²_pos` (m0 = wells per position, unbalanced: `(N − Σm²/N)/(P − 1)`);
- the rest is the **field (residual)** mean square, `σ²_field = MS_e`;
- between wells of a condition: `E[MS_well] = σ²_field + k0 σ²_well`,
  `k0 = (N − Σ n_j² / N)/(a − 1)` pooled over conditions (handles unbalanced wells).

Field position cancels in a treatment / reference ratio when every well is imaged at the same
positions (**standardised positions**): it is a plate-wide pattern. With random positions it adds to
field noise. So the field variance used for predictions is

    σ²_f = σ²_field                (standardised positions, the default when positions are "fixed")
    σ²_f = σ²_field + σ²_pos       (random positions)

Single wells: `σ²_well` is not estimable (reported as such); well variation is then inside the
biological part, so predictions for more wells are conservative. One field per well: field and well
cannot be separated; all of it is counted as well.

### Across n's: biological variance, pooled within classes

For reference `r` and a **class** of treatments, `d_nt = log(mean_t / mean_r)` per n. Classes come
from `[analysis] classes`, with the same rule as the statistics (`models.config`): the default `"auto"` = one
class per non-control treatment, references (`reference` + `normalise_to`) pooled as `controls`; a table
(class → names or glob patterns, e.g. `mirna = ["miR-*"]`) = those classes; `{}` = references vs everything
else (`treatments`). With `"auto"`, generic planning effects (`planning 0.7`, the cost analysis) use the
`default_class` = the treatment whose biological variance vs the main reference is the (upper) median, noted in
`notes`.

`d_nt = δ_t + u_n + ε_nt`: `u_n` is shared by all treatments of one n (the reference's own noise and
day effects on the response, e.g. n3 responding more strongly), `ε` is treatment-specific. The
per-n variance of a ratio is `V = Var(u) + Var(ε)`; it is estimated without bias by the pooled
within-treatment variance, and its uncertainty from the two mean squares of an additive
n + treatment fit (missing cells allowed):

    V = MS_n / t̄ + (1 − 1/t̄) MS_e      t̄ = cells per n,  df_n = N − 1,  df_e = cells − N − T + 1
    df(V) = V² / ((MS_n/t̄)²/df_n + ((1 − 1/t̄) MS_e)²/df_e)        (Satterthwaite)

Pooling treatments that share a reference this way is honest about their correlation: the shared
part has only N − 1 df. A class with one treatment uses its own variance (df = pairs − 1).
The biological variance is what is left after removing the measurement part of the actual design:

    σ²_bio = V − Σ (a_well σ²_well + a_field σ²_field)      a_well = Σ n_j² / N²,  a_field = 1/N per condition
                                                           (treatment and reference, averaged over cells)

### Uncertainty: generalized pivotal quantities

Each mean square has a pivot `MS · df / χ²_df`. 1000 joint draws of the pivots give draws of every
component (negative values set to 0). The 2.5–97.5% quantiles are the **95% CI** of a component;
the 10–90% quantiles of the predicted variance of a design give the **band**
(optimistic / pessimistic). A component estimated as 0 still has a positive upper limit, so the
pessimistic predictions stay honest when 3–4 n's happen to look similar.

### Predictions for a design

N n's, w wells × k fields per treatment, w_r wells per reference:

    V(w, k, w_r) = σ²_bio + σ²_well (1/w + 1/w_r) + σ²_f (1/(w k) + 1/(w_r k))
    power        = P(|T| > t_{1−α/2, N−1}),  T ~ noncentral t(df = N − 1, nc = |log effect| / sqrt(V/N))
    CI half-width (95%, on the ratio) = exp(t_{0.975, N−1} sqrt(V/N)) − 1      reported as ±x %
    MDE (ratio, target power)        = exp(−nc* sqrt(V/N)),  nc* solves power(nc*, N−1) = target

The ±x% CI means: upper limit = ratio × (1 + x/100), lower limit = ratio / (1 + x/100).
The paired test has df = N − 1 because each n contributes one log ratio; n's missing the reference
are left out of that reference's comparisons (as in the analysis).

### Cost / effort

`[planner] cost_per_n`, `cost_per_well`, `cost_per_field` in any unit (default minutes: 240 per n for
seeding, transfection and staining a plate; 5 per well; 1 per image). Total effort of a series:

    cost = N × (cost_per_n + (treatments × w + references × w_r) × (cost_per_well + k × cost_per_field))

- classic optimal fields per well: `k* = sqrt((cost_per_well / cost_per_field) × (σ²_f / σ²_well))`;
- wells per n: `w* = sqrt(B cost_per_n / (A D))` with `A = σ²_bio`, `B = 2(σ²_well + σ²_f/k*)`,
  `D = conditions × (cost_per_well + k* cost_per_field)` (minimises variance × cost);
- the Pareto frontier of (expected) power vs total effort, and the cheapest design reaching
  `target_power`, on the expected and on the pessimistic variance.

### Reference wells

Dunnett's rule: with k comparisons per reference, give the reference `√k` times the wells of a
treatment. The honest gain is smaller than the rule suggests, because only the reference's
measurement part shrinks (`σ²_bio` does not). `allocation.csv` lists every split
(treatment wells, reference wells) that fits `plate_format`, with its variance (averaged over all
treatment × reference comparisons), the n's needed for the focus effect, whether it fits the inner
wells, and the recommended split (lowest variance).

### How many more n's?

For every treatment × reference: draw σ² from a scaled-inverse-χ²(df(V), V) (the class-pooled
variance at the current design), the true mean log ratio from Normal(observed mean, σ²/n_now),
simulate k = 1..8 new n's at the current design, append them to the observed log ratios and run the
final paired t-test. The share of simulations with p < α is the predictive probability.

**Optional stopping.** Adding n's one at a time and stopping as soon as p < 0.05 inflates false
positives (the report shows by how much, from a null simulation: ~17% instead of 5% when testing after
every n from 3 to 11). Use these numbers to *choose* a final number of n's in advance, run them all,
then test once. They are planning information, not a stopping rule. The posterior is centred on the
observed mean of few n's, so an effect seen by chance can look promising: treat high probabilities for
small observed effects with caution.

## Outputs (`sampling/gN/`, next to the sampling files)

| file | content |
|---|---|
| `planner.json` | everything the browser planner needs (below) |
| `planner_components.csv` | field, position, well, biological variances with df, 95% and 80% intervals |
| `design_forward.csv` | per reference, class, effect (planning and observed) and design (n 2..`max_n`, wells 1..6, fields 1..max(3 × current, 12)): power and CI half-width low/mid/high, SD, cost |
| `mde.csv` | minimum detectable ratio (low/mid/high) and CI half-width per design |
| `cost_frontier.csv` | frontier designs per planning effect; `cheapest_reaching_target` |
| `allocation.csv` | treatment / reference well splits (see above) |
| `more_n.csv` | predicted probability of p < α after k more n's |
| `variance_components.png`, `power_vs_n.png`, `ci_width_vs_fields.png`, `cost_frontier.png`, `more_n.png` | plots (+ PDF) |
| `manifest.json` | `{"files": [{path, kind, title, description, pdf?}]}` for every file in the folder |
| `REPORT.txt` | now starts with **Recommendations for the next experiment** |

"low/mid/high" are in numeric order: `power_low` is the pessimistic power, `ci_halfwidth_pct_high`
the pessimistic width, `mde_ratio_low` the pessimistic (smallest) detectable ratio.

## Plate maps

    kerncount platemap <experiment> --plates 3 --format 24 [--wells 2] [--reference-wells 4|auto]
                       [--avoid-edge] [--seed 1] [--no-balance]

Conditions are the treatments of the project (counts.csv, else layout.csv) in `treatment_order`.
References (`reference` + `normalise_to`) get `--reference-wells` (default `[planner]
reference_wells`; `auto` = Dunnett, reduced until the plate fits). Writes
`analysis/<experiment>/plans/platemap_<timestamp>/`: `plate_K.csv/.png/.pdf` (A4 landscape),
`all_plates.csv`, `edge_usage.csv`, `README.txt` (with the command that recreates the same maps).

- Every plate is randomised separately; the seed is printed and stored.
- `--avoid-edge` (default from settings): inner wells first; edge wells only when needed, then shared
  out evenly.
- Balance (default on): each plate is the best of 400 random layouts for a score that rotates every
  condition through edge / inner wells, rows and columns across plates and keeps replicate wells
  apart.
- Unused wells stay empty (medium only). Too many wells for the plate is a clear error.

Python: `platemap.generate(conditions, plate_format, wells_per_condition, n_plates, avoid_edge, seed,
balance) -> list[DataFrame]` (columns `plate, plate_well, row, col, condition, replicate, edge`; empty
wells have condition `""`), `platemap.render(df, path)`, `platemap.edge_usage(plates)`,
`platemap.plan_project(project_dir, ...)`.

## Settings (`[planner]`)

| key | default | meaning |
|---|---|---|
| `effects` | `[0.8, 0.7, 0.5]` | ratios to plan for (0.7 = 30% fewer nuclei) |
| `target_power` | `0.8` | power that defines "enough" and the MDE |
| `alpha` | `0.05` | test level |
| `cost_per_n` | `240.0` | effort per n (plate / day), any unit |
| `cost_per_well` | `5.0` | effort per well |
| `cost_per_field` | `1.0` | effort per image |
| `max_n` | `12` | largest n in the tables |
| `plate_format` | `24` | 6, 12, 24, 48 or 96 |
| `avoid_edge` | `true` | plate maps: inner wells first |
| `reference_wells` | `"auto"` | wells per reference; `"auto"` = Dunnett √k |

`[analysis] classes` (owned by the statistics settings) is read if present.

## API for the browser planner

Pure functions, no file writes. `evaluate` takes well under 5 ms (tested < 50 ms).

```python
from kerncount import planner

plan = planner.load(project_dir, generation=None)       # reads sampling/<gN>/planner.json (latest if None)

r = planner.evaluate(plan, n=5, wells=2, fields=6, effect=0.7,   # or effect="miR-132-5p" (observed ratio)
                     reference_wells=None,               # None = wells; "auto" = Dunnett; or an int
                     alpha=0.05, target_power=0.8,
                     reference=None,                     # default: plan["design"]["main_reference"]
                     cls=None,                           # default: plan["default_class"]
                     positions=None)                     # "standardised" | "random"; default from [design]
# r = {n, wells, fields, reference_wells, effect_ratio, effect_label, reference, cls, positions,
#      power: {low, mid, high}, ci_halfwidth_pct: {low, mid, high}, mde_ratio: {low, mid, high},
#      sd_log_ratio: {low, mid, high}, measurement_share, cost, images, wells_per_n, t_critical,
#      notes: [plain-language strings]}

pts = planner.curve(plan, vary="n", fixed={"wells": 2, "fields": 6}, effect=0.7, values=None, **evaluate_kwargs)
# vary: "n" | "wells" | "fields" | "reference_wells"; returns [{"x": value, **evaluate(...)}, ...]
```

Useful parts of `planner.json`:

| key | content |
|---|---|
| `design` | `n_count, ns, conditions, condition_names, references, main_reference, wells, fields_per_well, images_per_condition, positions_meaningful, positions_default, unbalanced, single_well, missing_reference {ref: [n]}` |
| `components` | `field / position / well: {var, sd, ci95, ci80, df, note}` (`var` null = not estimable) |
| `bio` | `"<reference>|<class>": {var, ci95, ci80, V, V_ci95, df, kind, treatments, meas_actual, ...}` |
| `effects` | per reference × treatment: `ratio, log_ratio, n_pairs, n_list, values, sd_own, p, ci95_ratio, cls` |
| `classes`, `class_source`, `default_class` | class membership |
| `options`, `costs` | `[planner]` settings, cost unit |
| `summaries` | `cost` (frontier cheapest designs, optimum k*, w*), `allocation`, `more_n` (incl. optional-stopping false positive rate), `positions`, `missing_reference`, `focus_effect` |
| `draws` | 1000 joint GPQ draws of `field, position, well` and `bio["<ref>|<class>"]` (used for bands) |
| `notes` | plain-language caveats |

## Assumptions and limits

- Log-normal-ish counts, delta method (relative variance ≈ log-scale variance); fine for CVs up to
  ~70% per field because predictions concern means of many fields.
- Classes share one biological variance; a treatment with a much more variable response than its
  class is under-predicted. Per-treatment raw values are in `effects.csv` / `planner.json`.
- Planning effects are taken as true effects; observed effects from few n's are usually too large.
- The position pattern is assumed shared by all wells of a plate (true on the current data);
  standardised positions assume the same spots in treatment and reference wells.
- Future n's are assumed to have every reference, the same design, and the same biology.
- Plate positions (edge effects) are not part of the variance model; the plate maps randomise and
  balance them instead.
