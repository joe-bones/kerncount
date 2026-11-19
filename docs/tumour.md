# Tumour volume experiments (`[experiment] type = "tumour_volume"`)

mCherry tumours in brain sections, imaged as z-stacks. One **tumour** (hemisphere A or B of a section) is one n.
Vocabulary: `CONTEXT.md` (Tumour experiments). Decisions: ADR 0002 (measurement), 0003 (statistics), 0004 (types,
factor tables). Status: layout, primary segmentation, main-mass volume (ticket 02) and the pre-registered factorial
statistics with group, forest and interaction plots (ticket 04). Shape readouts, QC flags, Fiji legacy readout and
sensitivity models follow in tickets 05–08.

## Quick start
```
kerncount init tumour_TLR7 --type tumour_volume   # scans data/tumour_TLR7/, writes analysis/tumour_TLR7/
kerncount run tumour_TLR7
```

## Raw data and layout.csv
`data/<experiment>/` holds one sub-folder per **imaging session**, each with flat files as delivered (never modified):

| file | meaning |
|---|---|
| `030826 HERV -clodro n3 A.ome.tiff` (`.ome.tif`) | original z-stack = one tumour |
| `030826 HERV -clodro n3 Af.tif` | Fiji legacy binary mask (recorded in `legacy_mask`, not a tumour) |
| `Statistics for 030826 HERV -clodro n3 Af.csv` / `.txt` | Fiji 3D Objects Counter table (`legacy_stats`) |

`layout.csv` gets one row per original: `tumour_id`, `imaging_session`, `date` (from the 6-digit token,
`[project] date_order`), one column per declared factor, `section` (`n3`), `hemisphere` (`A`), `tumour` (`n3A`),
`file` (relative to the raw folder), `legacy_mask`, `legacy_stats`, `include`, `note`.
`tumour_id` = every factor level + section + hemisphere (`TLR7 KO / HERV / -clodro / n3A`), unique across the
experiment because `n3A` exists in several groups and both genotypes. A name the parser cannot read is set
`include = no` with a note (full validation: ticket 03).

### Factors (`[[tumour.factors]]`)
Each factor declares `name`, `levels`, `reference` (comparison baseline), `role` (where it goes in group plots:
`row`, `panel`, `x`, `colour`), `source` (`filename` token or imaging-`session` folder name) and `aliases` (other
spellings mapped onto a level, e.g. `KO7 = "TLR7 KO"`, the typo `-clodron = "-clodro"`; a filename alias is noted).
Defaults: genotype (WT ref, TLR7 KO; session; row), tumour_line (cntrl ref, HERV, N4; filename; x), clodronate
(-clodro ref, +clodro; filename; panel). Nothing downstream assumes three factors: adding a factor is a settings entry.

## Measurement (primary method, ADR 0002)
Per stack, identical settings for every tumour:
1. **Voxel size** from the file's OME metadata (`PhysicalSizeX/Y/Z` + units). The z-step differs between stacks
   (1.75 / 2.0 / 2.36 µm), so it is read per file. Missing calibration stops the run — never a pixel unit.
2. **Smoothing**: Gaussian, `smooth_sigma_um` (2 µm) in physical units (so per-axis in voxels).
3. **Threshold** per tumour, **no rolling-ball subtraction**: `intermedians` (below).
4. **3D hole filling**: background pockets not connected (6-connected) to the stack border become tumour.
5. **Connected components** (26-connected); objects < `min_object_um3` (1000 µm³) are noise.
   **Main mass** = largest object (label 1); the others are **satellites** (labels 2.., largest first).
6. Readouts: main-mass volume (µm³, mm³), threshold, background and tumour levels, number of satellites, saturated
   fraction of the stack and of the main mass (8-bit value 255).

### Threshold choice: `intermedians`
Iterative intermeans (Ridler–Calvard / IsoData) with **class medians**, started at Otsu:
`t = (median of voxels ≤ t + median of voxels > t) / 2`, repeated until stable, on a 1024-bin histogram of the
smoothed stack. Why:

- **Uniform brightness cancels.** t is proportional to the signal: a tumour recorded ×0.4 dimmer or ×2.5 brighter
  (gain, expression) gets the same mask up to 8-bit rounding and clipping. A fixed absolute threshold would turn a
  group-level brightness difference (seen in the data: N4 +clodro dim, N4 −clodro bright) into a volume difference.
- **Half-way between tissue background and tumour signal** is the half-maximum edge of a blurred boundary, i.e.
  the true edge, whatever share of the stack the tumour fills. Otsu weights its classes by size and drifts toward
  the background for small tumours in a large field.
- **Saturation-tolerant**: the tumour median ignores the spike at 255 while less than half of the tumour clips.
  Triangle (the Fiji macro's method) fails both ways in 8-bit data: on unsaturated stacks it lands on the
  background floor, and when the 255 spike becomes the histogram peak it jumps to 254.

Synthetic evidence (blurred ellipsoids, noise, 2.36 × 1.25 × 1.25 µm voxels; `tests/test_tumour.py`):

| | small (15×30×25 µm) | medium (30×70×60) | large (40×140×120) | ×0.4 / ×2.5 brightness (×2.5: 60 % of main mass at 255) |
|---|---|---|---|---|
| intermedians | +1.4 % | −0.9 % | −0.3 % | +0.1 % / +1.0 % |
| Otsu | +5.8 % | +1.2 % | +0.3 % | +0.0 % / +0.9 % |
| Triangle (exploration only) | +127 % | +53 % | +34 %, −20 % when saturated | — |
| Li (exploration only) | +44 % | +18 % | +11 % | — |

Tests: volume within **5 %** of the analytic volume across sizes and z-steps (1.75 vs 2.36 µm); ×0.4 and ×2.5
(clipped) change the volume by **< 3 %**; satellites 8–10 µm from the surface stay separate; an enclosed dark core is
filled. `threshold = "otsu"` remains selectable for comparison.

### What the real data show (one stack each, 2026-11-19)
| stack | main mass | Fiji largest object | ratio | runtime / peak RSS |
|---|---|---|---|---|
| cntrl +clodro n1 B (43 × 1024², z 2 µm) | 5.31 × 10⁶ µm³ | 1.11 × 10⁷ µm³ | 0.48 | 3 s / 0.8 GB |
| HERV +clodro n3 A (74 × 1946², z 1.75 µm, largest) | 1.72 × 10⁷ µm³ | 5.32 × 10⁷ µm³ | 0.32 | 14 s / 2.8 GB |

The xy outlines agree closely with Fiji's; the difference is mostly **z-extent**. The stacks have a strong axial
intensity profile (layer means ~3 → 70 → 4 grey levels): Fiji's rolling ball + Triangle keeps dim layers far above
and below the bright core, the half-way threshold does not (n1 B: substantial main-mass area in layers 10–32, Fiji 4–34). 2D per-layer hole
filling (the Fiji macro) would add ~8 % on n1 B. How much of the dim top/bottom layers is tumour versus out-of-focus
light is the open question for the threshold sweep and Fiji concordance (ticket 07) and truncation flags (ticket 06).

## Outputs
- `layout.csv` (init), `tumours.csv` — one row per included tumour: factor columns, imaging session, IDs, voxel
  size, `n_layers`, main-mass volume (µm³, mm³), threshold method / value / levels, satellites, saturation fractions.
- Statistics outputs (see **Statistics** below), per readout `<name>` (e.g. `main_mass_volume`):
  - `plots/<name>_absolute.png|pdf` — one row of panels per genotype present, −/+clodro panels, tumour line on x,
    every tumour a dot, model geometric mean with 95 % CI, n under each group, log y.
  - `plots/<name>_normalised.png|pdf` — the same as % of the model-estimated cntrl in the same genotype and
    clodronate level (dashed line = 100 %); bars = the F1 fold-change CI × 100, stars from F1 (Holm).
  - `plots/<name>_forest.png|pdf` — every contrast grouped by family, fold-change and 95 % CI on a log axis, with
    stars / labels (exploratory, descriptive, confounded, needs … data).
  - `plots/<name>_interaction.png|pdf` — model means (95 % CI): line on x, one line per clodronate level, one plot
    per genotype.
  - `prism/<name>_prism.csv` and `prism/<name>_pct_prism.csv` — Prism column tables, one column per group
    (`TLR7 KO / HERV / -clodro`), one tumour per row (absolute / normalised %).
- `stats_long.csv` — every estimated mean and contrast in kerncount's long stats format (`docs/statistics.md`):
  `level = all_n`, `metric` = readout name, `family` = F1…F7 / `estimated means`, `estimate_type = ratio`, Holm `p_adj`,
  `significance`, `note` (family role, labels, `lock <hash>`).
- `results.xlsx` — sheets About (analysis plan, lock hash, levels without data), Tumours (tumours.csv + `<name>_pct`
  columns), Statistics (= stats_long), Estimated means (absolute and % of reference with CIs), Families.
- `stats_report.md` — the locked plan (readout split, families) and one table per readout and family.
- Sensitivity (ticket 08): `plots/<name>_residuals.png|pdf` (residuals by group with each group's SD + normal QQ
  of the studentised residuals), sheet **Sensitivity** in `results.xlsx` (one row per contrast: Welch CI/df/p/Holm p,
  permutation p/Holm p, permutations, seed, `concordant`, `concordance`), stats_long rows with family
  `<id> (Welch-type)` / `<id> (permutation)` (never starred), and a `## Sensitivity` section in `stats_report.md`.
- Planning aid (ticket 09; confirmatory readouts only, never in stats_long): sheets **Planning aid** (one row per
  readout and family) and **Planning aid (contrasts)** in `results.xlsx`, section `## Planning aid` in the report.
- `_work/tumour/<slug>.json` + `<slug>.labels.tif` — cache per tumour (zlib uint8/uint16 labels, ~1 MB for the
  largest stack). Re-runs reuse it unless the file or the segmentation settings changed; `--resegment` forces it.

## Statistics (ADR 0003 + amendment)
**Model.** For each readout: one linear model on log(readout) with one mean per **cell** (combination of all
declared factor levels = full N-way interaction) and one residual variance pooled across all cells, fitted with
whatever cells are present (one genotype today). Tumours are independent n's. A contrast is a set of weights on cell
means; `exp` gives a fold-change (ratio of geometric means) or a ratio of fold-changes, with a 95 % t-CI on the pooled
residual df (tumours − cells). With every group n = 1 there is no residual df: estimates only, nothing tested.

**Readouts** (`[[tumour.readouts]]`: `name`, `column` in tumours.csv, `label`, `role` = `confirmatory` /
`exploratory`). Adding a readout (ticket 05) = add the column to the table and a settings entry; nothing in the
statistics knows readout names. Declared readouts missing from the table are reported, not fatal.

**Families** (`[[tumour.families]]`: `id`, `effect`, `factors`, `role`, `label`) are generic over factor names:
- `effect = "simple"`, one factor X: each non-reference level of X vs the reference, within every combination of the
  other factors (`HERV vs cntrl | TLR7 KO, -clodro`);
- `effect = "interaction"`, factors X1…Xk: the product of (non-reference − reference) over the listed factors, within
  every combination of the rest. Two factors = ratio of ratios, the last factor the inner ratio:
  `(HERV / cntrl in TLR7 KO) ÷ (HERV / cntrl in WT) | -clodro`.

| default | effect | role |
|---|---|---|
| F1 | simple, tumour_line (within genotype × clodronate) | confirmatory |
| F2 | simple, clodronate (within genotype × line) | confirmatory |
| F3 | interaction tumour_line × clodronate (within genotype) | exploratory |
| F4 | interaction genotype × tumour_line (within clodronate) | confirmatory |
| F5 | interaction genotype × clodronate (within line) | exploratory |
| F6 | simple, genotype (within line × clodronate) | descriptive |
| F7 | interaction genotype × line × clodronate | exploratory |

**Rules.** Holm within each family and readout, over the contrasts that are tested. Stars only for confirmatory
families of confirmatory readouts; exploratory families/readouts get p and Holm p but no stars and are labelled;
descriptive families show estimate + CI only. **Imaging session:** a contrast whose weights do not cancel within every
imaging session (a multiplicative day shift would move it) is labelled "confounded with imaging session" and never
tested — direct WT vs KO (F6) is; within-session ratios and ratios of them (F1–F5, F7) are not. **Missing levels:**
contrasts needing cells without tumours stay in every table, labelled `needs WT data` (whole level absent) or
`needs data for <group>`. A family naming an undeclared factor is skipped with a note.

**Normalised values** (`[tumour] normalise_by`, default the factor on x): % of the model-estimated reference level
(cntrl) in the same cell of all other factors. Reference tumours scatter around 100 % (their geometric mean is exactly
100 %); a group's normalised estimate and CI equal its F1 fold-change × 100.

**Sensitivity models** (ADR 0003; `tumour/sensitivity.py`). Every contrast the primary model tests is re-tested
twice; contrasts it does not test (confounded with imaging session, descriptive, needing absent groups) stay untested.
- *Welch-type*: same point estimate, SE² = Σ w²·s²ᵢ/nᵢ with each group's own variance, Satterthwaite df. Needs ≥ 2
  tumours in every group the contrast uses.
- *Permutation*, simple effects (two groups in one cell of the other factors): the focus-factor labels of the two
  compared groups are permuted. Other levels of the focus factor are left out: under the null of *this* contrast only
  these two groups are exchangeable. Every split is enumerated when there are at most `sensitivity_permutations`
  (6 + 6 tumours = 924 splits, smallest p 2/924 = 0.0022, still < 0.05/8 for Holm over the 8 contrasts of F1);
  otherwise random permutations, p = (1 + #extreme)/(B + 1) (B = 9999: smallest p 0.0001).
- *Permutation*, interactions: **Freedman–Lane** — fit the reduced model (cell means with this contrast constrained
  to 0: μᵢ − λ·wᵢ/nᵢ), permute its residuals across the tumours of the contrast's cells, add them back to the reduced
  fit, refit, compare |t|. Unlike raw label permutation this keeps the main effects in the null, so it tests only the
  ratio of ratios. For two groups it coincides with label permutation. The t of the permutation tests pools the
  variance over the contrast's own groups (not all groups), so it does not lean on equal variances elsewhere.
- Random permutations are seeded: `[tumour] sensitivity_seed` (default 20261119) combined with the contrast label,
  so reruns are identical and independent of order. Seed and permutation count are in every output.
- Holm within each family and readout, per method. **Concordant** = both methods reach the same Holm-adjusted
  conclusion at 0.05 as the primary model (direction always agrees: the point estimate is shared). A `no` means the
  conclusion depends on the pooled-variance assumption or on the t approximation — report it alongside.

**Planning aid** (power planner; ADR 0003: never feeds confirmatory results; `tumour/planner.py`). From the pooled
residual SD of the current data (log scale), per confirmatory readout and every family that is tested (descriptive
families and contrasts confounded with imaging session are not planned):
- **MDFC** — minimum detectable fold-change (up or down) of each contrast at the current group sizes with 80 % power
  (`planner_power`), two-sided α = 0.05 per family (`planner_alpha`): SE = s·√(Σ w²/n), noncentral t on the current
  residual df. The family row shows best – worst over its contrasts.
- **Multiplicity:** Holm is approximated by **Bonferroni within family** (α / m, m = contrasts the family tests in the
  full declared design, also those still needing data). Holm's first step is Bonferroni and later steps are less
  strict, so for any one contrast Holm power ≥ Bonferroni power: the planner is conservative.
- **Tumours per group** for `planner_fold_change` (default 1.5): *balanced* = every group of the design at k
  (df = groups × (k − 1)); *top-up* = raise only the hardest contrast's groups below k, others unchanged. The
  **limiting group** is the contrast's group with the largest share of its variance (the smallest; ties listed).
- **Absent levels:** a contrast needing groups without tumours (F4 genotype × line with only KO imaged) has no MDFC
  yet; its tumours-per-group are planned from the present data's variance with the note "assumes WT variance equals
  that of TLR7 KO" (wording built from the factor levels).
- Example (synthetic, 48 KO tumours in groups of 12/6 like the KO set, pooled SD 0.37 on log scale): F1 MDFC
  1.75–2.21, F2 1.95–2.16; for 1.5× F1 needs 22 and F2 21 tumours per group; F4 needs 37–38 per group
  in each genotype. Validated in `tests/test_tumour_planner.py` by simulation (power at the MDFC and at the planned
  n ≈ 80 %, Bonferroni within family) and against the normal approximation at large df.

**Lock.** The analysis plan (factors + levels + references, readouts + roles, families, normalise_by) is hashed
(`lock <12 hex>`); the hash is in every stats_long note, the workbook's About sheet and `stats_report.md`. Change the
plan only before results are viewed; a different hash means a different plan.

## Memory
One stack at a time: uint8 stack + float32 smoothed copy, then mask + int32 labels; peak ≈ 10 × the stack size
(2.8 GB for the largest, 280 MB stack). Hole filling labels the background instead of iterative dilation.
