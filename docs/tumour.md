# Tumour volume experiments (`[experiment] type = "tumour_volume"`)

mCherry tumours in brain sections, imaged as z-stacks. One **tumour** (hemisphere A or B of a section) is one n.
Vocabulary: `CONTEXT.md` (Tumour experiments). Decisions: ADR 0002 (measurement), 0003 (statistics), 0004 (types,
factor tables). Confirmatory readouts: main-mass **cross-sectional area and perimeter over a fixed-thickness
mid-section** (ADR 0002/0003 amendments 2026-11-20; see Readouts) — the 3D volume depends on the acquired depth.

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

`layout.csv` gets one row per original: `tumour_id`, `imaging_session` (the folder), `session_label` (the folder's
words that are no factor level, e.g. `5d`), `file_date` (6-digit token of the file name, `[project] date_order`;
optional, may differ between folders and files, **not** the imaging date), `acquisition_day` (from the original's OME
metadata; header only), one column per declared factor, `section` (`n3`), `hemisphere` (`A`), `tumour` (`n3A`),
`file` (relative to the raw folder), `legacy_mask`, `legacy_stats`, `include`, `note`.
`tumour_id` = every factor level + section + hemisphere (`TLR7 KO / HERV / -clodro / n3A`), unique across the
experiment because `n3A` exists in several groups and both genotypes. A name the parser cannot read is set
`include = no` with a note (full validation: ticket 03).

### Session folders (ticket 12)
- **Any number of session folders.** Every sub-folder of the experiment folder is an imaging session; hidden / system
  folders (`.…`, `_…`, `@…`, `$…`, `#…`, `~…`, `System Volume Information`) are ignored; flat files directly in the
  experiment folder are one session named after it.
- **Factor levels from the folder name**, word by word and case-insensitive, through levels and `aliases` (multi-word
  levels such as `TLR7 KO` are matched as a phrase): `5d K07`, `ready 5d KO7 030826`, `TLR7-KO 10d`, `5d TLR7KO` →
  TLR7 KO; `5d WT` → WT. Words that are no level stay in `session_label`; a `timepoint` factor with
  `source = "session"` and levels `["5d", "10d"]` would read them. **A new genotype** (e.g. a second KO line) is a new
  level plus its folder spelling as an alias — no code change. A folder whose genotype cannot be read blocks the run
  with a message naming the folder, the words not read and both fixes (alias in settings.toml, or the genotype
  column of its rows in layout.csv).
- **Repeat sessions** of one genotype (`5d WT`, `5d WT repeat`): when the same tumour ID comes from several session
  folders, the folder is appended to those IDs only (`WT / HERV / -clodro / n1A / 5d WT repeat`); other IDs stay
  plain. Two files giving one ID inside one folder still block the run, naming both files. (Adding a folder that
  creates a new clash changes the IDs of the clashing tumours; their caches and review codes are then new.)
- **Renaming or reorganising folders.** `layout.csv` rows are matched to the files by path; a row whose file is gone
  is matched to the one scanned file of the same name. Folder-derived cells (paths, session, session label, genotype
  from the folder) come from the new place; values you typed (cells differing from what the old place reads as),
  `include` and your part of `note` are carried over. `kerncount run` does this in memory (a note lists every moved
  file); `kerncount init` writes the updated layout.csv after backing the old one up to `_work/backups/` and lists
  every moved / added file. When nothing moved and no file or column was added, init keeps layout.csv and writes the
  fresh scan to `layout.new.csv` as before. Older layouts get the new columns (`date` becomes `file_date`).
- **Caches survive a rename**: every cache key (segmentation, QC, Fiji, sweep, shared threshold) uses the file
  name, size and modification time (not the folder) and the cache file is named after the tumour ID, which does not
  contain the folder; review codes are keyed by tumour ID. A renamed folder therefore re-segments nothing.

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

## Readouts (ticket 05; ADR 0002 / 0003 amendments 2026-11-20)
**Why area, not volume.** Every tumour is ~0.5–1 mm wide, but only ~80–150 µm of it is imaged. The lab's Fiji masks
were thresholded per layer, so Fiji's main object spans every layer of every stack and its volume ≈ footprint ×
acquired depth (a microscope setting). The primary method's z-extent in turn follows the axial brightness fall-off
(brighter tumours look deeper). Neither is a tumour volume. What *is* measured reliably is the tumour's cross-section
in the layers where it is fully in focus, so the **confirmatory readouts are 2D**:

| readout (column) | definition | role |
|---|---|---|
| `cross_sectional_area_mm2` (`…_um2`) | mean main-mass area per layer over the **mid-section** | **confirmatory** |
| `perimeter_mm` (`…_um`) | mean main-mass outline length per layer over the mid-section, at `outline_smooth_um` | **confirmatory** |
| `footprint_area_mm2` | area of the main mass's outline merged over all layers (projection) | exploratory |
| `slab_volume_mm3` | main-mass volume within ± `slab_um`/2 of the mid-section centre (`slab_fits`) | exploratory |
| `main_mass_volume_mm3` | 3D main-mass volume (depends on acquired depth / brightness fall-off) | exploratory |
| `main_mass_surface_mm2` | 3D mesh surface: isotropic `surface_grid_um` grid, smoothed at `outline_smooth_um`, marching cubes; includes the faces cut by the first/last layer | exploratory |
| `tumour_burden_mm3` | main mass + satellites | exploratory |
| `elongation` | √(λ₁/λ₂) of the section's second moments (major / minor axis), mid-section mean | exploratory |
| `solidity` | outline area / convex-hull area (lobes, branching), mid-section mean | exploratory |
| `edge_roughness` | outline length at `roughness_fine_um` (2 µm) / at `roughness_smooth_um` (20 µm) (ragged border), mid-section mean | exploratory |
| `circularity` | 4πA/P² (2D sphericity; confounds elongation and roughness — reference only) | exploratory |
| `invasion_fraction`, `n_satellites` | satellites' share of the burden; number of satellites | columns + shape panel only (zeros: no log model) |

**Mid-section.** Per tumour, the per-layer main-mass area (2D holes filled) is computed; the centre is the layer with
the largest area (middle one on ties); the mid-section is every layer whose centre lies within ± `mid_section_um`/2
of it — **20 µm by default = 9 layers at 2.36 µm, 11 at 2.0 or 1.75 µm**. Columns `mid_section_centre_layer`,
`mid_section_first_layer`, `mid_section_last_layer` (1-based, as in Fiji), `mid_section_layers`, `mid_section_fits`
(`no` = it runs off the stack and was clipped; ticket 06 flags these as truncated). Why 20 µm (2026-11-20, 4 real
stacks, main-mass area per layer): the main mass is within 10 % of its largest cross-section over 9–19 µm of depth and
within 50 % over 28–40 µm; the brightest layer sits ~20–80 µm below the first layer in all 48 stacks (data survey).
20 µm therefore averages ≥ 9 layers (noise) while staying on the near-flat top of the area profile and fitting into
the stack; a thicker slab would mix in the brightness-dependent fall-off that makes volumes unreliable.

**Section and area.** The section of the main mass in a layer is all its pieces in that layer with enclosed holes
filled (a dark core opening to the top or bottom of the stack is still inside the outline). Area = pixel count × pixel
area (unbiased).

**Perimeter estimator.** Counting pixel edges overestimates every curved outline by 4/π (+27 %) at any resolution, and
marching squares on the raw binary still follows the staircase (+5.5 %). kerncount Gaussian-smooths the binary
section at `outline_smooth_um` and traces its 0.5 iso-line with marching squares (sub-pixel). Synthetic digitised
shapes at the real 0.6239 µm pixel (`tests/test_tumour_shape.py`; error vs 2πr / Ramanujan):

| shape | raw marching squares | σ = 2 µm | σ = 4 µm | **σ = 8 µm (default)** |
|---|---|---|---|---|
| disc r = 20 µm | +5.6 % | −0.26 % | −1.8 % | −8.4 % |
| disc r = 50 / 80 µm | +5.4 % | −0.02 / +0.01 % | −0.27 / −0.09 % | −1.2 / −0.47 % |
| disc r = 150 / 300 / 500 µm | +5.5 % | +0.02–0.04 % | 0.00 % | −0.12 / −0.03 / −0.01 % |
| ellipses 400×150, 300×120 (rotated), 250×200 µm | | +0.01–0.03 % | −0.02–0.00 % | −0.06 – −0.14 % |

The smoothing shrinks an outline by ≈ σ²/2r; real mid-sections are r ≈ 250–400 µm (0.19–0.5 mm²), where 8 µm costs
< 0.1 %.

**Why an 8 µm outline scale** (2026-11-20, 3 real stacks, `segment.sweep` at threshold × 0.8 / 1.2). Real tumour
borders are ragged at the µm scale, so a fine-scale perimeter measures border texture (circularity 0.09–0.40 at 2 µm)
and moves with the threshold; area does not:

| stack | area, threshold ×0.8 / ×1.2 | perimeter σ = 2 µm | σ = 4 µm | σ = 8 µm | perimeter / Fiji's at 2 / 4 / 8 µm |
|---|---|---|---|---|---|
| cntrl +clodro n1 B | +5.5 / −4.8 % | −16.8 / +16.4 % | −13.8 / +13.1 % | −9.8 / +5.1 % | 0.85 / 0.95 / 1.01 |
| HERV −clodro n2 A | +4.4 / −5.2 % | −13.8 / +10.6 % | −11.0 / +9.7 % | −2.0 / +4.5 % | 0.68 / 0.69 / 0.63 |
| N4 +clodro n3 A (dim) | +15.8 / −10.8 % | +13.6 / −10.9 % | +11.8 / −8.2 % | +6.5 / −6.5 % | 0.49 / 0.58 / 0.69 |

At 8 µm the confirmatory perimeter is as threshold-stable as the area and describes the outline, not its texture; the
texture is kept as the exploratory edge roughness (outline length at `roughness_fine_um` = 2 µm ÷ at
`roughness_smooth_um` = 20 µm). Tests: synthetic ellipsoid tumours at the default scale — mid-section area −0.6 % and
perimeter −0.5 % vs the analytic mid-section mean (tolerance 2 %); ×0.4 / ×2.5 brightness (clipped) +0.005 / +0.2 %
(tolerance 3 %); z-step 1.75 vs 2.36 µm +0.4 % / +0.2 % (tolerance 1 %); the threshold sweep recomputes both (they
are computed in `segment._segment_at`, not from saved labels).

**Fiji counterparts** (`legacy` on the readout entries; concordance plots and side-by-side effects, ticket 07):
`fiji_cross_sectional_area_um2|mm2`, `fiji_perimeter_um|mm`, `fiji_mid_section_centre_layer`,
`fiji_mid_section_layers`, `fiji_mid_section_fits`, `fiji_section_source` — the same mid-section rule and estimator
applied to Fiji's largest 3D object on the lab's mask (or the re-implementation's mask when there is none and
`fiji_fallback = true`). Volume → `fiji_main_volume_mm3`, surface → `fiji_main_surface_mm2`, burden →
`fiji_total_volume_mm3` (exploratory).

Real stacks (default settings, Fiji = the lab's masks; slow test `tests/test_tumour_shape.py`, 2026-11-20):

| stack | centre layer (layers) | area mm² (Fiji, ratio) | perimeter mm (Fiji, ratio) | circularity | footprint mm² |
|---|---|---|---|---|---|
| cntrl +clodro n1 B | 17 (11) | 0.187 (0.203, 0.92) | 2.05 (2.03, 1.01) | 0.67 | 0.201 |
| HERV −clodro n2 A | 29 (9) | 0.417 (0.506, 0.82) | 3.08 (4.88, 0.63) | 0.56 | 0.448 |
| N4 +clodro n3 A (dim) | 11 (9) | 0.191 (0.385, 0.50) | 2.57 (3.70, 0.69) | 0.37 | 0.230 |

Fiji's centre layer is within one layer of ours. Fiji's per-layer Triangle keeps faint haze around dim tumours, so
its sections are larger and more ragged the dimmer the tumour (N4 +clodro: twice the area) — the brightness
dependence ADR 0002 warns about, now visible in 2D as well. Runtime of the readouts: +1–4 s per stack.

**Shape panel** `plots/shape_panel.png|pdf` (exploratory, never starred): area vs perimeter of every tumour with the
circle reference P = 2√(πA) (above it = elongated, lobed or ragged outline), and elongation, solidity, edge roughness,
circularity, invasion fraction and satellites by group (dots, group median). Colour = tumour line, marker =
clodronate.

**Settings** (`[tumour]`, all in the measurement cache key): `mid_section_um = 20`, `outline_smooth_um = 8`, `roughness_fine_um = 2`,
`roughness_smooth_um = 20`, `slab_um = 40`, `surface_grid_um = 2`. The per-tumour cache (`_work/tumour/<slug>.json`)
also carries a readout version: caches written before these readouts existed are re-measured, never reused without
the new columns (likewise the sweep and Fiji caches).

## Outputs
- `layout.csv` (init), `tumours.csv` — one row per included tumour: factor columns, imaging session, IDs, voxel
  size, `n_layers`, the readouts above (confirmatory first), threshold method / value / levels, saturation fractions.
- `plots/shape_panel.png|pdf` — see Readouts.
- Statistics outputs (see **Statistics** below), per readout `<name>` (e.g. `cross_sectional_area`):
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
- `lock_log.csv` — unlocks and changes of the analysis plan after results existed, and runs with a changed plan
  (copied into the About lines of the report and `results.xlsx`; see "Web interface and the analysis-plan lock").
- Sensitivity (ticket 08): `plots/<name>_residuals.png|pdf` (residuals by group with each group's SD + normal QQ
  of the studentised residuals), sheet **Sensitivity** in `results.xlsx` (one row per contrast: Welch CI/df/p/Holm p,
  permutation p/Holm p, permutations, seed, `concordant`, `concordance`), stats_long rows with family
  `<id> (Welch-type)` / `<id> (permutation)` (never starred), and a `## Sensitivity` section in `stats_report.md`.
- Planning aid (ticket 09; confirmatory readouts only, never in stats_long): sheets **Planning aid** (one row per
  readout and family) and **Planning aid (contrasts)** in `results.xlsx`, section `## Planning aid` in the report.
- `acquisition.csv` (+ report section) and, with `shared_threshold_check = true`, `shared_threshold.csv` — see
  **Fiji legacy readout and method concordance** below.
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
descriptive families show estimate + CI only. **Imaging batch:** a contrast whose weights do not cancel within every
imaging batch (a multiplicative day shift would move it) is labelled "confounded with imaging session" (or "…with
acquisition day") and never tested — with genotype read from the session folder, direct WT vs KO (F6) is; within-
session ratios and ratios of them (F1–F5, F7) are not. What a batch is, is `[tumour] batch`: `"session_folder"`
(default: one folder = one batch) or `"acquisition_day"` (one OME acquisition day = one batch; groups imaged on
different days are then not compared by a test, e.g. HERV vs cntrl when HERV was imaged on another day). tumours.csv
has the per-tumour `batch` column; the report's About lines state the definition and the batches, and § Acquisition
settings always shows the group × acquisition-day table. The setting enters the analysis-plan lock only when changed
from the default. With ≥ 3 genotype levels every family compares each non-reference level with the reference (F4 =
each KO line's ratio of line effects vs WT), group plots get one row per present level, and declared levels without
data are listed as `needs <level> data`. **Missing levels:**
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

## Fiji legacy readout and method concordance (ticket 07, ADR 0002)
**Legacy readout** (`fiji_*` columns of `tumours.csv`; exploratory, labelled legacy, never starred):
`fiji_main_volume_um3|mm3` (largest 3D object), `fiji_main_surface_um2|mm2` (Fiji's border-voxel surface),
`fiji_n_objects` (objects ≥ `fiji_min_voxels`), `fiji_total_volume_mm3`, `fiji_main_z_start` (1 = first layer),
`fiji_main_z_layers`, `fiji_main_full_depth` (main object spans every layer), `fiji_source`, `fiji_note`.
Source per tumour, in order: the lab's 3D Objects Counter table (`legacy_stats`), objects counted on the lab's mask
(`legacy_mask`), the macro re-run on the original stack (`[tumour] fiji_fallback = true`). Statistics files are read
with or without header (`Statistics for …f.csv`, `…f_3Dresults.csv`, the header-less `.txt`), with trailing `;;`, a
truncated header or a leading row-number column; anything else (e.g. Excel sum rows with decimal commas in
`Statistics for 030826 N4 -clodro n3 Bf.csv`) is skipped and listed in `fiji_note` — never guessed. All 48 real files
parse; that one has 3 skipped rows. Cache: `_work/tumour/<slug>.fiji.json`.

**Re-implementation** (`tumour/fiji.py`), written from the ImageJ / 3D Objects Counter sources, not textbook versions:
rolling ball 50 px per layer exactly as ImageJ does it (3 × 3 mean pre-smoothing, 4× min-shrink, trimmed 19 × 19 ball
patch, grey opening with that non-flat element, ImageJ's linear re-enlargement; scipy, no
skimage `rolling_ball`), median radius 1 (= 3 × 3), Triangle (ImageJ's variant), mask, 2D fill holes, 3 × 3 open per
layer, 26-connected objects ≥ 200 000 voxels. Two findings from matching the lab's masks:
- **The masks were thresholded per layer**: Convert to Mask re-computed Triangle on every layer ("Calculate threshold
  for each image" ticked), keeping values ≥ the layer's level — not one stack threshold as the macro text suggests.
  With the stack threshold the agreement is only Dice 0.80 / volume 0.69 (cntrl +clodro n1 B). Per-layer
  thresholds are why Fiji's largest object starts in layer 1 and spans every layer in **all 48 tumours**: even a
  nearly empty top or bottom layer gets its own threshold near its background, so Fiji's volume is close to
  area × acquired stack depth.
- 3D Objects Counter's surface is the sum of exposed faces of the surface voxels, accumulated in float32 (≈1 % below
  the exact face sum at 10⁶ µm²); reproduced bit-exactly.

Measured agreement (slow test `tests/test_tumour_fiji.py`, 2026-11-19): cntrl +clodro n1 B and N4 −clodro n3 B
(the malformed-statistics one) — re-implementation vs lab mask **Dice 1.00000, volume ratio 1.00000** (voxel
identical); objects on the lab masks reproduce Fiji's voxels, surface voxels, surface and B-depth exactly. Runtime
10 s / 33 s per stack; peak ~2 × the stack plus the int32 labels.

**Concordance** (confirmatory readouts only; readout-generic: the Fiji counterpart is `legacy = "<column>"` on a
`[[tumour.readouts]]` entry, else the default map volume → `fiji_main_volume_mm3`, surface → `fiji_main_surface_mm2`):
- `plots/<name>_vs_fiji.png|pdf` — per tumour, primary (y) vs Fiji (x), log axes, identity line, median ratio and
  Spearman ρ in the title.
- `concordance.csv` / sheets **Concordance** and **Concordance summary** — the same families, contrasts and Holm
  correction re-run on the Fiji column; per contrast `direction_agree` (same side of 1) and `significance_agree`
  (Holm p < 0.05 in both or neither, same direction); per readout × family counts and `conclusion_differs`.
  `stats_report.md` gets a `## Method concordance` section that names every contrast whose conclusion differs.

**Threshold sweep** (`[tumour] threshold_sweep = true`, `sweep_scales = [0.8, 0.9, 1.1, 1.2]`): every stack is read and
smoothed once and thresholded at the primary threshold × each scale (`segment.sweep`, identical to
`segment(..., threshold_scale=s)`); only readouts are kept (`_work/tumour/<slug>.sweep.json`), the per-scale tables
go through `table.build`. Outputs: `threshold_sweep.csv` / sheet **Threshold sweep** (every contrast at every scale),
`plots/<name>_threshold_sweep.png|pdf` (one panel per family, fold-change vs scale, 95 % CI band, filled = Holm
p < 0.05) and a `## Threshold sweep` report table (fold-change range, whether significance changes). Off by default
because the first run re-segments every stack (~20–60 s each, 1 stack in memory). Only readouts produced by
`segment.segment` / `table.build` are sweepable; others are listed as not available.

**Acquisition settings** (ticket 11, always): `acquisition.csv` — one row per tumour from the original's OME metadata
(no pixels read): microscope, objective, detector, **detector gain = PMT voltage** (`detector_voltage`), offset, zoom,
laser line / power / attenuation, pinhole, channel / fluor / ex / em, pixel type, voxel size, frame size, layers,
acquisition day and time. `stats_report.md` § Acquisition settings gives the verdict (every brightness-relevant
setting identical, or which differ; settings never recorded are named as unchecked), what varies harmlessly
(z-step, frame size), and a **group × acquisition-day table** with the groups imaged on a single day. Note:
`imaging_session` in layout.csv is the data folder, not the acquisition day; a folder can span several days (the
`acquisition_day` column has the day; `[tumour] batch` decides which of the two the confounding rule uses).

**Shared-threshold robustness check** (ADR 0002; `[tumour] shared_threshold_check = true`, off by default): only when
no checked acquisition setting differs. One experiment-wide grey value = **median of the per-tumour automatic
thresholds** (fixed rule); every stack is re-segmented at it (`segment.sweep` with scale = shared / own threshold, so
everything else is the primary method; cache `_work/tumour/<slug>.shared.json`), and the same statistics run on the
result. Outputs: `shared_threshold.csv` / sheet **Shared threshold** (per contrast primary vs shared estimate, CI,
Holm p, direction / significance agreement), `plots/<name>_vs_shared_threshold.png|pdf` (per tumour) and a report
section. An effect present only at the shared threshold is a brightness difference (a dimmer tumour shrinks at a
fixed grey value), not a size difference; never starred.

## QC, problem tumours and blinded review (ticket 06, ADR 0002 + amendment 2026-11-20)
**Rules, set before results are seen — never per-tumour judgement.** Computed per tumour from its (possibly
corrected) labels and the original stack (`tumour/qc.py`), columns in `tumours.csv`:

| column | meaning |
|---|---|
| `edge_contact` | `yes` when the main mass touches the xy image frame in any layer |
| `edge_outline_frac` | share of the main mass's outline (footprint, all layers merged) lying on the frame; sub-pixel contour lengths (a disc cut in half: 0.39) |
| `edge_outline_frac_mid` | the same per mid-section layer, averaged |
| `mid_section_fits` | `no` when the mid-section (same rule as the readouts, `shape.mid_section`) runs off the imaged layers |
| `saturated_frac_mid_section`, `saturation_warning` | share of the mid-section main mass at 255; warning above `qc_saturation_warn` (default 40 %: the intermedians threshold is protected up to ~50 %). A warning, not a flag |
| `truncated`, `flagged` | edge contact **or** a mid-section that does not fit. Extending beyond the imaged depth is **not** a flag (every tumour does) |
| `qc_flags` | the reasons in words |
| `main_mass_object` | which object is the main mass (1 = automatic; another number after a review correction) |

Flagged tumours stay in the primary analysis (excluding them would bias toward "smaller": large tumours touch the
edge more often), are drawn as **hollow markers** in every group plot, and are dropped in a **sensitivity analysis**:
the same model and families on the table without them, shown beside the primary in sheet **Without flagged**
(primary estimate / CI / Holm p / N beside the same without flagged, `conclusion_changed`) and as stats_long rows with
family `<id> (without flagged)` (never starred). Real data (2 stacks): HERV −clodro n2 A has 51 % of its mid-section
main mass saturated (warning).

**Flow table** per group (every factor): files → tumours → excluded (reasons from `excluded_tumours.csv`: layout.csv
exclusions and review exclusions) → analysed → flagged → analysed without flagged, plus a total row: `flow.csv`,
sheet **QC flow**, report section `## QC and problem tumours` (also lists flagged tumours, saturation warnings and
the review-corrections generation). Sheet **QC** has the QC columns per tumour.

**A–B hemisphere check** (exploratory, not a test): `plots/<name>_hemispheres_AB.png|pdf` per confirmatory readout —
hemisphere A vs B of the same section (same group, imaging session and section label), log-log with the identity line,
one colour per group, r of the logs. Each hemisphere stays its own n; this only shows whether sections share an effect.

### Blinded review sheets (`review/`)
- `review/<code>.jpg` per tumour (downsampled to `review_px` per panel, JPEG `review_jpeg_quality`; ~70–150 kB on the
  real stacks): left the maximum projection over all layers with every object's footprint, right the mid-section
  centre layer. Main mass (object 1) orange, satellites cyan; the 40 largest objects are numbered.
- `review/overview_<k>.jpg` — contact sheets, 12 tumours each, in code order.
- `review/corrections.csv` — one row per tumour in code order: `code`, `objects`, `qc_flags`, `manual_1..3`, `note`.
- **Blinding:** codes are random (5 characters without vowels) and fixed across runs; the key `review_key.csv` (code
  → tumour ID) is written in the project folder, **outside** `review/`. Nothing in `review/` — file names, image
  titles, the corrections file — carries a factor level, imaging session, date, section label or file name; order is
  by code, not by group; each image is contrast-stretched on its own (brightness differences between lines are not
  displayed as such). Hand the reviewer only the `review/` folder. (The tests search every file name and byte in
  `review/` for level, session and file-name strings.)

### Corrections (re-pick the main mass, or exclude — nothing else)
Type into `review/corrections.csv`, like `manual_1` in the cell-count `counts.csv`:
- a whole number — that object (number on the review image) becomes the main mass (`1` = back to automatic);
- `x: reason` — exclude the tumour; the reason goes to `excluded_tumours.csv` (`excluded at review (manual_1): …`)
  and the flow table;
- empty — no change. Each `manual_<k>` column is one round; later rounds override earlier ones; the chain stops at
  the first completely empty column. Generation **g1** = automatic, **g2** = + `manual_1`, … (reported in
  `stats_report.md`; `corrections_used.csv` in the project folder lists per tumour the code, ID, main-mass object,
  source and generation).

Then `kerncount run <experiment>` again. Nothing is re-segmented: the cached labels are relabelled (picked object
and object 1 swap numbers) and the main-mass readouts (area, perimeter, volume, shape, saturation, surface…) are
recomputed from them with the same functions as the segmentation (`corrections.main_mass_readouts` → `shape.readouts`,
`shape.surface_area`); QC is recomputed. An unknown code, an unreadable cell or an object number that does not exist
stops the run with every problem listed. **Per-tumour threshold changes are not possible** (ADR 0002: one rule for
every tumour). Unlike the cell-count type, earlier generations are not kept as separate folders: the outputs are
rebuilt for the current generation. Caches: `_work/tumour/<slug>.qc.json` (key = measurement key + picked object +
QC settings; the stack is only read again when one of them, the code or the review settings change).
`--no-review` skips drawing (QC is still computed).

Known limitation: the threshold sweep (ticket 07) keeps no labels, so for a re-picked tumour its swept values still
follow the automatic main mass.

## Web interface and the analysis-plan lock (ticket 10)
Everything above is available in `kerncount gui` (details: `docs/ui.md`, "Tumour experiments"): set up with
**experiment type = Tumour volume** (`init --type tumour_volume`); **Layout** shows `layout.csv` with the validation of
`layout.validate` (problems block the run) and edits factor values, section, hemisphere, include and note;
**Settings** shows the `[tumour]` settings with explanations and the **analysis plan** (factors with levels,
reference and display role; readouts confirmatory / exploratory; families) with its lock; **Review images** and
**Corrections** are blinded (codes only; they read nothing outside `review/`) and edit the `manual_k` rounds of
`review/corrections.csv`; **Results** shows the plots per readout, the headline statistics (sensitivity families
` (Welch-type)`, ` (permutation)`, ` (without flagged)` hidden unless asked for), tables and the `results.xlsx`
download; the home page runs the analysis (*skip review images*, *segment every stack again*, *threshold sweep*).

**Lock (`tumour/lock_log.py`, ADR 0003).** Once `stats_long.csv` exists its `lock <hash>` is the results' plan. The
web interface refuses any settings write that changes the plan hash unless the plan was unlocked with a reason.
`lock_log.csv` (project folder) records `unlock`, `plan changed` (old → new), `lock` (locked again) and `run` (a
rebuild with a different plan, noting the unlock — or "without a logged unlock" when settings.toml was edited
outside the web interface; a run also ends a pending unlock). `results.analyse` calls `lock_log.record_run` before
writing and appends `lock_log.report_lines` to the About lines, so the history (old → new lock, reasons) is printed in
`stats_report.md` and the About sheet of `results.xlsx`. Not part of the plan (free to change): display roles,
measurement, QC, sweep, sensitivity and planner settings.

## Memory
One stack at a time: uint8 stack + float32 smoothed copy, then mask + int32 labels; peak ≈ 10 × the stack size
(2.8 GB for the largest, 280 MB stack). Hole filling labels the background instead of iterative dilation.
