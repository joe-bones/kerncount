# kerncount

Automated nucleus counting (DAPI), per-cell intensity measurement, quality control, manual
correction and full statistics from single fields up to all experiments (n) combined.
Replaces `NucleiCount new.ijm`, the Mean-intensity macros, and the Prism bookkeeping.

## Install (once)

1. Mac/Linux: install **Python 3.12** from python.org. Windows: nothing; `setup.bat` gets Python 3.12 via uv into `%LOCALAPPDATA%\kc` (no admin, short path to avoid TensorFlow hitting the 260-char limit).
2. Double-click `setup.bat` (Windows) or run `./setup.command` (Mac).
   This creates a private environment in `.venv` and installs everything (~1 GB, mostly TensorFlow).

## Use

Work in a workspace folder with `data/<experiment>/` (raw images) and `analysis/<experiment>/` (outputs).
From the workspace folder:
```
kerncount init <experiment>        # or: kerncount init --all
kerncount run  <experiment>        # or: kerncount run --all
kerncount sampling <experiment>
```
Explicit paths also work: `kerncount init <raw folder> <project folder>`, `kerncount run <project folder>`.
Workspace lookup: the current folder or nearest parent containing `data/`, or `$KERNCOUNT_WORKSPACE`.
`settings.toml` stores `raw_dir` relative to the project (`../../data/<experiment>`), so workspaces are movable.

**init** writes `layout.csv` (folder → clean n / treatment names; check once) and `settings.toml` (all
parameters, commented). **run** does everything, is safe to repeat, and caches segmentation.

### Raw folder layout expected
```
RAW/<one folder per n, with date e.g. 05.06.26>/<one folder per treatment>/1.tif, 2.tif, ...
```
Any number of n's (1 is fine) and any number of images per folder. Images are numbered by
their file name. Folders named `output` and `*_mask.tif` files are ignored. Images may be 8- or
16-bit, greyscale, RGB or multi-channel, at any magnification. Nucleus size is measured per folder and
the detector is rescaled to it.

## Project folder after a run

```
project/
  settings.toml, layout.csv
  counts.csv          <- THE file to work from: one row per image, all n's
  qc_details.csv      count range, flags, focus/merge/debris numbers, Fiji count, intensities
  fiji_comparison.csv old Fiji macro count vs new count, per image
  review/<n>/<treatment>/NN.jpg       original | overlay, with count + QC in the header
  review/<n>/<treatment>_overview.jpg all fields of one treatment on one sheet
  generations/g1, g2, ...             analysis results (see below)
  _work/                              cache (safe to delete; it will just re-segment)
```

### counts.csv and generations

| n | treatment | image | auto_count | qc_ok | qc_level | manual_1 | manual_2 | ... |
|---|---|---|---|---|---|---|---|---|

- **g1** is always fully automatic (`auto_count`).
- To correct images, type the new count into `manual_1`, or type `x` to exclude an image.
  Leave every other cell empty, then `run` again. This creates **g2**: g1 with those corrections applied.
- Later rounds go into `manual_2` (which creates g3, built on top of g2), and so on. The chain stops
  at the first manual column that is completely empty. A fresh empty column is always added
  at the end.
- Existing generation folders are **never overwritten by default**. If you change a manual column
  after its generation was built, `run` warns you. To rebuild it on purpose:
  `kerncount run <experiment> --rebuild g2`, or `--rebuild all`.
- Before `counts.csv` is rewritten, the previous version is backed up to `_work/backups/`. The same holds for
  `settings.toml` and `layout.csv` (`init --force`, `--force-layout`, the page's reset / fresh-scan buttons) and
  `wells.csv` (`wells-template --force`): `_work/backups/<name>_<date-time>.<ext>`, never overwritten.
- A generation is built in `generations/gN.tmp/` and renamed into place only when complete. Files you put
  into a generation folder yourself are moved to `gN/_kept/` on a rebuild, never deleted. On Windows, a file of
  gN open in Excel stops the swap with a clear message (close it and run again).
- Bare experiment names are looked up in `data/` first (`kerncount run Ki67` uses `data/Ki67` even when the
  working directory has a `Ki67` folder); a symlinked `data/<experiment>` is kept as the link (settings store
  `../../data/<experiment>`), not its current target.
- If your Excel saves CSV with `;`, that is fine; it is detected automatically.

### Inside a generation
```
g1/
  README.txt, generation.json, counts_used.csv
  results.xlsx                 human-readable workbook (all tables below)
  <n>/<treatment>/fields.csv   per-field values + count_fields.png
  <n>/count_summary.csv        per treatment: fields, mean, SD, SEM, % of LyoVec / Control
  <n>/count_prism.csv          Prism Column table (fields as rows)
  <n>/count_stats.csv          ANOVA + Dunnett vs each control group present (Control, then LyoVec ->
                               count_vs_LyoVec.png) (fields = technical replicates)
  <n>/count.png/.pdf
  all_n/count_summary.csv      per-n means, mean of n, SD, SEM, % of controls
  all_n/count_prism.csv        Prism Column table: rows = n, columns = treatments
  all_n/count_prism_pct_of_LyoVec.csv   normalised within each n
  all_n/count_stats.csv        RM one-way ANOVA + Dunnett vs each control (pooled error; stars) + ratio paired
                               t-tests vs each control (Holm per control; conservative check)
  <n>/count_twoway*, all_n/count_twoway*   two-way ANOVA (per n: treatment x well; all n: treatment x n; fields =
                               replicates; Type II SS) + Tukey-Kramer all pairs matched within well / n
  all_n/*.png/.pdf             bar = mean of n, error = SEM, dot colour = n
  ... the same set for each intensity metric (intensity_<channel>_<region>)
```

## How counting works

1. **StarDist** (a deep-learning model for fluorescent nuclei; pretrained model downloaded on first run and saved
   next to the package, or in `~/.kerncount/models` when the package folder is read-only, e.g. a system-wide
   install or an app bundle; `KERNCOUNT_MODELS` overrides) segments every nucleus. The
   image is first rescaled so that nuclei are about 30 px across, which makes it independent of
   magnification and camera.
2. **Over-split repair** merges objects that touch when they are one lobed or fragmented nucleus:
   a small fragment, or a nucleus-sized pair with no dark seam between the parts.
3. Each object is classified:
   - **counted**
   - **possible merge** (counted as 1, but the upper bound is raised)
   - **debris/micronucleus** (smaller than 0.2 × a typical nucleus)
   - **out of focus**: blurred edges AND no chromatin texture, i.e. a floating or dead cell outside the
     focal plane. Texture is measured with camera noise subtracted. This stops bright, sharp nuclei
     that sit in glare from being thrown out. Saturated nuclei are never called out of focus.
   - **edge** (unbiased counting frame: nuclei cut by the bottom/right edge are dropped and those
     cut by the top/left are kept, so on average edge nuclei are counted exactly once)
4. **Error prediction** for each image:
   - a *plausible range* (count_low to count_high) built from every borderline decision
     (near the size or focus threshold, low detector confidence, possible merges)
   - agreement with an independent classical threshold-plus-watershed count
   - saturation

   These give `qc_level` = high / medium / low. They are heuristics until they can be calibrated
   against hand counts (ground truth).

## How many images are needed?

```
kerncount sampling "<project folder>"          (uses the latest generation; --generation g2 to choose)
```
Uses the design in `settings.toml` `[design]`: one plate per n, each condition in 2 wells × 6 fields, images 1–6 = well 1,
7–12 = well 2. Old key names (`plates`, `fields_per_plate`, `plate_layout`) are still read. Everything is computed once per
reference group (`reference` + `normalise_to`), matching the across-n statistics, so an n missing one reference still
counts against the others.
Writes `sampling/<gN>/`:

- `REPORT.txt` gives the plain-language summary.
- `variance_components.csv` splits the variation into field-to-field (within a well) and well-to-well (within a plate = n).
- `batch_by_condition.csv`: per n × treatment, the mean per well and log(well 2 / well 1) as a z-score against pooled field
  noise, SE = sqrt(cv2_field · (1/k1 + 1/k2)); |z| > 2.5 is flagged.
- `batch_by_n.csv`: per n, reference levels, missing references, geometric-mean well 2 / well 1, and the effect-size shift:
  this n's log ratio minus the across-n mean, averaged over non-control treatments. It is flagged when the shift is > 1.5× either
  way and at least 80% of treatments move in the same direction.
- `design_power.csv` gives, for every design from 1 well × 1 field to 2 × 6 and for each reference and treatment: the predicted noise of
  the across-n test relative to now, the power at the current number of n's, and the n's needed for 80% power.
- `subsampling.csv` comes from drawing real images many times for each design and re-running the across-n test. It records
  how often the result matches the conclusion from all 12 images, and how far each n's ratio moves.
- `images_needed.png/.pdf` plots the same results.

The key point: biological n-to-n variation is only reduced by more n's. Extra images help only while
field and well noise are still large compared with it. Well-to-well noise is only reduced by imaging more wells.
Plate/day batch effects are confounded with n (one plate per n): multiplicative plate-wide ones cancel in the ratio,
effect-size ones (e.g. transfection efficiency) remain in the n-to-n variance and are surfaced by the shift flag.

## Intensity (replaces the Mean-intensity macros)

For every counted nucleus, each channel is measured in the **nucleus**, the **ring** (cytoplasm band)
and the **cell** (nucleus + ring; the Fiji macro's "enlarge 6 px" region), with background
subtracted. Ring width is 0.25 × nucleus diameter by default, so it adapts to magnification.
Per-cell values are in `_work/cells.csv`. Per-image means flow into the same per-n / all-n
analysis as the counts. For two-channel images (C1 DAPI, C2 marker), set in `settings.toml`:
```
[channels]
names = ["DAPI", "marker"]
segment = 1
[intensity]
measure = [2]
```
For whole-cell mCherry segmentation instead of nuclei, set `segment` to the mCherry channel and
`backend = "classical"`.

## Significance on plots
Bars get stars only when their comparison has `p_adj` < 0.05 (no mark = not significant); the control is
labelled "ref" and a footer names the test. Default reference is Control; every group in `normalise_to`
(Control, LyoVec) gets its own test family and plots. Per-n plots: Dunnett vs each control. All-n plots (Prism standard): RM
one-way ANOVA (value ~ treatment + n on n means, least squares, so missing groups are allowed) then Dunnett vs the
plotted control with the pooled residual error (`analysis.rm_anova_dunnett`; adjusted p from the multivariate t with
the contrasts' correlation, checked against `scipy.stats.dunnett`). RM ANOVA p in the footer. Implemented in
`analysis.significance_marks` / `analysis._bar_plot`; covered by `test_significance_marks_only_where_significant`.

## Two-way ANOVA (`analysis.two_way`)
`value ~ treatment + B + treatment:B`, fields as replicates, B = well (per n) or n (all n). Sums of squares are
Type II by model comparison on dummy-coded least squares, so empty cells (a group missing in one n) only drop
their interaction column; checked against statsmodels `anova_lm(typ=2)`. Without replication the additive model
is used. Pairs: Tukey-Kramer on the mean of matched cell-mean differences (B levels containing both groups),
SE from the model MSE, p from the studentized range with k treatments. One B level = one-way ANOVA + Tukey,
matching `scipy.stats.tukey_hsd` (tested). Toggle with `[analysis] two_way_anova`.

## Experiment types (`experiment_types`, ADR 0004)
`[experiment] type` in `settings.toml` selects how an experiment is handled; omitted = `cell_count` (everything
above). A type is an `ExperimentType` subclass registered with `@register` and implements four seams:
`scan_layout(raw, settings)` (raw folder → `layout.csv` rows; `init`), `measure(...)` (raw measurements per unit),
`readout_table(...)` (one row per datapoint with its layout metadata, written to the project) and `analyse(...)`
(the type's statistical design). `kerncount init` calls `scan_layout` + `report_layout`; `kerncount run` calls
`run()`, which chains measure → readout table → analyse. `cell_count` overrides `run()` to call `pipeline.run`
unchanged (its stages are interleaved there). An unknown type stops with the list of valid types. Data system,
backups and GUI are shared; other commands (`regen`, `sampling`, `platemap`, …) are cell-count only.
`kerncount init <exp> --type <name>` writes that type's own settings template (`settings_template()`; default: the
cell-count template with the type line changed). `tumour_volume` lives in the `kerncount.tumour` package (imported at
the bottom of `experiment_types.py` so it is always registered); see [tumour.md](tumour.md).

## Tests
`pip install -e .[test]` then `pytest`.
