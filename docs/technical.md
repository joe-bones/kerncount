# cellcount

Automated nucleus counting (DAPI), per-cell intensity measurement, quality control, manual
correction and full statistics from single fields up to all experiments (n) combined.
Replaces `NucleiCount new.ijm`, the Mean-intensity macros, and the Prism bookkeeping.

## Install (once)

1. Install **Python 3.12** from python.org (Windows: tick "Add python.exe to PATH").
2. Double-click `setup.bat` (Windows) or run `./setup.command` (Mac).
   This creates a private environment in `.venv` and installs everything (~1 GB, mostly TensorFlow).

## Use

```
cellcount init  "<raw image folder>"  "<project folder>"
cellcount run   "<project folder>"
```
On Windows, `cellcount.bat` runs the same commands without activating anything:
`cellcount.bat run "D:\analysis"`.

**init** scans the raw folder and writes two files into the project folder:

- `layout.csv` lists one row per raw folder, mapping messy folder names to clean names (`n1_2026-06-05`,
  `miR-21-5p`, `LyoVec`, ...). Check it once; edit `treatment` or set `include` to `no` if needed.
- `settings.toml` holds all parameters, with comments. The defaults are sensible.

**run** does everything, and is safe to repeat at any time. Segmentation is cached, so a repeat run
after a manual correction takes seconds.

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
  `cellcount run <project> --rebuild g2`, or `--rebuild all`.
- Before `counts.csv` is rewritten, the previous version is backed up to `_work/backups/`.
- If your Excel saves CSV with `;`, that is fine; it is detected automatically.

### Inside a generation
```
g1/
  README.txt, generation.json, counts_used.csv
  results.xlsx                 human-readable workbook (all tables below)
  <n>/<treatment>/fields.csv   per-field values + count_fields.png
  <n>/count_summary.csv        per treatment: fields, mean, SD, SEM, % of LyoVec / Control
  <n>/count_prism.csv          Prism Column table (fields as rows)
  <n>/count_stats.csv          ANOVA + Dunnett vs reference (fields = technical replicates)
  <n>/count.png/.pdf
  all_n/count_summary.csv      per-n means, mean of n, SD, SEM, % of controls
  all_n/count_prism.csv        Prism Column table: rows = n, columns = treatments
  all_n/count_prism_pct_of_LyoVec.csv   normalised within each n
  all_n/count_stats.csv        RM one-way ANOVA + ratio paired t-tests vs reference (Holm)
  all_n/*.png/.pdf             bar = mean of n, error = SEM, dot colour = n
  ... the same set for each intensity metric (intensity_<channel>_<region>)
```

## How counting works

1. **StarDist** (a deep-learning model for fluorescent nuclei; pretrained model downloaded on first run) segments every nucleus. The
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
cellcount sampling "<project folder>"          (uses the latest generation; --generation g2 to choose)
```
Uses the plate design in `settings.toml` `[design]`: 2 plates × 6 fields, images 1–6 = plate 1, 7–12 = plate 2.
Writes `sampling/<gN>/`:

- `REPORT.txt` gives the plain-language summary.
- `variance_components.csv` splits the variation into field-to-field (within a plate) and plate-to-plate (within an n).
- `design_power.csv` gives, for every design from 1 plate × 1 field to 2 × 6 and for each treatment: the predicted noise of
  the across-n test relative to now, the power at the current number of n's, and the n's needed for 80% power.
- `subsampling.csv` comes from drawing real images many times for each design and re-running the across-n test. It records
  how often the result matches the conclusion from all 12 images, and how far each n's ratio moves.
- `images_needed.png/.pdf` plots the same results.

The key point: biological n-to-n variation is only reduced by more n's. Extra images help only while
field and plate noise are still large compared with it. Plate-to-plate noise is only reduced by imaging more plates.

## Intensity (replaces the Mean-intensity macros)

For every counted nucleus, each channel is measured in the **nucleus**, the **ring** (cytoplasm band)
and the **cell** (nucleus + ring; the Fiji macro's "enlarge 6 px" region), with background
subtracted. Ring width is 0.25 × nucleus diameter by default, so it adapts to magnification.
Per-cell values are in `_work/cells.csv`. Per-image means flow into the same per-n / all-n
analysis as the counts. For two-channel images (C1 DAPI, C2 p-eIF2α), set in `settings.toml`:
```
[channels]
names = ["DAPI", "pEIF2A"]
segment = 1
[intensity]
measure = [2]
```
For whole-cell mCherry segmentation instead of nuclei, set `segment` to the mCherry channel and
`backend = "classical"`.

## Tests
`pip install -e .[test]` then `pytest`.
