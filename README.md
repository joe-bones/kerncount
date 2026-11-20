# kerncount: automated nucleus counting and analysis

*[Deutsche Fassung](README.de.md)*

`kerncount` replaces a manual Fiji-macro + Prism workflow for fluorescence cell-count experiments with one
Python tool. It does the following:

1. Counts nuclei in every image (DAPI), with a deep-learning detector.
2. Estimates how reliable each count is and draws a review image for each field.
3. Lets you correct counts by hand in one CSV file, keeping every round of corrections as a "generation".
4. Produces tables, Prism-ready CSVs, plots and statistics for each treatment, each n, and all n's together.
5. Measures intensity per cell (any channel), as the Mean-intensity macros did.
6. Estimates how many images per condition you actually need.
7. Measures mCherry tumours in brain-section z-stacks (tumour experiments, section 9), with blinded review and pre-registered statistics.

```
kerncount/                 (this repository)
├── README.md              ← this guide (README.de.md = German version)
├── docs/technical.md      ← algorithms and implementation notes
├── src/kerncount/         ← the code
├── tests/                 ← unit tests (pytest)
├── launcher/              ← source of the one-click kerncount.exe / kerncount.app
└── setup.bat / setup.command / kerncount.bat / kerncount-gui.*   ← install and launch helpers
```
Your images and results live in a separate **workspace** folder, with input and output kept apart:
```
my_workspace/                    ← run all commands from this folder
├── data/                        ← INPUT: one folder per experiment (raw images; never modified)
│   ├── my_experiment/
│   └── another_experiment/
└── analysis/                    ← OUTPUT: one project per experiment, same name as in data/
    ├── my_experiment/
    └── another_experiment/
```
Nothing in `data/` is ever changed; everything produced for an experiment is in `analysis/<experiment>/`.

---

## 1. Input: how the raw data must be organised

Each experiment is one folder in `data/`:

```
data/<experiment>/                              e.g. data/my_experiment/
├── <one folder per n (experiment)>/            e.g. "05.06.26 DAPI count"
│   ├── <one folder per treatment>/             e.g. "hsamiRNA 21-5p", "LyoVec veh. cntrl", "cntrl"
│   │   ├── 1.tif
│   │   ├── 2.tif
│   │   └── ...
│   └── ...
└── ...
```

- **One folder per n**, with the date in the name (`05.06.26`, read as day.month.year). The n's are
  numbered by date: `n1_2026-06-05`, `n2_2026-07-06`, and so on. A single n works fine.
- **One folder per treatment** inside each n. Messy names are cleaned up automatically, and you can
  check the mapping in `layout.csv`:
  `hsamiRNA 21-5p` / `hsa miRNA 21-5p` / `21 5p` / `21` → `miR-21-5p`,
  `LyoVec veh. cntrl` / `lyo vec` → `LyoVec`, `cntrl` / `010125 cntrl` / `-PDL cntrl` → `Control`. (`LyoVec`, `CL075` and miRNA patterns are built-in rules; anything
  unrecognised is kept as-is and flagged.)
- **Images are named by number** (`1.tif`, `02.tif`; stray characters such as `` `10.tif `` are fine).
  Files without a number (e.g. `5 zoom.tif`) are skipped and listed in `layout.csv`.
  Any number of images per folder is fine.
- **Ignored:** sub-folders called `output`, files containing `mask` (old Fiji results), and other
  sub-folder levels (images are found at any depth).
- **Image format:** TIFF, 8- or 16-bit, greyscale, RGB or multi-channel, at any magnification and from any camera.
  Nucleus size is measured per folder and the detector rescales to it.
- **Wells:** one n = one plate. Each condition sits in 2 wells of that plate; by default images 1–6 = well 1
  and 7–12 = well 2 (only used by the sampling analysis; change it in `settings.toml`).

---

## 2. Install (once per computer)

### Easiest: one-click app (no admin rights, no Python, no commands)
Download from the [latest release](https://github.com/joe-bones/kerncount/releases/latest):
- **Windows:** `kerncount.exe`
- **Mac (Apple silicon):** `kerncount-mac.zip`. Double-click to unzip; this gives `kerncount.app`.

Put it in your workspace folder (next to `data`) and double-click it. The first start needs internet, takes
a few minutes and about 1 GB: it installs Python 3.12 and everything kerncount needs into your user folder.
After that, every start opens kerncount in your browser straight away. The first time, the system warns about an
unsigned app: on Windows click **More info → Run anyway**; on the Mac click **Done**, then
**System Settings → Privacy & Security → Open Anyway**. Details: [launcher/README.md](launcher/README.md).

### Alternative: install from this repository (for the command line)
1. **Mac/Linux:** install **Python 3.12** from python.org (3.11 or newer works; 3.12 is recommended for TensorFlow).
   **Windows:** nothing to install first, and no admin rights needed: `setup.bat` fetches Python 3.12 itself
   (via [uv](https://docs.astral.sh/uv/)) into your user folder (`%LOCALAPPDATA%\kc`).
2. Download this repository (green **Code** button → *Download ZIP*, or `git clone`).
3. Double-click `setup.bat` (Windows) or run `setup.command` (Mac/Linux).
   This takes a few minutes and about 1 GB. The StarDist model (~5 MB) downloads automatically on the
   first run, so that run needs internet.

Commands are run from your workspace folder. On Windows, every command below starts with the full path
to `kerncount.bat` in this folder instead of `kerncount` (e.g. `C:\tools\kerncount\kerncount.bat run my_experiment`).
On Mac/Linux it is `<this folder>/.venv/bin/kerncount`. Adding that folder to your PATH lets you type just `kerncount`.

---

## 3. How to use it

Run every command from the workspace folder (the one that contains `data/`). Experiments are referred
to by their folder name in `data/`.

### Without commands: the web interface
Double-click **`kerncount.exe`** / **`kerncount.app`** (section 2), or, with a repository install,
**`kerncount-gui.bat`** (Windows) / **`kerncount-gui.command`** (Mac). The same thing from a terminal: `kerncount gui`.
Everything below can be done in the browser, without editing files:

| Tab | What you do there |
|---|---|
| **Run** | Set up, run the analysis (with the command options below as tick boxes), sampling; live log, **Stop** |
| **Layout** | Fix n / treatment names, include or leave out folders |
| **Settings** | Every toggle from section 4 as a form, with explanations; or edit `settings.toml` as text |
| **Review images** | Overview sheet per treatment, then each image full size (← / → to browse); type a correction directly |
| **Corrections** | Edit a generation: type corrected counts or **x** for a round (`manual_1` → g2, …), see the value before and after, then **build** the generation (about a minute; detection is reused) |
| **Results** | All plots, tables and statistics per generation; download `results.xlsx` |
| **Sampling** | Run and read the "how many images?" analysis |
| **All files** | Browse and preview everything in the project folder |

Tumour experiments (section 9) use the same pages: Layout shows the tumours and what blocks a run, Settings the
analysis plan with its lock, Review images and Corrections show blinded codes only, Results the plots per readout.

Keep the black (Windows) or Terminal (Mac) window open while you work; close it or press **Quit** to stop.
The page only works on this computer, and nothing is sent anywhere else.

### Step 1: create a project (once per experiment)
```
kerncount init my_experiment
```
This scans `data/my_experiment/` and creates `analysis/my_experiment/` with two files
(`kerncount init --all` does this for every experiment in `data/`):

- **`layout.csv`** has one row per raw folder. Check the `n` and `treatment` columns. To fix a name, edit it;
  to leave a folder out, set `include` to `no`. The `note` column shows anything unusual
  (skipped files, inferred names, "72h").
- **`settings.toml`** holds all the toggles (section 4). The defaults are sensible.

### Step 2: run the analysis
```
kerncount run my_experiment
```
(`kerncount run --all` runs every initialised experiment.)
The first run takes about 2–3 s per image (roughly 15 min for 370 images). After that, results are
cached and a re-run takes under a minute. It is **always safe to re-run**.

### Step 3: review the counts
Open `analysis/my_experiment/review/<n>/<treatment>_overview.jpg`. There is one sheet per treatment, showing every field
with its count and a border coloured by quality: **green = high, orange = medium, red = low**.
For a closer look, open the single image `review/<n>/<treatment>/<NN>.jpg`
(left: original; right: what was counted). Start with the red ones.

### Step 4: correct counts by hand (optional)
Open **`analysis/my_experiment/counts.csv`** in Excel:

| n | treatment | image | auto_count | qc_ok | qc_level | manual_1 | manual_2 | manual_3 |
|---|---|---|---|---|---|---|---|---|
| n1_2026-06-05 | Control | 3 | 85 | yes | high | | | |
| n1_2026-06-05 | miR-21-5p | 6 | 140 | no | low | **121** | | |
| n1_2026-06-05 | miR-21-5p | 9 | 4 | no | low | **x** | | |

- In `manual_1`, type a corrected count, or **`x`** to exclude that image. Leave all other cells empty.
- Save, then `kerncount run my_experiment` again. This creates generation **g2** (see section 5).
- A later round of corrections goes into `manual_2`, which makes g3, and so on.

### Step 5: how many images do we need? (optional)
```
kerncount sampling my_experiment
```
Results are written to `analysis/my_experiment/sampling/g<N>/` (section 5.5).

### Measuring a second channel (e.g. a marker in C2)
For multi-channel TIFFs (C1 = DAPI, C2 = signal), set in `settings.toml`:
```
[channels]
names = ["DAPI", "marker"]
segment = 1            # find nuclei in DAPI

[intensity]
measure = [2]          # measure the marker
regions = ["nucleus", "ring", "cell"]
```
Then run as usual. Per-cell values go to `_work/cells.csv` (one row per detected object, with its
category). Per-image means go to `qc_details.csv` and through the same per-n and all-n statistics as the
counts (`intensity_marker_cell_*` files). To segment whole cells from an mCherry channel instead of nuclei,
set `segment` to that channel and `[segmentation] backend = "classical"`.

### Command options
| Command | What it does |
|---|---|
| `kerncount init <experiment>` | Create `analysis/<experiment>/` from `data/<experiment>/` |
| `kerncount init --all` | Same for every experiment in `data/` (existing projects are kept) |
| `kerncount init <experiment> --force` | Re-scan and overwrite `layout.csv` and `settings.toml` |
| `kerncount run <experiment>` / `--all` | Count, check and analyse one / every experiment |
| `kerncount run <experiment> --rebuild g2` | Overwrite an existing generation folder (never done by default). `--rebuild all` for all. |
| `kerncount run <experiment> --resegment` | Ignore the cache and re-detect every nucleus |
| `kerncount run <experiment> --no-review` | Skip drawing review images (faster) |
| `kerncount regen <experiment> [--rebuild g2 g3]` | Rebuild generations from `counts.csv` only (fast; after corrections) |
| `kerncount sampling <experiment> --generation g2` | Sampling analysis on a particular generation (default: latest) |

Paths work too, for folders outside a workspace: `kerncount init <raw folder> <project folder>`,
`kerncount run <project folder>`. Setting the environment variable `KERNCOUNT_WORKSPACE` lets you run
by name from anywhere.

---

## 4. Toggles (`settings.toml`)

Edit with any text editor. A setting you delete falls back to its default.

### `[project]`
| Setting | Default | Meaning |
|---|---|---|
| `raw_dir` | (set by init) | Where the raw images are. Stored relative to the project folder (`../../data/<experiment>`), so the whole workspace can be moved or copied. |
| `date_order` | `"DMY"` | How dates in n-folder names are written (`DMY`, `MDY`, `YMD`) |
| `csv_delimiter` | `","` | Use `";"` if your Excel (German/European locale) shows CSVs in one column. Reading auto-detects either. |

### `[channels]`
| Setting | Default | Meaning |
|---|---|---|
| `names` | `[]` | Channel names, e.g. `["DAPI", "marker"]`. Empty = C1, C2, ... |
| `segment` | `1` | Channel used to find nuclei |

### `[segmentation]`: finding nuclei
| Setting | Default | Meaning |
|---|---|---|
| `backend` | `"stardist"` | `stardist` = deep-learning detector (recommended); `classical` = Fiji-style threshold + watershed |
| `nucleus_diameter_px` | `"auto"` | `auto` measures nucleus size per folder (adapts to magnification). Or a fixed number. |
| `target_diameter_px` | `30` | Size the detector works best at; images are rescaled to it internally |
| `prob_thresh` | `"default"` | Detector confidence threshold (lower = more objects) |
| `pixel_size_um` | `""` | µm per pixel, if known (from a stage micrometer). Adds µm² and cells/mm². Changing it does not trigger re-detection. |

### `[counting]`: what counts as a nucleus
| Setting | Default | Meaning |
|---|---|---|
| `edge_rule` | `"half"` | `half` = unbiased counting frame (nuclei cut by the top/left edge count, bottom/right do not, so on average each edge nucleus counts once); `exclude` = drop all edge nuclei; `include` = count all |
| `min_area_frac` | `0.2` | Objects smaller than 0.2 × a typical nucleus = debris or micronucleus (not counted) |
| `exclude_out_of_focus` | `true` | Drop nuclei outside the focal plane (floating or dead cells) |
| `focus_frac` | `0.5` | Edge sharpness below 0.5 × typical ... |
| `texture_frac` | `0.8` | ... AND chromatin texture below 0.8 × typical = out of focus. Both must be low, so bright nuclei sitting in glare are not lost. |
| `merge_area_frac` | `2.0` | Objects larger than 2 × typical and non-convex = "possible merge" (counted as 1, flagged) |

### `[qc]`: how strict the quality levels are
| Setting | Default | Meaning |
|---|---|---|
| `high_max_uncertainty` | `0.10` | QC high requires the plausible range to be ≤ 10% of the count ... |
| `high_max_disagreement` | `0.15` | ... and the second method to agree within 15% ... |
| `max_saturated_frac` | `0.02` | ... and ≤ 2% of nuclear pixels saturated |
| `low_min_uncertainty` | `0.25` | QC low if the range is ≥ 25% of the count ... |
| `low_min_disagreement` | `0.35` | ... or the two methods differ by ≥ 35% |

For very small counts, the range is taken relative to at least 5, so an image with 2 cells
isn't automatically "low". An image with 0 nuclei is always "low" (check for a blank or defocused field).

### `[intensity]`: per-cell intensity (replaces the Mean-intensity macros)
| Setting | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Measure intensities |
| `measure` | `"all"` | Channels to measure, e.g. `[2]` for channel 2 only |
| `regions` | `["nucleus", "cell"]` | `nucleus`; `ring` (band around the nucleus, i.e. cytoplasm); `cell` (nucleus + ring, which is the Fiji macro's "enlarge 6 px") |
| `ring_px` | `"auto"` | Band width; `auto` = 0.25 × nucleus diameter |
| `subtract_background` | `true` | Subtract each image's background (median of empty areas) |

### `[analysis]`: statistics
| Setting | Default | Meaning |
|---|---|---|
| `reference` | `"Control"` | Main comparison group (headline plots) |
| `normalise_to` | `["Control", "LyoVec"]` | Groups to express results as % of (within each n), and to test against. Every group listed gets its own tests and plots, so you always get both *vs Control* and *vs LyoVec* |
| `two_way_anova` | `true` | Extra two-way ANOVA with Tukey all-pairs comparisons: per n treatment × well, all n treatment × n (5.4.3) |
| `treatment_order` | `["Control", "LyoVec", "CL075"]` | Order in tables and plots; other treatments follow alphabetically |
| `exclude_low_qc` | `false` | If `true`, images with QC low (and no manual value) are left out of the statistics |

### `[design]`: used by `kerncount sampling`
| Setting | Default | Meaning |
|---|---|---|
| `wells` | `2` | Wells imaged per condition per n (one plate per n) |
| `fields_per_well` | `6` | Fields (images) per well |
| `well_layout` | `"consecutive"` | `consecutive` = images 1–6 well 1, 7–12 well 2; `alternating` = odd/even |

Older settings files using `plates` / `fields_per_plate` / `plate_layout` still work; the report asks you to rename them.

### `[review]`
| Setting | Default | Meaning |
|---|---|---|
| `overlay_alpha` | `0.35` | Transparency of the coloured overlay |
| `jpeg_quality` | `88` | Review image quality |

---

## 5. Output: what comes out and how to read it

```
analysis/<experiment>/
├── settings.toml, layout.csv     your settings and folder mapping
├── counts.csv                    ★ THE file to work from: one row per image, all n's, manual columns
├── qc_details.csv                every per-image measurement (see 5.2)
├── fiji_comparison.csv           old Fiji count vs new count, per image
├── review/                       pictures to check the counting (5.1)
├── generations/g1, g2, ...       results (5.3, 5.4)
├── sampling/g1/                  how many images are needed (5.5)
└── _work/                        cache, cells.csv (per-cell measurements), backups/ of counts.csv
                                  (deleting _work forces re-detection)
```

### 5.1 Review images: `review/`
```
review/<n>/<treatment>_overview.jpg     all fields of one treatment; border = QC level
review/<n>/<treatment>/<NN>.jpg         one field: original | overlay
```
The header gives the count, the plausible range, the QC level and the reasons. Overlay colours:

| Colour | Meaning | Counted? |
|---|---|---|
| green | nucleus | yes |
| orange | possible merge of 2+ touching nuclei | yes, as 1 (raises the upper bound of the range) |
| magenta | debris, micronucleus or fragment | no |
| blue | out of focus (outside the focal plane) | no |
| yellow | cut by the bottom/right image edge | no (edge rule) |

### 5.2 Per-image numbers: `counts.csv` and `qc_details.csv`
- **`auto_count`**: the automatic count.
- **`count_low` – `count_high`** (*plausible range*): what the count would be if every borderline
  decision (near the size or focus threshold, low detector confidence, possible merges) went the
  other way. A narrow range means an unambiguous image.
- **`disagreement`**: difference from an independent classical (threshold + watershed) count.
  Large values mean the image is hard (dense, clumped, bright cytoplasm).
- **`qc_level`** (high / medium / low) combines the above. **Low means "look at this image", not
  "this count is wrong".** Low-QC images cluster in sparse, damaged wells (e.g. miRNA-treated), where
  counting really is harder.
- `n_debris`, `n_out_of_focus`, `n_edge`, `n_merged`: how many objects of each type were left out or flagged.
- `diameter_px`: typical nucleus size in that folder (useful as a magnification/biology check).
- `fiji_count`: what the old macro gave, for comparison.
- `intensity_<channel>_<region>`: mean background-subtracted intensity of the counted cells.
- `nuclei_per_mm2`: only present once `pixel_size_um` is set.
- `stardist_scale`, `saturated_frac`, `source_file`: detector rescaling factor, fraction of saturated
  nuclear pixels, and the raw file the row came from.

*About the QC levels:* they are rules of thumb until checked against hand counts (ground truth).
Hand-counting about 15 images across densities would turn them into real error estimates.

### 5.3 Generations: `generations/g1`, `g2`, ...
- **g1** = fully automatic counts. Always present.
- **g2** = g1 + the corrections typed into `manual_1`. **g3** = g2 + `manual_2`, and so on.
- The chain stops at the first manual column that is completely empty.
- An existing generation folder is **never overwritten** unless you ask (`--rebuild g2`). If counts.csv
  no longer matches an existing generation (you edited `manual_1` after g2 was built), `run` warns you.
- Each generation has `README.txt` (what it is built from), `counts_used.csv` (the value used for every image
  and where it came from: `auto`, `manual_1`, `manual_1: excluded`, ...) and `generation.json`.

### 5.4 Inside a generation: treatment → n → all n
```
g1/
├── results.xlsx                      ★ everything below in one readable workbook
├── counts_used.csv
├── <n>/                              one experiment
│   ├── <treatment>/fields.csv        per-field values (count + intensities)
│   ├── <treatment>/count_fields.png  bar per field (x = excluded)
│   ├── count_summary.csv             per treatment: fields, mean, SD, SEM, median, % of LyoVec, % of Control
│   ├── count_prism.csv               Prism Column table: columns = treatments, rows = fields
│   ├── count_stats.csv               one-way ANOVA + Dunnett vs each control (Control, LyoVec)
│   ├── count.png / .pdf              bar = mean, error = SD, dots = fields, stars = Dunnett vs Control (5.4.2)
│   ├── count_vs_LyoVec.png / .pdf    same, stars = Dunnett vs LyoVec (only if this n has LyoVec)
│   └── count_twoway.png / .pdf       two-way ANOVA treatment × well + Tukey all pairs (5.4.3)
│       (+ count_twoway_anova.csv, count_twoway_pairs.csv)
└── all_n/                            all experiments combined
    ├── count_summary.csv             per treatment: value per n, N, mean, SD, SEM, % of controls
    ├── count_prism.csv               Prism Column table: columns = treatments, rows = n
    ├── count_prism_pct_of_LyoVec.csv same, normalised to LyoVec within each n (also _Control)
    ├── count_stats.csv               RM one-way ANOVA + Dunnett vs each control (stars); ratio paired t-tests (check)
    ├── count.png / .pdf              bar = mean of n, error = SEM, dot colour = n (legend), stars (5.4.2)
    ├── count_pct_of_LyoVec.png/.pdf  normalised version (dashed line = 100%); also _Control
    ├── count_twoway.png / .pdf       two-way ANOVA treatment × n + Tukey all pairs (5.4.3)
    │   (+ count_twoway_anova.csv, count_twoway_pairs.csv)
    └── n_colours.csv
```
The same set of files exists for each intensity metric (`intensity_<channel>_<region>_*`).

**Prism:** open a `*_prism*.csv` and paste it into a Column table, or use File → Import.
The all-n files have the n label as a row title in the first column.

**`results.xlsx` sheets:** *About* (what this generation is) · *Counts used* · *Count - all n* ·
*Count - stats all n* · *Count - per n* · *Count - stats per n* · and, if intensity is on,
*Intensity - all n* · *Intensity - stats all n* · *Intensity - per n* · and the two-way results:
*Two-way anova - all n* · *Two-way pairs - all n* · *Two-way anova - per n* · *Two-way pairs - per n*.

#### 5.4.1 Columns of `*_stats.csv`
| Column | Meaning |
|---|---|
| `test` | Which test, and which control it compares against |
| `comparison` | `<treatment> vs <control>`, or `all` for an overall ANOVA |
| `N` | (all n) number of n's that contain both groups |
| `effect_type`, `effect`, `ci_low`, `ci_high` | (all n) Dunnett: mean difference vs the control and its 95% CI (Dunnett-adjusted); ratio t-test: geometric-mean fold change and its 95% CI |
| `statistic` | t, F or Dunnett statistic |
| `p` / `p_adj` | raw / multiple-comparison-adjusted p (Dunnett is adjusted by construction; t-tests by Holm) |
| `significance` | stars from `p_adj`: ns, \*, \*\*, \*\*\*, \*\*\*\* (< 0.05 / 0.01 / 0.001 / 0.0001) |

#### 5.4.2 How significance is shown on the plots
- **Stars** (\* to \*\*\*\*) above a bar appear only when that treatment's comparison is significant
  (`p_adj` < 0.05). **No mark = not significant.** The control used is labelled **"ref"** above its bar.
- **Per-n plots** use Dunnett vs the control: `count.png` vs `reference` (Control), `count_vs_<X>.png`
  vs each other control (LyoVec). If the reference is missing in that n, `count.png` uses the next control
  group, and the footer says which. These stars describe one experiment only.
- **All-n plots** use the Prism standard: repeated-measures one-way ANOVA (n = matched replicate), then
  Dunnett vs the control shown, using the pooled error from all treatments: the raw plot uses
  `reference` (Control); each `pct_of_<X>` plot uses X, so `pct_of_LyoVec` gives the stars vs LyoVec. The footer also gives the overall RM ANOVA p-value.
  A significant ANOVA with no starred bars means the treatments differ as a group, but no single
  comparison is significant yet with this many n's.
- The footer of every plot names the test, so a figure can be read without the tables.

**How to read the statistics:**
- **Within one n** (`<n>/count_stats.csv`): the fields are *technical* replicates, so these p-values only
  describe that experiment. Do not report them as the result. Dunnett is run against every control
  group present in that n (Control and LyoVec), each as its own family.
- **All n** (`all_n/count_stats.csv`): each n is one biological replicate (the mean of its fields).
  - *RM one-way ANOVA*: is there any difference at all? Each n is a matched replicate. It is fitted by
    least squares, so an n that misses one group (no LyoVec in n2) still counts for all its other groups.
    With nothing missing this is exactly Prism's RM one-way ANOVA; with missing values Prism would use its
    mixed-effects model, which gives very similar results.
  - *Dunnett vs Control / vs LyoVec* (**drives the stars**): each treatment vs the control, using the
    ANOVA's pooled residual error, adjusted for the number of comparisons. `effect` = mean difference
    with a 95% CI. Pooling gives far more degrees of freedom than a single paired test (e.g. 20 instead
    of 3 with 4 n's), so real effects reach significance with few n's. It assumes the n-to-n scatter is
    similar for all treatments.
  - *Ratio paired t-test vs LyoVec / vs Control* (**check**, no stars): for each treatment, the ratio to the control within
    each n, tested across n's. `ratio` = geometric-mean fold change (0.30 = 70% fewer cells), with a
    95% CI. `p_adj` is Holm-corrected within each control family. Stars: ns, \*, \*\*, \*\*\*, \*\*\*\*
    (0.05 / 0.01 / 0.001 / 0.0001).
  - A test only uses n's that contain *both* groups, so `N` can differ between comparisons when a
    group is missing from some n.
- **Normalised values** (`pct_of_...`) are generally the better way to compare across n's, because they
  cancel day-to-day differences in seeding, staining and lamp brightness.
- **Intensity values** are only comparable if exposure and gain were fixed (auto-exposure off) within each n.

#### 5.4.3 Two-way ANOVA with all-pairs comparisons (`*_twoway*`)
An extra analysis at both levels, like Prism's two-way ANOVA with "compare all pairs":
- **Per n:** treatment × well (wells from `[design]`); fields are the replicates.
- **All n:** treatment × n; fields are the replicates.

The model is `value ~ treatment + second factor + interaction` with Type II sums of squares, so a group
missing from one n (e.g. no LyoVec in n2) is handled correctly. `*_twoway_anova.csv` gives F, df, p and the
% of total variation for treatment, the second factor and the interaction. `*_twoway_pairs.csv` compares
every pair of treatments (Tukey–Kramer, `p_adj`) on cell means matched within each well / n, using only the
wells / n's that contain both groups (`matched_levels`), with the mean difference and its 95% CI. The plot
shows grouped bars (mean ± SEM of fields per well / n) next to a matrix of the adjusted p-values.

**Caution, all-n level:** treating fields as replicates makes these p-values more optimistic than the ratio
RM ANOVA + Dunnett, which uses one value per n. Report RM ANOVA + Dunnett as the result; use the two-way analysis to
see whether treatments differ from each other and whether the effect depends on the n (the interaction).

### 5.5 How many images are needed? `sampling/g1/`
| File | Contents |
|---|---|
| `REPORT.txt` | Plain-language summary: start here |
| `images_needed.png/.pdf` | Predicted noise, subsampling agreement and per-n error vs images per condition (vs the main reference; `images_needed_vs_<group>` for the others) |
| `batch_by_n.csv` | Per n: reference-group levels, missing references, systematic well difference, and whether this n's treatment effects are stronger/weaker than in the other n's |
| `batch_by_condition.csv` | Per n × treatment: mean per well, well-to-well ratio with a flag when it exceeds field noise, ratio vs each reference |
| `variance_components.csv` | Variation field-to-field (within a well) and well-to-well (within a plate), as CV² |
| `effects.csv` | Per reference and treatment: observed fold change, p with all images, biological (n-to-n) variance |
| `design_power.csv` | For every design (1–2 wells × 1–6 fields), per reference and treatment: noise vs now, power at your current number of n's, and n's needed for 80% power |
| `subsampling.csv` | Real images re-drawn many times per design: how often the result matches the conclusion from all images, and how far each n's ratio moves |

Every analysis is done once per reference group (`reference` + `normalise_to`), like the statistics in
section 5.4: an n without LyoVec is left out of the comparisons vs LyoVec but still counts vs Control.

**How to read it:** the across-n test is limited by two things. The first is *measurement noise*
(from fields and wells), which shrinks with more images. The second is *n-to-n variation*,
which shrinks only with more n's. More images help only while measurement noise is still comparable to the
n-to-n variation. Well-to-well noise shrinks only with more wells.

**Batch effects:** each n is one plate on one day, so plate and day effects are part of the n-to-n
variation and cannot be separated from it. Effects that scale a whole plate (seeding, staining, lamp) cancel
in the ratio to the reference. Effects that change how strongly the cells respond (e.g. transfection
efficiency) do not cancel; the report flags an n whose treatment effects are consistently stronger or weaker
than in the others. A flag is a prompt to check the lab notes for that day; it is not a reason to exclude the n.

**Typical pattern:** when cells clump, field-to-field variation is large. Imaging at lower
magnification (10x covers about 16× the area of 40x) then reduces measurement noise far more than adding
more high-magnification fields. With only a few n's, treat the results as guidance for planning, not as
exact values.

---

## 6. Typical workflows
**A new n for an existing experiment**
1. Put the new n folder into `data/<experiment>/` (same structure as the others).
2. `kerncount init <experiment>`. Your `layout.csv` is kept and a fresh scan is written to `layout.new.csv`;
   copy the new rows into `layout.csv` (check the names).
3. `kerncount run <experiment>`. Only the new images are processed; the rest come from the cache.
   Existing generations are not overwritten. Rebuild them on purpose with `--rebuild` once the new n is checked.

**A new experiment**
1. Create `data/<new experiment>/` with one folder per n.
2. `kerncount init <new experiment>`, check `layout.csv`, then `kerncount run <new experiment>`.

**Then, for either:** check the red and orange tiles in `review/*_overview.jpg`, correct counts in
`counts.csv` → `manual_1`, and run again. Report from the latest generation: `generations/g<N>/results.xlsx`,
the `all_n/` plots, and the Prism CSVs.

## 7. Troubleshooting
| Problem | Fix |
|---|---|
| "project not found" / "raw folder not found" | Run the command from the workspace folder (the one containing `data/`), check the experiment name, and run `kerncount init <experiment>` first. |
| First run fails while downloading the model | The StarDist model downloads once on first use; connect to the internet and run again. |
| `counts.csv` opens as one column in Excel | Set `csv_delimiter = ";"` in `settings.toml` (European Excel). Files Excel re-saves with `;` are read correctly either way. |
| `run` fails to write `counts.csv` (Windows) | Close the file in Excel first; Excel locks open files. |
| "cannot read manual value" | A manual cell contains something other than a whole number or `x`; the message gives the row and column. |
| "g2 already exists but counts.csv now gives different values" | You changed a manual column after g2 was built. Re-run with `--rebuild g2` to update it on purpose. |
| A treatment name is wrong or a folder should be left out | Edit `layout.csv` (`treatment`, `include`) and run again. |
| "has 1 channel(s) but settings ask for channel 2" | `[channels] segment` points to a channel the images do not have. (Channels in `[intensity] measure` that do not exist are skipped silently; check the `intensity_*` columns.) |
| Run is slow | First runs take ~2–3 s per image; later runs are cached. `--no-review` skips the review images. |
| Counts look wrong in a whole folder | Check `diameter_px` in `qc_details.csv`; if the size estimate is off (e.g. almost empty folder), set `nucleus_diameter_px` to a number. |

## 8. Notes and known limitations
- 8-bit camera images saved as 16-bit TIFF (every value × 257) are handled. Saturation is detected either way.
- Compared with a simple Fiji threshold macro, counts are usually lower. Threshold macros count
  background noise in blank or defocused fields, split lobed nuclei, and count every edge-touching nucleus.
- Very dense fields where cells have bright cytoplasm are the hardest case. They are flagged QC low/medium;
  check them in the review images.
- Without a known pixel size, results are per field rather than per mm².
- QC levels are heuristics until calibrated against hand counts.
- Technical details (algorithms, tests): `docs/technical.md`.

## 9. Tumour experiments (mCherry tumours in brain sections)

The same tool also analyses **tumour experiments**: mCherry-labelled tumours in brain sections, imaged as z-stacks
(`[experiment] type = "tumour_volume"`). One **tumour** (hemisphere A or B of a section) is one n; A and B are
never averaged. Words: `CONTEXT.md` (Tumour experiments); technical details: `docs/tumour.md`.

### 9.1 Input
One folder per experiment in `data/`, with **one sub-folder per imaging session** (one day at the microscope),
files flat and as delivered:
```
data/tumour_TLR7/
├── 030826 KO7/                                 ← session folder: its name gives the genotype (KO7 = TLR7 KO)
│   ├── 030826 HERV -clodro n3 A.ome.tiff        ← one original z-stack = one tumour
│   ├── 030826 HERV -clodro n3 Af.tif            ← Fiji mask (legacy, optional; read for comparison only)
│   └── Statistics for 030826 HERV -clodro n3 Af.csv   ← Fiji 3D Objects Counter table (legacy, optional)
└── 100926 WT/
    └── 100926 cntrl +clodro n1 B.ome.tif
```
- **Genotype comes from the session-folder name**: it must contain `WT` or `KO7` / `KO` / `TLR7 KO` (other spellings
  can be added as aliases in the settings). WT and KO imaged on different days = genotype confounded with the
  imaging session: ratios within a session are compared across genotypes, a direct WT-vs-KO difference is shown
  but never tested.
- **File names** give date, tumour line, clodronate, section and hemisphere. For **new acquisitions** use
  `DDMMYY <line> <±clodro> n<section> <A|B>.ome.tif`, e.g. `100926 HERV +clodro n3 A.ome.tif`: date as 6 digits,
  the tumour line and clodronate level exactly as in the settings (`cntrl`, `HERV`, `N4`; `-clodro`, `+clodro`),
  the section label `n<k>` (restarts per group; `s3` is not read as a section), hemisphere `A` (upper) or `B` (lower).
  Underscores instead of spaces (`100926_HERV_+clodro_n3_A.ome.tif`), `.ome.tiff` and known typos (`-clodron`) are
  read too and noted in `layout.csv`.
- Voxel size is read from each file's OME metadata (never typed in); z-steps may differ between files.

### 9.2 Set up, check the layout, run
In the web interface pick the experiment, choose **experiment type: Tumour volume (z-stacks)** under *Set up* and
press **Set up** (command line: `kerncount init tumour_TLR7 --type tumour_volume`). Then:
- **Layout**: one row per tumour with its factor values (genotype, tumour line, clodronate), section, hemisphere,
  tumour ID (`TLR7 KO / HERV / -clodro / n3A`), include and note. Problems that would make a file a wrong datapoint
  (an unreadable name, a duplicate tumour ID, a session without a genotype, a missing file) are listed at the top
  and **block the run** until you fix the cell or untick *include* and give the reason in *note*.
- **Settings** → *Analysis plan*: the factors (levels, reference level, where they go in the plots), which readouts
  are confirmatory, and the contrast families (below). *Tumour measurement*: smoothing, threshold rule,
  mid-section thickness, outline smoothing, threshold sweep, QC and review settings, each with an explanation.
- **Run** (`kerncount run tumour_TLR7`): options *skip review images*, *segment every stack again*, *threshold
  sweep*. The first run takes about 20–60 s per stack; segmentations are cached, so later runs are quick.

### 9.3 Readouts, and why area rather than volume
Only ~80–150 µm of each ~0.5–1 mm tumour is imaged, so a 3D volume mostly measures how deep the stack was acquired.
The **confirmatory readouts** are therefore 2D, on a fixed-thickness **mid-section** (20 µm around the layer where
the tumour's main mass is largest — the same physical slab for every tumour): the **cross-sectional area** (mm²) and
the **perimeter** (mm) of the main mass (the largest connected object). Everything else is **exploratory**
(effect sizes with CIs, labelled, no stars): footprint area, slab volume, main-mass volume and surface, tumour
burden (main mass + satellites), invasion fraction, elongation, solidity, edge roughness, circularity, and the Fiji
legacy values. One automatic threshold rule for every tumour, never adjusted by hand per tumour.

**Analysis plan lock.** Factors, readout roles and contrast families are fixed before results are viewed (ADR 0003).
Once results exist the plan is **locked** in the web interface; changing it needs **Unlock** with a reason. The
unlock and every change are written to `lock_log.csv`, and the next run prints the history (old → new plan lock) in
the report and the About sheet of `results.xlsx`. Every statistics row carries the lock it was built with.

### 9.4 Output (`analysis/<experiment>/`)
| file / folder | content |
|---|---|
| `results.xlsx` | everything in one workbook: tumours, statistics, estimated means, families, QC, flow, sensitivity analyses, concordance with Fiji, planning aid (download it on the **Results** page) |
| `stats_report.md` | the report: analysis plan and its lock (and lock history), results per readout and family, QC, sensitivity, Fiji concordance, planning aid |
| `tumours.csv` | one row per analysed tumour: factor values, every readout, QC flags |
| `excluded_tumours.csv`, `flow.csv` | excluded tumours with their reason (layout or review); files → tumours → analysed |
| `stats_long.csv` | every estimate in one long table (the **Statistics** page; sensitivity rows hidden unless asked for) |
| `plots/` | per readout: absolute, normalised (% of cntrl in the same genotype and clodronate), effects (forest), interaction, A vs B, diagnostics, and with Fiji data / the sweep: vs Fiji, threshold sweep; `shape_panel` |
| `prism/` | Prism-ready tables, one column per group (absolute and normalised) |
| `review/` | the **blinded** review images and `corrections.csv` (below) |

Statistics: one linear model on log(readout) over all factors; effects are fold-changes with 95 % CIs, Holm within
each family. Only confirmatory families of confirmatory readouts carry stars. Sensitivity analyses (Welch-type,
permutation, without QC-flagged tumours) are reported alongside, never instead.

### 9.5 Review and corrections (blinded)
`review/` holds one image per tumour, named by a random **code** (no group, session, date or file name anywhere;
the key `review_key.csv` lies outside `review/`): left the maximum projection, right the mid-section layer; the main
mass orange, satellites cyan, objects numbered. The **Review images** and **Corrections** pages show only the codes.
In a round (`manual_1` → g2, `manual_2` → g3, …) type:
- an **object number** from the image to make that object the main mass (`1` = back to automatic);
- **`x: reason`** to exclude the tumour (the reason goes to `excluded_tumours.csv`);
- nothing = no change.

Then run the analysis again: nothing is re-segmented; the picked object's readouts are recomputed. You can also edit
`review/corrections.csv` in Excel and hand only the `review/` folder to a colleague for a blinded check.

## 10. Development
```
pip install -e ".[test]"
pytest
```
Algorithm and implementation notes: `docs/technical.md`.

## Acknowledgements
Nucleus detection uses [StarDist](https://github.com/stardist/stardist) (Schmidt et al., MICCAI 2018;
Weigert et al., WACV 2020) with its pretrained `2D_versatile_fluo` model.
