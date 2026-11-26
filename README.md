# kerncount: automated nucleus counting and analysis

*[Deutsche Fassung](README.de.md)*

`kerncount` turns folders of fluorescence microscope images into checked counts, tables, plots and statistics.
It was written to replace a lab's manual workflow (ImageJ/Fiji macros to count nuclei and measure intensity, then
copying numbers into GraphPad Prism), and it keeps that workflow's conventions: one folder per experiment day,
Prism-ready tables, and comparisons with the old Fiji results where they exist. You do not need Fiji or Prism to
use it.

It does the following:

1. Counts nuclei in every image (DAPI or another nuclear stain), with a deep-learning detector (StarDist).
2. Estimates how reliable each count is and draws a review image for each field.
3. Lets you correct counts by hand in one CSV file (or in the browser), keeping every round of corrections as a
   "generation".
4. Produces tables, Prism-ready CSVs, plots and statistics for each treatment, each n, and all n's together, with a
   hierarchical (mixed) model as the main test.
5. Measures a marker channel per nucleus (e.g. Ki-67): % positive nuclei, intensities, flat-field correction, and a
   "layer comparison" that shows which analysis step changes the conclusion.
6. Adds extra per-image readouts from the detected nuclei (out-of-focus fraction, micronuclei, nuclear size, ...).
7. Plans the next experiment: how many n's, wells and images you need, and randomised plate maps.
8. Measures fluorescent tumours in brain-section z-stacks (tumour experiments, section 9), with blinded review, a
   locked analysis plan and pre-registered statistics.

Everything can be done from a point-and-click page in your browser, or from the command line.

## Contents

1. [Input: how the raw data must be organised](#1-input-how-the-raw-data-must-be-organised)
2. [Install](#2-install-once-per-computer)
3. [How to use it](#3-how-to-use-it)
4. [Settings (`settings.toml`)](#4-settings-settingstoml)
5. [Output: what comes out and how to read it](#5-output-what-comes-out-and-how-to-read-it)
6. [Typical workflows](#6-typical-workflows)
7. [Troubleshooting](#7-troubleshooting)
8. [Notes and known limitations](#8-notes-and-known-limitations)
9. [Tumour experiments](#9-tumour-experiments-fluorescent-tumours-in-brain-sections)
10. [Development](#10-development)

```
kerncount/                 (this repository)
├── README.md              ← this guide (README.de.md = German version)
├── docs/                  ← detailed documentation (see section 10)
├── src/kerncount/         ← the code
├── tests/                 ← tests (pytest)
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

Words used throughout: an **n** is one independent experiment (one plate, imaged on one day); a **treatment**
(condition) sits in one or more **wells** of that plate; each well is imaged in several **fields** (images).

---

## 1. Input: how the raw data must be organised

This section is for cell-count experiments (the default). Tumour experiments are organised differently (section 9).

Each experiment is one folder in `data/`:

```
data/<experiment>/                              e.g. data/my_experiment/
├── <one folder per n (experiment)>/            e.g. "12.01.26 DAPI count"
│   ├── <one folder per treatment>/             e.g. "Control", "LyoVec veh. cntrl", "miRNA 12-5p"
│   │   ├── 1.tif
│   │   ├── 2.tif
│   │   └── ...
│   └── ...
└── ...
```

- **One folder per n**, with the date in the name (`12.01.26`, read as day.month.year; see `date_order`). The n's
  are numbered by date: `n1_2026-01-12`, `n2_2026-02-09`, and so on. A single n works fine. A name with only day
  and month gets the year from the file dates (with a note to check it). A timepoint in the name (`72h`, `48 h`,
  `day 3`) is read into the `timepoint` column of `layout.csv`.
- **One folder per treatment** inside each n. Messy names are cleaned up automatically, and you can check the
  mapping in `layout.csv`. Built-in rules: `LyoVec` (e.g. `lyo vec`, `LyoVec veh. cntrl`), `CL075`, microRNAs
  (`hsa miRNA 12-5p` / `12 5p` → `miR-12-5p`; a bare `12` borrows its arm from other n's), and controls
  (`cntrl`, `ctrl`, `control` → `Control`). Anything unrecognised is kept as-is and flagged.
- **Images are named by number** (`1.tif`, `02.tif`; stray characters such as `` `10.tif `` are fine).
  Files without a number (e.g. `5 zoom.tif`) are skipped and listed in `layout.csv`.
  Any number of images per folder is fine.
- **Ignored:** sub-folders called `output`, files containing `mask` (old Fiji results), and other
  sub-folder levels (images are found at any depth). Old Fiji count tables in a treatment folder are read for
  comparison (`fiji_comparison.csv`).
- **Image format:** TIFF, 8- or 16-bit, greyscale, RGB or multi-channel, at any magnification and from any camera.
  Nucleus size is measured per folder and the detector rescales to it.
- **Channels saved as separate files** (e.g. `1.tif` = marker, `2.tif` = DAPI, alternating) are supported:
  set `[channels] files = "interleaved"` (section 3, *Measuring a marker channel*). `kerncount init` notices files
  that alternate dim / bright and suggests this.
- **Wells:** one n = one plate. By default each condition sits in 2 wells, and images 1–6 = well 1, 7–12 = well 2
  (worked out from the highest image number; change it in `[design]`, or list it per image in `wells.csv`).

---

## 2. Install (once per computer)

### Easiest: one-click app (no admin rights, no Python, no commands)
Download from the [latest release](https://github.com/joe-bones/kerncount/releases/latest):
- **Windows:** `kerncount.exe`
- **Mac (Apple silicon or Intel):** `kerncount-mac.zip`. Double-click to unzip; this gives `kerncount.app`. On Intel Macs it installs TensorFlow 2.16, the last version built for Intel Macs.

Put it in your workspace folder (next to `data`) and double-click it. The first start needs internet, takes
a few minutes and about 1 GB: it installs Python 3.12 and everything kerncount needs into your user folder.
After that, every start opens kerncount in your browser straight away. The first time, the system warns about an
unsigned app: on Windows click **More info → Run anyway**; on the Mac click **Done**, then
**System Settings → Privacy & Security → Open Anyway**. Updating = replacing the app with a newer one.
Details (where it installs, repair with `--reinstall`): [launcher/README.md](launcher/README.md).

A launcher contains the code of the release it was downloaded from. Check the release notes; for changes made after
the latest release, install from the repository (below).

### Alternative: install from this repository (for the command line)
1. **Mac/Linux:** install **Python 3.12** from python.org (`setup.command` looks for `python3.12`). kerncount itself
   needs Python 3.11 or newer; 3.12 is recommended because of TensorFlow.
   **Windows:** nothing to install first, and no admin rights needed: `setup.bat` fetches Python 3.12 itself
   (via [uv](https://docs.astral.sh/uv/)) into your user folder (`%LOCALAPPDATA%\kc`).
2. Download this repository (green **Code** button → *Download ZIP*, or `git clone`).
3. Double-click `setup.bat` (Windows) or run `setup.command` (Mac/Linux).
   This takes a few minutes and about 1 GB (mostly TensorFlow). The StarDist model (~5 MB) downloads automatically
   on the first run, so that run needs internet.

Or by hand, in any Python ≥ 3.11 environment: `pip install .` (or `pip install -e .` to edit the code).

Commands are run from your workspace folder. On Windows, every command below starts with the full path
to `kerncount.bat` in this folder instead of `kerncount` (e.g. `C:\tools\kerncount\kerncount.bat run my_experiment`).
On Mac/Linux it is `<this folder>/.venv/bin/kerncount`. Adding that folder to your PATH lets you type just `kerncount`.

If TensorFlow will not load (on Windows usually a missing Microsoft Visual C++ runtime; `setup.bat` says so),
kerncount still opens and the `classical` detector (`[segmentation] backend = "classical"`) works.

---

## 3. How to use it

Run every command from the workspace folder (the one that contains `data/`). Experiments are referred
to by their folder name in `data/`.

### Without commands: the web interface
Double-click **`kerncount.exe`** / **`kerncount.app`** (section 2), or, with a repository install,
**`kerncount-gui.bat`** (Windows) / **`kerncount-gui.command`** (Mac). The same thing from a terminal: `kerncount gui`
(options: a workspace folder, `--port` (default 8765), `--no-browser`). If the page is not started inside a
workspace, it lets you choose or create one.

The left sidebar lists the experiments; the open experiment expands into these pages:

| Page | What you do there |
|---|---|
| **Experiment home** | Set up the experiment (with its **experiment type**), run the analysis (the command options below as tick boxes), rebuild generations, run the sampling analysis, write a wells template; the **Attention** panel (duplicates, missing conditions, failed checks, settings changed since the last run, ...); the n's, the design summary and the generations. Live log and **Stop** at the bottom. |
| **Results explorer** | All plots per generation, by level (all n / per n / per treatment), metric and view (absolute, relative, model, comparisons, robustness, diagnostics, design, tables); download `results.xlsx` |
| **Statistics table** | Every test result (`stats_long.csv`) as a filterable, sortable table; export the filtered rows |
| **Layer comparison** | (marker channels only) Which analysis step changes the conclusion (section 5.4) |
| **Intensity & cut-off** | (marker channels only) Positivity cut-offs, flat-field maps, distributions, channel checks, and a **blinded tuning view** to choose and lock a cut-off |
| **Review images** | Overview sheet per treatment, then each image full size (← / → to browse); type a correction directly |
| **Corrections** | Edit a generation: type corrected counts or **x** for a round (`manual_1` → g2, …), see the value before and after, then **build** the generation (about a minute; detection is reused) |
| **Planner** | Plan the next experiment interactively: power, CI width and effort for any number of n's, wells and fields |
| **Plate maps** | Randomised plate maps for the next experiment (PNG/PDF/CSV) |
| **Layout** | Fix n / treatment names, dates and timepoints; include or leave out folders |
| **Settings** | Every setting from section 4 as a form, with explanations; or edit `settings.toml` as text. After saving, the page says what needs re-running, with a Start button. |
| **All files** | Browse and preview everything in the project folder |

Every view has its own link (in the address bar), so a bookmark reopens exactly that view. Tumour experiments
(section 9) use the same pages with their own content; Planner and Plate maps are not shown for them.

Keep the black (Windows) or Terminal (Mac) window open while you work; close it or press **Quit** to stop.
The page only works on this computer, and nothing is sent anywhere else. Anything that overwrites a file you may
have edited asks first, and the old file is copied to `_work/backups/`.

### Step 1: create a project (once per experiment)
```
kerncount init my_experiment
```
This scans `data/my_experiment/` and creates `analysis/my_experiment/` with two files
(`kerncount init --all` does this for every experiment in `data/`):

- **`layout.csv`** has one row per raw folder. Check the `n`, `date`, `timepoint` and `treatment` columns. To fix a
  name, edit it; to leave a folder out, set `include` to `no`. The `note` column shows anything unusual
  (skipped files, inferred names, unpaired channel files). You can add columns named like the `[experiment]`
  settings (e.g. `cell_line`) to override them per n or per treatment.
- **`settings.toml`** holds all the settings (section 4). The defaults are sensible; fill in the `[experiment]`
  description (cell line, magnification, stain, ...) if you like, and set `reference` / `normalise_to` to your
  control groups' names.

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
(left: original; right: what was counted). Start with the red ones. Also look at `design_check.csv` and
`duplicates.csv` (if present) or the **Attention** panel: they flag missing conditions, odd well sizes and
byte-identical images in different conditions.

### Step 4: correct counts by hand (optional)
Open **`analysis/my_experiment/counts.csv`** in Excel (or use the **Corrections** page):

| n | treatment | image | auto_count | qc_ok | qc_level | manual_1 | manual_2 | manual_3 |
|---|---|---|---|---|---|---|---|---|
| n1_2026-01-12 | Control | 3 | 85 | yes | high | | | |
| n1_2026-01-12 | miR-12-5p | 6 | 140 | no | low | **121** | | |
| n1_2026-01-12 | miR-12-5p | 9 | 4 | no | low | **x** | | |

- In `manual_1`, type a corrected count, or **`x`** to exclude that image. Leave all other cells empty.
- Save, then `kerncount run my_experiment` again (or the faster `kerncount regen my_experiment`). This creates
  generation **g2** (see section 5.3).
- A later round of corrections goes into `manual_2`, which makes g3, and so on.

### Step 5: how many images do we need, and what next? (optional)
```
kerncount sampling my_experiment
```
Results are written to `analysis/my_experiment/sampling/g<N>/` (section 5.5): a look back ("could we have used
fewer images?") and a plan for the next experiment (how many n's, wells and fields), which the **Planner** page
shows interactively. It needs at least 2 n's with the reference group. Then, for randomised plate maps:
```
kerncount platemap my_experiment --plates 3 --format 24
```
This writes `analysis/my_experiment/plans/platemap_<date-time>/` with one map per plate (PNG, PDF, CSV) and a
`README.txt` with the seed and the command that recreates the same maps. Details: [docs/planner.md](docs/planner.md).

### Measuring a marker channel (e.g. Ki-67 in C2)
For multi-channel TIFFs (C1 = DAPI, C2 = marker), set in `settings.toml`:
```
[channels]
names = ["DAPI", "Ki67"]
segment = "DAPI"       # find nuclei in DAPI (a name or a number)
```
If each channel is saved as its own file, alternating (1.tif = Ki-67, 2.tif = DAPI, 3.tif = Ki-67, ...):
```
[channels]
files = "interleaved"
file_order = ["Ki67", "DAPI"]
segment = "DAPI"
```
then rewrite the layout with `kerncount init my_experiment --force-layout`. Each pair of files is checked for
shifts and mis-pairing (`intensity/alignment.csv`).

Then run as usual. Every channel other than the nuclear one is measured **per nucleus** (`[intensity]
per_nucleus`), corrected for uneven illumination (flat field) and local background, and classified as positive or
negative (`[positivity]`). The headline readout is **% positive nuclei** per image; mean intensities are secondary.
The cut-off (per n and channel) comes, in this order, from negative-control wells in that n (`[experiment]
negative_control`, e.g. secondary-only), a value you set (`manual_cutoff`) or locked in the browser's blinded
tuning view, or a limit-of-detection rule. Without a negative control, "positive" means "detectably stained above
background". Outputs go to `intensity/` and through the same statistics as the counts. Details:
[docs/intensity.md](docs/intensity.md).

The older region-based means still exist: `[intensity] regions = ["nucleus", "ring", "cell", "fiji_cell"]` gives
`intensity_<channel>_<region>` columns (`fiji_cell` = nucleus enlarged by 6 px, as the old Fiji macro). To segment
whole cells from a cytoplasmic channel instead of nuclei, set `segment` to that channel and
`[segmentation] backend = "classical"`.

### Wells and field positions (optional)
kerncount needs to know which well each image comes from, because well-to-well and field-to-field variation are
treated separately. The default rule (images 1–6 = well 1, 7–12 = well 2) is set in `[design]`. For anything the
rules cannot describe:
```
kerncount wells-template my_experiment
```
writes `wells.csv` (one row per image: `well`, `position`; or `wells.new.csv` if one exists). Edit it and set
`[design] well_layout = "file"`. Details: [docs/design.md](docs/design.md).

### Command options
| Command | What it does |
|---|---|
| `kerncount init <experiment>` | Create `analysis/<experiment>/` from `data/<experiment>/`. An existing `layout.csv` is kept and a fresh scan written to `layout.new.csv`. |
| `kerncount init --all` | Same for every experiment in `data/` (existing projects are kept) |
| `kerncount init <experiment> --force` | Re-scan and overwrite `layout.csv` and `settings.toml` (old copies go to `_work/backups/`) |
| `kerncount init <experiment> --force-layout` | Keep `settings.toml`, rewrite `layout.csv` (e.g. after changing `[channels] files`) |
| `kerncount init <experiment> --type tumour_volume` | Create a project of another experiment type (default `cell_count`). With `--all`, experiments that already have another type are reported and skipped. |
| `kerncount run <experiment>` / `--all` | Count, check and analyse one / every experiment |
| `kerncount run <experiment> --rebuild g2` | Overwrite an existing generation folder (never done by default). `--rebuild all` for all. |
| `kerncount run <experiment> --resegment` | Ignore the cache and re-detect every nucleus |
| `kerncount run <experiment> --no-review` | Skip drawing review images (faster) |
| `kerncount run <experiment> --unlock "reason"` | Tumour experiments: accept a changed analysis plan, with a logged reason (section 9.3) |
| `kerncount regen <experiment> [--rebuild g2 g3]` | Rebuild generations from `counts.csv` only (fast; after corrections or a statistics setting change) |
| `kerncount sampling <experiment> [--generation g2] [--iterations 1000]` | Sampling analysis and planner on a generation (default: latest) |
| `kerncount wells-template <experiment> [--force]` | Write `wells.csv` from the current well assignment (`--force` overwrites an existing one) |
| `kerncount platemap <experiment>` | Randomised plate maps. Options: `--plates` (default 3), `--format` 6/12/24/48/96, `--wells` per treatment (default 2), `--reference-wells` (a number or `auto`), `--avoid-edge`, `--seed`, `--no-balance`, `--negative-control NAME`, `--control-wells` 1/2 |
| `kerncount gui [workspace] [--port N] [--no-browser]` | Open the web interface |

`regen`, `sampling`, `wells-template` and `platemap` (and `--all` on each) are for cell-count experiments.
Paths work too, for folders outside a workspace: `kerncount init <raw folder> <project folder>`,
`kerncount run <project folder>`. Setting the environment variable `KERNCOUNT_WORKSPACE` lets you run
by name from anywhere. `kerncount <command> --help` lists every option.

---

## 4. Settings (`settings.toml`)

Edit with any text editor or on the **Settings** page. A setting you delete falls back to its default. The tables
below are for cell-count experiments; tumour experiments have their own settings (section 9).

### `[project]`
| Setting | Default | Meaning |
|---|---|---|
| `raw_dir` | (set by init) | Where the raw images are. Stored relative to the project folder (`../../data/<experiment>`), so the whole workspace can be moved or copied. |
| `date_order` | `"DMY"` | How dates in n-folder names are written (`DMY`, `MDY`, `YMD`) |
| `csv_delimiter` | `","` | Use `";"` if your Excel (German/European locale) shows CSVs in one column. Reading auto-detects either. |

### `[experiment]`: what was imaged
| Setting | Default | Meaning |
|---|---|---|
| `type` | `"cell_count"` | Experiment type: `cell_count` (this guide) or `tumour_volume` (section 9). Set by `init --type`. |
| `cell_line` | `""` | e.g. `"LN229"`. n's of different cell lines are analysed separately, never pooled. |
| `magnification`, `camera_pixel_um` | `""` | Objective (e.g. `20`) and camera pixel size in µm; together they give the pixel size |
| `seeding_density`, `seeding_unit` | `""`, `"cells/well"` | Recorded in the outputs |
| `stain`, `notes` | `""` | Free text, recorded in the outputs |
| `negative_control` | `""` | Name of a negative-control treatment (e.g. secondary-only wells). Used only to set positivity cut-offs; left out of every statistic. |
| `negative_control_type` | `"secondary_only"` | `secondary_only`, `isotype` or `unstained` |

Any of these can be overridden per n or per treatment with a same-named column in `layout.csv`.

### `[channels]`
| Setting | Default | Meaning |
|---|---|---|
| `names` | `[]` | Channel names, e.g. `["DAPI", "Ki67"]`. Empty = C1, C2, ... |
| `segment` | `1` | Channel used to find nuclei: a number (1 = first) or a name |
| `files` | `"single"` | `single` = all channels in one file; `interleaved` = one file per channel, alternating |
| `file_order` | `[]` | `interleaved`: channel names in file order, e.g. `["Ki67", "DAPI"]` |
| `align` | `"check"` | `interleaved`: `check` (report shifts between a field's files), `correct` (also shift small offsets back), `off` |
| `align_warn_px` | `2.0` | Warn when a field's channels are shifted by more than this |
| `align_max_correct_px` | `10.0` | `correct` only fixes shifts up to this; larger ones are flagged |

### `[segmentation]`: finding nuclei
| Setting | Default | Meaning |
|---|---|---|
| `backend` | `"stardist"` | `stardist` = deep-learning detector (recommended); `classical` = Fiji-style threshold + watershed |
| `nucleus_diameter_px` | `"auto"` | `auto` measures nucleus size per folder (adapts to magnification). Or a fixed number. |
| `target_diameter_px` | `30` | Size the detector works best at; images are rescaled to it internally |
| `prob_thresh` | `"default"` | Detector confidence threshold (lower = more objects) |
| `pixel_size_um` | `""` | µm per pixel, if known (from a stage micrometer). Adds µm² and cells/mm². Empty = `camera_pixel_um / magnification` if both are set. Changing it does not trigger re-detection. |

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

### `[intensity]`: intensity measurement
| Setting | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Measure intensities |
| `measure` | `"all"` | Channels for the region means, e.g. `[2]` for channel 2 only |
| `regions` | `["nucleus", "cell"]` | `nucleus`; `ring` (band around the nucleus, i.e. cytoplasm); `cell` (nucleus + ring); `fiji_cell` (nucleus + `fiji_ring_px`, the old Fiji macro) |
| `ring_px` | `"auto"` | Band width; `auto` = 0.25 × nucleus diameter |
| `subtract_background` | `true` | Subtract each image's background (median of empty areas) |
| `per_nucleus` | `"auto"` | Marker channels measured per nucleus: `auto` = every channel except the nuclear one; a list; or `"off"` |
| `erode_px` | `"auto"` | Shrink the nucleus outline before measuring (`auto` = 1 px, 2 px for nuclei ≥ 16 px across) |
| `bg_gap_px`, `bg_ring_px`, `bg_min_px` | `3`, `"auto"`, `40` | Local background ring per nucleus: gap to the nucleus, width (`auto` = 0.5 × diameter), minimum pixels before falling back to a wider ring / the image background |
| `fiji_ring_px` | `6` | Enlargement for the `fiji_cell` region |
| `flatfield` | `"auto"` | Uneven-illumination correction: `auto` (estimated per n and channel from that n's images), `off`, or a folder of reference images |
| `flatfield_dark` | `true` | Remove the additive background pattern (camera offset, autofluorescence) before estimating the flat field |

### `[readouts]`: extra per-image readouts
| Setting | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Add `readout_*` columns computed from the detected nuclei (no re-detection) |
| `include` | `"all"` | `all`, or a list such as `["frac_out_of_focus", "coverage"]` |

The readouts are `frac_out_of_focus`, `micronuclei_per_nucleus`, `nuclear_area_median` (and `_um2`),
`nuclear_area_cv`, `frac_irregular`, `coverage` and `nn_distance_median`; they help tell cell death, damage and
reduced growth apart. Details: [docs/readouts.md](docs/readouts.md).

### `[positivity]`: % positive nuclei per marker channel
| Setting | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Classify nuclei as positive / negative |
| `channels` | `"auto"` | `auto` = every per-nucleus channel, or a list of names |
| `control_percentile` | `99.0` | With negative-control wells: cut-off = this percentile of their nuclei |
| `manual_cutoff` | `""` | A number (all channels) or e.g. `{ Ki67 = 1200 }`. Empty = a cut-off locked in the browser, else automatic. Negative-control wells, if present, take precedence. |
| `lod_k` | `3.0` | Automatic (limit-of-detection) rule: positive when the signal is above k × the background noise |
| `min_nuclei` | `5` | Images with fewer measurable nuclei get no % positive |

### `[analysis]`: statistics
| Setting | Default | Meaning |
|---|---|---|
| `reference` | `"Control"` | Main comparison group (headline plots) |
| `normalise_to` | `["Control", "LyoVec"]` | Control groups to express results as % of (within each n), and to test against. Every group listed gets its own tests and plots; a listed group that is not in your data gets none. |
| `treatment_order` | `["Control", "LyoVec", "CL075"]` | Order in tables and plots; other treatments follow alphabetically |
| `exclude_low_qc` | `false` | If `true`, images with QC low (and no manual value) are left out of the statistics |
| `two_way_anova` | `true` | Extra two-way ANOVA with Tukey all-pairs comparisons: per n treatment × well, all n treatment × n (5.4.3) |
| `primary_test` | `"mixed"` | What draws the stars on the all-n plots: `mixed` (hierarchical model, log scale) or `rm_anova` (Prism's repeated-measures ANOVA + Dunnett) |
| `primary_family` | `""` | The control whose comparisons are the main claim; empty = `reference` |
| `vehicle_check` | `"auto"` | Equivalence test for a vehicle control, e.g. `"LyoVec vs Control"`; `auto` = every other control vs `reference`; `""` = off |
| `equivalence_margin` | `0.2` | Vehicle check: the ratio must lie within 1/1.2 – 1.2 to be called equivalent |
| `equivalence_margin_pp` | `5.0` | Vehicle check for % positive: difference within ± 5 percentage points |
| `positive_control` | `""` | e.g. `"CL075"`: every non-control treatment is also compared with it |
| `model_unit` | `"well"` | Unit of the mixed model: `well` (well means; ratios of mean nuclei per field) or `field` |
| `classes` | `"auto"` | Groups of treatments that get their own n-to-n variance: `auto` = one per treatment; `{}` = controls vs everything else; or e.g. `{ mirna = ["miR-*"] }` |
| `density_covariate` | `"readout_coverage"` | For % metrics, also give an estimate adjusted for cell density; `""` = off |
| `among_classes` | `[]` | Classes whose members are also compared with each other (all pairs), e.g. `["mirna"]` |

### `[acceptance]`: checks per n
Checked on the counts of every n and written to `acceptance.csv` in each generation. By default they are only
reported.

| Setting | Default | Meaning |
|---|---|---|
| `vehicle_within` | `0.3` | Vehicle / control in that n must lie within 1 ± 0.3 (`""` = off) |
| `positive_control_max_ratio` | `""` | e.g. `0.8`: the positive control / reference must reach this ratio (`""` = off) |
| `min_images_per_condition` | `3` | Every condition needs at least this many usable images |
| `max_low_qc_fraction` | `0.5` | At most this share of the n's images may be QC low |
| `no_duplicate_images` | `true` | An n fails if one of its raw images is byte-identical to an image of another condition |
| `enforce` | `false` | `true`: n's that fail are left out of the all-n statistics (and named in the outputs) |

### `[design]`: wells and fields
| Setting | Default | Meaning |
|---|---|---|
| `wells` | `2` | Wells imaged per condition per n (one plate per n); `"auto"` = from `wells.csv` if present, else 2 |
| `fields_per_well` | `"auto"` | `auto` = per condition, highest image number ÷ wells (rounded up), so 18 images = 9 + 9; or a fixed number |
| `well_layout` | `"consecutive"` | `consecutive` = images 1–6 well 1, 7–12 well 2; `alternating` = odd/even; `file` = read `wells.csv` |
| `field_positions` | `"fixed"` | `fixed` = field k is the same spot in every well (position effects are modelled); `random` = position means nothing |

Older settings files using `plates` / `fields_per_plate` / `plate_layout` still work; the report asks you to rename them.

### `[review]`
| Setting | Default | Meaning |
|---|---|---|
| `overlay_alpha` | `0.35` | Transparency of the coloured overlay |
| `jpeg_quality` | `88` | Review image quality |

### `[planner]`: planning the next experiment
| Setting | Default | Meaning |
|---|---|---|
| `effects` | `[0.8, 0.7, 0.5]` | Treatment / reference ratios to plan for (0.7 = 30% fewer nuclei) |
| `target_power` | `0.8` | Chance of a significant result if the effect is real |
| `alpha` | `0.05` | Significance level |
| `cost_per_n`, `cost_per_well`, `cost_per_field` | `240`, `5`, `1` | Effort per n, per extra well and per image (any unit, e.g. minutes); used to find the cheapest design |
| `max_n` | `12` | Largest number of n's in the design tables |
| `plate_format` | `24` | 6, 12, 24, 48 or 96 wells |
| `avoid_edge` | `true` | Plate maps: inner wells first, edge wells only when needed |
| `reference_wells` | `"auto"` | Wells per reference condition; `auto` = Dunnett's rule (wells × √comparisons) |
| `adjust_for_family` | `true` | Power = chance of a star on the headline plot (adjusted for all comparisons); `false` = one unadjusted test |

After a settings change, the **Settings** page tells you what needs re-running: detection settings re-detect every
nucleus (slow); most others need a run with `--rebuild all`; `[analysis]` / `[acceptance]` only need
`regen --rebuild all`; `[planner]` only the sampling analysis.

---

## 5. Output: what comes out and how to read it

```
analysis/<experiment>/
├── settings.toml, layout.csv     your settings and folder mapping (layout.new.csv = a fresh scan to compare)
├── counts.csv                    ★ THE file to work from: one row per image, all n's, manual columns
├── qc_details.csv                every per-image measurement (see 5.2)
├── design_check.csv              one row per n × treatment: wells, fields, missing images, flags
├── duplicates.csv                raw images whose content is identical (if any)
├── fiji_comparison.csv           old Fiji count vs new count, per image (only if old Fiji results were found)
├── wells.csv                     (optional) which well / position each image is
├── review/                       pictures to check the counting (5.1)
├── intensity/                    marker channels: cut-offs, flat field, distributions, channel checks
├── generations/g1, g2, ...       results (5.3, 5.4)
├── sampling/g1/                  how many images are needed, and the planner (5.5)
├── plans/                        plate maps
└── _work/                        cache, cells.csv and nuclei_intensity.csv (per-cell / per-nucleus values),
                                  backups/ of every file kerncount overwrote (deleting _work forces re-detection)
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
  "this count is wrong".** Low-QC images tend to cluster in sparse or damaged wells, where counting really is
  harder. `legacy` = no raw image, only an old Fiji count.
- `n_debris`, `n_out_of_focus`, `n_edge`, `n_merged`: how many objects of each type were left out or flagged.
- `diameter_px`: typical nucleus size in that folder (useful as a magnification/biology check).
- `fiji_count`: what the old macro gave, for comparison (when available).
- `intensity_<channel>_<region>`: mean background-subtracted intensity of the counted cells.
- `readout_*`: the extra readouts (section 4, `[readouts]`); `readout_frac_positive_<channel>`: % positive nuclei.
- `channel_check`: result of the per-field channel checks (separate channel files).
- `nuclei_per_mm2`: only present once the pixel size is known.
- `stardist_scale`, `saturated_frac`, `source_file`: detector rescaling factor, fraction of saturated
  nuclear pixels, and the raw file(s) the row came from.

*About the QC levels:* they are rules of thumb until checked against hand counts (ground truth).
Hand-counting about 15 images across densities would turn them into real error estimates.

### 5.3 Generations: `generations/g1`, `g2`, ...
- **g1** = fully automatic counts. Always present.
- **g2** = g1 + the corrections typed into `manual_1`. **g3** = g2 + `manual_2`, and so on.
- The chain stops at the first manual column that is completely empty.
- If counts.csv no longer matches an existing generation (you edited `manual_1` after g2 was built), the folder
  is **not overwritten**: `run` warns you, and `--rebuild g2` updates it on purpose. If only the intensity
  results, the design or the experiment description changed (counts unchanged), the generation is rebuilt
  automatically and the log says why.
- A generation is built in `gN.tmp/` and swapped in only when complete. Files you put into a generation folder
  yourself are moved to `gN/_kept/` on a rebuild, never deleted.
- Each generation has `README.txt` (what it is built from), `counts_used.csv` (the value used for every image
  and where it came from: `auto`, `manual_1`, `manual_1: excluded`, ...; plus well, position and the experiment
  description), `generation.json`, `acceptance.csv` (5.4.4), `stats_long.csv` (every test result in one table) and
  `manifest.json` (a list of every plot and table, which the browser reads).

### 5.4 Inside a generation: treatment → n → all n
```
g1/
├── results.xlsx                      ★ everything below in one readable workbook
├── counts_used.csv, stats_long.csv, acceptance.csv, manifest.json
├── <n>/                              one experiment
│   ├── <treatment>/fields.csv        per-field values (count + intensities + readouts)
│   ├── <treatment>/count_fields.png  bar per field (x = excluded)
│   ├── count_summary.csv             per treatment: fields, mean, SD, SEM, median, % of each control
│   ├── count_prism.csv               Prism Column table: columns = treatments, rows = fields
│   ├── count_stats.csv               one-way ANOVA + Dunnett vs each control
│   ├── count.png / .pdf              bar = mean, error = SD, dots = fields, stars = Dunnett vs reference (5.4.2)
│   ├── count_vs_<control>.png / .pdf same, vs each other control
│   ├── count_fields_by_well.png, count_position_effect.png   design checks
│   └── count_twoway.png / .pdf       two-way ANOVA treatment × well + Tukey all pairs (5.4.3)
└── all_n/                            all experiments combined
    ├── count_summary.csv             per treatment: value per n, N, mean, SD, SEM, % of controls
    ├── count_prism.csv               Prism Column table: columns = treatments, rows = n
    ├── count_prism_pct_of_<control>.csv   same, normalised to each control within each n
    ├── count_stats.csv               RM one-way ANOVA + Dunnett vs each control; ratio paired t-tests (check)
    ├── count.png / .pdf              bar = mean of n, error = SEM, dot colour = n (legend), stars (5.4.2)
    ├── count_pct_of_<control>.png    normalised version (dashed line = 100%)
    ├── count_model_absolute.png, count_model_means.csv      mixed model: estimated mean per treatment
    ├── count_forest_<family>.png, count_contrasts.csv       mixed model: every comparison with its CI
    ├── count_variance_components.*, count_diagnostics.png   where the variation comes from; model checks
    ├── count_loo*.csv, count_per_n_ratio*                   robustness: leave one n out; the effect in every n
    ├── count_twoway.png / .pdf       two-way ANOVA treatment × n + Tukey all pairs (5.4.3)
    ├── layers_<channel>_*            layer comparison (marker channels only, 5.4.5)
    └── n_colours.csv
```
The same set of files exists for every other metric: each intensity (`intensity_<channel>_*`), each readout
(`readout_*`) and % positive. With several cell lines, the all-n analysis is done once per line in
`all_n/<cell line>/`. Details of every file: [docs/statistics.md](docs/statistics.md).

**Prism:** open a `*_prism*.csv` and paste it into a Column table, or use File → Import.
The all-n files have the n label as a row title in the first column.

**`results.xlsx`** has an *About* sheet (what this generation is, which test drew the stars, acceptance, duplicates,
notes), *Counts used*, and sheets for the counts, intensities and other metrics at the all-n and per-n levels, the
model (means, ratios, comparison families, variance components), robustness, two-way ANOVA and acceptance.

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

The mixed-model results are in `*_contrasts.csv` and `stats_long.csv` (same idea: estimate, CI, p, adjusted p,
stars, plus `estimand`, which says exactly what the number is).

#### 5.4.2 How significance is shown on the plots
- **Stars** (\* to \*\*\*\*) above a bar appear only when that treatment's comparison is significant
  (adjusted p < 0.05). **No mark = not significant.** The control used is labelled **"ref"** above its bar.
- **Per-n plots** use Dunnett vs the control: `count.png` vs `reference`, `count_vs_<X>.png`
  vs each other control. If the reference is missing in that n, `count.png` uses the next control
  group, and the footer says which. These stars describe one experiment only.
- **All-n plots** take their stars from `[analysis] primary_test`:
  - `mixed` (default): a hierarchical model on the log scale with fields within wells within n's, and a separate
    n-to-n variance for each treatment class. Comparisons are ratios (0.70 = 30% fewer nuclei) with 95% CIs,
    adjusted for all comparisons against the same control. If the model cannot be fitted for a metric, that
    metric falls back to the RM ANOVA and the footer says so.
  - `rm_anova`: the Prism standard, a repeated-measures one-way ANOVA (n = matched replicate) then Dunnett vs the
    control, using the pooled error from all treatments.
- The footer of every plot names the test, so a figure can be read without the tables.

**How to read the statistics:**
- **Within one n** (`<n>/count_stats.csv`): the fields are *technical* replicates, so these p-values only
  describe that experiment. Do not report them as the result.
- **All n** (`all_n/`): each n is one biological replicate. Several tests are always given side by side:
  - *Mixed model* (`*_contrasts.csv`, the default for the stars): uses every field and well, but treats the n as the
    unit of replication. Its degrees of freedom come from the number of n's, so the p-values stay honest with 3–6
    n's. It does not assume that every treatment varies equally between n's.
  - *RM one-way ANOVA + Dunnett* (`*_stats.csv`): the Prism standard. It pools one error for all treatments, which
    gives many degrees of freedom but is overconfident when some treatments' effects vary much more between n's
    than the controls do. An n that misses one group still counts for its other groups.
  - *Ratio paired t-test* (check, no stars): for each treatment, the ratio to the control within each n, tested
    across n's, Holm-corrected. The most cautious; with 3 n's it rarely reaches p < 0.05.
  - *Vehicle check* (equivalence test): "not significant" does not mean "no vehicle effect". Only **equivalent**
    supports that; **inconclusive** means more n's are needed.
  - *Robustness*: leave-one-n-out (`stable = yes` when no single n changes the conclusion) and the effect in every n.
  - A test only uses n's that contain *both* groups, so `N` can differ between comparisons when a
    group is missing from some n. With one n, estimates are given but no across-n test ("needs ≥ 2 n").
- **Percentages** (% positive, `readout_frac_*`): read the difference in percentage points and the absolute
  percentages; ratios of percentages are secondary.
- **Normalised values** (`pct_of_...`) are generally the better way to compare across n's, because they
  cancel day-to-day differences in seeding, staining and lamp brightness.
- **Intensity values** are only comparable if exposure and gain were fixed (auto-exposure off) within each n.

#### 5.4.3 Two-way ANOVA with all-pairs comparisons (`*_twoway*`)
An extra analysis at both levels, like Prism's two-way ANOVA with "compare all pairs":
- **Per n:** treatment × well (wells from `[design]`); fields are the replicates.
- **All n:** treatment × n; fields are the replicates.

The model is `value ~ treatment + second factor + interaction` with Type II sums of squares, so a group
missing from one n is handled correctly. `*_twoway_anova.csv` gives F, df, p and the
% of total variation for treatment, the second factor and the interaction. `*_twoway_pairs.csv` compares
every pair of treatments (Tukey–Kramer, `p_adj`) on cell means matched within each well / n, using only the
wells / n's that contain both groups (`matched_levels`), with the mean difference and its 95% CI. The plot
shows grouped bars (mean ± SEM of fields per well / n) next to a matrix of the adjusted p-values.

**Caution, all-n level:** treating fields as replicates makes these p-values more optimistic than the tests that
use the n as the unit. Report the primary test as the result; use the two-way analysis to see whether treatments
differ from each other and whether the effect depends on the n (the interaction).

#### 5.4.4 Acceptance and data checks
`acceptance.csv` lists, for every n, each criterion of `[acceptance]` with its value, limit and pass / fail.
`design_check.csv` (project folder) has one row per n × treatment with the number of images, wells and fields,
missing image numbers and plain-language flags; a row with 0 images is a condition the n does not have (for example
a missing reference). `duplicates.csv` lists raw files with identical content: results for the affected conditions
are not interpretable until the real images are found. All of these also appear in the **Attention** panel.
Unusual designs (one well, a missing reference, an empty folder, a single n, ...) are described in
[docs/design.md](docs/design.md).

#### 5.4.5 Marker channels: `intensity/` and the layer comparison
`intensity/` (rewritten on every run) holds the cut-off used per n and channel (`cutoffs.csv`), % positive per
condition with and without flat-field correction, a cut-off sweep, per-nucleus distributions, the flat-field maps
and the channel checks. In each generation, every marker channel gets a **layer comparison**
(`layers_<channel>_*`): the same comparison at five levels of processing (1 whole field raw, 2 whole field
corrected, 3 per nucleus raw, 4 per nucleus corrected, 5 % positive), and a table that names the step at which the
conclusion changes. Details: [docs/intensity.md](docs/intensity.md).

### 5.5 How many images are needed, and planning: `sampling/g1/`
| File | Contents |
|---|---|
| `REPORT.txt` | Plain-language summary, starting with **recommendations for the next experiment**: start here |
| `images_needed.png/.pdf` | Predicted noise, subsampling agreement and per-n error vs images per condition (vs the main reference; `images_needed_vs_<group>` for the others) |
| `batch_by_n.csv` | Per n: reference-group levels, missing references, systematic well difference, and whether this n's treatment effects are stronger/weaker than in the other n's |
| `batch_by_condition.csv` | Per n × treatment: mean per well, well-to-well ratio with a flag when it exceeds field noise, ratio vs each reference |
| `variance_components.csv` | Variation field-to-field (within a well) and well-to-well (within a plate), as CV² |
| `effects.csv` | Per reference and treatment: observed fold change, p with all images, biological (n-to-n) variance |
| `design_power.csv`, `summary_by_design.csv` | For every design (wells × fields), per reference and treatment: noise vs now, power at your current number of n's, and n's needed for 80% power |
| `subsampling.csv` | Real images re-drawn many times per design: how often the result matches the conclusion from all images, and how far each n's ratio moves |
| `planner.json`, `planner_components.csv` | The forward planner (field, position, well and n-to-n variation with intervals), used by the **Planner** page |
| `design_forward.csv`, `mde.csv`, `cost_frontier.csv`, `allocation.csv`, `more_n.csv` | Power and CI width for every design, the smallest detectable effect, the cheapest designs that reach the target power, how to split wells between treatments and references, and the chance that more n's settle a comparison (+ plots) |

Every analysis is done once per reference group (`reference` + `normalise_to`), like the statistics in
section 5.4: an n without one of the controls is left out of the comparisons vs that control but still counts vs
the others. Predictions come as pessimistic / expected / optimistic, because with 3–4 n's every estimate of
variation is uncertain. Cell lines are planned separately (`sampling/g1/<cell line>/`). Details:
[docs/planner.md](docs/planner.md).

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
   copy the new rows into `layout.csv` (check the names), or use the fresh scan on the **Layout** page.
3. `kerncount run <experiment>`. Only the new images are processed; the rest come from the cache.
   Generations with manual corrections are not overwritten. Rebuild them on purpose with `--rebuild` once the new n
   is checked.

**A new experiment**
1. Create `data/<new experiment>/` with one folder per n.
2. `kerncount init <new experiment>`, check `layout.csv` and the `[experiment]` / `[analysis]` settings, then
   `kerncount run <new experiment>`.

**Then, for either:** check the **Attention** panel and the red and orange tiles in `review/*_overview.jpg`,
correct counts in `counts.csv` → `manual_1`, and run again. Report from the latest generation:
`generations/g<N>/results.xlsx`, the `all_n/` plots, and the Prism CSVs.

**Planning the next experiment:** `kerncount sampling <experiment>`, then the **Planner** page (or
`sampling/g<N>/REPORT.txt`), then `kerncount platemap <experiment>` for the plate layouts.

## 7. Troubleshooting
| Problem | Fix |
|---|---|
| "project not found" / "raw folder not found" | Run the command from the workspace folder (the one containing `data/`), check the experiment name, and run `kerncount init <experiment>` first. |
| First run fails while downloading the model | The StarDist model downloads once on first use; connect to the internet and run again. |
| TensorFlow will not load (Windows) | Install the Microsoft Visual C++ Redistributable, or use `[segmentation] backend = "classical"`. |
| `counts.csv` opens as one column in Excel | Set `csv_delimiter = ";"` in `settings.toml` (European Excel). Files Excel re-saves with `;` are read correctly either way. |
| `run` fails to write `counts.csv`, or "cannot replace gN" (Windows) | Close the file in Excel first; Excel locks open files. The new results wait in `gN.tmp` until you run again. |
| "cannot read manual value" | A manual cell contains something other than a whole number or `x`; the message gives the row and column. |
| "g2 already exists but counts.csv now gives different values" | You changed a manual column after g2 was built. Re-run with `--rebuild g2` to update it on purpose. |
| A treatment name is wrong or a folder should be left out | Edit `layout.csv` (`treatment`, `include`) and run again. |
| "has 1 channel(s) but settings ask for channel 2" | `[channels] segment` points to a channel the images do not have. (Channels in `[intensity] measure` that do not exist are skipped silently; check the `intensity_*` columns.) |
| Marker readouts empty for a whole folder | The channel check failed (files not paired as `file_order` says, or a file missing so every pair after it is shifted). See `intensity/alignment.csv` and the Attention panel. |
| Stars are missing, "needs ≥ 2 n" | An across-n test needs at least 2 n's that contain both groups. |
| Run is slow | First runs take ~2–3 s per image; later runs are cached. `--no-review` skips the review images. |
| Counts look wrong in a whole folder | Check `diameter_px` in `qc_details.csv`; if the size estimate is off (e.g. almost empty folder), set `nucleus_diameter_px` to a number. |
| Tumour run "blocked" | See the list at the top of the **Layout** page (or the run output): fix the cell, or untick *include* and give the reason in *note* (section 9.2). |
| Tumour run refuses a changed analysis plan | The plan is locked once results exist; unlock with a reason (section 9.3). |

## 8. Notes and known limitations
- 8-bit camera images saved as 16-bit TIFF (every value × 257) are handled. Saturation is detected either way.
- Compared with a simple Fiji threshold macro, counts are usually lower. Threshold macros count
  background noise in blank or defocused fields, split lobed nuclei, and count every edge-touching nucleus.
- Very dense fields where cells have bright cytoplasm are the hardest case. They are flagged QC low/medium;
  check them in the review images.
- Without a known pixel size, results are per field rather than per mm².
- QC levels are heuristics until calibrated against hand counts.
- With 3–4 n's every between-n estimate (variances, power) is rough; read the confidence intervals, not only the
  stars.
- The built-in treatment-name rules and the default control names (`Control`, `LyoVec`, `CL075`) come from the lab
  the tool was written for; other names work, but set `reference`, `normalise_to` and `treatment_order` to yours.
- Algorithms and implementation notes: [docs/technical.md](docs/technical.md).

## 9. Tumour experiments (fluorescent tumours in brain sections)

The same tool also analyses **tumour experiments**: fluorescently labelled (e.g. mCherry) tumours in brain
sections, imaged as confocal z-stacks (`[experiment] type = "tumour_volume"`). One **tumour** (hemisphere A or B of a
section) is one n; A and B are never averaged. Technical details: [docs/tumour.md](docs/tumour.md).

### 9.1 Input
One folder per experiment in `data/`, with **one sub-folder per imaging session folder as delivered** (e.g.
`5d KO`, `5d WT`), files as delivered:
```
data/my_tumours/
├── 5d KO/                                      ← session folder: its name gives the genotype
│   ├── 010126 lineA +drug n3 A.ome.tiff        ← one original z-stack = one tumour
│   ├── 010126 lineA +drug n3 Af.tif            ← Fiji mask (legacy, optional; read for comparison only)
│   └── Statistics for 010126 lineA +drug n3 Af.csv   ← Fiji 3D Objects Counter table (legacy, optional)
└── 5d WT/
    └── 150126 control -drug n1 B.ome.tif
```
- **Factors.** What the groups are is set in `settings.toml` as `[[tumour.factors]]`: each factor has its levels, a
  reference level, other spellings (`aliases`), where it goes in the plots, and where it is read from: the
  **session-folder name** or a word in the **file name**. The settings template ships with an example design of
  three factors: a genotype (from the session-folder name) and two factors read from the file name (a tumour line
  and a ± treatment). Edit the levels to match your experiment; a new level or a new factor is a settings entry, no
  code change. The example above assumes levels named `lineA` / `control` and `+drug` / `-drug`.
- **Session folders** are read word by word (upper/lower case does not matter); words that are not a level (`5d`) are
  kept as the session label. Any number of session folders is fine, also several of one genotype (repeat sessions).
  A folder whose genotype cannot be read stops the run with a message saying how to fix it.
- **Sub-folders do not matter**: below a session folder the files may sit in any folder structure; everything except
  the session-folder factors is read from the file name. `_s1` endings, `1n` for `n1` and a missing date are
  handled; a `.czi` next to the `.ome.tiff` of the same tumour is recorded, a tumour that exists only as `.czi` is
  reported (export it as OME-TIFF). A session folder may be a link to an external drive; if the drive is not
  connected, the run says so and existing results stay viewable.
- **Renaming a session folder** loses nothing: `kerncount run` finds the files by name, keeps your edits in
  layout.csv and re-uses every cached measurement; `kerncount init` then updates layout.csv (old copy backed up).
- **File names** give the date, the file-name factors, section and hemisphere. For **new acquisitions** use
  `DDMMYY <factor levels> n<section> <A|B>.ome.tif`: date as 6 digits, each factor level exactly as in the
  settings, the section label `n<k>` (restarts per group; `s3` is not read as a section), hemisphere `A` (upper) or
  `B` (lower). Underscores instead of spaces, `.ome.tiff`, and aliases (e.g. known typos) are read too and noted in
  `layout.csv`.
- **Batches.** What counts as a batch is the setting `[tumour] batch`: by default the **culture batch** (the
  6-digit date at the start of the file names = the culture preparation), or the *session folder*, or the
  *acquisition day* (from each stack's metadata). Several batches that share a group enter the model as a fixed
  block effect. A comparison that cannot be separated from the batch (e.g. two genotypes that come from different
  batches) is shown, labelled "confounded", but never tested; ratios within a genotype can still be compared across
  genotypes. The report says which definition was used.
- Voxel size and acquisition settings are read from each file's OME metadata (never typed in); z-steps may differ
  between files.

### 9.2 Set up, check the layout, run
In the web interface pick the experiment, choose **experiment type: Tumour volume (z-stacks)** under *Set up* and
press **Set up** (command line: `kerncount init my_tumours --type tumour_volume`). Then:
- **Layout**: one row per tumour with its factor values, section, hemisphere, tumour ID (factor levels + section and
  hemisphere, e.g. `KO / lineA / +drug / n3A`), include and note. Problems that would make a file a wrong datapoint
  (an unreadable name, a duplicate tumour ID, a session without a genotype, a missing file, a stack in the raw
  folder that layout.csv does not list yet, an exclusion without a reason) are listed at the top and **block the
  run** until you fix the cell or untick *include* and give the reason in *note*. Stacks added to the raw folder
  later: run **Set up** again (`kerncount init my_tumours`); it adds them and keeps your edits.
- **Settings** → *Analysis plan*: the factors (levels, reference level, where they go in the plots), which readouts
  are confirmatory, and the contrast families (below). *Tumour measurement* (`[tumour]`): smoothing, threshold rule,
  minimum object size, mid-section thickness, outline smoothing, threshold sweep, batch definition, QC, review and
  planning-aid settings, each with an explanation.
- **Run** (`kerncount run my_tumours`): options *skip review images* (`--no-review`), *segment every stack again*
  (`--resegment`), *threshold sweep*. The first run takes about 20–60 s per stack, plus about as much again for the
  **threshold sweep** (on by default: every confirmatory effect re-computed at the threshold × 0.8–1.2; untick it,
  or set `threshold_sweep = false`, to save that time). Segmentations are cached, so later runs are quick.

### 9.3 Readouts, and why area rather than volume
Often only a slab (e.g. ~80–150 µm) of a much larger tumour is imaged, so a 3D volume mostly measures how deep the
stack was acquired. The **confirmatory readouts** are therefore 2D, on a fixed-thickness **mid-section**
(`mid_section_um`, default 20 µm, around the layer where the tumour's main mass is largest; exactly the same
physical thickness for every tumour, whatever the z-step; the thickness used is in `tumours.csv`): the
**cross-sectional area** (mm²) and the **perimeter** (mm) of the main mass (the largest connected object).
Everything else is **exploratory** (effect sizes with CIs, labelled, no stars): footprint area, slab volume,
main-mass volume and surface, tumour burden (main mass + satellites), elongation, solidity, edge roughness,
circularity, and the Fiji legacy values. One automatic threshold rule for every tumour, never adjusted by hand per
tumour.

**Analysis plan lock.** Factors, readout roles, contrast families, the batch rule and the measurement settings
behind the confirmatory readouts (smoothing, threshold rule, minimum object size, mid-section thickness, outline
smoothing) are fixed before results are viewed. Once results exist the plan is **locked**; changing it needs
**Unlock** with a reason — in the web interface, or on the command line
`kerncount run my_tumours --unlock "reason"` (`kerncount run` refuses a changed plan without it). The unlock and
every change are written to `lock_log.csv`, and the next run prints the history (old → new plan lock) in the
report and the About sheet of `results.xlsx`. Every statistics row carries the lock it was built with. QC, sweep,
sensitivity, review and planning-aid settings are not part of the plan and can change at any time.

### 9.4 Output (`analysis/<experiment>/`)
| file / folder | content |
|---|---|
| `results.xlsx` | everything in one workbook: tumours, statistics, estimated means, families, QC, flow, sensitivity analyses, threshold sweep, concordance with Fiji, planning aid (download it on the **Results** page) |
| `stats_report.md` | the report: analysis plan and its lock (and lock history), results per readout and family, QC, acquisition settings, sensitivity, Fiji concordance, planning aid |
| `tumours.csv` | one row per analysed tumour: factor values, batch, every readout, QC flags |
| `excluded_tumours.csv`, `flow.csv` | excluded tumours with their reason (layout or review); files → tumours → analysed, per group |
| `stats_long.csv` | every estimate in one long table (the **Statistics** page; sensitivity rows hidden unless asked for) |
| `acquisition.csv` | per tumour, the microscope settings from the file metadata (detector gain, laser, pinhole, voxel size, acquisition day, ...) |
| `concordance.csv`, `threshold_sweep.csv` | new method vs Fiji per contrast; every confirmatory effect at each threshold scale |
| `lock_log.csv`, `corrections_used.csv` | plan unlocks and changes; the review corrections applied |
| `plots/` | per readout: absolute, normalised (% of the reference level within the other factors), effects (forest), interaction, A vs B, residuals, and with Fiji data / the sweep: vs Fiji, threshold sweep; `shape_panel` |
| `prism/` | Prism-ready tables, one column per group (absolute and normalised) |
| `review/` | the **blinded** review images and `corrections.csv` (below) |
| `review_key.csv` | code → tumour (outside `review/`, so the review folder stays blind) |

**Statistics:** one linear model on log(readout) over all factors; effects are fold-changes with 95 % CIs, Holm
within each family. Only confirmatory families of confirmatory readouts carry stars. **Sensitivity analyses**
(Welch-type, permutation, without QC-flagged tumours) are reported alongside, never instead. The **planning aid**
gives the smallest detectable fold-change at the current group sizes and the tumours per group needed for
`planner_fold_change` (default 1.5×); it never feeds the results. The **acquisition check** says whether every
brightness-relevant microscope setting was identical, and which groups were imaged on a single day; with
`shared_threshold_check = true` (and identical settings) the analysis is repeated at one shared threshold. The report
also says that the Fiji legacy volume mostly reflects the acquired stack depth, and prints the section thickness and
how the z-range was chosen at the microscope (settings `section_thickness_um`, `z_range_rule`; "not recorded" while
empty).

**QC:** a tumour is **flagged** when its main mass touches the image frame or its mid-section does not fit in the
stack; a **saturation warning** is given above `qc_saturation_warn`. Flagged tumours stay in the main analysis (drawn
hollow) and are dropped in a sensitivity analysis.

### 9.5 Review and corrections (blinded)
`review/` holds one image per tumour, named by a random **code** (no group, session, date or file name anywhere):
left the maximum projection, right the mid-section layer; the main mass orange, satellites cyan, objects numbered.
The **Review images** and **Corrections** pages show only the codes. In a round (`manual_1` → g2, `manual_2` → g3,
…) type:
- an **object number** from the image to make that object the main mass (`1` = back to automatic);
- **`x: reason`** to exclude the tumour (the reason goes to `excluded_tumours.csv`);
- nothing = no change.

Then run the analysis again: nothing is re-segmented; the picked object's readouts are recomputed (the threshold
sweep follows the picked object too). You can also edit `review/corrections.csv` in Excel and hand only the
`review/` folder to a colleague for a blinded check.

## 10. Development
```
pip install -e ".[test]"   # Python 3.11 or newer
pytest                     # fast test suite
pytest -m slow             # heavy simulations and end-to-end runs (skipped by default)
```
`pytest -m "slow or not slow"` runs everything. The browser interface is plain HTML/JavaScript in
`src/kerncount/web/` (no build step); the one-click launchers are built from `launcher/` (Go; see
[launcher/README.md](launcher/README.md)).

Detailed documentation:

| Document | Contents |
|---|---|
| [docs/technical.md](docs/technical.md) | How counting works, algorithms, experiment types |
| [docs/statistics.md](docs/statistics.md) | Every statistical output, the mixed model, `stats_long.csv`, `manifest.json`, when to trust what |
| [docs/intensity.md](docs/intensity.md) | Marker channels: separate channel files, per-nucleus measurement, flat field, positivity, layer comparison |
| [docs/readouts.md](docs/readouts.md) | The extra per-image readouts |
| [docs/design.md](docs/design.md) | Wells, field positions, `wells.csv`, timepoints, unusual designs |
| [docs/planner.md](docs/planner.md) | The planner and plate maps |
| [docs/ui.md](docs/ui.md) | The browser interface |
| [docs/tumour.md](docs/tumour.md) | Tumour experiments in full |

Some of these documents refer to a glossary (`CONTEXT.md`) and design decisions (ADRs) that are kept outside this
repository; the documents are complete without them.

## Acknowledgements
Nucleus detection uses [StarDist](https://github.com/stardist/stardist) (Schmidt et al., MICCAI 2018;
Weigert et al., WACV 2020) with its pretrained `2D_versatile_fluo` model.

## License
MIT: see [LICENSE](LICENSE).
