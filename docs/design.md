# Plate design: wells, field positions, timepoint

kerncount needs to know which **well** each image (field) comes from and where in the well it was
taken (its **position**), because both are sources of variation:

- well-to-well differences only shrink with more wells, field-to-field ones with more fields;
- on the real data, field position matters a lot (fields 2/5/6 count up to ~2x more than 1/3/4;
  position explains ~35% of the field-to-field variance).

Nesting (see `CONTEXT.md`): field ⊂ well ⊂ condition ⊂ n (one n = one plate).

All of this lives in `src/kerncount/design.py`. No other module works out wells by itself.

## Settings (`[design]` in settings.toml)

| key | default | values |
|---|---|---|
| `wells` | `2` | a number of wells per condition per n, or `"auto"`: use `wells.csv` if the project has one, else 2 |
| `fields_per_well` | `"auto"` | `"auto"`: worked out per condition (below), or a fixed number |
| `well_layout` | `"consecutive"` | `"consecutive"`, `"alternating"` or `"file"` |
| `field_positions` | `"fixed"` | `"fixed"` or `"random"` |

The old names `plates`, `fields_per_plate` and `plate_layout` are still read; a note asks you to rename them.

### How `auto` works

For every condition (n × treatment) separately:

    fields_per_well = ceil(highest image number / wells)

| images present | wells | fields per well | result |
|---|---|---|---|
| 1-12 | 2 | 6 | 6 + 6 |
| 1-12 without 8 | 2 | 6 | 6 + 5 (unbalanced, image 8 missing) |
| 1-11 | 2 | 6 | 6 + 5 (unbalanced, image 12 missing) |
| 1-18 | 2 | 9 | 9 + 9 |
| 1-18 | 3 | 6 | 6 + 6 + 6 |
| 1-6 | 2 | 3 | 3 + 3 |

The highest image number is used, not the number of images, so a gap does not shift the images after it
into the wrong well. With a fixed number, images above `wells × fields_per_well` stay in the last well
(consecutive) or get positions above `fields_per_well` (alternating), and are flagged.

### Layouts

- **consecutive**: images 1..k are well 1, k+1..2k well 2, ...; position = place within the well
  (image 9 with 6 fields per well = well 2, position 3).
- **alternating**: round robin; with 2 wells, odd images are well 1 and even images well 2;
  position = (image - 1) // wells + 1.
- **file**: read `wells.csv` in the project folder (see below). `wells` and `fields_per_well` are then
  ignored: the file defines them.

### Field positions

- **fixed** (default): field k of every well is the same physical spot (e.g. the microscope's stage
  positions are saved and reused), so `position` can explain differences between fields.
- **random**: fields were picked by hand anywhere in the well. `position` is still filled in (it is the
  order in which fields were taken) but carries no meaning, and models must not use it.

## wells.csv

For designs the rules cannot describe. Columns:

| column | required | meaning |
|---|---|---|
| `image` | yes | image number |
| `well` | yes | well number (1-based) |
| `position` | no | field position in the well (1-based); blank = number the well's fields in image order |
| `n` | no | blank = every n; `n1` matches the full label `n1_2026-06-05` |
| `treatment` | no | blank = every treatment |

For every image the most specific row wins: n + treatment, then n only, then treatment only, then the
generic row (neither). Images not covered by any row fall back to the consecutive rule and are flagged.

To start from the current assignment, write a template:

    kerncount wells-template <experiment>        # writes wells.csv (or wells.new.csv if one exists)

It holds one generic row per image number plus n + treatment rows only where a condition differs.
Edit it, then set `well_layout = "file"` (or `wells = "auto"`).

## Timepoint

`layout.csv` has a `timepoint` column. `kerncount init` fills it from the n-folder name and normalises
it: `72h`, `48 h`, `24hr`, `72 hrs` → `72h`, `48h`, `24h`, `72h`; `day 3`, `d3`, `3 days` → `d3`. It is
blank when the name has no timepoint. Edit it freely. Older `layout.csv` files without the column still
work: the timepoint is then read from the `n_folder` names on the fly.

## Outputs

- `design_check.csv` in the project folder (written by `run` and `regen`): one row per n × treatment
  (columns below). Check it whenever the summary line reports something unusual.
- `generations/gN/counts_used.csv` and the "Counts used" sheet of `results.xlsx` have the columns, in this
  order: `n`, `treatment`, `image`, `count`, `source`, `excluded`, `qc_level`, `well`, `position`, `timepoint`,
  `cell_line`, `meta_magnification`, `meta_camera_pixel_um`, `meta_seeding_density`, `meta_seeding_unit`,
  `meta_stain`, `meta_negative_control`, `meta_negative_control_type` (from `design.output_columns`; metadata are
  strings, `""` = not set), then the field's source files: `source_file` (all files of the field, `; `-separated)
  and, for interleaved channel files, `source_file_<channel>`. Per-treatment `fields.csv` also lists the source files.
- A change of the design or the metadata (e.g. a new `cell_line`) rebuilds existing generations whose counts are
  unchanged (see docs/statistics.md, "When a generation is rebuilt").
- `run` / `regen` log one line, e.g.

      Design: 4 n, 8 conditions, 2 wells x 6 fields (auto), timepoint 72h (some n unknown);
      1 condition unbalanced (n1 / CL075: 11 images); see design_check.csv

### design_check.csv columns

| column | meaning |
|---|---|
| `n`, `treatment`, `timepoint` | the condition |
| `images` | number of distinct image numbers |
| `wells_expected`, `wells_present` | wells in the design / wells that have at least one image |
| `fields_per_well`, `fields_source` | design size of a well; `auto`, `fixed` or `file` |
| `fields_min`, `fields_max` | fewest / most images in any present well |
| `max_image` | highest image number |
| `missing_images` | numbers absent from 1..max(max_image, wells × fields_per_well), `; `-separated |
| `beyond_design` | image numbers above wells × fields_per_well (fixed size only) |
| `unbalanced` | wells have different numbers of images, or an expected well is empty |
| `single_well` | all images are in one well |
| `flags` | plain-language list (`; `-separated); empty = nothing unusual |

Rows with `images` = 0 are conditions an n does **not** have (`design.absent`): a condition present in another n of
the same cell line, or an included `layout.csv` folder without images. Their well/field columns are empty.

## Unusual designs

Every case below runs end to end (detection, design check, statistics, readout ladder, duplicates, planner, plate
maps, browser) without a traceback. What drops out is said in `design_check.csv`, the run log, the generation notes
(`manifest.json` `notes`, README, About sheet) and the **Attention** panel; a headline never turns into a silent
empty value. Tests: `tests/test_unusual_designs.py` (rules: default suite; every case end to end: `-m slow`;
synthetic projects in `tests/unusual.py`).

| case | behaviour |
|---|---|
| **a. Interleaved channel files**: gaps in the numbering | ignored: files pair up in file order (`discover.group_files`); the field's files are listed in `source_file_<channel>`. |
| odd file count (unpaired last file) | the leftover file is skipped and named in the log and in `layout.csv` `note`; Attention: "folder(s) with files that were not used". |
| a folder holding one channel only (or pairs shifted by a missing file) | the channel check fails in most fields (mis-paired, shifted by more than `align_max_correct_px`, or the nuclear channel not the brightest): at ≥ 50% of a folder's fields (and ≥ 2) the whole folder is *suspect*: its marker readouts are left empty (counts stay, from every second file), the log warns, `channel_check` says "folder fails the channel check", and Attention names the folder. |
| **b. An n without a reference** (no LyoVec, or no Control) | `design_check.csv` row with `images` = 0, "reference LyoVec absent from this n (present in n1, n3)"; log warning; Attention "condition(s) missing from an n". Comparisons vs that reference use the n's that have it (`N` in the contrasts); the vehicle check and `vehicle_within` are "not applicable" in that n. |
| an n without a treatment the others have | the same flag ("condition absent from this n"); its comparisons use the n's that have it; with fewer than 2 such n's: estimate only, "needs ≥ 2 n". |
| **c. One well instead of 2** (images 1-6 only) / last images missing (1, 7, 8, 9) | `auto` reads the design from the highest image number, so it sees a smaller well; when most of the n has a larger size, the row is flagged "fields per well read as 3 ... but 6 in the other conditions of this n" (fix it in `wells.csv`). |
| wells with 1, 3 or 9 fields; 18 images | `auto` gives 9 + 9; unequal wells are flagged "unbalanced"; the well-mean model handles unequal fields. |
| a well with 0 usable fields after QC (`exclude_low_qc`, manual `x`) | generation note "n3 / LyoVec: well 2 has no usable field (all 6 excluded: QC low); this condition rests on 1 well in n3"; a condition with no usable field at all: "left out of every statistic". |
| **d. An empty condition** (folder without images) | `design_check.csv` row with `images` = 0 ("no images: folder ... is empty"); Attention lists it twice (unused folder, missing condition). |
| all-zero condition (0 nuclei everywhere) | counts: ratio rests on the +0.5 offset (note "read them as 'close to 0'"); means and their CIs are never below 0. Readouts without nuclei (% positive, intensities): note "has no value in any field of any n ...; not testable for this readout". |
| every fraction 0 or 1 | analysed (empirical logit keeps it finite); note "0% / 100% in every field ... read the percentage-point difference" (the odds ratio is arbitrary there). |
| **e. One condition + a reference** | ordinary analysis with one comparison per reference. |
| a single n | estimates only, "needs ≥ 2 n" everywhere at all n; the planner page says it needs ≥ 2 n. |
| three cell lines, one with 1 n | never pooled: one all-n block per line (`all_n/<line>/`), the 1-n line gets "needs ≥ 2 n" and no plan; a condition missing in another line is not "missing". |
| **f. A whole condition duplicates another** (the real n3 case) | `duplicates.csv`, Attention error, acceptance `no_duplicate_images` fails. Combined with a missing reference elsewhere, the note says how many n's are left with a reference that is neither missing nor a copy, e.g. "comparisons vs LyoVec: only 1 n (n1) has a LyoVec that is neither missing nor a duplicate - needs ≥ 2 n for a trustworthy test". |

## Python interface (for developers)

```python
from kerncount import design

out = design.assign(table, settings, project_dir=None, layout=None)
```

`table` needs `n`, `treatment`, `image` (int). Returns a copy with `well` (int, 1-based), `position`
(int, 1-based) and `timepoint` (str, `""` if unknown) added (replaced if present), followed by the experiment
metadata: `cell_line` (str, `""` = unknown) and `meta_magnification`, `meta_camera_pixel_um`,
`meta_seeding_density`, `meta_seeding_unit`, `meta_stain`, `meta_negative_control`, `meta_negative_control_type`
(all str) from `[experiment]`, overridden per n / condition by same-named `layout.csv` columns (blank = experiment
value; see `docs/intensity.md`). `design.metadata_columns()` lists them, `design.output_columns(assigned)` gives
every design + metadata column to merge into per-field tables. Rows are never dropped
or reordered; the index is kept. `project_dir` is where `wells.csv` and `layout.csv` are read; `layout`
(a layout DataFrame) overrides `layout.csv` for timepoints. Without either, an existing `timepoint`
column is kept. Raises `ValueError` for invalid settings or a missing/broken `wells.csv`.
`out.attrs["design"]` holds the resolved per-condition design used by `describe`.

| function | returns |
|---|---|
| `design.describe(assigned)` | the design_check table (one row per n × treatment, columns above) |
| `design.absent(assigned, settings, layout=None)` | rows (`images` = 0) for conditions an n lacks: present in another n of the same cell line, or an empty included `layout.csv` folder; appended to design_check.csv by `pipeline.check_design` |
| `design.summary_line(desc, settings)` | the one-line log summary (lists `images` = 0 rows as "missing from an n") |
| `design.spec(settings)` | `Spec(wells, fields_per_well, well_layout, field_positions, notes)`; `"auto"` stays a string; old names resolved |
| `design.positions_meaningful(settings)` | `True` if `field_positions = "fixed"`. **Check this before using `position` in any model.** |
| `design.nominal(desc)` | `(wells, fields_per_well)` of the largest condition: the "full design" for the sampling analysis |
| `design.wells_template(assigned)` / `design.write_wells_template(project_dir, table, settings, overwrite=False)` | wells.csv contents / write it |
| `discover.parse_timepoint(name)` | normalised timepoint from a folder name (`""` if none) |

Where it is wired in:

- `pipeline.build_generations` → `pipeline.check_design` assigns once per run (from the images in
  counts.csv), writes `design_check.csv`, logs the summary, and passes the table to `generations.build`
  (`design_table=`), which merges `well / position / timepoint` into every row it analyses.
- `analysis.per_n` uses the `well` column for the treatment × well two-way ANOVA (it falls back to
  `design.assign` when the column is missing) and skips that ANOVA when an n has only one well.
- `sampling.run` takes wells from `design.assign` (or a `well` column) and simulates designs up to
  `design.nominal`. `sampling.design()` and `sampling.well_of()` remain as backward-compatible wrappers.
