# Fluorescence intensity and % positive nuclei (e.g. Ki-67)

Design decision: `docs/adr/0001-ki67-intensity-analysis.md` (in the workspace). Words: `CONTEXT.md`.

kerncount finds nuclei on the nuclear channel (DAPI) as for counting, then measures every **marker channel**
(e.g. Ki-67) inside each nucleus. The primary readout is **% marker-positive nuclei** per image; the per-nucleus
intensities, the mean intensity of the positive nuclei and the old Fiji "nucleus + 6 px" mean are secondary.
Nothing in this page changes detection or the counts.

Code: `intensity.py` (measurement, project stage), `flatfield.py`, `positivity.py`, `intensity_report.py`
(outputs), `tuning.py` (browser API). Tests: `tests/test_intensity.py`.

## 1. Channels saved as separate files

```toml
[channels]
files = "interleaved"          # default "single": all channels inside one file
file_order = ["Ki67", "DAPI"]  # channel names in file order
segment = "DAPI"               # a name or a 1-based number
```

- In each treatment folder the numbered files are sorted by number and grouped in `len(file_order)`: (1, 2),
  (3, 4), ... Gaps in the numbering are ignored (a folder numbered 2-25 or one that skips 3 still pairs
  correctly). The **field number** is the pair index 1..k; `[design]` assigns wells from it as usual
  (12 fields = 2 wells x 6).
- A leftover file (odd count) is reported (`layout.csv` note, `run` log) and skipped.
- Greyscale saved as RGB (three identical channels) is read as one channel. A file with several *different*
  channels is an error in this mode.
- `kerncount init` notices single-channel files that alternate dim / bright and suggests these settings; it does
  not change them. After changing `files`, rewrite the layout with `kerncount init <exp> --force-layout`.
- A folder name without a full date (e.g. `Ki67 0611 ...`) gets day/month from the name and the year from the
  file dates, with a note to check it. Edit `n` and `date` in `layout.csv` freely.

**Checks per field** (`intensity/alignment.csv`, column `channel_check` in `qc_details.csv`):
- *alignment*: the shift between the marker and the nuclear image by phase correlation (band-passed, 0.1 px).
  Shifts above `channels.align_warn_px` (2 px) are flagged. `align = "correct"` shifts the marker back by whole
  pixels when the shift is between 1 px and `align_max_correct_px`; the shift is always recorded.
- *pairing*: the correlation of the two (band-passed) images after the shift; below 0.2 the files do not show the
  same cells (mis-paired, duplicated or wrong file). The field keeps its count, but its marker readouts are left
  empty (`excluded` column in the per-nucleus table). In the first Ki-67 set this caught `cntrl/16.tif` (looks
  like the same field as `cntrl/18.tif`, not as its DAPI `17.tif`) and `lyovec/19.tif` vs `20.tif`.
- *channel order*: the nuclear channel should be the brightest on its nuclei; otherwise the field is flagged
  (`file_order` reversed?).

## 2. Experiment metadata

`[experiment]` (right after `[project]`): `cell_line`, `magnification`, `camera_pixel_um`, `seeding_density`,
`seeding_unit`, `stain`, `notes`, `negative_control` (a treatment name, "" = none), `negative_control_type`
(`secondary_only` | `isotype` | `unstained`). Pixel size = `camera_pixel_um / magnification` when
`[segmentation] pixel_size_um` is empty.

Overrides: a `layout.csv` column with the same name (`cell_line`, `seeding_density`, ...); a blank cell keeps the
experiment value. A layout row is one n x treatment folder (fill every row of an n to override the whole n).

`design.assign` adds `cell_line` (str, "" = unknown) and `meta_<key>` for `magnification`, `camera_pixel_um`,
`seeding_density`, `seeding_unit`, `stain`, `negative_control`, `negative_control_type` (all strings) after
`well`, `position`, `timepoint`. `design.metadata_columns()` lists them; `design.output_columns(assigned)` gives all
design + metadata columns to merge into per-field tables (`counts_used.csv`).

## 3. Measurement per nucleus

Per counted object and marker channel (`_work/nuclei_intensity.csv`, one row per nucleus and channel):

| column | meaning |
|---|---|
| `area` | nucleus area (px) |
| `erode_px`, `measured_px` | outline shrunk by this much before measuring; pixels left (N) |
| `mean`, `p90` | raw mean / 90th percentile inside the shrunk outline |
| `total` | integrated intensity over the whole nucleus minus area x local background |
| `bg`, `bg_px`, `bg_source` | local background (median of a ring), its pixels (M), `ring` / `wide ring` / `image` |
| `fiji_cell_mean` | legacy: raw mean of the nucleus enlarged by `fiji_ring_px` (6 px, the old macro) |
| `saturated` | a nucleus pixel at the camera maximum in this channel: **excluded from intensity readouts, still counted** |
| `noise_sd`, `noise_sd_p90`, `p90_null` | background-noise SD expected for (mean - bg), (p90 - bg); p90 of pure noise |
| `ff` | flat field at the nucleus centre |
| `mean_bgsub_ff`, `p90_bgsub_ff`, `total_bgsub_ff` | (value - bg) / ff: the **signal** used everywhere below |
| `counted`, `measurable` | counted nucleus; counted and not saturated |
| `thr_lod`, `thr90_lod` | this nucleus' automatic (limit-of-detection) cut-off, corrected units |
| `thr`, `thr_p90`, `positive`, `positive_p90` | cut-off in force and the call |

- **Shrinking rule** (`intensity.erode_px = "auto"`): 2 px for nuclei at least 16 px across (equivalent circle
  diameter), else 1 px; lowered again where fewer than max(9 px, 25% of the area) would remain. It leaves out the
  blurred rim, which mixes in background and neighbours' light. Contacts with a touching nucleus count as rim.
- **Local background**: a ring from `bg_gap_px` (3) to `bg_gap_px + bg_ring_px` (auto = 0.5 x nucleus diameter)
  outside the nucleus, split between neighbours (each pixel belongs to its nearest nucleus), keeping more than
  `bg_gap_px` away from every detected object (nuclei, debris, out-of-focus ones). Its median is the background.
  Fewer than `bg_min_px` (40) pixels: a 3x wider ring; still too few: the image background (median of all pixels
  away from objects).
- **Noise** (for the automatic cut-off): the robust (MAD) SD sigma of background pixels after removing their smooth
  trend, and the noise correlation area a_c (from 4 x 4 and 8 x 8 block means; camera and optics make
  neighbouring pixels share noise, so N pixels carry only N / a_c independent values). Then
  SD(mean - bg) = sigma x sqrt(a_c / N + (pi / 2) a_c / M) and SD(p90 - bg) = sigma x sqrt(2.92 a_c / N + (pi / 2) a_c / M).
- Measurements are cached per field (`_work/<n>/<treatment>/NN_intensity.csv/.json`) and redone only when the
  detection, the files or the measurement settings change.

## 4. Flat-field correction

`intensity.flatfield = "auto"` (default) | `"off"` | a folder of reference images.

- **auto**: BaSiC-style, per n and channel from all of that n's images: observed = signal x F + D. Each image is
  reduced to 16 x 16 blocks; D = median over images of the blocks' 20th percentile (background), F = mean over
  images of (block mean - D) / image mean. Both are fitted by a smooth 2-D quadratic (robust for D) and F is
  scaled to mean 1. Cells sit at random positions, so averaging many fields leaves the illumination. Needs at
  least 10 images (else no correction, with a note). No extra dependency (`basicpy` would add a large JAX stack;
  this is enough for smooth vignetting).
- **reference folder**: TIFFs whose name contains the channel name (else all), averaged; files with `dark` in the
  name are dark frames subtracted first.
- Applied to background-subtracted values only: signal = (mean - bg) / F(nucleus centre). Counts, detection and
  the legacy intensity columns are not corrected.
- The automatic cut-off is unaffected by F by design (signal and its noise scale together); control-based and
  manual cut-offs are.
- Reported (`intensity/flatfield.csv`, `flatfield/<n>__<channel>.png/.npy`): edges vs centre
  (`edge_dim_pct`), darkest side, range, a **check from the nuclei** (robust fit of log signal on position: edge
  dimming before and after correction; after should be near 0) and % positive with and without correction.

## 5. Positivity

Per marker channel (`[positivity] channels = "auto"`), per n. **Cut-off precedence**:

1. **negative control** - the `[experiment] negative_control` condition has >= 20 measurable nuclei in that n:
   cut-off = `control_percentile` (99) of their signal. The condition is measured but excluded from the
   statistics (the statistics exclude `experiment.negative_control`).
2. **manual** - `[positivity] manual_cutoff` (a number, or `{ Ki67 = 1200 }`), else a lock in `cutoffs.json`
   written by the tuning view (`tuning.lock`; value, time, fields tuned on; an n-specific lock beats an all-n
   lock). Lock it before looking at the statistics.
3. **limit of detection** (automatic) - positive when signal > `lod_k` (3) x SD(mean - bg) of that nucleus,
   i.e. the noise SD scaled to the nucleus area. On pure background this calls about 0.1-1% positive (tested).
   The cut-off differs per nucleus; its median is reported.

The 90th-percentile variant (dotted early-G1 pattern) uses the same rule on `p90_bgsub_ff`; under the automatic
rule its null level is 1.28 sigma + k x SD(p90 - bg). A manual lock is per measure (`measure="p90"`); without one
the p90 variant falls back to the automatic rule (stated in `cutoffs.csv`).

**Without negative-control wells, "positive" means "detectably stained above background", not "above
non-specific staining"**: secondary-antibody background or autofluorescence that is brighter in nuclei counts as
signal. Add control wells (`kerncount platemap --negative-control "secondary only"`).

### Per-image readouts (the registry contract)

Columns in `qc_details.csv` (one row per image) that flow into every statistic through `metrics.registry`
(labels, kinds and units from `classify.READOUTS.get(column)` / `classify.channel_metric_info`):

| column | kind | meaning |
|---|---|---|
| `readout_frac_positive_<ch>` | fraction | positive / measurable nuclei (mean rule). Empty below `min_nuclei` (5). |
| `readout_frac_positive_p90_<ch>` | fraction | the same with the 90th percentile |
| `intensity_<ch>_nucleus_mean_bgsub` | positive | median over the image's measurable nuclei of `mean_bgsub_ff` |
| `intensity_<ch>_positive_mean_bgsub` | positive | mean `mean_bgsub_ff` of the positive nuclei |
| `intensity_<ch>_fiji_cell` | positive | legacy, only with region `fiji_cell`: mean of (nucleus + 6 px) minus image background over counted nuclei |

Not registered (helpers): `n_measured_<ch>` (denominator of the fractions; use it as N for the logit weights),
`n_saturated_<ch>`, `cutoff_rule_<ch>`, `channel_check`. `<ch>` is the channel name with unsafe characters
replaced by `_`. The existing `intensity_<ch>_nucleus` / `_cell` / `_ring` columns are unchanged.

## 6. Outputs: `<project>/intensity/`

| file | content |
|---|---|
| `manifest.json` | `{schema: "kerncount.intensity.manifest/1", channels, nuclear_channel, cutoffs, flatfield, files: [{path, kind, title, description, n, channel}]}`; kinds: `flatfield`, `flatfield_map`, `flatfield_table`, `cutoffs`, `sweep`, `sweep_table`, `distribution`, `density`, `summary`, `alignment`, `alignment_table` |
| `cutoffs.csv` | n, channel, measure (mean / p90), rule, value, nuclei, % positive (with / without flat field), details |
| `positivity_sweep.csv`, `sweep_<n>__<ch>.png` | % positive per condition when the cut-off is scaled 1/4 to 4 x the one used |
| `distribution_<n>__<ch>.png` | per-nucleus signal per condition (violins + quartiles, % positive) and log histograms |
| `positive_by_condition.csv` | per n, channel, condition: nuclei, saturated left out, % positive (pooled and mean of images), p90 variant, median signal |
| `density_<n>__<ch>.png` | % positive vs nuclear coverage, one dot per image (descriptive; adjusted models are in the statistics) |
| `flatfield.csv`, `flatfield/<n>__<ch>.png/.npy` | the illumination maps and their effect |
| `alignment.csv`, `alignment_<n>.png` | shift, correlation and flags per field |

`_work/nuclei_intensity.csv` is the per-nucleus table (section 3).

## 7. Tuning API (for the browser)

`kerncount.tuning`, pure functions on a project that has been run. Blinded: fields have opaque ids (salted hash,
salt in `_work/tuning_salt`); no function returns a condition name, a file path or a per-condition result.
Cached in memory per project (rebuilt when `nuclei_intensity.csv`, `layout.csv` or `settings.toml` change).
Measured on the 1024 x 1280 Ki-67 fields: `render` about 30-60 ms per call after the first one per field,
`preview` a few ms.

```python
tuning.channels(project_dir) -> ["Ki67"]
tuning.sample(project_dir, k=12, seed=None, channel=None, n=None)
    -> {"seed", "channel", "n", "fields": [{"field_id", "nuclei"}]}    # round robin over conditions, shuffled
tuning.nuclei(project_dir, field_id, channel=None, contours=True)
    -> {"field_id", "channel", "shape": [H, W], "nuclei": [{"label", "x", "y", "area", "signal", "signal_p90",
        "lod_cutoff", "measurable", "saturated", "contour": [[x, y], ...]}]}
tuning.render(project_dir, field_id, channel=None, cutoff=None, measure="mean", outline=True) -> PNG bytes
    # flat-field corrected, contrast-stretched (0.5-99.8 percentile, gamma 0.5); outlines green = positive,
    # grey = negative, orange = saturated (excluded); cutoff=None = the cut-off now in force
tuning.preview(project_dir, cutoff, channel=None, n=None, measure="mean", field_ids=None)
    -> {"cutoff", "channel", "measure", "nuclei", "positive", "frac_positive", "lod_frac_positive"}
tuning.suggest_range(project_dir, channel=None, measure="mean") -> {"min", "max", "lod_median", ...}
tuning.lock(project_dir, channel, cutoff, field_ids, measure="mean", n=None) -> the stored lock
```

Errors: `tuning.TuningError` (unknown field / channel, project not run). Cut-off units: signal = mean inside the
shrunk nucleus minus local background, divided by the flat field. A lock takes effect at the next `kerncount run`
(cached detections and measurements are reused). Note: an existing generation is only rebuilt when the counts
change, so after a new lock run `kerncount run <exp> --rebuild gN` to refresh its statistics. Suggested flow: `sample` (per n), show `render` with a slider over
`suggest_range`, `preview` for the overall %, `lock` with the sampled ids, then run.

## 8. Plate maps

`kerncount platemap <exp> --negative-control "secondary only" --control-wells 2` (default: `[experiment]
negative_control`) reserves 1-2 wells per plate for the control; `platemap.generate(...,
negative_control=..., negative_control_wells=...)`.

## 9. Caveats

- Without control wells the cut-off is a detection limit; the % positive then depends on staining brightness and
  exposure, so compare conditions within an n, and treat absolute values with care.
- The automatic rule assumes the local background ring represents the background under the nucleus. Cytoplasmic
  marker signal raises the ring and lowers the nuclear signal.
- Saturated nuclei are excluded; if many are saturated (see `n_saturated_<ch>`), % positive is biased low.
- The flat field is estimated from the data: with few or very sparse fields it is noisy (check the nuclei-based
  before / after numbers). Reference images are better when available.
- The cut-off sweep shows how robust condition differences are to the cut-off; curves that cross mean the
  ranking depends on it.
- Dense fields often proliferate less; look at the density scatter before reading % positive differences.
