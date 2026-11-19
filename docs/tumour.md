# Tumour volume experiments (`[experiment] type = "tumour_volume"`)

mCherry tumours in brain sections, imaged as z-stacks. One **tumour** (hemisphere A or B of a section) is one n.
Vocabulary: `CONTEXT.md` (Tumour experiments). Decisions: ADR 0002 (measurement), 0003 (statistics), 0004 (types,
factor tables). Status: tracer bullet (ticket 02) — layout, primary segmentation, main-mass volume, one group plot.
Statistics, shape readouts, QC flags, Fiji legacy readout and validation follow in tickets 03–07.

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
- `plots/main_mass_volume_absolute.png|pdf` — one row of panels per genotype present, −/+clodro panels, tumour
  line on x, every tumour a dot, geometric-mean line, n under each group, log y. No statistics yet (ticket 04).
- `_work/tumour/<slug>.json` + `<slug>.labels.tif` — cache per tumour (zlib uint8/uint16 labels, ~1 MB for the
  largest stack). Re-runs reuse it unless the file or the segmentation settings changed; `--resegment` forces it.

## Memory
One stack at a time: uint8 stack + float32 smoothed copy, then mask + int32 labels; peak ≈ 10 × the stack size
(2.8 GB for the largest, 280 MB stack). Hole filling labels the background instead of iterative dilation.
