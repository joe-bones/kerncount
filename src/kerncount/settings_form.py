"""Settings as a form: form items (one per setting) for the web interface, and a comment-preserving settings.toml writer.

Standard library only (the web interface must start without the scientific stack).
"""
from __future__ import annotations

import json
import re
import tomllib

from . import experiment_types
from .settings import DEFAULTS

# short names of classify.READOUTS (without "readout_"); kept here because this module must not import numpy/pandas
READOUT_NAMES = ["frac_out_of_focus", "micronuclei_per_nucleus", "nuclear_area_median", "nuclear_area_median_um2",
                 "nuclear_area_cv", "frac_irregular", "coverage", "nn_distance_median"]

# FORM_ITEMS: one form item per setting (not to be confused with a microscope field). kind: bool | int | float | choice | multi | list | num_or (a number, or one of `special`) | measure | names ("all" or some of `options`) | text
# | table (TOML inline table of name = [strings], e.g. { mirna = ["miR-*"] }; shown and edited as text)
# reseg = changing it re-detects every nucleus (the slow step, ~2-3 s per image)
FORM_ITEMS: dict[str, dict[str, dict]] = {
    "project": {
        "raw_dir": dict(kind="text", label="Raw images folder",
                        help="Where the raw images are. Relative paths are relative to the project folder."),
        "date_order": dict(kind="choice", options=["DMY", "MDY", "YMD"], label="Date order in n-folder names",
                           help="DMY: 05.06.26 = 5 June 2026."),
        "csv_delimiter": dict(kind="choice", options=[",", ";"], label="CSV delimiter",
                              help='Use ";" if Excel (German/European locale) shows CSVs in one column.'),
    },
    "experiment": {
        "type": dict(kind="choice", options=sorted(experiment_types.TYPES), label="Experiment type",
                     help="How the raw data are found, measured and analysed. cell_count = nuclei and marker intensity "
                          "on plates; tumour_volume = mCherry tumour z-stacks (one tumour = one n, docs/tumour.md). "
                          "After changing it, scan the raw folder again (init with a fresh layout)."),
        "cell_line": dict(kind="text", label="Cell line",
                          help="e.g. LN229. Experiments (n's) of different cell lines are never pooled. Override per n "
                               "or condition with a cell_line column in layout.csv."),
        "magnification": dict(kind="num_or", special=[""], label="Objective magnification",
                              help="e.g. 20 for a 20x objective. With the camera pixel size it gives the pixel size."),
        "camera_pixel_um": dict(kind="num_or", special=[""], label="Camera pixel size (µm)",
                                help="Physical size of one camera pixel, e.g. 6.5. Pixel size in the image = this ÷ "
                                     "magnification (used when 'Pixel size (µm)' under Finding nuclei is empty)."),
        "seeding_density": dict(kind="num_or", special=[""], label="Seeding density",
                                help="Cells seeded, e.g. 30000. Override per n or condition in layout.csv."),
        "seeding_unit": dict(kind="text", label="Seeding unit", help='e.g. "cells/well" or "cells/cm²".'),
        "stain": dict(kind="text", label="Stain / antibodies",
                      help="e.g. Ki-67 (rabbit) + anti-rabbit AF488, DAPI. For the record only."),
        "notes": dict(kind="text", label="Notes", help="Anything worth knowing about this experiment."),
        "negative_control": dict(kind="text", no_dash=True, label="Negative-control condition",
                                 help="Name of the condition with negative-control wells (stained without the specific "
                                      "signal). If present in an n, its nuclei set that n's positivity cut-off and it is "
                                      "left out of the statistics. Empty = no control wells."),
        "negative_control_type": dict(kind="choice", options=["secondary_only", "isotype", "unstained"],
                                      label="Negative-control type",
                                      help="secondary_only = secondary antibody without the primary (indirect "
                                           "staining); isotype = non-binding antibody of the same isotype; unstained "
                                           "= no antibody (only autofluorescence)."),
    },
    "channels": {
        "names": dict(kind="list", label="Channel names",
                      help='e.g. DAPI, marker. Empty = C1, C2, ... (separate files: empty = the file order below).',
                      reseg=True),
        "segment": dict(kind="text", label="Channel used to find nuclei",
                        help='A number (1 = first channel) or a channel name, e.g. DAPI.', reseg=True),
        "files": dict(kind="choice", options=["single", "interleaved"], label="How channels are saved", reseg=True,
                      help="single = all channels inside one file per image; interleaved = one file per channel, "
                           "alternating in file order (e.g. 1 = Ki67, 2 = DAPI, 3 = Ki67, 4 = DAPI, ...). Gaps in the "
                           "numbering are ignored. Run init again (--force-layout) after changing it."),
        "file_order": dict(kind="list", label="Channel order of the files",
                           help="Interleaved only: channel names in file order, e.g. Ki67, DAPI.", reseg=True),
        "align": dict(kind="choice", options=["check", "correct", "off"], label="Channel alignment",
                      help="Interleaved only: check = measure the shift between the files of a field and flag big "
                           "ones; correct = also shift the marker back (whole pixels, up to the limit below); off."),
        "align_warn_px": dict(kind="float", min=0, label="Alignment warning above (px)",
                              help="Fields whose channels are shifted by more than this are flagged."),
        "align_max_correct_px": dict(kind="float", min=0, label="Correct shifts up to (px)",
                                     help="Bigger shifts are not corrected (they suggest mis-paired files) but flagged."),
    },
    "segmentation": {
        "backend": dict(kind="choice", options=["stardist", "classical"], label="Detector", reseg=True,
                        help="stardist = deep-learning detector (recommended); classical = Fiji-style threshold + watershed."),
        "nucleus_diameter_px": dict(kind="num_or", special=["auto"], label="Nucleus diameter (px)", reseg=True,
                                    help="auto measures nucleus size per folder (adapts to magnification), or a fixed number."),
        "target_diameter_px": dict(kind="float", label="Detector working diameter (px)", reseg=True,
                                   help="Images are rescaled internally so nuclei are this size."),
        "prob_thresh": dict(kind="num_or", special=["default"], label="Detector confidence threshold", reseg=True,
                            help="Lower = more objects. default = the model's own value."),
        "pixel_size_um": dict(kind="num_or", special=[""], label="Pixel size (µm)",
                              help="If known: adds µm² areas and cells/mm². Leave empty if unknown. Does not re-detect."),
    },
    "counting": {
        "edge_rule": dict(kind="choice", options=["half", "exclude", "include"], label="Edge nuclei",
                          help="half = unbiased frame (top/left-touching count, bottom/right don't); exclude = drop all; include = count all."),
        "min_area_frac": dict(kind="float", label="Debris cut-off (× typical area)",
                              help="Smaller objects are debris / micronuclei and are not counted."),
        "exclude_out_of_focus": dict(kind="bool", label="Drop out-of-focus nuclei",
                                     help="Nuclei outside the focal plane (floating or dead cells) are not counted."),
        "focus_frac": dict(kind="float", label="Focus: edge sharpness (× typical)",
                           help="Out of focus = sharpness below this ..."),
        "texture_frac": dict(kind="float", label="Focus: chromatin texture (× typical)",
                             help="... AND texture below this. Both must be low."),
        "merge_area_frac": dict(kind="float", label="Possible merge (× typical area)",
                                help="Larger and non-convex = flagged as a possible merge (counted as 1)."),
    },
    "qc": {
        "high_max_uncertainty": dict(kind="float", label="High: max plausible range", help="Fraction of the count, e.g. 0.10 = 10%."),
        "high_max_disagreement": dict(kind="float", label="High: max disagreement with 2nd method", help="e.g. 0.15 = 15%."),
        "max_saturated_frac": dict(kind="float", label="High: max saturated nuclear pixels", help="e.g. 0.02 = 2%."),
        "low_min_uncertainty": dict(kind="float", label="Low: plausible range at least", help="e.g. 0.25 = 25%."),
        "low_min_disagreement": dict(kind="float", label="Low: disagreement at least", help="e.g. 0.35 = 35%."),
    },
    "intensity": {
        "enabled": dict(kind="bool", label="Measure intensities", reseg=True, help="Per-cell mean intensity (replaces the Mean-intensity macros)."),
        "measure": dict(kind="measure", label="Channels to measure", reseg=True, help='"all", or channel numbers, e.g. 2 or 2, 3.'),
        "regions": dict(kind="multi", options=["nucleus", "ring", "cell", "fiji_cell"], label="Regions", reseg=True,
                        help="ring = band around the nucleus (cytoplasm); cell = nucleus + ring; fiji_cell = nucleus "
                             "enlarged by the Fiji ring width below (the old macro's 'Enlarge 6 px'). Changing only "
                             "fiji_cell does not re-detect."),
        "ring_px": dict(kind="num_or", special=["auto"], label="Ring width (px)", reseg=True, help="auto = 0.25 × nucleus diameter."),
        "subtract_background": dict(kind="bool", label="Subtract background", reseg=True, help="Median of empty areas, per image."),
        "per_nucleus": dict(kind="text", label="Marker channels measured per nucleus",
                            help='auto = every channel except the nuclear one; or names, e.g. Ki67; off = none. '
                                 "Per nucleus: mean, total, 90th percentile, local background. Does not re-detect."),
        "erode_px": dict(kind="num_or", special=["auto"], label="Shrink nucleus outline by (px)",
                         help="Measured inside the outline shrunk by this much, so the blurry rim and neighbours' "
                              "light are left out. auto = 1 px, 2 px for nuclei at least 16 px across."),
        "bg_gap_px": dict(kind="int", min=0, label="Background: gap around nuclei (px)",
                          help="The local background ring starts this far outside the nucleus and keeps this far "
                               "from every other nucleus."),
        "bg_ring_px": dict(kind="num_or", special=["auto"], label="Background ring width (px)",
                           help="auto = half a nucleus diameter. The median of the ring is the nucleus' background."),
        "bg_min_px": dict(kind="int", min=1, label="Background: minimum ring pixels",
                          help="Crowded nuclei with fewer ring pixels get a wider ring, then the image background."),
        "fiji_ring_px": dict(kind="int", min=1, label="Fiji ring width (px)",
                             help="Region fiji_cell = nucleus enlarged by this (the old macro used 6 px)."),
        "flatfield": dict(kind="text", label="Flat-field correction",
                          help="Uneven illumination (e.g. darker edges). auto = estimated per n and channel from all "
                               "of that n's images; off = none; or a folder of reference images (uniform slide; "
                               "files named with the channel; 'dark' in the name = dark frames). Intensities only, "
                               "never counts."),
        "flatfield_dark": dict(kind="bool", label="Flat field: remove background pattern first",
                               help="Estimate and remove the additive background (camera offset, autofluorescence) "
                                    "pattern before estimating the illumination."),
    },
    "readouts": {
        "enabled": dict(kind="bool", label="Extra readouts",
                        help="Per-image readouts besides the count (out-of-focus fraction, micronuclei, nuclear size and "
                             "shape, crowding), computed from the detected nuclei. Fast: does not re-detect."),
        "include": dict(kind="names", options=READOUT_NAMES, label="Readouts to include",
                        help='"all", or names separated by commas: ' + ", ".join(READOUT_NAMES) + "."),
    },
    "positivity": {
        "enabled": dict(kind="bool", label="% positive nuclei",
                        help="Classify each nucleus as positive or negative for the marker; per image the % positive "
                             "(primary readout) and the same with the 90th percentile."),
        "channels": dict(kind="text", label="Channels", help="auto = every per-nucleus channel; or names, e.g. Ki67."),
        "control_percentile": dict(kind="float", min=50, max=100, label="Negative control: percentile",
                                   help="With negative-control wells in an n: cut-off = this percentile of their nuclei "
                                        "(99 = 1% of control nuclei would be called positive)."),
        "manual_cutoff": dict(kind="text", label="Manual cut-off",
                              help="A fixed cut-off (corrected intensity units), used when there are no control wells. "
                                   "Empty = the value locked in the tuning view (cutoffs.json), else the automatic "
                                   "rule. Fix it before looking at the statistics."),
        "lod_k": dict(kind="float", min=0.5, max=20, label="Automatic rule: k x noise",
                      help="Limit of detection: a nucleus is positive when its signal exceeds k x the background "
                           "noise expected for its size. 3 = about 1 false positive in 700 unstained nuclei."),
        "min_nuclei": dict(kind="int", min=1, label="Minimum nuclei per image",
                           help="Images with fewer measurable nuclei get no % positive."),
    },
    "analysis": {
        "reference": dict(kind="text", label="Main comparison group", help="Used for the headline plots."),
        "normalise_to": dict(kind="list", label="Control groups",
                             help="Each gets its own % plots and tests, e.g. Control, LyoVec."),
        "two_way_anova": dict(kind="bool", label="Two-way ANOVA (extra)",
                              help="Per n: treatment × well; all n: treatment × n; Tukey all pairs."),
        "treatment_order": dict(kind="list", label="Treatment order", help="Order in tables and plots; others follow alphabetically."),
        "exclude_low_qc": dict(kind="bool", label="Leave out QC-low images",
                               help="Images with QC low (and no manual value) are left out of the statistics."),
        "primary_test": dict(kind="choice", options=["mixed", "rm_anova"], label="Test behind the stars",
                             help="mixed = hierarchical model on the log scale (fields in wells in n; each treatment "
                                  "class gets its own n-to-n variability; results as ratios). rm_anova = Prism-style "
                                  "repeated-measures ANOVA + Dunnett on the per-n means. Both are always computed; this "
                                  "only decides which one draws the stars. If the model cannot be fitted, rm_anova is "
                                  "used and the plot says so."),
        "primary_family": dict(kind="text", label="Main claim: compare with",
                               help="The control whose comparisons are the main result (headline stars), e.g. LyoVec. "
                                    "Empty = the main comparison group above."),
        "vehicle_check": dict(kind="text", label="Vehicle check (equivalence test)",
                              help='e.g. "LyoVec vs Control": tests whether the vehicle has NO meaningful effect '
                                   "(equivalence, not just 'not significant'). auto = every other control group vs the "
                                   "main comparison group; empty = off. Several: separate with ;"),
        "equivalence_margin": dict(kind="float", min=0.01, max=1, label="Equivalence margin",
                                   help="0.2 = the vehicle counts as equivalent if its ratio to the control (90% CI; "
                                        "with several vehicle pairs the Holm-matched CI) lies entirely within 1/1.2 "
                                        "to 1.2 (about ±20%)."),
        "equivalence_margin_pp": dict(kind="float", min=0.1, max=50, label="Equivalence margin, fractions (points)",
                                      help="For % positive and other fractions: equivalent if the difference vehicle "
                                           "minus control lies entirely within ± this many percentage points "
                                           "(default 5). Odds ratios are not used for this."),
        "positive_control": dict(kind="text", label="Positive control",
                                 help="e.g. CL075. Every treatment that is not a control is also compared with it. "
                                      "Empty = none."),
        "model_unit": dict(kind="choice", options=["well", "field"], label="Model unit",
                           help="well = each well's mean per field (pooled % for fractions) is one value: results are "
                                "ratios of mean nuclei per field, as biologists read them; missing fields are "
                                "position-adjusted. field = every field is one value (older model; compares geometric "
                                "means of fields, which uneven, sparse fields pull down)."),
        "classes": dict(kind="table", special=["auto"], label="Treatment classes",
                        help="Groups of treatments whose effects may vary differently between n's, each with its own "
                             "n-to-n variance in the model. auto (recommended) = every treatment, each control included, "
                             "gets its own class (conservative: each comparison rests on the n-to-n variability of "
                             "its own two groups). Or a TOML inline table of name = [patterns], * = wildcard, e.g. "
                             '{ mirna = ["miR-*"] }. Empty {} = controls form one class, everything else another.'),
        "density_covariate": dict(kind="text", label="Density covariate (fractions)",
                                  help="For every fraction metric (% out of focus, % Ki-67-positive) a second model "
                                       "adjusts for local crowding with this per-field column (default "
                                       "readout_coverage). The raw estimate stays the headline. Empty = off."),
        "among_classes": dict(kind="list", label="Compare within classes",
                              help="Class names whose members are compared with each other (all pairs), e.g. mirna."),
    },
    "acceptance": {
        "vehicle_within": dict(kind="num_or", special=[""], label="Vehicle within ± (ratio)",
                               help="Per n, the vehicle / control ratio of nuclei per field must lie within 1 ± this "
                                    "(0.3 = 0.7 to 1.3). Uses the vehicle check pair. Empty = off."),
        "positive_control_max_ratio": dict(kind="num_or", special=[""], label="Positive control must reach (ratio)",
                                           help="Per n, positive control / main comparison group must be at most this "
                                                "ratio (e.g. 0.8 = at least 20% fewer nuclei); a value above 1 means "
                                                "at least this ratio. Needs a positive control. Empty = off."),
        "min_images_per_condition": dict(kind="int", min=0, label="Minimum images per condition",
                                         help="Every condition in that n needs at least this many usable images."),
        "max_low_qc_fraction": dict(kind="float", min=0, max=1, label="Maximum share of QC-low images",
                                    help="e.g. 0.5 = at most half of that n's images may have QC low."),
        "no_duplicate_images": dict(kind="bool", label="No duplicate images",
                                    help="An n fails if one of its raw image files is byte-identical to an image of "
                                         "another condition (a copied folder). Listed in duplicates.csv."),
        "enforce": dict(kind="bool", label="Leave out n's that fail",
                        help="Off = only report (acceptance.csv). On = n's that fail any criterion are left out of "
                             "every all-n statistic; they are listed in the About sheet and the manifest."),
    },
    "design": {
        "wells": dict(kind="num_or", special=["auto"], label="Wells per condition per n",
                      help="Usually 2. auto = take the wells from wells.csv in the project folder if it exists, else 2."),
        "fields_per_well": dict(kind="num_or", special=["auto"], label="Fields (images) per well",
                                help="auto = per condition, highest image number ÷ wells, rounded up (11 images = 6 + 5, "
                                     "18 images = 9 + 9). A number fixes it; images beyond wells × fields are flagged."),
        "well_layout": dict(kind="choice", options=["consecutive", "alternating", "file"], label="Well layout",
                            help="consecutive = images 1-6 well 1, 7-12 well 2; alternating = odd/even; "
                                 "file = wells.csv in the project folder (columns image, well, position; optional n, treatment)."),
        "field_positions": dict(kind="choice", options=["fixed", "random"], label="Field positions",
                                help="fixed = field 1 of every well is the same spot (e.g. centre), so position can "
                                     "explain differences between fields; random = fields were taken anywhere, "
                                     "position is ignored."),
    },
    "review": {
        "overlay_alpha": dict(kind="float", label="Overlay transparency", help="0 = invisible, 1 = opaque."),
        "jpeg_quality": dict(kind="int", min=1, max=100, label="Review image JPEG quality"),
    },
    "planner": {
        "effects": dict(kind="list", label="Effects to plan for (ratios)",
                        help="Treatment ÷ reference ratios, e.g. 0.8, 0.7, 0.5 (0.7 = 30% fewer nuclei)."),
        "target_power": dict(kind="float", min=0.05, max=0.999, label="Target power",
                             help="Chance of finding a real effect (p < alpha); 0.8 is usual."),
        "alpha": dict(kind="float", min=0.0001, max=0.5, label="Significance level (alpha)", help="Usually 0.05."),
        "cost_per_n": dict(kind="float", min=0, label="Effort per n (plate / day)",
                           help="Seeding, transfection, staining for one independent experiment. Any unit, e.g. minutes."),
        "cost_per_well": dict(kind="float", min=0, label="Effort per well", help="Same unit; pipetting and handling one more well."),
        "cost_per_field": dict(kind="float", min=0, label="Effort per field (image)",
                               help="Same unit; taking and checking one more image."),
        "max_n": dict(kind="int", min=3, max=100, label="Largest number of n's in the tables"),
        "plate_format": dict(kind="int", min=6, max=96, label="Plate format (wells)", help="6, 12, 24, 48 or 96."),
        "avoid_edge": dict(kind="bool", label="Plate maps: avoid edge wells",
                           help="Use inner wells first; edge wells (evaporation) only when the inner ones are not enough."),
        "reference_wells": dict(kind="num_or", special=["auto"], label="Wells per reference condition",
                                help="auto = Dunnett's rule: wells × √(number of comparisons), as far as the plate allows."),
        "adjust_for_family": dict(kind="bool", label="Power for a star on the headline plot (family-adjusted)",
                                  help="On: power and smallest detectable effect use the multiple-comparison adjustment "
                                       "of the statistics (all comparisons vs the reference). Off: one unadjusted test."),
    },
    # tumour experiments ([experiment] type = "tumour_volume"; docs/tumour.md). reseg = the next run measures every
    # stack again (segmentation + readouts, ~20-60 s per stack); the cache key notices the change by itself.
    # The analysis plan (factors, readouts, families: arrays of tables) has its own editor (gui_tumour), not items here.
    "tumour": {
        "smooth_sigma_um": dict(kind="float", min=0, label="Smoothing before the threshold (µm)", reseg=True,
                                help="Gaussian smoothing at a physical scale, identical for every tumour, so noise and "
                                     "single bright pixels do not change the outline. 2 µm = about 3 pixels in xy."),
        "threshold": dict(kind="choice", options=["intermedians", "otsu"], label="Threshold rule", reseg=True,
                          help="One automatic rule for every tumour (ADR 0002; never set by hand per tumour). "
                               "intermedians = half-way between the background and the tumour level (medians), so "
                               "it does not depend on tumour size and ignores saturated pixels; otsu = classical "
                               "two-class split (shifts with tumour size). Check its robustness with the threshold "
                               "sweep."),
        "min_object_um3": dict(kind="float", min=0, label="Smallest object (µm³)", reseg=True,
                               help="Connected objects smaller than this are noise. Larger ones besides the main "
                                    "mass are satellites (tumour burden, invasion fraction)."),
        "mid_section_um": dict(kind="float", min=0, label="Mid-section thickness (µm)", reseg=True,
                               help="The confirmatory area and perimeter are averaged over the layers within this "
                                    "thickness, centred on the layer where the main mass is largest - the same "
                                    "physical slab for every tumour, whatever the imaged depth. 20 µm = 9-11 layers."),
        "outline_smooth_um": dict(kind="float", min=0, label="Outline smoothing for perimeters (µm)", reseg=True,
                                  help="The outline is smoothed at this scale before its length is measured: at 8 µm "
                                       "the perimeter measures the tumour's shape, not pixel noise or single cells "
                                       "at the border (a finer scale moves it with the threshold)."),
        "roughness_fine_um": dict(kind="float", min=0, label="Edge roughness: fine scale (µm)", reseg=True,
                                  help="Edge roughness = outline length at this smoothing ÷ at the coarse scale "
                                       "below (1 = smooth border; higher = ragged, infiltrating). Exploratory."),
        "roughness_smooth_um": dict(kind="float", min=0, label="Edge roughness: coarse scale (µm)", reseg=True,
                                    help="The heavily smoothed reference outline of the edge roughness."),
        "slab_um": dict(kind="float", min=0, label="Slab depth for the slab volume (µm)", reseg=True,
                        help="Exploratory slab volume = main-mass volume within this depth around the mid-section."),
        "surface_grid_um": dict(kind="float", min=0.1, label="3D surface grid (µm)", reseg=True,
                                help="Isotropic grid of the exploratory 3D mesh surface area."),
        "normalise_by": dict(kind="text", label="Normalise to the reference of factor",
                             help="Normalised values = % of the model-estimated reference level of this factor in "
                                  "the same cell of the other factors (e.g. tumour_line: % of cntrl with the same "
                                  "genotype and clodronate). Empty = the factor on the x axis. Part of the locked "
                                  "analysis plan."),
        "threshold_sweep": dict(kind="bool", label="Threshold sweep",
                                help="Recompute every confirmatory effect with the threshold × each scale below and "
                                     "report whether a conclusion depends on the threshold (ADR 0002). On by "
                                     "default. Slow on the first run (every stack segmented once more; cached "
                                     "afterwards); switch off to save that time."),
        "sweep_scales": dict(kind="list", numbers=True, label="Threshold sweep scales",
                             help="Factors applied to each tumour's automatic threshold, e.g. 0.8, 0.9, 1.1, 1.2."),
        "shared_threshold_check": dict(kind="bool", label="Shared-threshold robustness check",
                                       help="Re-segment every tumour at ONE threshold (median of the per-tumour "
                                            "thresholds) and report whether a confirmatory effect depends on the "
                                            "per-tumour adaptation, i.e. on brightness (ADR 0002). Only run when the "
                                            "acquisition settings are identical. Slow on the first run."),
        "batch": dict(kind="choice", options=["culture_batch", "session_folder", "acquisition_day"], label="Batch",
                      help="What counts as one batch in the statistics. culture_batch (default) = the culture "
                           "preparation, the date at the start of the file names (tissue fixed before imaging, so the "
                           "imaging day is not a batch); session_folder = each session folder; acquisition_day = each "
                           "day in the stacks' OME metadata. Batches that share a group enter the model as a fixed "
                           "block effect; a contrast between groups of unconnected batches is shown but never tested. "
                           "Part of the locked analysis plan; decide before results are viewed."),
        "qc_saturation_warn": dict(kind="float", min=0, max=1, label="Saturation warning above (share)",
                                   help="A tumour gets a warning (never a flag) when more of its mid-section main mass "
                                        "is saturated (255) than this; the threshold rule is protected up to ~0.5."),
        "review_px": dict(kind="int", min=100, max=4000, label="Review image size (px per panel)",
                          help="Longest side of each panel of the blinded review images (review/<code>.jpg)."),
        "review_jpeg_quality": dict(kind="int", min=1, max=100, label="Review image JPEG quality"),
        "fiji_legacy": dict(kind="bool", label="Fiji legacy readout",
                            help="Read the lab's Fiji outputs (masks, 3D Objects Counter tables) as fiji_* columns "
                                 "and report the new method side by side with Fiji (method concordance). Legacy / "
                                 "exploratory: never starred."),
        "fiji_fallback": dict(kind="bool", label="Fiji: re-run the macro when outputs are missing",
                              help="A tumour without Fiji statistics or mask gets the macro re-implementation run on "
                                   "its original stack (slow)."),
        "fiji_min_voxels": dict(kind="int", min=1, label="Fiji: smallest object (voxels)",
                                help="Minimum object size of the macro's 3D Objects Counter (200000 in the macro)."),
        "sensitivity_permutations": dict(kind="int", min=99, label="Sensitivity: permutations per contrast",
                                         help="Random permutations of the permutation test (enumerated exactly when "
                                              "there are fewer possible splits). Sensitivity only."),
        "sensitivity_seed": dict(kind="int", min=0, label="Sensitivity: random seed",
                                 help="Same seed = same permutation p-values; printed in the report."),
        "planner_fold_change": dict(kind="float", min=1.01, max=100, label="Planning aid: fold-change to detect",
                                    help="The power planner reports how many tumours per group are needed to detect "
                                         "this fold-change (e.g. 1.5 = 50 % larger or 1/1.5 smaller). Planning only."),
        "planner_power": dict(kind="float", min=0.05, max=0.999, label="Planning aid: target power"),
        "planner_alpha": dict(kind="float", min=0.0001, max=0.5, label="Planning aid: alpha per family"),
    },
}

# Sections (and, per section, the form items) shown for each experiment type. Types not listed show every section
# except those that belong to another type only.
TYPE_SECTIONS: dict[str, dict[str, list[str] | None]] = {
    "tumour_volume": {"project": None, "experiment": ["type", "notes"], "tumour": None},
}
TYPE_ONLY = {"tumour": "tumour_volume"}  # sections of one experiment type only

SECTION_TITLES = {
    "project": "Project", "experiment": "Experiment (cell line, imaging, controls)", "channels": "Channels", "segmentation": "Finding nuclei", "counting": "What counts as a nucleus",
    "qc": "Quality levels", "intensity": "Intensity per cell / nucleus",
    "positivity": "% positive nuclei (e.g. Ki-67)", "readouts": "Extra readouts (beyond the count)",
    "analysis": "Statistics", "acceptance": "Acceptance criteria per n",
    "design": "Plate design (wells, fields, positions)", "review": "Review images",
    "planner": "Planning the next experiment",
    "tumour": "Tumour measurement, QC and robustness",
}


def _defaults(section: str) -> dict:
    if section == "tumour":
        from .tumour.config import DEFAULTS as TUMOUR_DEFAULTS

        return TUMOUR_DEFAULTS
    return DEFAULTS[section]


def schema(etype: str = "cell_count") -> list[dict]:
    """Form sections of an experiment type (default cell_count: every section it has always had)."""
    shown = TYPE_SECTIONS.get(etype)
    out = []
    for s, keys in FORM_ITEMS.items():
        if shown is not None and s not in shown or shown is None and TYPE_ONLY.get(s, etype) != etype:
            continue
        only = shown.get(s) if shown is not None else None
        out.append(dict(section=s, title=SECTION_TITLES.get(s, s),
                        items=[dict(key=k, default=_defaults(s)[k], **meta) for k, meta in keys.items()
                               if only is None or k in only]))
    return out


def _num(v, as_int=False):
    if isinstance(v, bool):
        raise ValueError
    f = float(str(v).strip().replace(",", "."))
    if as_int:
        if f != int(f):
            raise ValueError
        return int(f)
    return int(f) if f == int(f) and "." not in str(v) else f


def coerce(values: dict) -> dict:
    """Turn form values (strings, booleans, lists) into typed settings; raise ValueError naming the bad form item."""
    out: dict[str, dict] = {}
    for sec, keys in values.items():
        if sec not in FORM_ITEMS:
            continue
        for key, v in keys.items():
            meta = FORM_ITEMS[sec].get(key)
            if meta is None:
                continue
            k, where = meta["kind"], f"{SECTION_TITLES[sec]} → {meta.get('label', key)}"
            if meta.get("no_dash") and str(v).strip().startswith("-"):
                raise ValueError(f"{where}: may not start with - (it would be read as a command-line option)")
            try:
                if k == "bool":
                    val = bool(v)
                elif k in ("int", "float"):
                    val = _num(v, as_int=k == "int")
                    if "min" in meta and val < meta["min"] or "max" in meta and val > meta["max"]:
                        raise ValueError
                elif k == "choice":
                    if v not in meta["options"]:
                        raise ValueError
                    val = v
                elif k == "multi":
                    val = [o for o in meta["options"] if o in v]
                    if not val:
                        raise ValueError
                elif k == "list":
                    items = v if isinstance(v, list) else str(v).split(",")
                    strs = [str(i).strip() for i in items if str(i).strip()]
                    val = [_num(i) for i in strs] if meta.get("numbers") else strs
                elif k == "num_or":
                    s = str(v).strip()
                    val = s if s in meta["special"] else _num(s)
                elif k == "table":
                    sv = v.strip() if isinstance(v, str) else None
                    val = sv.lower() if sv is not None and sv.lower() in meta.get("special", []) else parse_table(v)
                elif k == "measure":
                    s = str(v).strip().lower()
                    val = "all" if s in ("", "all") else [_num(c, as_int=True) for c in s.replace(";", ",").split(",") if c.strip()]
                elif k == "names":  # "all", or a subset of `options` (the readout_ prefix is optional)
                    items = v if isinstance(v, list) else str(v).replace(";", ",").split(",")
                    items = [str(i).strip().removeprefix("readout_") for i in items if str(i).strip()]
                    if not items or "all" in items:
                        val = "all"
                    else:
                        if any(i not in meta["options"] for i in items):
                            raise ValueError
                        val = [o for o in meta["options"] if o in items]
                else:
                    val = str(v).strip()
            except (ValueError, TypeError):
                raise ValueError(f"{where}: '{v}' is not a valid value") from None
            out.setdefault(sec, {})[key] = val
    return out


def parse_table(v) -> dict:
    """kind "table": a TOML inline table of name = [strings] (dict or its text form, e.g. '{ mirna = ["miR-*"] }')."""
    if isinstance(v, dict):
        t = v
    else:
        s = str(v).strip()
        if not s:
            return {}
        if not s.startswith("{"):
            s = "{ " + s + " }"
        t = tomllib.loads(f"x = {s}")["x"]
    if not isinstance(t, dict):
        raise ValueError
    out = {}
    for name, pats in t.items():
        if isinstance(pats, str):
            pats = [pats]
        if not isinstance(pats, list) or not all(isinstance(p, str) for p in pats):
            raise ValueError
        out[str(name)] = [p.strip() for p in pats if p.strip()]
    return out


def table_text(v: dict) -> str:
    """Display form of a kind "table" value ("auto" stays a word)."""
    if isinstance(v, str):
        return v
    return toml_value(v or {})


def toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(i) for i in v) + "]"
    if isinstance(v, dict):
        if not v:
            return "{}"
        keys = [k if re.fullmatch(r"[A-Za-z0-9_\-]+", k) else json.dumps(k) for k in v]
        return "{ " + ", ".join(f"{k} = {toml_value(x)}" for k, x in zip(keys, v.values())) + " }"
    return json.dumps(str(v), ensure_ascii=False)  # a JSON string is a valid TOML basic string


_SECTION = re.compile(r"^\s*\[([A-Za-z0-9_\-]+)\]\s*(#.*)?$")
_ARRAY = re.compile(r"^\s*\[\[\s*([A-Za-z0-9_.\-]+)\s*\]\]\s*(#.*)?$")  # [[a.b]]: array-of-tables entry
_HEADER = re.compile(r"^\s*\[")
_KEY = re.compile(r"^(\s*)([A-Za-z0-9_\-]+)(\s*=\s*)(.*)$")


def _split_comment(rest: str) -> tuple[str, str]:
    """'"a#b"   # note' -> ('"a#b"', '   # note'): find a # outside quotes."""
    quote = None
    i = 0
    while i < len(rest):
        c = rest[i]
        if quote:
            if c == "\\" and quote == '"':
                i += 1
            elif c == quote:
                quote = None
        elif c in "\"'":
            quote = c
        elif c == "#":
            value = rest[:i].rstrip()
            return value, rest[len(value):]
        i += 1
    return rest.rstrip(), ""


def update_text(text: str, typed: dict) -> str:
    """Write `typed` values into settings.toml text, keeping every comment and every key not being changed."""
    lines = text.splitlines()
    todo = {s: dict(kv) for s, kv in typed.items()}
    section, last_line_of = None, {}
    for i, line in enumerate(lines):
        m = _SECTION.match(line)
        if m:
            section = m.group(1)
            last_line_of[section] = i
            continue
        if _ARRAY.match(line) or _HEADER.match(line):  # [[tumour.factors]] etc.: keys below belong to that entry
            section = None
            continue
        m = _KEY.match(line)
        if m and section:
            last_line_of[section] = i
            key = m.group(2)
            if key in todo.get(section, {}):
                value, comment = _split_comment(m.group(4))
                new = toml_value(todo[section].pop(key))
                if comment and len(value) > len(new):  # keep comments roughly aligned
                    new = new.ljust(len(value))
                lines[i] = f"{m.group(1)}{key}{m.group(3)}{new}{comment}"
    # keys or sections not in the file yet: append (insert from the bottom so line numbers stay valid)
    for sec in sorted((s for s in todo if todo[s] and s in last_line_of), key=lambda s: -last_line_of[s]):
        at = last_line_of[sec] + 1
        lines[at:at] = [f"{k} = {toml_value(v)}" for k, v in todo[sec].items()]
    for sec, kv in todo.items():
        if kv and sec not in last_line_of:
            lines += ["", f"[{sec}]"] + [f"{k} = {toml_value(v)}" for k, v in kv.items()]
    out = "\n".join(lines) + "\n"
    parsed = tomllib.loads(out)  # never write a file that cannot be read back
    for sec, kv in typed.items():
        for k, v in kv.items():
            if parsed.get(sec, {}).get(k) != v:
                raise ValueError(f"could not write {sec}.{k} (please edit settings.toml directly)")
    return out


def replace_array(text: str, name: str, rows: list[dict]) -> str:
    """Replace every `[[name]]` entry (e.g. tumour.factors) in settings.toml text by `rows`, written where the first
    entry was (appended at the end when there was none). Comments above the first entry and everything outside the
    entries are kept; comments inside the entries are not."""
    lines = text.splitlines()
    starts = [i for i, ln in enumerate(lines) if (m := _ARRAY.match(ln)) and m.group(1) == name]
    keep: list[str] = []
    at: int | None = None
    i = 0
    while i < len(lines):
        if i in starts:
            if at is None:
                at = len(keep)
            j = i + 1
            while j < len(lines) and not _HEADER.match(lines[j]):
                j += 1
            k = j  # comment / blank lines right above the next header belong to it
            while k > i + 1 and (not lines[k - 1].strip() or lines[k - 1].lstrip().startswith("#")):
                k -= 1
            keep += lines[k:j]
            i = j
            continue
        keep.append(lines[i])
        i += 1
    block: list[str] = []
    for r in rows:
        block += [f"[[{name}]]"] + [f"{k} = {toml_value(v)}" for k, v in r.items()] + [""]
    if at is None:
        keep += ([""] if keep and keep[-1].strip() else []) + block
    else:
        while at < len(keep) and not keep[at].strip():  # one blank line after the entries, not several
            del keep[at]
        keep[at:at] = block
    out = "\n".join(keep).rstrip("\n") + "\n"
    parsed: object = tomllib.loads(out)
    for part in name.split("."):
        parsed = parsed.get(part, {}) if isinstance(parsed, dict) else {}
    if (parsed or []) != rows:
        raise ValueError(f"could not write [[{name}]] (please edit settings.toml directly)")
    return out
