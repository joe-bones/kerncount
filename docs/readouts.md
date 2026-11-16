# Extra readouts (beyond the count)

A lower nucleus count can mean that cells stopped dividing, died, or detached. The detector already sorts every
object it finds (counted nucleus, possible merge, out of focus, debris, edge) and measures its size, shape and
position. The **readouts** turn those measurements into a few numbers per image that help tell death, damage and
growth apart, without new imaging.

- Columns `readout_*` in `qc_details.csv`, one value per image. The statistics treat them like the count.
- Computed from the classified objects only: **no re-detection**, no pixel intensities. Turning them on or off, or
  changing which ones are included, never invalidates the detection cache.
- Settings: `[readouts]` in `settings.toml`: `enabled = true`, `include = "all"` or a list such as
  `["frac_out_of_focus", "coverage"]` (the `readout_` prefix is optional).
- Units are relative wherever possible (fractions, multiples of the nucleus diameter), so images taken at different
  magnifications stay comparable. The only exception is the median area in px².
- A readout needs a minimum number of nuclei in the image (*min nuclei* below). Below that it is left empty (NaN)
  rather than giving a wild value from 2-3 nuclei. In the test data (371 images) this empties about 6 % of images
  for the median area and nearest-neighbour distance (≥ 5 nuclei) and about 15 % for the others (≥ 10), almost
  all of them sparse miRNA fields.
- The readouts inherit the image's quality. A defocused or glare-filled field (QC low) gives unreliable shapes and
  sizes as well as an unreliable count. Check the review image before trusting an extreme value.
- Two per-object columns are added to `_work/cells.csv`: `nn_distance_px` (counted nuclei: distance to the
  nearest in-focus nucleus) and `micronucleus_candidate` (debris objects that pass the micronucleus rule).

"Whole nuclei" below means objects classified *nucleus* that touch no image border. Nuclei cut by the border and
suspected merges are left out of size and shape statistics, because their outline is not a real nucleus outline.

## readout_frac_out_of_focus: out-of-focus fraction

- **Measures:** the share of nuclei outside the focal plane:
  out of focus / (counted + out of focus). *min nuclei:* 10 (counted + out of focus, so a field where
  every cell floats still scores 1.0).
- **Biology:** adherent U251 nuclei lie flat in one plane. Cells that round up (mitosis, dying cells) or float
  (detached or dead) sit above it and look blurred. A treatment that kills cells or makes them detach raises this
  fraction before the cells are washed away.
- **Caveats:** it uses the classifier's focus rule (blurred edges *and* no chromatin texture), so it moves with
  `counting.focus_frac` / `texture_frac`. A field that is partly defocused (tilted plate, wrong focus) also
  raises it, and is usually flagged QC low. It is 0 in most control images, so compare it as a *difference* from
  the reference, not as a ratio (a ratio to 0 is undefined).

## readout_micronuclei_per_nucleus: micronucleus candidates per nucleus

- **Measures:** debris objects that are (a) at most 1/9 of the image's median nucleus area (diameter ≤ 1/3 of the
  nucleus, the usual micronucleus size criterion) and (b) centred within one nucleus diameter of a nucleus centre
  (about half a diameter from its rim, i.e. inside the same cell). The number of such objects is divided by the
  counted nuclei. *min nuclei:* 10.
- **Biology:** micronuclei are chromosome fragments or whole chromosomes left outside the nucleus after a faulty
  division. They indicate DNA damage, chromosomal instability, or mitotic errors. Fragments of a nucleus breaking
  up during apoptosis (karyorrhexis) also pass the rule, so read it as "small DNA bodies next to nuclei" (damage or
  death).
- **How debris is separated from micronuclei:** the classifier itself only knows "debris" (smaller than
  `counting.min_area_frac` × typical nucleus). Dust and background specks are spread at random, while micronuclei
  hug a nucleus, so the proximity rule removes most dust. With the proximity rule, randomly placed dust adds
  roughly a constant amount *per nucleus*, so this readout is not inflated in sparse fields. (Debris per nucleus
  without the proximity rule *was*: see "Dropped".) In the test data about 40 % of debris objects lay within one diameter
  of a nucleus centre, compared with 25 % of random positions in the same images. So many are real
  nucleus-associated bodies, but not all of them.
- **Caveats:** this is a **lower bound**. The detector rarely finds objects smaller than about 40 px², so small
  micronuclei (visible by eye as tiny bright dots) are often missed, and a micronucleus touching its nucleus is
  usually merged into it. At lower magnification even fewer are found, so compare it within an n, not across
  magnifications. The values are sparse (about 0.5 candidates per image), so a single image says little; use it at
  the condition level.

## readout_nuclear_area_median (px²) and readout_nuclear_area_median_um2

- **Measures:** median area of whole nuclei. *min nuclei:* 5. The µm² column only appears when the pixel size is
  known (`segmentation.pixel_size_um`, or read from the image file), the same rule as `nuclei_per_mm2`.
- **Biology:** nuclei grow during S/G2 and in cells arrested in G2/M or becoming polyploid. They shrink and condense
  in apoptosis (pyknosis) and in senescent or quiescent states in some lines.
- **Caveats:** px² depends on magnification and camera binning. Compare it as a ratio to the reference *within one
  n*. Nuclei also flatten and spread when the cells are sparse, and they are compressed in dense fields, so
  density partly drives this readout (in the test data, miRNA nuclei were about 10-15 % smaller than control
  nuclei at the *same* density).

## readout_nuclear_area_cv: nuclear size heterogeneity

- **Measures:** coefficient of variation (SD / mean) of the area of whole nuclei. Dimensionless. *min nuclei:* 10.
- **Biology:** a uniform, cycling population has a CV of about 0.4-0.5 (the cell cycle alone doubles DNA content).
  It rises when the population splits into small condensed nuclei (dying) and large ones (arrested or polyploid).
- **Caveats:** sensitive to a few outliers (an undetected merge counts as one huge nucleus; suspected merges are
  already excluded). Noisy with fewer than about 20 nuclei.

## readout_frac_irregular: irregular-nucleus fraction

- **Measures:** share of whole nuclei with solidity < 0.9 (area / convex-hull area; a smooth oval is about 0.97).
  The same convexity cut-off the classifier uses for merges. *min nuclei:* 10.
- **Biology:** blebbed, lobed, kidney-shaped or fragmenting nuclei: nuclear envelope damage, mitotic catastrophe,
  early apoptosis.
- **Caveats:** two touching nuclei that were not split also look irregular, more so in crowded fields. In the test
  data the miRNA effect did *not* come from density (it was there at matched density). Shape measurements get
  coarser below about 20 px nucleus diameter.

## readout_coverage: nuclear coverage

- **Measures:** total area of all in-focus nuclei (including border-cut nuclei and suspected merges) / image area.
  Dimensionless. *min nuclei:* 0 (an empty field is truly 0).
- **Biology:** a confluence proxy, roughly count × nuclear size. It does not depend on the counting rules (edge rule,
  merges), so it is a cross-check of the count. If coverage falls less than the count, the remaining nuclei are
  larger.
- **Caveats:** it measures nuclei, not cytoplasm, so it is not true confluence. In the test data it follows the
  count closely (Spearman 0.96). Analysed like the count (ratio to reference).

## readout_nn_distance_median: nearest-neighbour distance

- **Measures:** for each counted nucleus, the distance to the nearest in-focus nucleus (border-cut nuclei count as
  neighbours). The median of these distances is divided by the median nucleus diameter of the image (equivalent
  circle diameter of the median whole-nucleus area). Unit: nucleus diameters. *min nuclei:* 5.
- **Biology:** small = crowded or growing in clusters/colonies, large = sparse. At the same count, a smaller distance
  means clustering (clonal growth, or survivors in islands), and a larger one means cells spread evenly.
- **Caveats:** partly redundant with the count (Spearman −0.5). Nuclei near the border can have their true
  neighbour outside the image, which inflates their distance; the median damps this. Because it is scaled by the
  nucleus diameter, a treatment that shrinks nuclei raises it slightly at unchanged spacing.

## Dropped (tested on the real data, not kept)

- **Merged fraction:** only 8 suspected merges in about 19,000 objects; it was 0 in nearly every image. StarDist
  separates touching nuclei well, so this says nothing about clumping. Crowding is covered by coverage and
  nearest-neighbour distance.
- **Median solidity:** it varied only between 0.967 and 0.974 across conditions. The irregular *fraction* captures
  the same information with a clear signal.
- **Median eccentricity:** it mainly tracked density (nuclei elongate in dense fields), so it would be read as a
  shape change when it is really a density change.
- **All debris per nucleus:** dust and background specks dominate. Divided by the count, it rises in sparse fields
  even when the debris itself is unchanged, so it would make toxic treatments look "dirty". The micronucleus rule
  above keeps only the debris that hugs nuclei.

## Marker channels (% positive, e.g. Ki-67)

Readouts from a second (marker) channel, `readout_frac_positive_<channel>` and friends, are described in
`docs/intensity.md`. They are computed from per-nucleus intensities, not from the shapes above.

## Not yet: readouts that need intensity calibration

Intensity calibration is still pending, so no readout uses brightness, absolute or relative to the other nuclei in
the image. Candidates for later: pyknosis (small *and* bright, condensed nuclei), DNA content per nucleus
(integrated DAPI, giving a G1/G2/polyploid cell-cycle profile), chromatin condensation or texture contrast, and
confirming micronuclei by DAPI brightness similar to the main nucleus.
