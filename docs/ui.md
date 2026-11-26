# The browser interface (`kerncount gui`)

Everything a biologist needs is possible in the browser, without code. The page talks to a small local server
(`src/kerncount/gui.py`, standard library only, so it starts instantly) and is plain HTML/CSS/JavaScript in
`src/kerncount/web/` (no build step, no external libraries: it works offline).

Words follow `CONTEXT.md` (n, condition / treatment, well, field, form item). Result contracts:
`docs/statistics.md` (generation `manifest.json`, `stats_long.csv`), `docs/intensity.md` (`intensity/manifest.json`,
tuning API), `docs/planner.md` (`planner.json`, planner API).

## Pages

The left sidebar lists the experiments; the active experiment expands into its sections. Every view has a deep link
in the URL hash, `#<experiment>/<section>?key=value&...`, so a link (or a bookmark, or a headless screenshot)
reopens exactly that view.

| section | what it shows | hash parameters |
|---|---|---|
| **Experiment home** (`home`) | header card with the experiment metadata (cell line, magnification, seeding density, stain, channel layout, negative control, references); the **Attention** panel; the n's with date, timepoint, conditions, images, QC-low share and acceptance (dot colour = the n's colour in every plot, from `n_colours.csv`); the design summary (`design_check.csv`); the generations and their status; the **Run** actions (set up, run, rebuild generations, sampling, wells template, plate maps) with their options. The live log is docked at the bottom on every page. | `anchor` |
| **Results explorer** (`results`) | driven by `generations/gN/manifest.json`. Choose the generation, the level (all n, per cell line when several; per n; per treatment), the metric (grouped: Counts, Readouts, Intensity per channel; labels, units, kinds and estimands from the manifest) and the view (Absolute, Relative, Model, Comparisons, Robustness, Diagnostics, Design, Tables). Each plot shows its title, one-sentence description, PDF/PNG download, and enlarges on click. **Relative, comparison and robustness views always show the absolute plots of the same metric next to them.** With one n the "needs ≥ 2 n" state is stated above the plots. Only the selected plots are loaded (lazy images). | `gen, level, folder, n, t, metric, view` |
| **Statistics table** (`stats`) | `stats_long.csv` as an interactive table: filters (level, n, cell line, metric, family, estimate type, significance incl. "needs ≥ 2 n", test source: mixed model / classical / descriptive, comparison text), sorting by any column, readable numbers (ratios `0.29×`, percentage points `+1.2 pp`, p `<0.0001`), the estimand on every row, significant rows highlighted, "needs ≥ 2 n" badges, export of the filtered view to CSV (in the browser), and a click on a row opens the plot it belongs to. | `gen, level, n, cl, metric, fam, et, sig, src, q, sort, dir, limit` |
| **Layer comparison** (`layers`, fluorescence only) | the readout ladder per channel: the five layers and what each step changes; the ladder plot per reference (all n and per n) with the absolute version below it; the concordance table with a coloured call per layer and its plain-language summary ("the call changes at step 2→3 …"; red rows flip direction); the density scatters (generation and `intensity/`) and a link to the counts. | `gen, ch, lvl, ref` |
| **Intensity & cut-off** (`intensity`, fluorescence only) | tabs: *Cut-offs & flat field* (cut-off rule and value per n and measure, % positive with / without flat field, locks and their history; flat-field maps with their magnitude; cut-off sweep), *Distributions & density*, *Channel checks* (alignment plot, flagged and mis-paired fields with their files), and the **blinded tuning view** (below). | `tab`; tuning: `ch, measure, scope, k, seed, field, cutoff` |
| **Review images** (`review`) | review sheets per n and treatment; every field shows its source file(s) (`source_file`, `source_file_<channel>`), duplicates are marked; the lightbox steps through fields with arrow keys and takes corrections. Fluorescence experiments: **Show Ki67 positivity** draws the marker image with positive / negative / saturated outlines at the cut-off in force, with the condition shown (this is the non-blinded view). | `n, t, field, overlay, q` |
| **Corrections** (`corrections`) | the manual-count editor that builds generations (g1 + round 1 → g2 …), unchanged; field numbers carry their source files as a tooltip. | |
| **Planner** (`planner`) | from `sampling/gN/planner.json` (several cell lines: one plan each, chosen with the **cell-line** pills; n's left out by enforced acceptance are named) via `planner.load / evaluate / curve`: sliders and selects for n, wells, fields, effect (a ratio or the observed effect of a treatment), reference, reference wells, alpha, target power and field positions, and the **"Power for a star"** checkbox (on by default = `[planner] adjust_for_family`: power and smallest detectable ratio at the family-adjusted level of the stars; off = one unadjusted test; the test, critical value and df are printed under the notes); live power, CI half-width, smallest detectable ratio and effort, each with its pessimistic–optimistic band; a live SVG chart of the curve (vary n / wells / fields; power or CI width, with the band, target line and "your design" marker); the plain-language notes; "Recommendations for the next experiment" from `REPORT.txt`; the static planner plots and tables; and the sampling run. Without ≥ 2 n with the reference the page says so and shows the report only. | `gen, cl, n, wells, fields, effect, ref, rw, alpha, power, pos, adj, vary, show` |
| **Plate maps** (`plates`) | a form (plates, format 6/12/24/48/96, wells per condition, reference wells or auto, avoid edge, balance, seed, negative-control name and wells) → `platemap.plan_project` in-process → the plate PNGs, edge-use table, README (seed and the command that recreates the maps), PDF/CSV downloads; earlier plans in `plans/`. | `plan` |
| **Layout**, **Settings**, **All files** | layout.csv editor (n, date, timepoint, treatment, include); every setting as a form (all form-item kinds: bool, int, float, choice, multi, list, num_or, measure, names, table, text, with help text and defaults; settings that re-detect nuclei are marked); after saving, the page lists the steps that need re-running with a Start button each; the project folder browser. | `sec` (jump to a settings section) |

Visual language: absolute = blue, relative = orange, model = purple, comparisons = teal, robustness = green,
diagnostics = grey, design = brown, layers = magenta (`--k-*` in `app.css`); n's always use the plot colours.
Confirmations are in-page (`confirmBox`); the browser's own dialogs are never used. Anything that overwrites a
file you may have edited asks first and says what is reset: **Set up → "start over"** (settings.toml, i.e. every
setting *including the experiment metadata* — cell line, magnification, pixel size, seeding density, stain,
negative control — and layout.csv), "rewrite the layout only", **Settings → reset**, **use the fresh scan**, and
**overwrite wells.csv**. The old file is copied to `_work/backups/<name>_<date-time>` first (also from the
command line: `init --force`, `wells-template --force`); the message names the copy.

### The Attention panel

`gui.attention_items` gathers every issue, most serious first (error → warn → info), each with a plain-language
explanation and a link to the evidence (a table, or the page that shows it): duplicate images across conditions
(`duplicates.csv`), acceptance failures (`acceptance.csv` of the latest generation), mis-paired and shifted channel
files (`intensity/alignment.csv`; whole folders that fail the channel check are named first), folders with files
that were not used (an unpaired channel file, non-numbered files, an empty folder: `layout.csv` notes), conditions
missing from an n (`design_check.csv` rows with 0 images: a missing reference or treatment), other design flags
(`design_check.csv`; docs/design.md, "Unusual designs"), a high QC-low share per n, model fallbacks
(`manifest.model_status`), the one-n notice, an automatic cut-off that calls almost every nucleus positive, a
missing negative control, a locked cut-off not yet applied, settings changed after the latest generation, and a
fresh folder scan (`layout.new.csv`).

### What needs re-running after a settings change

`gui.rerun_advice(before, after)`: a changed form item marked `reseg` → a run that re-detects every nucleus (slow);
another change in project / experiment / channels / segmentation / counting / qc / intensity / positivity /
readouts / design / review → a run (detections reused) with `--rebuild all`; analysis / acceptance → `regen
--rebuild all` (statistics only); planner → the sampling analysis. The page shows these steps with Start buttons.

## Tumour experiments (`[experiment] type = "tumour_volume"`)

The type comes from `settings.toml` (`gui.experiment_type`; `/api/state` lists it per experiment with `results` =
`stats_long.csv` exists, and the registered `types` for the set-up form). `web/tumour.js` draws the tumour views and
wraps `RENDER.home / layout / review / corrections / results / stats` (a tumour project gets its own view, every other
project the unchanged cell-count one) and `RENDER.settings` (the analysis-plan card is added on top). Planner and
plate maps are hidden for tumour projects. Server side: `src/kerncount/gui_tumour.py` (standard library at import;
the tumour modules that need pandas are imported on first use).

| section | tumour view |
|---|---|
| **Experiment home** | factors and levels (with where they are read from), attention items (layout problems that block the run, layout notes, an unreadable or unlocked plan, a plan changed after the results, settings changed after the run, QC-flagged and saturation-warning tumours, a fresh scan), included tumours per group, results summary (tumours analysed / flagged / excluded, review generation, `results.xlsx`, report), the plan lock, and the **Run** steps: re-scan (start over / layout only), run (*skip review images* = `--no-review`, *segment every stack again* = `--resegment`, *threshold sweep* = `[tumour] threshold_sweep`, saved before the run starts), review and corrections. A project that is not set up yet shows the generic set-up step with an **experiment type** select (`init --type <name>`). |
| **Layout** | `layout.csv` with the validation of `tumour.layout.validate`: problems (block the run) and warnings above the table; per row file, session folder, one select per factor (its levels; a value that is not a level is kept and marked "?"), section, hemisphere, the rebuilt tumour ID, include, note; other files in the raw folder. Saving writes the edited cells by `file` and answers with the new validation. |
| **Settings** | the form sections `project`, `experiment` (type, notes) and `tumour` (`settings_form.schema("tumour_volume")`; items marked `reseg` = the next run measures every stack again, shown as *re-measures*), plus the **analysis plan** card: factors (levels, reference level, display role), readouts (confirmatory / exploratory), families (confirmatory / exploratory / descriptive), the lock and its history. Adding or removing factors, readouts or families: the settings.toml text editor. Re-run advice: one step, a run (`gui_tumour.rerun_advice`; slow when a `reseg` item changed). Reset to defaults uses the type's own template. |
| **Review images** | **blinded**: the review sheets by code (`review/<code>.jpg`), QC flags, object count, overview sheets; a lightbox per code takes a correction (object number or `x: reason`) and a note for the chosen round. Filter: QC-flagged only (`q=flagged`). |
| **Corrections** | **blinded**: `review/corrections.csv` as a table by code: earlier rounds, the round being edited (revise an earlier round or start the next), note, and what the next run applies (`auto` / `object k` / `excluded: reason`); Save, or Save & run. A new `manual_k` is added when the last one is used; values are checked like the run does (`tumour.corrections` rules; an object number must exist on the image). |
| **Results explorer** | one pill per readout (confirmatory first): its plots from `plots/<readout>_<kind>.png` (absolute, normalised, effects, interaction, A vs B, vs Fiji, threshold sweep, diagnostics; PDF/PNG), its headline statistics, the shape panel, tables (`tumours.csv`, `excluded_tumours.csv`, `flow.csv`, `concordance.csv`, `threshold_sweep.csv`, `corrections_used.csv`, `lock_log.csv`), Prism CSVs, **`results.xlsx` download** and the report. |
| **Statistics table** | `stats_long.csv`: readout, family, significance and text filters. **Headline view without the sensitivity rows**: families ending in ` (Welch-type)`, ` (permutation)` or ` (without flagged)` (`gui_tumour.is_sensitivity`) are hidden (the count is shown) unless *show sensitivity rows* is ticked (`sens=1`). |

Hash parameters: results `ro` (readout); stats `metric, fam, sig, sens, q`; review / corrections `q`.

### Blinding of the tumour review

`/api/tumour/review` and `/api/tumour/corrections` read and write only the `review/` folder: never
`review_key.csv` (code → tumour ID), `corrections_used.csv`, `layout.csv` or any result. Their answers carry the code,
object count, QC flags, the `manual_*` values, the note and the image path — no factor level, session, date, section
label or file name (codes are checked against `[A-Za-z0-9_-]` before an image path is built). The Review and
Corrections views call no other endpoint (`tests/test_tumour_gui.py` checks both the answers and the script). The key
and `corrections_used.csv` stay reachable through **All files** and the results tables: unblinding is a deliberate
step, not something these two pages do.

### Analysis-plan lock (ADR 0003)

`tumour/lock_log.py`. The plan = factors (levels, reference), readouts (column, role), families and the focus factor
of normalised values; `tumour.results.lock` hashes it and every `stats_long.csv` row carries `lock <hash>`. Results
exist = `stats_long.csv` in the project folder; its hash is the **results lock**. Every settings write of a tumour
project (form, text editor, plan editor, reset) goes through `gui_tumour.guard`: when results exist and the plan
hash would change, it is refused (`PlanLocked`) unless the plan is **unlocked**. **Unlock** needs a reason; it is
appended to `lock_log.csv` (`at, action, results_lock, settings_lock, reason`), as is every saved plan change
(`plan changed`, old → new) and **Lock again** (`lock`). The run (`results.analyse` → `lock_log.record_run`) appends a
`run` entry when it rebuilds the results with a different plan — noting the unlock, or "without a logged unlock" when
settings.toml was edited outside the web interface — and ends a pending unlock; then it prints the whole log in the
About lines of `stats_report.md` and `results.xlsx`. Display roles, QC, sweep, sensitivity and planner settings are
not part of the plan and can change at any time.

## Endpoints

All answers come from `127.0.0.1` only. Every request must be addressed to `127.0.0.1:<port>` or
`localhost:<port>` (Host header; an Origin header, if present, must match too), which blocks other websites and DNS
rebinding. Every API call and file needs the per-session token (`X-Token` header, or `?t=` for GETs such as images);
`/` puts it in a `<meta name="kerncount-token">` tag. File paths are resolved inside the project folder
(`gui.inside`) and generation names must match `g<number>`.

**Content-Security-Policy.** Every response carries `Content-Security-Policy: default-src 'self'; img-src 'self'
data: blob:; style-src 'self' 'unsafe-inline'; script-src 'self'; object-src 'none'; base-uri 'none'; form-action
'self'; frame-ancestors 'none'`, `X-Frame-Options: DENY`, `X-Content-Type-Options: nosniff`, `Referrer-Policy:
no-referrer`. So the page runs no inline script at all: there is no `<script>` block and no `on*=` attribute.
Templates attach handlers with `on("click", fn)` (app.js), which stores the function and writes
`data-on-click="<id>"`; one delegated listener per event type (click, change, input, dblclick; `<img>` errors in
the capture phase) calls it with `this` = the element, like an inline handler (returning false prevents the
default). The static buttons of `index.html` are wired in app.js. A test fails if any `on…=` attribute or inline
script comes back.

**Project files are never active content** (`/file`): only passive types are shown in the page (images, CSV/TXT/
TOML/MD/LOG as `text/plain`, JSON); HTML, HTM, XHTML, SVG, XML, XSL, JS/MJS, PDF and anything unknown are sent as
a download (`Content-Disposition: attachment`, `application/octet-stream` / PDF), and every file response carries
its own `Content-Security-Policy: default-src 'none'; …; sandbox`. File names in `Content-Disposition` follow
RFC 5987: an ASCII fallback (`filename="dose_5.csv"`; quotes, backslashes, `;`, control and non-ASCII characters
replaced) plus the exact name `filename*=UTF-8''dose%E2%89%A55.csv`; all headers are built (and checked to be
encodable) before the status line is sent, so an odd name gives a clean error, never a half-sent response.

**Internal files are never served**: `_work/` (detection caches, the tuning audit trail, backups) is refused by
`/file`, `/api/text`, `/api/table`, `/api/open` and left out of `/api/files` (also `_WORK/…` spellings, which
open the same folder on macOS / Windows). `cutoffs.json` is listed, but `/file` and `/api/text` return a
sanitized summary (locks and history with whitelisted keys, `fields` as a count), never the raw file.

**Writes** of `cutoffs.json`, `settings.toml`, `layout.csv` and `counts.csv` from the page take a per-file lock
(read-modify-write from two tabs cannot interleave) and are atomic (temporary file in the same folder, then
`os.replace`; permissions kept): a reader never sees half a file (`kerncount.fileio`).

**Job command lines** put the options first and the experiment after `--` (`run --rebuild g2 -- <experiment>`),
so a name can never be read as an option; experiment names and `negative_control` may not start with `-`
(refused by the page and the settings form).

| method | path | answer |
|---|---|---|
| GET | `/`, `/static/<file>` | the page (token embedded), the scripts and styles (only files directly in `web/`) |
| GET | `/file?exp&path[&dl]` | a file inside the project (passive types inline, everything else as an attachment; `_work/` refused; `cutoffs.json` sanitized) |
| GET | `/api/state`, `/api/log?since` | experiments (with `intensity` flag), running job; job log lines |
| GET | `/api/project?exp` | settings, form schema (`items`), generations, n colours, sampling / planner generations |
| GET | `/api/overview?exp` | experiment home: metadata, n's, design summary, generations, attention items, cut-offs |
| GET | `/api/manifest?exp&gen` | `{manifest, colours}` (cached by file mtime) |
| GET | `/api/stats?exp&gen` | `stats_long.csv` as `{columns, rows}` ("nan" → "") |
| GET | `/api/review`, `/api/counts`, `/api/table`, `/api/text`, `/api/files` | review index (with `sources`, `marker`, `duplicate` per field), counts.csv, any CSV / text, folder listing |
| GET | `/api/intensity?exp` | `intensity/manifest.json`, locks (without file names) and their history |
| GET | `/api/tuning/info?exp` | marker channels, lock scopes (n's), locks, history |
| GET | `/api/tuning/range?exp&channel&measure` | slider range: `min, max, lod_median` |
| POST | `/api/tuning/sample` `{experiment, channel, k, seed, n}` | `{seed, channel, fields: [{field_id}]}` |
| GET | `/api/tuning/render?exp&field&channel&cutoff&measure` | PNG (token-checked) |
| GET | `/api/tuning/preview?exp&cutoff&channel&measure&n&fields` | overall % positive in the scope and in the sampled fields (only for ≥ 3 fields); the automatic rule's % |
| POST | `/api/tuning/lock`, `/api/tuning/unlock` | lock (via `tuning.lock`) / remove a lock; both append to `history` in `cutoffs.json` |
| GET | `/api/field_overlay?exp&n&treatment&image&channel` | PNG of a known field with outlines at the cut-off in force (non-blinded; review page only) |
| GET | `/api/planner?exp&gen&cl` | design, options (incl. `adjust_for_family`), `family_size`, `cell_lines`, `cell_line`, `base` (folder of the plan's files), `excluded_ns`, effects, components, notes, recommendations, sampling manifest; `ok` / `error` |
| POST | `/api/planner/evaluate` `{experiment, gen, cell_line, n, wells, fields, effect, reference_wells, alpha, target_power, reference, positions, adjust_for_family, vary}` | `{result: planner.evaluate(...), curve: planner.curve(...), vary}` |
| GET | `/api/plans?exp` | earlier plate maps (plates, files, README, edge use) and form defaults |
| POST | `/api/platemap` | `platemap.plan_project(...)` in-process (validated options); returns the new folder |
| POST | `/api/start` `{action: init / run / regen / sampling / wells-template / platemap, options}`, `/api/stop`, `/api/quit` | jobs (`python -m kerncount.cli ...` subprocess, one at a time) |
| GET | `/api/tumour/overview?exp`, `/api/tumour/layout?exp`, `/api/tumour/plan?exp` | tumour home; layout + validation (`problems`, `warnings`, `ids`, `groups`, `unmatched`); analysis plan + `lock` status and log |
| GET | `/api/tumour/review?exp` | **blinded** review: `rows` (code, objects, qc_flags, manual, note, image), `manual` columns, `overviews`, `mtime` |
| GET | `/api/tumour/results?exp`, `/api/tumour/stats?exp&sensitivity` | plots per readout, tables, Prism files, `results.xlsx`, report, lock; `stats_long.csv` headline rows (`hidden` = sensitivity rows left out) |
| POST | `/api/tumour/layout` `{rows}`, `/api/tumour/plan` `{factors, readouts, families}`, `/api/tumour/unlock` `{reason}`, `/api/tumour/relock`, `/api/tumour/corrections` `{changes: [{code, col, value}], mtime}` | layout edits (by `file`); plan edits (lock checked); unlock / lock again (logged); blinded corrections |
| POST | `/api/settings`, `/api/settings_reset`, `/api/layout`, `/api/layout_use_new`, `/api/counts`, `/api/open`, `/api/workspace`, `/api/browse`, `/api/new_experiment` | small edits, with the command line's validation; `/api/settings` also returns `{changed, steps}` (what to re-run) |

The scientific modules (`planner`, `tuning`, `platemap`) are imported on first use and warmed up in a background
thread when the server starts; manifests, `stats_long.csv` and plans are cached and re-read when the file changes.
On the review workspace a page loads in well under a second; the tuning slider re-renders a field in about 30-90 ms.

## How blinding is guaranteed (cut-off tuning)

1. **Opaque, per-session ids only.** `/api/tuning/sample` returns `{seed, channel, fields: [{field_id}]}`: salted
   hash ids of `kerncount.tuning` in random order (round robin over conditions inside the module). The salt lives
   only in the server's memory and is new for every session, so an id cannot be matched to a field across sessions
   and nothing on disk maps ids to fields (an older `_work/tuning_salt` file is deleted). The server rebuilds the
   answer field by field, so the requested `n` (the lock scope) is not echoed. **No nucleus count per field**: it
   tracks the treatment (treatments that kill cells leave sparse fields), so it is neither returned nor shown.
2. **No per-condition result.** `/api/tuning/preview` returns only overall numbers: the % positive of the whole scope
   (every n, or the one n chosen as the lock scope), of the sampled fields together (only for 3 or more fields, so it
   cannot be used to count the nuclei of one field), and of the automatic rule.
3. **Images carry no label.** `/api/tuning/render` returns PNG bytes of the marker image with outlines; the field
   is addressed only by its opaque id.
4. **Locks without file names.** `cutoffs.json` stores the value, time, scope and the *number* of fields. The
   fields' files (the audit trail) go to `_work/cutoff_fields.json`, which the server never serves; file names an
   older version left in `cutoffs.json` are moved there on the next tuning visit. Every answer the page receives
   (`/api/tuning/lock`, `/api/tuning/info`, `/api/intensity`, and `/file` or `/api/text` of `cutoffs.json`) is a
   whitelist of keys with `fields` as a count.
5. **The page shows nothing else.** The tuning tab shows the field index, the opaque id, the slider,
   and the overall %; it never calls the review, overlay or results endpoints and displays nothing else. The only n
   names on the tab are the lock-scope choices the user makes ("every n" or one n) and the scope of earlier locks,
   which say nothing about the field on screen.
6. **Tested.** `tests/test_gui.py::test_tuning_endpoints_are_blinded` runs the pipeline on a small synthetic
   interleaved project and checks that no tuning answer contains a condition name, a folder name, a file name or
   an n label; `test_sample_lock_read_files_never_leaks_conditions` follows the whole path sample → lock → read
   every file the page can reach (`cutoffs.json` incl. a legacy one with field ids and files, `_work/…` through
   `/file`, `/api/text`, `/api/table`, the recursive file list) and checks that ids change with the session.

The non-blinded per-field overlay of the review page uses a separate endpoint (`/api/field_overlay`, addressed by
n / condition / field, via `tuning.field_id`), so the two views cannot be mixed up. The lock confirmation explains
that the value must be fixed before looking at the statistics and that a run applies it, and then offers to start
the run.

## Tests

`tests/test_tumour_gui.py` (tumour projects, synthetic, nothing segmented): the plan lock (no results → free; locked
→ refused through the plan editor, the form, the text editor; unlock needs a reason; change and run are logged; a
rebuild with a changed plan shows old → new lock in the report and About sheet), headline statistics without
sensitivity rows, blinded review / corrections answers and script, corrections validation and rounds, layout
validation and edits, settings form sections and re-run advice, `init --type`, results / overview endpoints.


`tests/test_gui.py` (default suite): token, Host and Origin checks on every kind of endpoint; the CSP / frame
headers and that no inline script or handler is left; active project files served as attachments and RFC 5987
names (`dose≥5.csv`, a name with `"`); backups when the page resets or replaces a file; path confinement
(`/file`, `/api/table`, `/static`, generation names, unknown experiments); the scripts parse (`node --check` when
node is installed) and never use the browser's dialogs; manifest / stats_long serving; the experiment home and its
attention items; the planner endpoints (bands ordered, curve, input validation); the plate-map endpoint and job
arguments; a settings round-trip for every form-item kind; the re-run advice; the blinded tuning endpoints and the
review overlay on a real (synthetic) run.

Visual check: headless Chrome (`--headless=new --window-size=1440,1000 --virtual-time-budget=8000 --screenshot`)
on deep links such as `#Ki67/intensity?tab=tune&cutoff=300&seed=7&field=1` or
`#DAPI_output/planner?n=6&effect=0.7`.
