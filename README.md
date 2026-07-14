# ADC Calculator

A centralized tool for antibody–drug conjugate (ADC) calculations: drug-to-antibody
ratio (DAR) determination, conjugation reaction design, and yield / formulation
planning. The scientific math lives in a **pure, tested calculation core** that is
implemented identically in Python and JavaScript; the web UI is a thin layer that
only calls core functions — no formulas are ever inlined in event handlers.

Every formula is the standard Beer–Lambert / DAR relationship from the
bioconjugation literature, restated once in `adc_spec.md` and covered by the test
suite. The methods implemented:

- UV/Vis DAR + embedded dye/payload library, with an extinction-coefficient
  calculator for new dyes and payloads.
- Conjugation routes and analytical methods: Ellman free-thiol assay, HIC,
  LC-MS, and concentration of unknown-payload ADCs.

---

## Files

| File | Role |
|------|------|
| `index.html` | The app. Single self-contained file (styles + markup + core + UI). Open it in any browser — no server, no build step, no network. |
| `adc_core.py` | Pure Python calculation core. Zero dependencies (stdlib only). |
| `adc_core.js` | Pure JavaScript ES-module port of the core. Identical signatures (camelCase) and identical numeric results. |
| `test_adc_core.py` | pytest suite — 21 tests pinning every formula to a golden value. |
| `adc_spec.md` | Formula specification and the full table of golden values. |
| `dye_library.json` | The 22-entry dye / payload extinction-coefficient library. |
| `shot_*.png` | Screenshots of each app panel. |

### Architecture: pure core / UI separation

```
          adc_spec.md   (formulas + golden values — the source of truth)
                │
      ┌─────────┴─────────┐
      ▼                   ▼
  adc_core.py         adc_core.js
      │                   │
 test_adc_core.py    selfCheck()     ← both pinned to the SAME golden values
                         │
                     index.html
                  (UI calls ADC.* only)
```

The same golden values (in `adc_spec.md`) are asserted by `test_adc_core.py`
(Python) and by `selfCheck()` (JavaScript, run automatically on page load — the
footer badge reads **"core self-check: 17/17 golden values ✓"**). Because both
language cores are pinned to one set of reference numbers, there is zero drift
between the Python analysis code and the numbers a scientist sees in the browser.

---

## The calculators (and the core functions behind them)

The UI is organized into **seven tabs**: **DAR**, Conjugation Designer, Ellman,
Yield & Formulation, In Vitro Dosing, Dye / Payload Library, and ADC Registry. The three DAR
methods (UV/Vis, HIC, LC-MS) live together under the single **DAR** tab and are
chosen with a segmented control — the same core functions back each one. Each
calculator maps to core functions; the UI passes inputs in and formats results out.

The UV/Vis calculator follows an "everyday path" layout: the three run-to-run
inputs (payload, A280, Aλmax) are shown directly, and the six constants
(antibody MW, ε values, dilution, ελmax) sit behind an **Advanced · constants**
disclosure. Selecting a payload auto-fills its extinction coefficients into
read-only (muted, dashed) fields so the derived numbers are visible but clearly
not hand-edited.

### 1. DAR · UV/Vis — `dar_uv`, `adc_properties`, `conc_unknown_payload_eps`
Two-wavelength Beer–Lambert deconvolution. Measures A280 and absorbance at the
payload λmax, then solves the coupled equations for DAR and ADC concentration.

- `R = Aλmax / A280`
- `DAR = (ελmax,mAb − R·ε280,mAb) / (R·ε280,LP − ελmax,LP)`
- `MW_ADC = MW_mAb + MW_LP · DAR`
- `ε280,ADC(molar) = ε280,mAb + DAR · ε280,LP`; `ε280,ADC(mass) = ε_molar / MW_ADC`
- `[ADC] mg/mL = A280 · DF / ε280,ADC(mass)`; `[ADC] µM = mg/mL / MW_ADC · 1e6`
- Unknown-payload fallback: `[ADC] = (A280 − Aλmax·ε280,LP/ελmax,LP) · DF / ε280,mAb(mass)`
> **Note on the circular-ε case.** If ελmax,mAb is defined as R·ε280,mAb it becomes
> circular and algebraically forces DAR → 0. The app exposes ελmax,mAb as an
> independent input (default 0 for a bare IgG in the visible range); the
> "Circular-ε case" preset restores that definition and reproduces DAR = 0 exactly,
> retained as a regression fixture (**Case A**). The "AF647 payload" preset (**Case B**)
> is the realistic case giving DAR ≈ 2.414.

### 2. Conjugation Designer — `plan_conjugation`, `moles_from_mass`, `reagent_volume_uL`
Plans a conjugation reaction from antibody mass and reagent equivalents.

- `n(mAb) = mass_mg / 1000 / MW`; `V_reagent(µL) = moles / (C_stock_mM·1e-3) · 1e6`
- `V_rxn = mass / [mAb]`; cosolvent (DMSO) topped to the target fraction.
- Equivalents can be specified *per mAb* or *per free thiol*. IgG1 full reduction
  gives 8 free thiols (4 interchain disulfides × 2).
- Route presets carry protocol defaults: `lysine`, `cys_native`, `cys_rebridge`,
  `thiomab`, `mtgase`. Source: PDF p.3–6.

### 3. Ellman (free thiols) — `ellman_free_thiols`
Free thiols per antibody from a DTNB assay.
- `free_thiols = (A412/ε412) / (A280/ε280,mAb)`, ε412 = 14150 M⁻¹cm⁻¹. Source: PDF p.10.

### 4. DAR · HIC — `dar_hic`, `dar_species_fractions`
Peak-area–weighted mean DAR from hydrophobic-interaction chromatography, plus the
% species distribution.
- `DAR = Σ(AUCᵢ · i) / Σ AUCᵢ`. Source: PDF p.14–15.

### 5. DAR · LC-MS — `dar_lcms_reduced`, `dar_lcms_intact`
DAR from deconvoluted mass-spec peak intensities.
- Reduced (denaturing + reducing, separate LC/HC): `DAR = 2·mean(LC) + 2·mean(HC)`
- Intact / native: `DAR = Σ(k·Iₖ) / Σ Iₖ`. Source: PDF p.17–18.

### 6. Yield & Formulation — `yield_and_formulation`
Recovered mass, moles, and the buffer volume needed to hit a target formulation
concentration, from starting mass and measured post-conjugation concentration.

### 6b. In Vitro Dosing — `plan_serial_dilution`, `plate_well_volume`, `plate_map`, `convert_concentration`, `series_shape_from_selection`, `assign_selection`, `aggregate_dosing_plans`
Serial-dilution planner for cell-based assays. Enter an ADC stock (µM), a plate
format (6- to 1536-well, with per-format default working volumes), an ideal top
in-well concentration, a constant fold-dilution, the number of concentration
points, and replicates. It returns the per-tube recipe (stock and assay-medium
volumes prepared to one uniform tube volume), the volume to add to each well,
the total stock and medium required, and feasibility warnings (stock too dilute
for the requested top concentration; sub-microliter pipetting). Two modes:
**spike** (add a small volume of an *f×*-concentrated dosing solution on top of
existing medium, default 10×) and **replace** (exchange the medium for dosing
solution already at the final concentration).

A **plate map** (`plate_map`, `plate_map_to_csv`) then lays the series onto the
physical wells of the chosen format. Three layouts: *dose across columns*
(default — each concentration is a column, replicates run down the rows), *dose
down rows*, and *sequential* row-major fill. The UI renders the grid with row
labels (A–H … up to A–AF for 1536-well) and column numbers, shades each well by
dose, annotates it with the concentration and replicate number, and exports the
assignment as CSV. Layouts that do not fit the grid fall back to sequential with
a warning; series exceeding plate capacity are flagged and the overflow wells
marked off-plate rather than dropped.

**Paint-first Plate Designer (UI workflow).** Rather than typing point/replicate
counts, you *paint* a rectangular block of wells directly on the plate; the
block geometry drives the series (`series_shape_from_selection`). By default the
dose **decreases top→bottom** (`by_row`: rows = concentration points, columns =
replicates); switching orientation puts the gradient across columns instead. On
release, a group editor captures the cell line, ADC label, ADC stock (µM), top
dose, and fold, and the block is committed as a **named group** in its own
colour. Multiple groups can share one plate — each with a different cell
line/ADC — and `assign_selection` maps every painted well to its point,
replicate, and final concentration (offset blocks handled). Doses are entered in
a selectable unit (**µM / nM / pM**, default nM) with live conversion against the
µM stock (`convert_concentration`); the ADC stock itself always stays µM.
`aggregate_dosing_plans` rolls the per-group plans into plate totals (wells,
stock, medium), and the CSV export emits one row per well across all groups with
`group,cell_line,adc` columns and the concentration in the selected unit.

### 7. Dye / Payload Library — `get_dye`, `DYE_LIBRARY`, `extinction_coefficient`
22 built-in reference entries (Alexa Fluor AF350–AF790 + 3 pHrodo dyes) with MW,
ελmax, the CF@280 correction factor, ε280 (= ελmax·CF280), and λmax. Selecting a
row auto-fills the UV/Vis panel.

You can add your own dyes and payloads: the **extinction-coefficient calculator**
fits a Beer–Lambert dilution series (least squares, ε = slope / path length) to
give ε280 and ελmax with R², and the **Add to library** form saves the result.
User dyes persist in the browser (`localStorage`, key `adc_user_dyes_v1`), are
merged into the table and the UV/Vis selector, marked as yours, and are deletable;
the 22 built-ins stay read-only so the core self-check is unaffected.

### 8. ADC Registry — `make_registry_record`, `registry_to_csv`, `summarize_registry`
A running log so information accumulates as more ADCs are made and tested. Every
calculator panel has a **⊞ Log to registry** button that captures its current
inputs and *core-computed* results (DAR, concentration, free thiols, conjugation
route + auto-note) against an ADC ID, date, operator, and free-text notes.

- **Persistence:** records are stored in the browser's `localStorage`, so the log
  survives closing the tab and reopening the file on the same machine. No server.
- **Portability:** *Export CSV* downloads the whole log for sharing, backup, or
  analysis in a spreadsheet. Nothing is locked in.
- **Aggregation:** a live Summary panel shows record count, distinct ADC IDs,
  counts by assay, and DAR statistics (n / mean / sample-SD / range) across every
  logged prep.
- The registry *logic* is pure and tested like every other module:
  - `make_registry_record(fields)` → canonical record (numeric coercion, blank→None,
    stripped text; requires a non-empty `id`).
  - `registry_to_csv(records)` → deterministic CSV over `REGISTRY_FIELDS`
    (RFC-style quote-escaping). **Byte-identical** between Python and JS.
  - `summarize_registry(records)` → `{n_total, n_ids, by_assay, by_id, dar_n,
    dar_mean, dar_sd, dar_min, dar_max}` (sample SD, ddof=1).
  The DAR mean/SD golden values (3.271333…, 1.673698…) are asserted by both
  `pytest` and the in-browser `selfCheck()`.

---

## How to run the tests

**Python core** (needs `pytest`):

```bash
python -m pytest test_adc_core.py -q
# 15 passed
```

**JavaScript core** (needs Node ≥ 18):

```bash
node --input-type=module -e "import {selfCheck} from './adc_core.js'; console.log(selfCheck());"
# 17/17 PASS
```

The JS `selfCheck()` also runs automatically every time the app loads; the result
is shown in the footer badge.

## How to open the app

Use the hosted version at **https://gitmcginty.github.io/Antibody-Drug-Conjugate-ADC-Calculator/**, or
open `index.html` in any modern browser (double-click, or `File → Open`). It is
fully self-contained — no server, no build, no internet connection required. The
footer self-check badge confirms the calculation core is
intact.

---

## Validation summary

| Check | Result |
|-------|--------|
| Python `pytest` | 26 / 26 pass |
| JavaScript `selfCheck()` | 16 / 16 golden values |
| Live DOM (headless browser, all panels) | UV DAR 2.414 · [ADC] 33.689 µM · HIC 2.2 · LC-MS 5.2 · Thiomab plan · ε-calculator ε280 50000 · 22-dye library + user dyes · registry log/persist/CSV — all match |
| Registry CSV, Python vs JS | byte-identical output on the same records |

Golden values (rel. tolerance 1e-6) are listed in full in `adc_spec.md §9`.
