# ADC Calculator — Calculation Core Specification

This document defines every formula implemented by the calculation core
(`adc_core.py` / `adc_core.js`). It is the single source of truth: both
implementations and the test suite are pinned to the **golden values** in §9,
which keeps the two cores and the embedded HTML copy from drifting apart. The
formulas are the standard Beer–Lambert / DAR relationships from the
bioconjugation literature.

All functions are **pure** (no UI, no I/O). Units are explicit in every
signature.

---

## 1. Symbols and units

| Symbol | Meaning | Unit |
|---|---|---|
| `MW_mAb` | antibody molar mass | g/mol (Da) |
| `MW_LP` | linker–payload molar mass | g/mol |
| `ε280_mAb` | antibody molar ε at 280 nm | M⁻¹·cm⁻¹ |
| `ε280_mAb_mgml` | antibody mass ε at 280 nm | (mg/mL)⁻¹·cm⁻¹ |
| `ελmax_mAb` | antibody molar ε at payload λmax | M⁻¹·cm⁻¹ (≈0 for bare IgG in visible) |
| `ε280_LP` | linker–payload molar ε at 280 nm | M⁻¹·cm⁻¹ |
| `ελmax_LP` | linker–payload molar ε at its λmax | M⁻¹·cm⁻¹ |
| `A280` | measured absorbance at 280 nm | AU |
| `Aλmax` | measured absorbance at payload λmax | AU |
| `DF` | nanodrop dilution factor | – |
| `DAR` | drug-to-antibody ratio | mol/mol |

**Extinction-coefficient conversion** (used throughout):
`ε_molar (M⁻¹cm⁻¹) = ε_mass (mg/mL⁻¹cm⁻¹) × MW (g/mol)` and the inverse.

---

## 2. DAR by UV/Vis (spreadsheet cell I16)

Two-equation Beer–Lambert deconvolution at 280 nm and payload λmax. Let
`R = Aλmax / A280` (the measured absorbance ratio).

```
DAR_UV = (ελmax_mAb − R · ε280_mAb) / (R · ε280_LP − ελmax_LP)
```

- This is the exact rearrangement of the spreadsheet formula
  `=(E9−(E8·E18/E17))/((E18/E17)·I7−I9)`.
- **Template caveat:** the spreadsheet defines `ελmax_mAb` (E9) circularly as
  `R · ε280_mAb`, which forces `DAR_UV = 0`. In the app, `ελmax_mAb` is an
  **independent input** and defaults to `0` (bare antibodies have negligible
  absorbance at typical payload λmax in the visible range). The DAR=0 case is
  retained as a regression fixture (golden case A).

---

## 3. ADC molecular properties (spreadsheet I17–I21; PDF p.12/15/17/18)

Given a DAR (from any method) and the payload/antibody properties:

```
MW_ADC          = MW_mAb + MW_LP · DAR                         # g/mol
ε280_ADC_molar  = ε280_mAb + DAR · ε280_LP                     # M⁻¹cm⁻¹
ε280_ADC_mass   = ε280_ADC_molar / MW_ADC                      # (mg/mL)⁻¹cm⁻¹
conc_mgml       = A280 · DF / ε280_ADC_mass                    # mg/mL
conc_uM         = conc_mgml / MW_ADC · 1e6                     # µM
```

Note the spreadsheet writes `[ADC] mg/mL = (A280·MW_ADC)/ε280_ADC_molar · DF`,
which is algebraically identical to `A280·DF / ε280_ADC_mass`.

**Alternative concentration when payload ε is unknown** (PDF p.12): DAR is
not obtainable, but the mass concentration is:
```
conc_mgml = ((A280 − Aλmax_LP_contribution) / A280_scale) ...
# PDF form: C = (A280 − A_LP280) where A_LP280 = (Aλmax/1)·(ε_LP,280/ε_LP,λmax)
conc_mgml = (A280 − Aλmax · ε280_LP/ελmax_LP) / ε280_mAb_mgml
```

---

### 3b. Scattering / turbidity baseline correction (`scatter_corrected_absorbance`)

Turbid or aggregated samples scatter light, adding a chromophore-independent
baseline that biases A280 (and therefore concentration and DAR) high. Read a
reference absorbance `A_ref` at a wavelength where neither antibody nor payload
absorbs (typically 320–340 nm) and extrapolate a Rayleigh-type baseline:

```
A_scatter(λ)   = A_ref · (λ_ref / λ)^p               # p = 4 (Rayleigh) default, exposed
A280_corrected = max(A280 − A_scatter(280), 0)
Aλmax_corrected= max(Aλmax − A_scatter(λmax), 0)
```

Feed the corrected absorbances into `dar_uv` / `adc_properties`. A zero
reference read is the identity (no correction); a reference large enough to
drive a corrected read negative is floored at 0 (physically meaningless).

---

## 4. Ellman's free-thiol assay (PDF p.10)

```
free_thiols = (A412 / ε412) / (A280 / ε280_mAb)      with ε412 = 14150 M⁻¹cm⁻¹
```
`A280` and `A412` are the absorbances of the DTNB-reacted protein aliquot;
`ε280_mAb` is the molar extinction coefficient of the protein.

### 4b. Reagent-blank / background subtraction

DTNB and its TNB²⁻ product, together with any residual reductant (TCEP/DTT),
absorb at **both** 412 and 280 nm. A matched reagent blank (all reagents, no
protein) measures that background; subtracting it before forming the ratio
removes the bias:
```
free_thiols = (A412 − A412_blank) / ε412  ÷  (A280 − A280_blank) / ε280_mAb
```
`A412_blank` and `A280_blank` default to 0, so the call reduces exactly to the
raw ratio when no blank is supplied. The net (blank-subtracted) A280 must be
positive.

---

## 5. DAR by analytical HIC (PDF p.14–15)

AUC-weighted average over resolved DAR species (peak area at 280 nm ∝ AUC):
```
DAR_HIC = Σ_i (AUC_i · i) / Σ_i AUC_i           i = 0,1,2,…,n
```

### 5b. Drug-load-corrected HIC DAR (`dar_hic_corrected`)

The uncorrected form above treats 280 nm peak area% as **molar** abundance%.
That is only true when every species has the same 280 nm extinction. The
payload absorbs at 280 nm, so species with k drugs have
`ε_k = ε280_mab + k·ε280_payload`, and measured area is `A_k = c_k·ε_k`. To
recover molar abundance, divide each area by its species extinction before
averaging:
```
c_k  ∝  A_k / (ε280_mab + k·ε280_payload)
DAR_HIC,corr = Σ_k (k · c_k) / Σ_k c_k
species%_k   = 100 · c_k / Σ_j c_j
```
When `ε280_payload = 0` this reduces exactly to the uncorrected `DAR_HIC`.
Because high-DAR peaks carry more 280 nm signal per mole, the correction
**lowers** the reported DAR relative to the naive area-weighted mean.

---

## 6. DAR by LC-MS (PDF p.17–18)

Let `I(DAR=k)` be the deconvoluted peak intensity for species with k drugs.

- **Denaturing + reducing** (light + heavy chains, LC max 2, HC max 4 drugs):
  ```
  DAR = 2·[Σ_k k·I_LC(k) / Σ_k I_LC(k)]  +  2·[Σ_l l·I_HC(l) / Σ_l I_HC(l)]
  ```
- **Denaturing + non-reducing** and **native** (intact, up to ~9 drugs):
  ```
  DAR = Σ_k k·I(k) / Σ_k I(k)
  ```

### 6b. Ionization-response correction (`dar_lcms_*_corrected`)

The forms above assume MS peak intensity is proportional to molar abundance.
In electrospray this fails: adding a hydrophobic payload changes a species'
ionization efficiency, so highly loaded species are typically **under-counted**
and the naive DAR is biased low. Given a relative response factor `r_k`
(signal per mole, normalised to any convenient reference, e.g. `r_0 = 1`) for
each drug-load species, molar abundance is recovered before averaging:
```
n_k  ∝  I(k) / r_k
DAR  = Σ_k k·n_k / Σ_k n_k          (intact)
DAR  = 2·[Σ k·n_LC / Σ n_LC] + 2·[Σ l·n_HC / Σ n_HC]   (reduced)
```
With no response factors supplied (or a missing key → `r_k = 1`) these reduce
exactly to the uncorrected LC-MS DAR. Response factors must be positive.

---

## 7. Conjugation designer (PDF p.3–6)

Given antibody amount, computes reagent volumes from equivalents.
`n_mAb (mol) = mass_mAb (g) / MW_mAb`. For a reagent at stock concentration
`C_stock` and `E` equivalents (relative to mAb, or per free thiol):

```
n_reagent = E · n_mAb            (or E · n_thiol, n_thiol = n_mAb · thiols_per_mAb)
V_reagent = n_reagent / C_stock
```

Route defaults (from protocols):
- **Lysine:** activated ester 5–20 equiv (screen); mAb 5 g/L; +10% v/v 1.0 M
  NaHCO₃; final 10% v/v DMSO; 2 h/25 °C; quench L-Lys 100 equiv.
- **Cys native / rebridge:** TCEP 6.0 equiv for full IgG1 reduction
  (1 equiv per disulfide; partial 1.0–3.0 equiv); maleimide 1.5 equiv per free
  thiol; rebridging equiv by screen; +5 mM EDTA; 10% v/v DMSO; quench L-Cys 100 equiv.
- **Thiomab:** TCEP 10 equiv full reduction → desalt → dhAA 30 equiv reoxidation
  (4 h/25 °C) → maleimide 1.5 equiv/thiol.
- **MTGase:** PNGase F 10 µL/mg; tag (e.g. azido-PEG3-amine) 80 equiv; MTGase
  5.5 U/mg; click reagent (DBCO) 3.0 equiv; site-specific DAR 2.0.

Cosolvent guardrails: ≤15% v/v organic (IPA/DMSO/DMF/DMA) generally; up to 30%
v/v propylene glycol for hydrophobic payloads. Final DMSO target 10% v/v.

**Free-thiol count** per full reduction: IgG1 = 8 (4 interchain disulfides →
implemented as `disulfides_reduced × 2`).

---

## 7b. Predicted DAR distribution (site-occupancy model, chemistry-keyed)

Drug load per antibody follows a site-occupancy model in which each of `n`
equivalent, independent **sites** reacts with probability `p` and, when it
reacts, adds `d = drugs_per_site` drugs. The occupied-site count is
**Binomial(n, p)**, so with `d` drugs per reacted site:
```
P(DAR = d·j) = C(n,j) · p^j · (1−p)^(n−j)     j = 0,1,…,n
mean DAR     = d · n · p
variance     = d² · n · p · (1−p)      SD = sqrt(variance)
```

The `drugs_per_site` parameter encodes the conjugation chemistry — this is the
whole point of the model:

- **Cysteine / TCEP interchain conjugation (`d = 2`, default).** An IgG1 has
  four interchain disulfides (`n = 4`); reducing each exposes a **pair** of
  thiols conjugated together, so a reacted site contributes 2 drugs. The DAR
  ladder is therefore **0, 2, 4, 6, 8** — the even-only support seen for real
  thiol-linked ADCs (e.g. brentuximab vedotin). Odd DAR species cannot occur.
- **Stochastic amine / lysine conjugation (`d = 1`).** Each surface lysine
  reacts independently carrying a single drug, giving the classic smooth
  binomial over **0, 1, 2, …, n**.

The per-site probability may be given directly, or derived from the molar feed
ratio and conjugation efficiency, spread over the total drug-carrying capacity
`n·d`:
```
p = feed_ratio · efficiency / (n · d)     ⟹   mean DAR = feed_ratio · efficiency
```
Implemented as `predict_dar_distribution(n_sites, p_site=…|feed_ratio=…, efficiency=…, drugs_per_site=2)`;
moments are taken through `distribution_dispersion` so definitions never drift.
The distribution is keyed by **DAR** (multiples of `d`). The SD is the intrinsic
drug-load **heterogeneity** predicted by the model, not a measurement error.
Errors if `n_sites<1`, `drugs_per_site<1`, neither `p_site` nor `feed_ratio`
given, or the implied `p` falls outside `[0,1]`.

---

## 8. Yield & formulation

```
recovered_mass_mg = conc_mgml · volume_mL
yield_pct         = recovered_mass_mg / starting_mass_mg · 100
molar_amount_nmol = recovered_mass_mg / MW_ADC · 1e6
# formulate to a target concentration:
V_final_mL        = recovered_mass_mg / target_conc_mgml
V_diluent_to_add  = V_final_mL − current_volume_mL   (if concentrating: negative → concentrate)
```

---

## 8b. Extinction-coefficient determination (Beer–Lambert regression)

Given a dilution series of known molar concentrations `c_i` (mol/L) and
measured absorbances `A_i` at one wavelength and path length `L` (cm),
Beer–Lambert gives `A = ε·L·c`. An ordinary least-squares fit of `A` vs `c`
yields:

```
slope     = Σ(c_i − c̄)(A_i − Ā) / Σ(c_i − c̄)²
intercept = Ā − slope · c̄
R²        = 1 − Σ(A_i − (slope·c_i + intercept))² / Σ(A_i − Ā)²
ε         = slope / L                       # M⁻¹ cm⁻¹
```

`linear_regression(xs, ys)` → `{slope, intercept, r_squared, n}` (needs ≥2
points and ≥2 distinct x). `extinction_coefficient(concentrations_M,
absorbances, path_length_cm=1.0)` → `{eps, slope, intercept, r_squared, n}`.
Run the fit separately on the A280 column and the Aλmax column to obtain
ε280 and ελmax for a new dye/payload.

**Golden fixture:** `c = [1,2,3,4,5] µM`, `A280 = 0.01 + 50000·c`,
`Aλmax = 120000·c`.
- `linear_regression(c, A280)` → slope **50000.0**, intercept **0.01**, R² **1.0**
- `extinction_coefficient(c, A280, 1.0).eps` → **50000.0**
- `extinction_coefficient(c, A280, 0.5).eps` → **100000.0**
- `extinction_coefficient(c, Aλmax, 1.0).eps` → **120000.0**

---

## 8c. In vitro dosing / serial dilution

Design a plate-dosing series from an ADC stock. Inputs: stock concentration
`C_stock` (µM), the ideal top in-well concentration `C_top` (µM), a constant
`fold` dilution, number of concentration points `n`, `replicates` per point,
and the per-well working volume `V_well` (µL) — set from the plate format
(defaults below) or overridden.

**Final (in-well) concentrations** for `i = 0 … n−1`:

```
C_final,i = C_top / fold^i
```

**Dosing modes.** In *spike* mode a small volume of a `f×`-concentrated
dosing solution is added on top of the existing medium; in *replace* mode the
medium is exchanged for dosing solution already at the final concentration
(`f = 1`).

```
working_conc,i = f · C_final,i           # concentration prepared in each tube
spike:   V_add = V_well / (f − 1)        # so post-mix well conc = C_final
         V_final_well = V_well + V_add
replace: V_add = V_well ,  V_final_well = V_well ,  f = 1
```

**Serial-dilution volumes.** Each tube is prepared to one uniform total so the
transfer drawn into the next tube plus the volume dispensed to the wells
balance the fold ratio. With `overage` (default 1.1) and `extra_dead_uL`:

```
V_dispense = V_add · replicates · overage + extra_dead_uL   # drawn from each tube for its wells
V_transfer = V_dispense / (fold − 1)                        # tube → next tube
V_tube      = V_dispense + V_transfer                        # uniform prep volume per tube
tube 1: stock  = V_tube · working_conc,0 / C_stock ,  medium = V_tube − stock
tube i>1: transfer = V_transfer , medium = V_dispense
```

`plate_well_volume(plate_type)` returns the default working volume (6→2000,
12→1000, 24→500, 48→250, 96→100, 384→50, 1536→10 µL).
`plan_serial_dilution(...)` → `DosingPlan` with per-tube recipe, `total_stock_uL`,
`total_diluent_uL`, `total_wells`, and feasibility `warnings` (top working conc
exceeding stock; sub-microliter stock draw).

**Golden fixture:** `C_stock=1000`, `C_top=10` µM, `fold=4`, `n=4`,
`replicates=3`, `V_well=100` µL, spike `f=10`, `overage=1.0`.
- `C_final` → **[10, 2.5, 0.625, 0.15625]** µM; `working[0]` → **100** µM
- `V_add` (spike) → **11.1111** µL; `V_final_well` → **111.1111** µL
- `V_dispense` → **33.3333** µL; `V_transfer` → **11.1111** µL; `V_tube` → **44.4444** µL
- tube-1 stock → **4.44444** µL; `total_wells` → **12**
- replace mode (`C_stock=500, C_top=50, fold=2, n=6, reps=4, V_well=200`):
  `working = C_final`, `V_add` → **200** µL

---

## 8d. Plate map / layout

Assign the dose series (concentration points × replicates) to physical wells of
a standard microplate. Supported geometries (rows × columns):

| Plate | Grid  | Rows          | Cols |
|-------|-------|---------------|------|
| 6     | 2×3   | A–B           | 1–3  |
| 12    | 3×4   | A–C           | 1–4  |
| 24    | 4×6   | A–D           | 1–6  |
| 48    | 6×8   | A–F           | 1–8  |
| 96    | 8×12  | A–H           | 1–12 |
| 384   | 16×24 | A–P           | 1–24 |
| 1536  | 32×48 | A–AF          | 1–48 |

Row labels are Excel-style (`0→A … 25→Z, 26→AA …`).

**Orientations.**

```
by_column : each concentration point occupies one column (dose increases
            left→right); replicates fill the rows of that column.
by_row    : each point occupies one row; replicates fill columns.
sequential: row-major fill, replicates grouped consecutively.
```

If the requested grid layout cannot hold the series (`n > cols` for `by_column`,
or `replicates > rows`; symmetric for `by_row`), the planner falls back to
`sequential` and records a warning. A second warning is added when
`n · replicates > rows · cols` (plate capacity); over-capacity wells are
returned with `in_bounds = False` rather than dropped.

`plate_map(final_concs_uM, replicates, plate_type, orientation="by_column")`
→ `PlateMap` with `rows`, `cols`, `orientation` (actual), `capacity`, `used`,
per-well `PlateWell` list (`row`, `col`, `row_label`, `col_label`, `point`,
`replicate`, `final_conc_uM`, `in_bounds`) and `warnings`.
`plate_map_to_csv(pmap)` → deterministic CSV, one row per assigned well
(`well,row,column,point,replicate,final_conc_uM`).

**Golden fixture:** `final_concs=[10, 3.333, 1.111, 0.37]`, `replicates=3`,
`plate="96"`, `orientation="by_column"`.
- `orientation` (actual) → **by_column**; `used` → **12**; `capacity` → **96**
- well[0] → **A1** (point 0, replicate 1); well[1] → **B1** (point 0, replicate 2)
- well[3] → **A2** (point 1, replicate 1); `_row_label(26)` → **AA**
- fallback: 14 points × 3 replicates, 96-well, `by_column` → `orientation` =
  **sequential**, non-empty `warnings`, `used` = **42**

---

## 8e. Dose units & paint-first selection helpers

The In Vitro Dosing UI is **paint-first**: the user marquee-selects a
rectangular block of wells on the plate, and the block geometry drives the
dilution series. The ADC stock is always entered in **µM**; per-group dose
concentrations are entered in a selectable unit (µM / nM / pM) and converted.

`CONC_UNIT_TO_UM = {uM: 1.0, nM: 1e-3, pM: 1e-6}` — multiplicative factors to
µM. `convert_concentration(value, from_unit, to_unit)` /
`convertConcentration` converts via µM: `value * CONC_UNIT_TO_UM[from] /
CONC_UNIT_TO_UM[to]`.

`series_shape_from_selection(n_rows, n_cols, orientation="by_column")` /
`seriesShapeFromSelection` maps a painted block's dimensions to
`(n_points, replicates)`: `by_column` → `(n_cols, n_rows)` (dose varies across
columns, replicates down rows); `by_row` → `(n_rows, n_cols)` (dose decreases
top→bottom, replicates across columns). `by_row` is the app default so a series
reads high-at-top.

`assign_selection(cells, final_concs_uM, orientation)` / `assignSelection`
takes a list of `(row, col)` wells and the computed final concentrations and
returns `SelectionAssignment` with per-well `SelectionWell`
(`row, col, point, replicate, final_conc_uM`), the block `origin (r0, c0)`,
`rows`, `cols`, `rectangular`, `orientation`, and `warnings`. Point and
replicate are read from each well's position **within the selection bounding
box** (handles offset/non-origin blocks). Concentration is `None` for wells
past the series length; a warning is emitted when the series length differs
from the dose-axis span.

`aggregate_dosing_plans(plans)` / `aggregateDosingPlans` sums a list of
per-group `DosingPlan`s → `{n_groups, total_stock_uL, total_diluent_uL,
total_wells}` (JS: `nGroups, totalStockUL, totalDiluentUL, totalWells`).

**Golden fixtures:**
- `convert_concentration(5, "uM", "nM")` → **5000**; `(250, "nM", "uM")` →
  **0.25**.
- `series_shape_from_selection(4, 3, "by_row")` → **(4, 3)**;
  `(4, 3, "by_column")` → **(3, 4)**.
- `assign_selection` on an offset 4×3 block (rows 2–5 × cols 5–8),
  `final_concs_uM = [10, 2.5, 0.625, 0.15625]`, `orientation="by_row"`:
  well `(2,5)` → point 0, replicate 1, conc **10**; `(2,6)` → point 0,
  replicate 2; `(5,5)` → point 3, conc **0.15625**; `rectangular` → **True**.

---

## 9. Golden values (regression fixtures)

From the spreadsheet inputs: `MW_mAb=145000`, `ε280_mAb_mgml=1.4`
(⇒ `ε280_mAb=203000`), `A280=7.4315`, `Aλmax=0.7827` (⇒ `R=0.105322`),
`MW_LP=1559.62`, `ε280_LP=7287`, `ελmax_LP=9624`, `DF=1`.

| Quantity | Case A (ελmax_mAb = R·ε280, circular) | Case B (ελmax_mAb = 0, realistic) |
|---|---|---|
| DAR_UV | **0.0** | **2.41408** |
| MW_ADC | 145000.0 | 148765.049 |
| ε280_ADC_molar | 203000.0 | 220591.408 |
| ε280_ADC_mass | 1.4 | 1.482817 |
| conc_mgml | **5.308214** | 5.011743 |
| conc_uM | **36.608374** | 33.688982 |

Other fixtures:
- Ellman `(A412=0.5, A280=0.6, ε280=203000)` → **11.95524** thiols
- Ellman blank-subtracted `(A412=0.5, A280=0.6, ε280=203000, ε412=14150, A412_blank=0.02, A280_blank=0.01)` → **11.671558** (lower than raw); zero blanks reduce to raw
- Scatter `(A_ref=0.05, λ_ref=320, λ=280, p=4)` → **0.085298**; corrected `(A280=1.0, Aλmax=0.5, A_ref=0.05, 320, 495, 4)` → A280→**0.914702**, Aλmax→**0.491267**
- HIC `{0:5,1:15,2:45,3:25,4:10}` → **2.2**
- HIC corrected `{0:10,2:50,4:30,6:10}, ε280_mab=203000, ε280_payload=5000` → **2.741596** (uncorrected 2.8); species% DAR2 → **50.870418**
- LC-MS reduced `LC{0:10,1:90}, HC{0:5,1:20,2:75}` → **5.2**
- LC-MS intact response-corrected `I{0:10,2:30,4:60}, r{0:1.0,2:0.8,4:0.6}` → **3.220339** (raw 3.0; rises as under-ionized high-load species are up-weighted); no response factors → raw
- LC-MS reduced response-corrected `LC{0:10,1:90}, HC{0:5,1:20,2:75}, r_LC{0:1.0,1:0.9}, r_HC{0:1.0,1:0.9,2:0.8}` → **5.285461**; no response factors → raw
- DAR distribution (cysteine, `d=2`) `n=4, feed_ratio=2.5, efficiency=0.8` → p_site=**0.25**, mean=**2.0**, variance=**3.0**, support **{0,2,4,6,8}**, P(DAR=2)=**0.421875**, P(DAR=0)=**0.31640625**
- DAR distribution (cysteine, `d=2`) `n=4, p_site=0.5` → mean=**4.0**, variance=**4.0**, P(DAR=4)=**0.375**, P(DAR=0)=P(DAR=8)=**0.0625**
- DAR distribution (lysine, `d=1`) `n=8, p_site=0.5` → mean=**4.0**, variance=**2.0**, support **{0..8}**, P(DAR=3)=**0.21875**
- Measured DAR distribution `measured_dar_distribution({0:5, 2:35, 4:40, 6:15, 8:5})` → mean_dar=**3.6**, variance=**3.44**, sd=**1.854723699099141**, total=**100.0**, fractions[2]=**0.35**, fractions[4]=**0.40**; negative abundance or zero total → `ValueError`
- DAR-UV uncertainty (Case B, σA280=σAλmax=0.01) → σ_DAR=**0.033701**; with ε terms (σε280_mAb=2030, σε280_LP=364.35, σελmax_LP=481.2) → σ_DAR=**0.137956** (ε-coefficient uncertainty dominates)
- Physical bounds `check_physical_bounds`: clean `(dar=2.5, r=0.1, ε280=203000, ελmax=210000, conc=5)` → **[]** (no warnings); `dar=-1.2` → **[dar_negative]**; `dar=20` → **[dar_high]** (but `dar=16` → []); `r=-0.5` → **[r_negative]**; `r=60` → **[r_high]** (but `r=50` → []); `ε280=203000, ελmax=5000` → **[eps_inconsistent]** (but ε280==ελmax → []); `conc=-3` → **[conc_negative]**; collect-all: `(dar=-1, r=-0.5, ε280=1000, ελmax=100, conc={c1:-2,c2:5})` → **4 warnings** `{conc_negative, dar_negative, eps_inconsistent, r_negative}`
- Payload library `PAYLOAD_LIBRARY`: **9 entries**; **3** with `sourced=True`
  (vc-MMAE, DM1, DM4). `get_payload("vc-MMAE")` → `eps_lmax_conj=15900`,
  `eps280_conj=1500`, `lambda_max=248`, `mw_conj=1316.63`, `cls="auristatin"`;
  `get_payload("DM1")` → `eps_lmax_free=26790`, `eps280_free=0`, `lambda_max=252`;
  `get_payload("DM4")` → `eps_lmax_free=28044`, `eps280_free=5700`;
  `get_payload("mc-MMAF")` → `sourced=False`, `eps_lmax_conj=None`;
  `get_payload("nonexistent")` → `None`.
- SEC purity `sec_purity(1900, 62, 38)` → %monomer=**95.0**, %HMW=**3.1**, %LMW=**1.9**, total=**2000** (sum = 100); two-peak `sec_purity(4900, 100)` → %monomer=**98.0**, %HMW=**2.0**, %LMW=**0.0**; negative area or zero total → `ValueError`
- Free-drug `free_drug_percent(1.5, 598.5)` → %free=**0.25**, %conjugated=**99.75**, total=**600**; negative amount or zero total → `ValueError`

Tolerance for all: relative 1e-6.

---

## 9c. Payload–linker reference library (`PAYLOAD_LIBRARY`, `get_payload`)

Extinction coefficients for the common ADC cytotoxic payloads, for DAR-by-UV
(the drug's `ελmax` plus its `ε280`, the latter used to correct the mAb A280
contribution). `ε` in M⁻¹cm⁻¹, `MW` in g/mol, `λmax` in nm.

Each entry carries: `name`, `cls` (payload class), `mw_free`, `mw_conj`,
`lambda_max`, `eps_lmax_free`, `eps280_free`, `eps_lmax_conj`, `eps280_conj`,
`sourced` (bool), `source` (citation string).

**FREE vs CONJUGATED are kept distinct.** Published ε are reported for either
the free payload or the drug-linker construct, and the two are NOT
interchangeable; a field that is not firmly citable is `null`/`None`. When an
entry's needed ε is null the UI must show "ε not in library — measure or enter
manually" and must NOT substitute a guessed value into a logged DAR.
`sourced=True` means at least one ε pair (free or conjugated) is a cited value;
`sourced=False` is a reference entry (MW/λmax/class only, all ε null).

`get_payload(name)` is a case-insensitive exact-name lookup returning the entry
dict or `None`.

Firmly-sourced ε (as of this build):
- **vc-MMAE** (mc-vc-PAB-MMAE, auristatin, λmax 248): conjugated ελmax=15,900,
  ε280=1,500. Cruz & Kayser DAR-by-UV/Vis, Cancers 2009 11(6):870.
- **DM1** (maytansinoid, λmax 252): free ελmax=26,790, ε280=0 (280 contribution
  treated as negligible; standard mAb correction ratio εab252/εab280≈0.378).
  US 7,666,414 Ex.17.
- **DM4** (SPDB-DM4, maytansinoid, λmax 252): free ελmax=28,044, ε280=5,700.
  US 10,780,179.

Reference-only entries (`sourced=False`, ε null): mc-MMAF, deruxtecan (DXd),
SN-38, tesirine (SG3199), calicheamicin γ1, α-amanitin.

---

## 9d. Product-quality aggregates: SEC purity & free-drug (`sec_purity`, `free_drug_percent`)

Two area-/amount-normalised critical quality attributes reported for every ADC lot.

**SEC purity (`sec_purity`).** Size-exclusion chromatography separates by
hydrodynamic size into three classes: **monomer** (the intact one-antibody
product), **HMW** (high-molecular-weight dimers/aggregates, eluting earlier
because they are larger — a stability- and immunogenicity-related CQA), and
**LMW** (low-molecular-weight fragments/clips, eluting later). Percentages are
area-normalised, the standard SEC reporting convention:

    %class_i = AUC_i / (AUC_monomer + AUC_HMW + AUC_LMW) × 100

`auc_lmw` defaults to 0 for two-peak (monomer/HMW-only) integration. All areas
must be non-negative and the total positive, else `ValueError`. Returns
`{pct_monomer, pct_hmw, pct_lmw, total}`.

**Free-drug (`free_drug_percent`).** The fraction of payload present as
unconjugated (free) small molecule — a safety-critical CQA, typically measured
by reversed-phase HPLC or an SEC/UPLC free-drug method. Inputs are amounts in
matching units (molar concentration, peak area, or mass); the percentage is
unit-independent:

    %free = free / (free + conjugated) × 100

Both inputs non-negative, sum positive, else `ValueError`. Returns
`{pct_free, pct_conjugated, total}`.

Neither aggregate models a distribution — they are direct normalisations of
measured integrals, so no binomial/site-occupancy assumptions apply.

## 10. Embedded dye/payload library

22 entries (Alexa Fluor series + pHrodo) parsed from the spreadsheet
(`dye_library.json`), each with `name, mw, e_lmax (M⁻¹cm⁻¹), cf280`
(correction factor: fraction of λmax ε contributed at 280 nm),
`e280 = e_lmax·cf280`, `lmax (nm)`, and `comment`.

## 10b. Physical-bounds guards (`check_physical_bounds`)

A pure plausibility check that flags computed quantities which are physically
impossible or almost certainly the result of a swapped/mis-entered input. It
**never raises** and never alters a calculation — it returns a list of
`{code, message}` warnings so the UI can surface them inline while still
showing the (suspect) number. Callers pass only the quantities they have;
`None`/absent inputs are skipped.

Checks and their `code`s:

| Condition | `code` | Rationale |
|---|---|---|
| DAR < 0 | `dar_negative` | A below-blank absorbance drives DAR negative; usually swapped Aλmax/A280 or ε. |
| DAR > `DAR_MAX_PLAUSIBLE` (16) | `dar_high` | Standard/site-specific ADCs are 2–8; > 16 signals bad ε or absorbances. |
| R = Aλmax/A280 < 0 | `r_negative` | A measured absorbance is below its blank. |
| R > `R_MAX_PLAUSIBLE` (50) | `r_high` | A280 and Aλmax likely swapped. |
| ε280 > ελmax (with ελmax > 0) | `eps_inconsistent` | λmax is the absorbance maximum, so ε280 ≤ ελmax; violation ⇒ swapped/mislabeled coefficients. |
| any concentration < 0 | `conc_negative` | Concentrations cannot be below zero. |

Boundaries are inclusive-below: exactly `DAR_MAX_PLAUSIBLE`, exactly
`R_MAX_PLAUSIBLE`, and ε280 == ελmax all pass without warning. Multiple
conditions are collected in one call (collect-all, not fail-fast).

Constants: `DAR_MAX_PLAUSIBLE = 16.0`, `R_MAX_PLAUSIBLE = 50.0`.
