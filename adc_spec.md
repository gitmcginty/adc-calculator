# ADC Calculator — Calculation Core Specification

This document defines every formula implemented by the calculation core
(`adc_core.py` / `adc_core.js`). It is the single source of truth: both
implementations and the test suite are validated against the **golden values**
in §9. Sources: `Calculation_DAR_conjugation.xlsx` (UV/Vis DAR sheet + dye
library) and `250811 Bioconjugation protocols and analyses.pdf`.

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

## 4. Ellman's free-thiol assay (PDF p.10)

```
free_thiols = (A412 / ε412) / (A280 / ε280_mAb)      with ε412 = 14150 M⁻¹cm⁻¹
```
`A280` and `A412` are the absorbances of the DTNB-reacted protein aliquot;
`ε280_mAb` is the molar extinction coefficient of the protein.

---

## 5. DAR by analytical HIC (PDF p.14–15)

AUC-weighted average over resolved DAR species (peak area at 280 nm ∝ AUC):
```
DAR_HIC = Σ_i (AUC_i · i) / Σ_i AUC_i           i = 0,1,2,…,n
```

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
- HIC `{0:5,1:15,2:45,3:25,4:10}` → **2.2**
- LC-MS reduced `LC{0:10,1:90}, HC{0:5,1:20,2:75}` → **5.2**

Tolerance for all: relative 1e-6.

---

## 10. Embedded dye/payload library

22 entries (Alexa Fluor series + pHrodo) parsed from the spreadsheet
(`dye_library.json`), each with `name, mw, e_lmax (M⁻¹cm⁻¹), cf280`
(correction factor: fraction of λmax ε contributed at 280 nm),
`e280 = e_lmax·cf280`, `lmax (nm)`, and `comment`.
