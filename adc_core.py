"""adc_core — pure calculation core for the ADC Calculator.

Zero UI, zero I/O (the dye library is embedded as a literal). Every function
maps directly to a formula in ``adc_spec.md`` and is validated by
``test_adc_core.py`` against the golden values in the spreadsheet
``Calculation_DAR_conjugation.xlsx`` and the protocols PDF.

Units are stated in every signature. Molar extinction coefficients are in
M^-1 cm^-1; mass extinction coefficients in (mg/mL)^-1 cm^-1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

# --------------------------------------------------------------------------
# Constants (from protocols PDF)
# --------------------------------------------------------------------------
ELLMAN_E412: float = 14150.0          # M^-1 cm^-1, DTNB product at 412 nm (p.10)
IGG1_INTERCHAIN_DISULFIDES: int = 4   # -> 8 free thiols on full reduction

# ==========================================================================
# 1. Extinction-coefficient conversions  (spec §1)
# ==========================================================================
def eps_mass_to_molar(eps_mass: float, mw: float) -> float:
    """(mg/mL)^-1 cm^-1  ->  M^-1 cm^-1.  eps_molar = eps_mass * MW."""
    return eps_mass * mw


def eps_molar_to_mass(eps_molar: float, mw: float) -> float:
    """M^-1 cm^-1  ->  (mg/mL)^-1 cm^-1.  eps_mass = eps_molar / MW."""
    if mw == 0:
        raise ValueError("MW must be non-zero")
    return eps_molar / mw


# ==========================================================================
# 2. DAR by UV/Vis  (spec §2; spreadsheet I16)
# ==========================================================================
def dar_uv(
    a280: float,
    a_lmax: float,
    eps280_mab: float,
    eps_lmax_mab: float,
    eps280_lp: float,
    eps_lmax_lp: float,
) -> float:
    """Drug-to-antibody ratio from a two-wavelength UV/Vis measurement.

    DAR = (eps_lmax_mab - R*eps280_mab) / (R*eps280_lp - eps_lmax_lp),
    with R = a_lmax / a280.

    eps_lmax_mab is the antibody's molar extinction at the payload lambda-max
    (~0 for a bare IgG in the visible; the spreadsheet's circular definition
    R*eps280_mab reproduces DAR = 0).
    """
    if a280 == 0:
        raise ValueError("A280 must be non-zero")
    r = a_lmax / a280
    denom = r * eps280_lp - eps_lmax_lp
    if denom == 0:
        raise ZeroDivisionError("Degenerate optics: R*eps280_LP == eps_lmax_LP")
    return (eps_lmax_mab - r * eps280_mab) / denom


# ==========================================================================
# 3. ADC molecular properties  (spec §3)
# ==========================================================================
@dataclass(frozen=True)
class AdcProperties:
    dar: float
    mw_adc: float            # g/mol
    eps280_adc_molar: float  # M^-1 cm^-1
    eps280_adc_mass: float   # (mg/mL)^-1 cm^-1
    conc_mgml: float | None  # mg/mL  (None if A280 not supplied)
    conc_uM: float | None    # micromolar


def adc_properties(
    dar: float,
    mw_mab: float,
    mw_lp: float,
    eps280_mab: float,
    eps280_lp: float,
    a280: float | None = None,
    dilution_factor: float = 1.0,
) -> AdcProperties:
    """MW, corrected epsilon280 (molar & mass) and, if A280 given, concentration."""
    mw_adc = mw_mab + mw_lp * dar
    eps_molar = eps280_mab + dar * eps280_lp
    eps_mass = eps_molar / mw_adc
    conc_mgml = conc_uM = None
    if a280 is not None:
        conc_mgml = a280 * dilution_factor / eps_mass
        conc_uM = conc_mgml / mw_adc * 1e6
    return AdcProperties(dar, mw_adc, eps_molar, eps_mass, conc_mgml, conc_uM)


def conc_unknown_payload_eps(
    a280: float,
    a_lmax: float,
    eps280_lp: float,
    eps_lmax_lp: float,
    eps280_mab_mass: float,
    dilution_factor: float = 1.0,
) -> float:
    """Mass concentration (mg/mL) when the payload eps is unknown (spec §3, PDF p.12).

    Subtracts the payload's estimated 280 nm contribution, then applies the
    naked-mAb mass extinction coefficient. DAR is not obtainable this way.
    """
    a_lp_280 = a_lmax * eps280_lp / eps_lmax_lp
    return (a280 - a_lp_280) * dilution_factor / eps280_mab_mass


# ==========================================================================
# 4. Ellman's free-thiol assay  (spec §4; PDF p.10)
# ==========================================================================
def ellman_free_thiols(
    a412: float,
    a280: float,
    eps280_mab: float,
    eps412: float = ELLMAN_E412,
) -> float:
    """Free sulfhydryl-to-protein ratio = (A412/eps412) / (A280/eps280_mab)."""
    if a280 == 0:
        raise ValueError("A280 must be non-zero")
    return (a412 / eps412) / (a280 / eps280_mab)


# ==========================================================================
# 5. DAR by analytical HIC  (spec §5; PDF p.14-15)
# ==========================================================================
def dar_hic(auc_by_species: Mapping[int, float]) -> float:
    """AUC(280 nm)-weighted average DAR over resolved species.

    auc_by_species maps DAR value (0,1,2,...) -> peak area.
    """
    total = sum(auc_by_species.values())
    if total == 0:
        raise ValueError("Total AUC must be non-zero")
    return sum(k * a for k, a in auc_by_species.items()) / total


def dar_species_fractions(auc_by_species: Mapping[int, float]) -> dict[int, float]:
    """Relative percentage of each DAR species (AUC_i / sum AUC * 100)."""
    total = sum(auc_by_species.values())
    if total == 0:
        raise ValueError("Total AUC must be non-zero")
    return {k: a / total * 100.0 for k, a in auc_by_species.items()}


# ==========================================================================
# 6. DAR by LC-MS  (spec §6; PDF p.17-18)
# ==========================================================================
def _weighted_mean(intensity_by_k: Mapping[int, float]) -> float:
    total = sum(intensity_by_k.values())
    if total == 0:
        raise ValueError("Total intensity must be non-zero")
    return sum(k * i for k, i in intensity_by_k.items()) / total


def dar_lcms_reduced(
    light_chain: Mapping[int, float],
    heavy_chain: Mapping[int, float],
) -> float:
    """Denaturing + reducing DAR = 2*mean(LC) + 2*mean(HC)  (PDF p.17)."""
    return 2 * _weighted_mean(light_chain) + 2 * _weighted_mean(heavy_chain)


def dar_lcms_intact(intensity_by_k: Mapping[int, float]) -> float:
    """Denaturing non-reducing / native DAR = intensity-weighted mean (PDF p.17-18)."""
    return _weighted_mean(intensity_by_k)


def _chain_stats(intensity_by_k: Mapping[int, float]) -> dict:
    """Per-chain intermediates for the reduced LC-MS DAR.

    Returns sum_i (ΣI), sum_ki (Σk·I), mean (Σk·I/ΣI), and contribution
    (2·mean, the chain's additive share of DAR since an IgG has two of each
    chain).  Surfacing these makes a mass-vs-intensity mix-up or a dropped
    peak visible instead of hidden inside a single DAR number.
    """
    sum_i = sum(intensity_by_k.values())
    if sum_i == 0:
        raise ValueError("Total intensity must be non-zero")
    sum_ki = sum(k * i for k, i in intensity_by_k.items())
    mean = sum_ki / sum_i
    return {"sum_i": sum_i, "sum_ki": sum_ki, "mean": mean, "contribution": 2 * mean}


def dar_lcms_reduced_breakdown(
    light_chain: Mapping[int, float],
    heavy_chain: Mapping[int, float],
) -> dict:
    """Full reduced-DAR breakdown: per-chain stats + total DAR.

    {'light': {sum_i, sum_ki, mean, contribution},
     'heavy': {...}, 'dar': <light.contribution + heavy.contribution>}.
    dar is identical to dar_lcms_reduced(light_chain, heavy_chain).
    """
    light = _chain_stats(light_chain)
    heavy = _chain_stats(heavy_chain)
    return {"light": light, "heavy": heavy,
            "dar": light["contribution"] + heavy["contribution"]}


# ==========================================================================
# 7. Conjugation designer  (spec §7; PDF p.3-6)
# ==========================================================================
@dataclass(frozen=True)
class ReagentAddition:
    name: str
    equivalents: float
    basis: str            # "per mAb" or "per free thiol"
    moles: float          # mol
    stock_conc_mM: float
    volume_uL: float


@dataclass(frozen=True)
class ConjugationPlan:
    route: str
    mass_mab_mg: float
    n_mab_nmol: float
    free_thiols_per_mab: float
    additions: list[ReagentAddition]
    reaction_volume_uL: float
    dmso_volume_uL: float
    notes: list[str] = field(default_factory=list)


def moles_from_mass(mass_mg: float, mw: float) -> float:
    """mol from mass(mg) and MW(g/mol).  mol = (mass_mg/1000)/MW."""
    return (mass_mg / 1000.0) / mw


def reagent_volume_uL(moles: float, stock_conc_mM: float) -> float:
    """Volume (uL) of a stock at C(mM) delivering `moles` mol.

    V(L) = mol / C(mol/L); mM = mmol/L, so V(uL) = mol/ (C_mM*1e-3) * 1e6.
    """
    if stock_conc_mM == 0:
        raise ValueError("Stock concentration must be non-zero")
    liters = moles / (stock_conc_mM * 1e-3)
    return liters * 1e6


def plan_conjugation(
    route: str,
    mass_mab_mg: float,
    mw_mab: float = 145000.0,
    mab_conc_g_per_L: float = 5.0,
    reagents: list[dict] | None = None,
    free_thiols_per_mab: float | None = None,
    dmso_fraction: float = 0.10,
) -> ConjugationPlan:
    """Build a conjugation plan: reagent volumes from equivalents + top-up volumes.

    `reagents` is a list of dicts: {name, equiv, basis, stock_mM} where basis is
    "mAb" or "thiol". If free_thiols_per_mab is None it defaults per route
    (cys/thiomab full reduction -> IGG1_INTERCHAIN_DISULFIDES*2 = 8).
    """
    n_mab = moles_from_mass(mass_mab_mg, mw_mab)          # mol
    if free_thiols_per_mab is None:
        free_thiols_per_mab = (
            IGG1_INTERCHAIN_DISULFIDES * 2
            if route in ("cys_native", "cys_rebridge", "thiomab")
            else 0.0
        )
    n_thiol = n_mab * free_thiols_per_mab

    base_vol_uL = mass_mab_mg / mab_conc_g_per_L * 1000.0  # mg / (g/L) = mL -> uL

    additions: list[ReagentAddition] = []
    for rg in (reagents or []):
        basis = rg.get("basis", "mAb")
        n_ref = n_thiol if basis == "thiol" else n_mab
        moles = rg["equiv"] * n_ref
        vol = reagent_volume_uL(moles, rg["stock_mM"])
        additions.append(
            ReagentAddition(
                name=rg["name"],
                equivalents=rg["equiv"],
                basis="per free thiol" if basis == "thiol" else "per mAb",
                moles=moles,
                stock_conc_mM=rg["stock_mM"],
                volume_uL=vol,
            )
        )

    # Final reaction volume so DMSO ends at dmso_fraction, given DMSO-borne reagents.
    reagent_vol = sum(a.volume_uL for a in additions)
    # total volume V such that dmso_volume / V = dmso_fraction, with dmso_volume
    # at least the reagent volume carried in DMSO. We size DMSO to hit the target.
    reaction_volume = max(base_vol_uL / (1 - dmso_fraction), base_vol_uL + reagent_vol)
    dmso_volume = reaction_volume * dmso_fraction

    notes = _route_notes(route)
    return ConjugationPlan(
        route=route,
        mass_mab_mg=mass_mab_mg,
        n_mab_nmol=n_mab * 1e9,
        free_thiols_per_mab=free_thiols_per_mab,
        additions=additions,
        reaction_volume_uL=reaction_volume,
        dmso_volume_uL=dmso_volume,
        notes=notes,
    )


def _route_notes(route: str) -> list[str]:
    return {
        "lysine": [
            "Activated ester 5-20 equiv (optimize by screening).",
            "mAb at 5 g/L; add 10% v/v 1.0 M NaHCO3; final 10% v/v DMSO.",
            "2 h at 25 C; quench with L-lysine 100 equiv, 20 min.",
        ],
        "cys_native": [
            "1 equiv TCEP per disulfide; 6.0 equiv for full IgG1 reduction.",
            "Maleimide 1.5 equiv per free thiol; +5 mM EDTA; 10% v/v DMSO.",
            "Maleimide 1 h / 25 C; quench L-cysteine 100 equiv, 10 min.",
        ],
        "cys_rebridge": [
            "Full reduction: 6.0 equiv TCEP, 1 h / 37 C.",
            "Rebridging reagent equivalents by screening; 16 h / 25 C.",
            "Purify by preparative HIC for defined DAR.",
        ],
        "thiomab": [
            "Full reduction 10 equiv TCEP -> desalt to remove capping cysteines.",
            "Reoxidize with dhAA 30 equiv, 4 h / 25 C.",
            "Maleimide 1.5 equiv per free thiol, 1 h / 25 C.",
        ],
        "mtgase": [
            "PNGase F 10 uL/mg to deglycosylate; tag 80 equiv (e.g. azido-PEG3-amine).",
            "MTGase 5.5 U/mg, 16 h / 25 C -> site-specific DAR 2.0.",
            "Click reagent (e.g. DBCO-LP) 3.0 equiv, 6 h / 25 C.",
        ],
    }.get(route, [])


# ==========================================================================
# 8. Yield & formulation  (spec §8)
# ==========================================================================
@dataclass(frozen=True)
class YieldResult:
    recovered_mass_mg: float
    yield_pct: float | None
    molar_amount_nmol: float
    v_final_mL: float | None       # to reach target concentration
    v_change_mL: float | None      # +add diluent / -concentrate


def yield_and_formulation(
    conc_mgml: float,
    volume_mL: float,
    mw_adc: float,
    starting_mass_mg: float | None = None,
    target_conc_mgml: float | None = None,
) -> YieldResult:
    """Recovered mass, % yield, molar amount, and volume to a target concentration."""
    recovered = conc_mgml * volume_mL
    yld = (recovered / starting_mass_mg * 100.0) if starting_mass_mg else None
    nmol = recovered / mw_adc * 1e6
    v_final = v_change = None
    if target_conc_mgml:
        v_final = recovered / target_conc_mgml
        v_change = v_final - volume_mL
    return YieldResult(recovered, yld, nmol, v_final, v_change)


# ==========================================================================
# 10. Embedded dye / payload library  (spreadsheet cols L-R)
# ==========================================================================
# name, mw (g/mol), e_lmax (M^-1 cm^-1), cf280 (280 nm correction factor),
# e280 (= e_lmax*cf280), lmax (nm), comment
DYE_LIBRARY: list[dict] = [
    {"name": "AF350", "mw": 410, "e_lmax": 19000, "cf280": 0.19, "e280": 3610, "lmax": 350, "comment": None},
    {"name": "AF405", "mw": 1028, "e_lmax": 34000, "cf280": 0.70, "e280": 23800, "lmax": 405, "comment": None},
    {"name": "AF430", "mw": 702, "e_lmax": 16000, "cf280": 0.28, "e280": 4480, "lmax": 430, "comment": None},
    {"name": "AF488", "mw": 643, "e_lmax": 71000, "cf280": 0.11, "e280": 7810, "lmax": 488, "comment": None},
    {"name": "AF500", "mw": 700, "e_lmax": 71000, "cf280": 0.18, "e280": 12780, "lmax": 500, "comment": None},
    {"name": "AF514", "mw": 714, "e_lmax": 80000, "cf280": 0.18, "e280": 14400, "lmax": 514, "comment": None},
    {"name": "AF532", "mw": 721, "e_lmax": 81000, "cf280": 0.09, "e280": 7290, "lmax": 532, "comment": None},
    {"name": "AF546", "mw": 1260, "e_lmax": 112000, "cf280": 0.12, "e280": 13440, "lmax": 546, "comment": None},
    {"name": "AF555", "mw": 1250, "e_lmax": 150000, "cf280": 0.08, "e280": 12000, "lmax": 555, "comment": None},
    {"name": "AF568", "mw": 792, "e_lmax": 91300, "cf280": 0.46, "e280": 41998, "lmax": 568, "comment": None},
    {"name": "AF594", "mw": 820, "e_lmax": 90000, "cf280": 0.56, "e280": 50400, "lmax": 594, "comment": None},
    {"name": "AF610X", "mw": 1285, "e_lmax": 132000, "cf280": 0.44, "e280": 58080, "lmax": 610, "comment": None},
    {"name": "AF633", "mw": 1200, "e_lmax": 100000, "cf280": 0.55, "e280": 55000, "lmax": 633, "comment": None},
    {"name": "AF647", "mw": 1300, "e_lmax": 239000, "cf280": 0.03, "e280": 7170, "lmax": 647, "comment": None},
    {"name": "AF660", "mw": 1100, "e_lmax": 132000, "cf280": 0.10, "e280": 13200, "lmax": 660, "comment": None},
    {"name": "AF680", "mw": 1150, "e_lmax": 184000, "cf280": 0.05, "e280": 9200, "lmax": 680, "comment": None},
    {"name": "AF700", "mw": 1400, "e_lmax": 192000, "cf280": 0.07, "e280": 13440, "lmax": 700, "comment": None},
    {"name": "AF750", "mw": 1300, "e_lmax": 240000, "cf280": 0.04, "e280": 9600, "lmax": 750, "comment": None},
    {"name": "AF790", "mw": 1750, "e_lmax": 260000, "cf280": 0.08, "e280": 20800, "lmax": 790, "comment": None},
    {"name": "pHRodo deep red", "mw": 1300, "e_lmax": 140000, "cf280": 0.33, "e280": 46200, "lmax": None, "comment": "Dilute 1:1 with H2SO4 1M"},
    {"name": "pHRodo iFL Red", "mw": 1000, "e_lmax": 65000, "cf280": 0.12, "e280": 7800, "lmax": 560, "comment": "Dilute 1:1 with H2SO4 1M"},
    {"name": "pHRodo iFL Green", "mw": 1000, "e_lmax": 74500, "cf280": 0.20, "e280": 14900, "lmax": None, "comment": "Dilute 1:1 with H2SO4 1M"},
]


def get_dye(name: str) -> dict | None:
    """Case-insensitive lookup in the embedded dye library."""
    key = name.strip().lower()
    for d in DYE_LIBRARY:
        if d["name"].lower() == key:
            return d
    return None


# ============================================================================
# ADC registry: log each preparation / assay and aggregate across batches.
#
# These are pure functions on plain dicts/lists so the same normalization and
# summary math can be exercised by the test suite and mirrored exactly in JS.
# Persistence itself (localStorage / files) lives in the UI layer.
# ============================================================================

# Canonical, ordered set of columns for a registry record (also the CSV header).
REGISTRY_FIELDS = [
    "id", "date", "operator", "assay", "route",
    "dar", "conc_mgml", "conc_uM", "free_thiols", "notes",
]
# Fields that are stored as numbers (coerced; blank -> None).
_REGISTRY_NUMERIC = {"dar", "conc_mgml", "conc_uM", "free_thiols"}


def _coerce_num(value):
    """Blank/None -> None; otherwise float(value) or None if not parseable."""
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    s = str(value).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def make_registry_record(fields: Mapping) -> dict:
    """Normalize a raw entry into a canonical registry record.

    Every canonical field is present in the output (missing -> None / ""),
    numeric fields are coerced to float or None, text fields are stripped.
    Raises ValueError if no ``id`` is supplied (a record must be identifiable).
    """
    rec = {}
    for key in REGISTRY_FIELDS:
        val = fields.get(key) if hasattr(fields, "get") else None
        if key in _REGISTRY_NUMERIC:
            rec[key] = _coerce_num(val)
        else:
            rec[key] = "" if val is None else str(val).strip()
    if rec["id"] == "":
        raise ValueError("registry record requires a non-empty 'id'")
    return rec


def _fmt_cell(value) -> str:
    """Format one CSV cell identically in Python and JS (shortest round-trip;
    integer-valued floats print without a trailing '.0')."""
    if value is None:
        return ""
    if isinstance(value, float):
        if value == int(value):
            return str(int(value))
        return repr(value)
    return str(value)


def _csv_escape(cell: str) -> str:
    if any(c in cell for c in [',', '"', '\n', '\r']):
        return '"' + cell.replace('"', '""') + '"'
    return cell


def registry_to_csv(records) -> str:
    """Deterministic CSV (header + one line per record) over REGISTRY_FIELDS."""
    lines = [",".join(REGISTRY_FIELDS)]
    for r in records:
        lines.append(",".join(_csv_escape(_fmt_cell(r.get(k))) for k in REGISTRY_FIELDS))
    return "\n".join(lines)


def summarize_registry(records) -> dict:
    """Aggregate stats across logged records.

    Returns total count, counts by assay and by ADC id, and DAR statistics
    (mean, sample SD [ddof=1], min, max) over records with a numeric DAR.
    """
    n_total = len(records)
    by_assay: dict[str, int] = {}
    by_id: dict[str, int] = {}
    dars = []
    for r in records:
        a = r.get("assay") or "(none)"
        by_assay[a] = by_assay.get(a, 0) + 1
        i = r.get("id") or "(none)"
        by_id[i] = by_id.get(i, 0) + 1
        d = r.get("dar")
        if isinstance(d, (int, float)) and not isinstance(d, bool):
            dars.append(float(d))
    out = {
        "n_total": n_total,
        "n_ids": len(by_id),
        "by_assay": by_assay,
        "by_id": by_id,
        "dar_n": len(dars),
        "dar_mean": None, "dar_sd": None, "dar_min": None, "dar_max": None,
    }
    if dars:
        mean = sum(dars) / len(dars)
        out["dar_mean"] = mean
        out["dar_min"] = min(dars)
        out["dar_max"] = max(dars)
        if len(dars) > 1:
            var = sum((x - mean) ** 2 for x in dars) / (len(dars) - 1)
            out["dar_sd"] = var ** 0.5
        else:
            out["dar_sd"] = 0.0
    return out
