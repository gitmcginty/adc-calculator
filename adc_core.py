"""adc_core — pure calculation core for the ADC Calculator.

Zero UI, zero I/O (the dye library is embedded as a literal). Every function
maps directly to a formula in ``adc_spec.md`` and is pinned by
``test_adc_core.py`` to fixed golden values, so the Python core, the JS core
and the embedded HTML copy cannot drift from one another.

Units are stated in every signature. Molar extinction coefficients are in
M^-1 cm^-1; mass extinction coefficients in (mg/mL)^-1 cm^-1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import comb
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


def scatter_absorbance(
    a_ref: float,
    ref_wavelength: float,
    target_wavelength: float,
    exponent: float = 4.0,
) -> float:
    """Estimate the light-scattering baseline at ``target_wavelength``.

    Turbid or aggregated samples scatter light with a strong wavelength
    dependence: Rayleigh-type scattering falls off as ``lambda^-4``, so a
    baseline read at a wavelength where neither antibody nor payload absorbs
    (typically 320-340 nm) can be extrapolated to the analytical wavelengths.

        A_scatter(lambda) = A_ref * (ref_wavelength / lambda) ^ exponent

    ``exponent`` defaults to 4 (Rayleigh); empirically fitted exponents for
    protein aggregates are often 2-4, so it is exposed as a parameter.
    """
    if ref_wavelength <= 0 or target_wavelength <= 0:
        raise ValueError("wavelengths must be positive")
    return a_ref * (ref_wavelength / target_wavelength) ** exponent


def scatter_corrected_absorbance(
    a280: float,
    a_lmax: float,
    a_ref: float,
    ref_wavelength: float = 320.0,
    lmax_wavelength: float = 495.0,
    exponent: float = 4.0,
) -> dict:
    """Subtract the extrapolated scattering baseline from both UV/Vis reads.

    Given a scattering reference read ``a_ref`` at ``ref_wavelength`` (a region
    free of chromophore absorbance), returns the scatter-corrected A280 and
    A_lmax that should be fed to :func:`dar_uv` / :func:`adc_properties` so the
    DAR and concentration are not biased high by turbidity.

    Returns a dict with the corrected absorbances and the scatter estimates.
    Corrected absorbances are floored at 0 (a negative value means the
    reference read over-estimated the baseline and is physically meaningless).
    """
    s280 = scatter_absorbance(a_ref, ref_wavelength, 280.0, exponent)
    slmax = scatter_absorbance(a_ref, ref_wavelength, lmax_wavelength, exponent)
    return {
        "a280_corrected": max(a280 - s280, 0.0),
        "a_lmax_corrected": max(a_lmax - slmax, 0.0),
        "scatter_280": s280,
        "scatter_lmax": slmax,
    }


# ==========================================================================
# 4. Ellman's free-thiol assay  (spec §4; PDF p.10)
# ==========================================================================
def ellman_free_thiols(
    a412: float,
    a280: float,
    eps280_mab: float,
    eps412: float = ELLMAN_E412,
    a412_blank: float = 0.0,
    a280_blank: float = 0.0,
) -> float:
    """Free sulfhydryl-to-protein ratio = (A412/eps412) / (A280/eps280_mab).

    DTNB and its TNB(2-) product, plus any residual reductant (TCEP/DTT),
    absorb at both 412 and 280 nm.  A matched reagent blank (all reagents,
    no protein) measures that background; subtracting it before the ratio
    removes the bias.  ``a412_blank`` and ``a280_blank`` default to 0.0, so
    the call reduces exactly to the raw ratio when no blank is supplied.
    The net (blank-subtracted) A280 must be positive.
    """
    net_a280 = a280 - a280_blank
    if net_a280 <= 0:
        raise ValueError("net A280 (A280 - blank) must be positive")
    return ((a412 - a412_blank) / eps412) / (net_a280 / eps280_mab)


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


def _species_molar_abundance(
    auc_by_species: Mapping[int, float],
    eps280_mab: float,
    eps280_payload: float,
) -> dict[int, float]:
    """Convert 280 nm peak areas to *molar* abundance per DAR species.

    In HIC-UV280 the measured area of a DAR-k species is
    ``A_k = c_k * eps_k`` with ``eps_k = eps280_mab + k*eps280_payload`` (the
    payload absorbs at 280 nm, so each additional drug adds to the molar
    extinction). Molar abundance is therefore ``c_k proportional to A_k/eps_k``.
    Dividing area by eps_k removes the bias that otherwise over-weights
    high-DAR peaks and inflates the average DAR.
    """
    mol = {}
    for k, a in auc_by_species.items():
        eps_k = eps280_mab + k * eps280_payload
        if eps_k <= 0:
            raise ValueError("per-species eps280 must be > 0")
        mol[k] = a / eps_k
    return mol


def dar_hic_corrected(
    auc_by_species: Mapping[int, float],
    eps280_mab: float,
    eps280_payload: float,
) -> float:
    """Drug-load-corrected average DAR from HIC-UV280 peak areas.

    Divides each peak area by that species' 280 nm molar extinction
    ``eps280_mab + k*eps280_payload`` before taking the abundance-weighted
    mean, so area% is converted to molar% first. Reduces to :func:`dar_hic`
    only when ``eps280_payload == 0`` (payload transparent at 280 nm).
    """
    mol = _species_molar_abundance(auc_by_species, eps280_mab, eps280_payload)
    total = sum(mol.values())
    if total == 0:
        raise ValueError("Total corrected abundance must be non-zero")
    return sum(k * m for k, m in mol.items()) / total


def dar_species_fractions_corrected(
    auc_by_species: Mapping[int, float],
    eps280_mab: float,
    eps280_payload: float,
) -> dict[int, float]:
    """Molar percentage of each DAR species after 280 nm eps correction."""
    mol = _species_molar_abundance(auc_by_species, eps280_mab, eps280_payload)
    total = sum(mol.values())
    if total == 0:
        raise ValueError("Total corrected abundance must be non-zero")
    return {k: m / total * 100.0 for k, m in mol.items()}


# ==========================================================================
# 6. DAR by LC-MS  (spec §6; PDF p.17-18)
# ==========================================================================
def _weighted_mean(intensity_by_k: Mapping[int, float]) -> float:
    total = sum(intensity_by_k.values())
    if total == 0:
        raise ValueError("Total intensity must be non-zero")
    return sum(k * i for k, i in intensity_by_k.items()) / total


def _response_corrected_mean(
    intensity_by_k: Mapping[int, float],
    response_by_k: Mapping[int, float] | None = None,
) -> float:
    """Intensity-weighted mean drug load with an optional ionization-response
    correction.

    Raw LC-MS DAR assumes MS peak intensity is proportional to molar
    abundance.  In electrospray this fails: adding a hydrophobic payload
    changes a species' ionization efficiency, so highly loaded species are
    typically under-counted.  ``response_by_k`` gives the relative response
    factor r_k of each drug-load species (r_k = signal per mole, normalised
    to any convenient reference, e.g. r_0 = 1).  Molar abundance is recovered
    as I_k / r_k before the weighted mean.  With ``response_by_k=None`` (or a
    missing key, which defaults to r_k = 1) the result is identical to the
    uncorrected intensity-weighted mean.
    """
    if response_by_k is None:
        molar = dict(intensity_by_k)
    else:
        molar = {}
        for k, i in intensity_by_k.items():
            r = response_by_k.get(k, 1.0)
            if r <= 0:
                raise ValueError("response factor must be positive")
            molar[k] = i / r
    total = sum(molar.values())
    if total == 0:
        raise ValueError("Total corrected abundance must be non-zero")
    return sum(k * m for k, m in molar.items()) / total


def dar_lcms_reduced(
    light_chain: Mapping[int, float],
    heavy_chain: Mapping[int, float],
) -> float:
    """Denaturing + reducing DAR = 2*mean(LC) + 2*mean(HC)  (PDF p.17)."""
    return 2 * _weighted_mean(light_chain) + 2 * _weighted_mean(heavy_chain)


def dar_lcms_intact(intensity_by_k: Mapping[int, float]) -> float:
    """Denaturing non-reducing / native DAR = intensity-weighted mean (PDF p.17-18)."""
    return _weighted_mean(intensity_by_k)


def dar_lcms_intact_corrected(
    intensity_by_k: Mapping[int, float],
    response_by_k: Mapping[int, float] | None = None,
) -> float:
    """Intact LC-MS DAR with an ionization-response correction (A2).

    See ``_response_corrected_mean``.  Reduces exactly to
    ``dar_lcms_intact`` when no response factors are supplied.
    """
    return _response_corrected_mean(intensity_by_k, response_by_k)


def dar_lcms_reduced_corrected(
    light_chain: Mapping[int, float],
    heavy_chain: Mapping[int, float],
    light_response: Mapping[int, float] | None = None,
    heavy_response: Mapping[int, float] | None = None,
) -> float:
    """Reduced LC-MS DAR with per-chain ionization-response correction (A2).

    = 2*corrected_mean(LC) + 2*corrected_mean(HC).  Reduces exactly to
    ``dar_lcms_reduced`` when no response factors are supplied.
    """
    return (
        2 * _response_corrected_mean(light_chain, light_response)
        + 2 * _response_corrected_mean(heavy_chain, heavy_response)
    )


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
# 8b. In vitro dosing / serial dilution  (spec §8b)
# ==========================================================================
# Default per-well working (assay) volumes in microliters, keyed by plate
# format (number of wells). These are typical adherent-cell working volumes;
# the UI lets the user override the well volume for any format.
PLATE_WELL_VOLUME_UL: dict[str, float] = {
    "6": 2000.0, "12": 1000.0, "24": 500.0, "48": 250.0,
    "96": 100.0, "384": 50.0, "1536": 10.0,
}


def plate_well_volume(plate_type) -> float:
    """Default working volume (uL) for a plate format; raises on unknown format."""
    key = str(plate_type).strip()
    if key not in PLATE_WELL_VOLUME_UL:
        raise ValueError(f"unknown plate type {plate_type!r}")
    return PLATE_WELL_VOLUME_UL[key]


@dataclass(frozen=True)
class DoseTube:
    index: int                 # 1-based tube position in the series
    final_conc_uM: float       # target concentration in the well
    working_conc_uM: float     # concentration of the dosing solution in the tube
    source: str                # where the concentrated aliquot comes from
    transfer_in_uL: float      # volume drawn from source (stock or previous tube)
    diluent_uL: float          # assay medium added to this tube
    total_uL: float            # prepared volume in this tube
    dispense_to_wells_uL: float  # volume dispensed from this tube to its wells


@dataclass(frozen=True)
class DosingPlan:
    final_concs_uM: list[float]
    working_concs_uM: list[float]
    dose_mode: str
    dose_factor: float
    v_add_per_well_uL: float
    final_well_volume_uL: float
    v_use_per_tube_uL: float
    v_transfer_uL: float
    v_tube_total_uL: float
    tubes: list[DoseTube]
    total_stock_uL: float
    total_diluent_uL: float
    total_wells: int
    warnings: list[str]


def plan_serial_dilution(
    stock_uM: float,
    top_uM: float,
    fold: float,
    n_points: int,
    replicates: int,
    well_volume_uL: float,
    dose_mode: str = "spike",
    dose_factor: float = 10.0,
    dose_volume_uL: float | None = None,
    overage: float = 1.1,
    extra_dead_uL: float = 0.0,
) -> DosingPlan:
    """Design a serial dilution + plate-dosing scheme for an in vitro assay.

    From an ADC stock (`stock_uM`), an ideal top in-well concentration
    (`top_uM`), a constant `fold` dilution, a number of concentration points
    (`n_points`), and `replicates` per point, compute the per-tube pipetting
    recipe and the volume to add to each well.

    Two dosing modes:

    * ``"spike"`` — add a small spike of a concentrated dosing solution on top
      of the existing medium. The dosing solution is prepared at
      ``dose_factor x`` the final concentration, and the spike volume is
      ``well_volume / (dose_factor - 1)`` so that after mixing into the well
      the concentration is the intended final value. By default
      ``dose_factor = 10`` (a 1:10 spike).
    * ``"replace"`` — aspirate the medium and replace it with dosing solution
      already at the final concentration (``dose_factor`` forced to 1). The
      volume added equals the full well volume.

    Serial-dilution volumes are self-consistent: each tube is prepared to a
    uniform total so that the transfer drawn into the next tube plus the
    volume dispensed to the wells exactly balance the fold ratio. ``overage``
    (default 1.1 = 10% extra) and ``extra_dead_uL`` add margin for pipetting
    losses.

    Returns a :class:`DosingPlan`. All volumes are in microliters, all
    concentrations in micromolar.
    """
    if stock_uM <= 0:
        raise ValueError("stock concentration must be > 0")
    if top_uM <= 0:
        raise ValueError("top concentration must be > 0")
    if fold <= 1:
        raise ValueError("fold-dilution must be > 1")
    if n_points < 1:
        raise ValueError("need at least 1 concentration point")
    if replicates < 1:
        raise ValueError("need at least 1 replicate")
    if well_volume_uL <= 0:
        raise ValueError("well volume must be > 0")
    if dose_mode not in ("spike", "replace"):
        raise ValueError("dose_mode must be 'spike' or 'replace'")
    if dose_mode == "spike" and dose_factor <= 1:
        raise ValueError("dose_factor must be > 1 for spike mode")

    n = int(n_points)
    final_concs = [top_uM / (fold ** i) for i in range(n)]
    df = 1.0 if dose_mode == "replace" else float(dose_factor)
    working_concs = [c * df for c in final_concs]

    if dose_mode == "replace":
        v_add = well_volume_uL
        final_well_vol = well_volume_uL
    else:
        v_add = dose_volume_uL if dose_volume_uL else well_volume_uL / (df - 1.0)
        final_well_vol = well_volume_uL + v_add

    v_use = v_add * replicates * overage + extra_dead_uL
    v_xfer = v_use / (fold - 1.0)
    v_tube_total = v_use + v_xfer

    stock_vol_tube1 = v_tube_total * working_concs[0] / stock_uM
    tubes: list[DoseTube] = []
    for i in range(n):
        if i == 0:
            source = "ADC stock"
            transfer_in = stock_vol_tube1
            diluent = v_tube_total - stock_vol_tube1
        else:
            source = f"tube {i} (transfer)"
            transfer_in = v_xfer
            diluent = v_use
        tubes.append(DoseTube(
            index=i + 1,
            final_conc_uM=final_concs[i],
            working_conc_uM=working_concs[i],
            source=source,
            transfer_in_uL=transfer_in,
            diluent_uL=diluent,
            total_uL=transfer_in + diluent,
            dispense_to_wells_uL=v_use,
        ))

    total_diluent = sum(t.diluent_uL for t in tubes)
    warnings: list[str] = []
    if working_concs[0] > stock_uM:
        warnings.append(
            f"Top working conc {working_concs[0]:.4g} uM exceeds stock "
            f"{stock_uM:.4g} uM - stock too dilute for this design."
        )
    if stock_vol_tube1 < 1.0:
        warnings.append(
            f"Tube-1 stock volume {stock_vol_tube1:.3g} uL < 1 uL - hard to "
            f"pipette accurately; use an intermediate dilution or a larger "
            f"prep volume."
        )

    return DosingPlan(
        final_concs_uM=final_concs,
        working_concs_uM=working_concs,
        dose_mode=dose_mode,
        dose_factor=df,
        v_add_per_well_uL=v_add,
        final_well_volume_uL=final_well_vol,
        v_use_per_tube_uL=v_use,
        v_transfer_uL=v_xfer,
        v_tube_total_uL=v_tube_total,
        tubes=tubes,
        total_stock_uL=stock_vol_tube1,
        total_diluent_uL=total_diluent,
        total_wells=n * replicates,
        warnings=warnings,
    )


# ==========================================================================
# 8c. Plate map  (spec §8c)
# ==========================================================================
# Physical row x column geometry for each supported plate format.
PLATE_GEOMETRY: dict[str, tuple[int, int]] = {
    "6": (2, 3), "12": (3, 4), "24": (4, 6), "48": (6, 8),
    "96": (8, 12), "384": (16, 24), "1536": (32, 48),
}


def _row_label(i: int) -> str:
    """Excel-style row label: 0->A, 25->Z, 26->AA, ..."""
    s = ""
    i += 1
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


@dataclass(frozen=True)
class PlateWell:
    row: int                 # 0-based grid row
    col: int                 # 0-based grid column
    row_label: str | None    # "A".."AF" (None if off-plate)
    col_label: int           # 1-based column number
    point: int               # 0-based concentration-point index
    replicate: int           # 1-based replicate number
    final_conc_uM: float
    in_bounds: bool          # False if the assignment falls outside the grid


@dataclass(frozen=True)
class PlateMap:
    plate_type: str
    rows: int
    cols: int
    orientation: str         # actual layout used
    capacity: int
    used: int
    wells: list[PlateWell]
    warnings: list[str]


def plate_map(
    final_concs_uM,
    replicates: int,
    plate_type,
    orientation: str = "by_column",
) -> PlateMap:
    """Assign a dose series (concentration points x replicates) to plate wells.

    ``orientation``:

    * ``"by_column"`` (default) — each concentration point occupies one column
      (dose increases left to right), replicates fill the rows of that column.
    * ``"by_row"`` — each point occupies one row, replicates fill columns.
    * ``"sequential"`` — row-major fill, replicates grouped consecutively.

    If the requested grid layout does not fit the plate (too many points for
    the columns, or too many replicates for the rows), the function falls back
    to sequential fill and records a warning. A second warning is added when
    the total well count exceeds the plate capacity.
    """
    key = str(plate_type).strip()
    if key not in PLATE_GEOMETRY:
        raise ValueError(f"unknown plate type {plate_type!r}")
    if replicates < 1:
        raise ValueError("need at least 1 replicate")
    concs = list(final_concs_uM)
    n = len(concs)
    if n < 1:
        raise ValueError("need at least 1 concentration point")

    rows, cols = PLATE_GEOMETRY[key]
    capacity = rows * cols
    used = n * replicates
    warnings: list[str] = []

    if orientation == "by_column" and n <= cols and replicates <= rows:
        used_orient = "by_column"
    elif orientation == "by_row" and n <= rows and replicates <= cols:
        used_orient = "by_row"
    else:
        used_orient = "sequential"
        if orientation in ("by_column", "by_row"):
            warnings.append(
                f"{n} points x {replicates} replicates do not fit the "
                f"{rows}x{cols} grid in {orientation} layout; using "
                f"sequential fill."
            )

    placed: list[tuple[int, int, int, int]] = []  # (row, col, point, rep0)
    if used_orient == "by_column":
        for i in range(n):
            for r in range(replicates):
                placed.append((r, i, i, r))
    elif used_orient == "by_row":
        for i in range(n):
            for r in range(replicates):
                placed.append((i, r, i, r))
    else:
        pos = 0
        for i in range(n):
            for r in range(replicates):
                placed.append((pos // cols, pos % cols, i, r))
                pos += 1

    if used > capacity:
        warnings.append(
            f"{used} wells needed exceeds {key}-well capacity ({capacity})."
        )

    wells: list[PlateWell] = []
    for (rr, cc, pi, rep0) in placed:
        inside = rr < rows and cc < cols
        wells.append(PlateWell(
            row=rr, col=cc,
            row_label=_row_label(rr) if inside else None,
            col_label=cc + 1,
            point=pi, replicate=rep0 + 1,
            final_conc_uM=concs[pi], in_bounds=inside,
        ))

    return PlateMap(
        plate_type=key, rows=rows, cols=cols, orientation=used_orient,
        capacity=capacity, used=used, wells=wells, warnings=warnings,
    )


def plate_map_to_csv(pmap: PlateMap) -> str:
    """Deterministic CSV of a plate map: one row per assigned well."""
    lines = ["well,row,column,point,replicate,final_conc_uM"]
    for w in pmap.wells:
        well = f"{w.row_label}{w.col_label}" if w.row_label else "(off-plate)"
        lines.append(
            f"{well},{w.row_label or ''},{w.col_label},{w.point + 1},"
            f"{w.replicate},{w.final_conc_uM:.6g}"
        )
    return "\n".join(lines)


# ==========================================================================
# 8e. Concentration units + selection-driven dosing  (spec §8e)
# ==========================================================================
# Multiplicative factor to convert a value in the given unit to micromolar.
CONC_UNIT_TO_UM: dict[str, float] = {"uM": 1.0, "nM": 1e-3, "pM": 1e-6}


def convert_concentration(value: float, from_unit: str, to_unit: str) -> float:
    """Convert a concentration between uM / nM / pM."""
    if from_unit not in CONC_UNIT_TO_UM:
        raise ValueError(f"unknown unit {from_unit!r}")
    if to_unit not in CONC_UNIT_TO_UM:
        raise ValueError(f"unknown unit {to_unit!r}")
    return value * CONC_UNIT_TO_UM[from_unit] / CONC_UNIT_TO_UM[to_unit]


def series_shape_from_selection(n_rows: int, n_cols: int,
                                orientation: str = "by_column") -> tuple[int, int]:
    """Geometry contract: a painted ``n_rows`` x ``n_cols`` rectangle becomes a
    ``(n_points, replicates)`` design.

    * ``by_column`` (default) — each column is a concentration point (dose
      decreases left to right), rows are replicates.
    * ``by_row`` — each row is a concentration point, columns are replicates.
    """
    if n_rows < 1 or n_cols < 1:
        raise ValueError("selection must be at least 1x1")
    if orientation == "by_row":
        return n_rows, n_cols
    return n_cols, n_rows


@dataclass(frozen=True)
class SelectionWell:
    row: int                 # 0-based grid row (absolute plate coordinate)
    col: int                 # 0-based grid column
    point: int               # 0-based concentration-point index within the group
    replicate: int           # 1-based replicate number
    final_conc_uM: float | None   # None if the well lies past the series length


@dataclass(frozen=True)
class SelectionAssignment:
    wells: list[SelectionWell]
    rows: int                # bounding-box height of the selection
    cols: int                # bounding-box width
    origin: tuple[int, int]  # (row, col) of the top-left selected well
    rectangular: bool        # True if the selection fills its bounding box
    orientation: str
    warnings: list[str]


def assign_selection(cells, final_concs_uM,
                     orientation: str = "by_column") -> SelectionAssignment:
    """Map a group's concentration series onto its selected wells.

    ``cells`` is an iterable of ``(row, col)`` 0-based plate coordinates (any
    offset on the plate). Concentration point and replicate index are read from
    each well's position within the selection's bounding box, per
    ``orientation``. Wells whose point index falls past the series length get a
    ``final_conc_uM`` of ``None``; a warning is recorded when the series length
    does not match the selection's span along the dose axis.
    """
    cleaned = sorted({(int(r), int(c)) for r, c in cells})
    if not cleaned:
        raise ValueError("no cells selected")
    concs = list(final_concs_uM)
    n = len(concs)

    rows = [r for r, _ in cleaned]
    cols = [c for _, c in cleaned]
    r0, c0 = min(rows), min(cols)
    h = max(rows) - r0 + 1
    w = max(cols) - c0 + 1
    rectangular = (len(cleaned) == h * w)

    dose_span = h if orientation == "by_row" else w
    warnings: list[str] = []
    if n != dose_span:
        warnings.append(
            f"{n} concentration points but selection spans {dose_span} "
            f"well(s) along the dose axis."
        )

    wells: list[SelectionWell] = []
    for (r, c) in cleaned:
        if orientation == "by_row":
            pt = r - r0
            rep = c - c0 + 1
        else:
            pt = c - c0
            rep = r - r0 + 1
        conc = concs[pt] if 0 <= pt < n else None
        wells.append(SelectionWell(row=r, col=c, point=pt,
                                   replicate=rep, final_conc_uM=conc))

    return SelectionAssignment(
        wells=wells, rows=h, cols=w, origin=(r0, c0),
        rectangular=rectangular, orientation=orientation, warnings=warnings,
    )


def aggregate_dosing_plans(plans) -> dict:
    """Sum reagent totals across several group :class:`DosingPlan` objects so a
    multi-group plate reports one combined stock / medium / well tally.

    (Stock volumes are summed as micromolar-equivalent volumes only when the
    groups share a stock; callers that mix stocks should read per-group totals.)
    """
    plans = list(plans)
    return {
        "n_groups": len(plans),
        "total_stock_uL": sum(p.total_stock_uL for p in plans),
        "total_diluent_uL": sum(p.total_diluent_uL for p in plans),
        "total_wells": sum(p.total_wells for p in plans),
    }


# ==========================================================================
# 9b. Extinction-coefficient determination  (spec §9b)
# ==========================================================================
# Ordinary least-squares fit of a Beer-Lambert dilution series. For a set of
# known concentrations c (mol/L) and measured absorbances A at a fixed
# wavelength and path length L (cm), Beer-Lambert gives A = eps * L * c, so a
# straight-line fit A vs c has slope = eps * L, hence eps = slope / L.
def linear_regression(xs, ys) -> dict:
    """Ordinary least-squares line y = slope*x + intercept.

    Returns slope, intercept, r_squared and n. Requires >= 2 points with at
    least two distinct x values.
    """
    xs = [float(x) for x in xs]
    ys = [float(y) for y in ys]
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs and ys must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        raise ValueError("need at least two distinct x values")
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    if ss_tot == 0.0:
        r_squared = 1.0
    else:
        ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
        r_squared = 1.0 - ss_res / ss_tot
    return {"slope": slope, "intercept": intercept, "r_squared": r_squared, "n": n}


def extinction_coefficient(
    concentrations_M, absorbances, path_length_cm: float = 1.0
) -> dict:
    """Molar extinction coefficient (M^-1 cm^-1) from a dilution series.

    Fits absorbance vs concentration by least squares; eps = slope / path length.
    Returns eps plus the underlying slope, intercept and r_squared.
    """
    if path_length_cm <= 0:
        raise ValueError("path length must be positive")
    fit = linear_regression(concentrations_M, absorbances)
    return {
        "eps": fit["slope"] / path_length_cm,
        "slope": fit["slope"],
        "intercept": fit["intercept"],
        "r_squared": fit["r_squared"],
        "n": fit["n"],
    }


# ==========================================================================
# 9b. Uncertainty / error propagation  (spec §2, §5, §6)
# ==========================================================================
def dar_uv_uncertainty(
    a280: float,
    a_lmax: float,
    eps280_mab: float,
    eps_lmax_mab: float,
    eps280_lp: float,
    eps_lmax_lp: float,
    sigma_a280: float = 0.0,
    sigma_a_lmax: float = 0.0,
    sigma_eps280_mab: float = 0.0,
    sigma_eps_lmax_mab: float = 0.0,
    sigma_eps280_lp: float = 0.0,
    sigma_eps_lmax_lp: float = 0.0,
) -> dict:
    """UV/Vis DAR with a 1-sigma error propagated from every input quantity.

    DAR = (eps_lmax_mab - R*eps280_mab) / (R*eps280_lp - eps_lmax_lp),
    R = a_lmax / a280. First-order (Gaussian) propagation sums independent
    contributions in quadrature:

      * absorbance reads enter only through R, so
        sigma_R  = |R| * sqrt((sigma_a280/a280)^2 + (sigma_a_lmax/a_lmax)^2)
        and their contribution is |dDAR/dR| * sigma_R;
      * the four extinction coefficients contribute via their analytic
        partials of DAR (denom = R*eps280_lp - eps_lmax_lp, numer = the DAR
        numerator):
            dDAR/d(eps_lmax_mab) =  1/denom
            dDAR/d(eps280_mab)   = -R/denom
            dDAR/d(eps280_lp)    = -numer*R/denom^2
            dDAR/d(eps_lmax_lp)  =  numer/denom^2

    A symmetric 95% interval uses 1.96*sigma_DAR. Zero input sigmas give
    sigma_DAR = 0 and a degenerate interval equal to the point estimate, so
    the extra terms are fully backward-compatible with the read-only form.
    """
    if a280 == 0:
        raise ValueError("A280 must be non-zero")
    if a_lmax == 0 and sigma_a_lmax != 0.0:
        raise ValueError("A_lmax must be non-zero to propagate its error")
    r = a_lmax / a280
    denom = r * eps280_lp - eps_lmax_lp
    if denom == 0:
        raise ZeroDivisionError("Degenerate optics: R*eps280_LP == eps_lmax_LP")
    dar = (eps_lmax_mab - r * eps280_mab) / denom
    # dDAR/dR via quotient rule on (eps_lmax_mab - R*eps280_mab)/(R*eps280_lp - eps_lmax_lp)
    numer = eps_lmax_mab - r * eps280_mab
    d_dar_d_r = (-eps280_mab * denom - numer * eps280_lp) / denom ** 2
    rel_var = 0.0
    if sigma_a280:
        rel_var += (sigma_a280 / a280) ** 2
    if sigma_a_lmax:
        rel_var += (sigma_a_lmax / a_lmax) ** 2
    sigma_r = abs(r) * (rel_var ** 0.5)
    var = (d_dar_d_r * sigma_r) ** 2
    # extinction-coefficient partials, added in quadrature
    if sigma_eps_lmax_mab:
        var += ((1.0 / denom) * sigma_eps_lmax_mab) ** 2
    if sigma_eps280_mab:
        var += ((-r / denom) * sigma_eps280_mab) ** 2
    if sigma_eps280_lp:
        var += ((-numer * r / denom ** 2) * sigma_eps280_lp) ** 2
    if sigma_eps_lmax_lp:
        var += ((numer / denom ** 2) * sigma_eps_lmax_lp) ** 2
    sigma_dar = var ** 0.5
    return {
        "dar": dar,
        "sigma_dar": sigma_dar,
        "ci95_low": dar - 1.96 * sigma_dar,
        "ci95_high": dar + 1.96 * sigma_dar,
        "sigma_r": sigma_r,
        "rel_sigma": (sigma_dar / dar) if dar != 0 else None,
    }


def distribution_dispersion(weights: Mapping[int, float]) -> dict:
    """Population mean, variance and SD of a weighted drug-load distribution.

    `weights` maps drug count k -> abundance (HIC peak area or LC-MS intensity).
    The SD is the *heterogeneity* of the sample: the spread of drug load across
    species, not a measurement error. Population moments (divide by total
    weight, not weight-1) since the peaks are the whole measured population.
    """
    total = sum(weights.values())
    if total == 0:
        raise ValueError("Total weight must be non-zero")
    mean = sum(k * w for k, w in weights.items()) / total
    variance = sum(w * (k - mean) ** 2 for k, w in weights.items()) / total
    return {"mean": mean, "variance": variance, "sd": variance ** 0.5}


def predict_dar_distribution(
    n_sites: int,
    p_site: float | None = None,
    feed_ratio: float | None = None,
    efficiency: float = 1.0,
    drugs_per_site: int = 2,
) -> dict:
    """Predict the full DAR distribution from a site-occupancy model, keyed to
    the conjugation chemistry.

    Real conjugation chemistries do not populate every integer DAR equally.
    The controlling variable is how many drugs each independent *site* carries
    once it reacts (``drugs_per_site``):

    * **Cysteine / TCEP interchain conjugation (drugs_per_site=2, default).**
      An IgG1 has four interchain disulfides; reducing each one exposes a
      *pair* of thiols that are conjugated together. A site is therefore
      either unreacted or contributes **two** drugs, so the achievable DAR
      ladder is 0, 2, 4, 6, 8 — the even values seen for real thiol-linked
      ADCs (e.g. brentuximab vedotin). Use ``n_sites=4`` for a standard IgG1.
    * **Stochastic amine / lysine conjugation (drugs_per_site=1).**
      Each of many surface lysines reacts independently and carries a single
      drug, giving the classic smooth binomial over 0, 1, 2, ... drugs.

    Each of ``n_sites`` sites reacts independently with probability ``p_site``;
    a reacted site adds ``drugs_per_site`` drugs. The occupied-site count is
    Binomial(n_sites, p_site), so with ``d = drugs_per_site``:

        P(DAR = d*j) = C(n, j) p^j (1-p)^(n-j),  j = 0..n
        mean DAR     = d * n * p
        variance     = d^2 * n * p * (1-p)

    Two ways to set the occupancy probability:

    * ``p_site`` directly (0..1), or
    * ``feed_ratio`` (drug:mAb molar feed) times ``efficiency`` (fraction of
      offered drug that conjugates), spread over the total drug-carrying
      capacity (``n_sites * drugs_per_site``):
      ``p_site = feed_ratio * efficiency / (n_sites * drugs_per_site)``.
      With this definition mean DAR = feed_ratio * efficiency exactly.

    Returns the per-species probabilities (keyed by DAR, i.e. multiples of
    ``drugs_per_site``) plus mean/variance/SD. The SD is the intrinsic
    drug-load *heterogeneity* the model predicts, not a measurement error.
    """
    if n_sites < 1:
        raise ValueError("n_sites must be >= 1")
    if drugs_per_site < 1:
        raise ValueError("drugs_per_site must be >= 1")
    max_dar = n_sites * drugs_per_site
    if p_site is None:
        if feed_ratio is None:
            raise ValueError("provide either p_site or feed_ratio")
        if feed_ratio < 0 or efficiency < 0:
            raise ValueError("feed_ratio and efficiency must be >= 0")
        p_site = (feed_ratio * efficiency) / max_dar
    if not (0.0 <= p_site <= 1.0):
        raise ValueError(
            "implied per-site probability outside [0,1]; "
            "feed_ratio*efficiency exceeds the drug-carrying capacity "
            "(n_sites * drugs_per_site)"
        )
    dist = {
        drugs_per_site * j: comb(n_sites, j) * p_site ** j
        * (1.0 - p_site) ** (n_sites - j)
        for j in range(n_sites + 1)
    }
    disp = distribution_dispersion(dist)
    return {
        "p_site": p_site,
        "distribution": dist,
        "mean_dar": disp["mean"],
        "variance": disp["variance"],
        "sd": disp["sd"],
        "drugs_per_site": drugs_per_site,
        "max_dar": max_dar,
    }


def dar_lcms_reduced_uncertainty(
    light_chain: Mapping[int, float],
    heavy_chain: Mapping[int, float],
) -> dict:
    """Reduced LC-MS DAR with an SD from per-chain load heterogeneity.

    DAR = 2*mean(LC) + 2*mean(HC). Treating the two light chains and two heavy
    chains as independent draws from their measured load distributions, the
    variance of the sum is 2*var(LC) + 2*var(HC); sigma_dar is its square root.
    Also returns each chain's SD so an asymmetric contributor is visible.
    """
    lc = distribution_dispersion(light_chain)
    hc = distribution_dispersion(heavy_chain)
    dar = 2 * lc["mean"] + 2 * hc["mean"]
    variance = 2 * lc["variance"] + 2 * hc["variance"]
    sigma_dar = variance ** 0.5
    return {
        "dar": dar,
        "sigma_dar": sigma_dar,
        "light_sd": lc["sd"],
        "heavy_sd": hc["sd"],
        "ci95_low": dar - 1.96 * sigma_dar,
        "ci95_high": dar + 1.96 * sigma_dar,
    }


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
