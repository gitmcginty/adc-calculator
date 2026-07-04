"""Regression tests for adc_core against the golden values in adc_spec.md §9.

All golden values are reproduced from Calculation_DAR_conjugation.xlsx and the
protocols PDF. Run:  pytest -q test_adc_core.py
"""
import math
import pytest

import adc_core as core

REL = 1e-6

# Spreadsheet inputs (spec §9)
MW_MAB = 145000.0
EPS280_MAB_MASS = 1.4
EPS280_MAB = EPS280_MAB_MASS * MW_MAB      # 203000
A280 = 7.4315
A_LMAX = 0.7827
MW_LP = 1559.62
EPS280_LP = 7287.0
EPS_LMAX_LP = 9624.0
R = A_LMAX / A280


def approx(x):
    return pytest.approx(x, rel=REL)


# ---- Extinction-coefficient conversions -------------------------------------
def test_eps_conversion_roundtrip():
    assert core.eps_mass_to_molar(1.4, MW_MAB) == approx(203000.0)
    assert core.eps_molar_to_mass(203000.0, MW_MAB) == approx(1.4)


# ---- DAR by UV/Vis ----------------------------------------------------------
def test_dar_uv_case_A_circular_is_zero():
    """Spreadsheet's circular eps_lmax_mab = R*eps280_mab forces DAR = 0."""
    eps_lmax_mab = R * EPS280_MAB
    dar = core.dar_uv(A280, A_LMAX, EPS280_MAB, eps_lmax_mab, EPS280_LP, EPS_LMAX_LP)
    assert dar == pytest.approx(0.0, abs=1e-9)


def test_dar_uv_case_B_realistic():
    """Bare mAb, eps_lmax_mab = 0 -> DAR 2.41408 (spec §9)."""
    dar = core.dar_uv(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP)
    assert dar == approx(2.4140809554780795)


# ---- ADC properties ---------------------------------------------------------
def test_adc_properties_case_A():
    p = core.adc_properties(0.0, MW_MAB, MW_LP, EPS280_MAB, EPS280_LP,
                            a280=A280, dilution_factor=1.0)
    assert p.mw_adc == approx(145000.0)
    assert p.eps280_adc_molar == approx(203000.0)
    assert p.eps280_adc_mass == approx(1.4)
    assert p.conc_mgml == approx(5.3082142857142856)
    assert p.conc_uM == approx(36.608374384236456)


def test_adc_properties_case_B():
    dar = 2.4140809554780795
    p = core.adc_properties(dar, MW_MAB, MW_LP, EPS280_MAB, EPS280_LP,
                            a280=A280, dilution_factor=1.0)
    assert p.mw_adc == approx(148765.04893978272)
    assert p.eps280_adc_molar == approx(220591.40792256876)
    assert p.eps280_adc_mass == approx(1.482817432553395)
    assert p.conc_mgml == approx(5.011743075614535)
    assert p.conc_uM == approx(33.68898213210816)


def test_conc_unknown_payload_eps():
    # payload 280 contribution subtracted, naked-mAb mass eps applied
    c = core.conc_unknown_payload_eps(A280, A_LMAX, EPS280_LP, EPS_LMAX_LP,
                                      EPS280_MAB_MASS)
    a_lp_280 = A_LMAX * EPS280_LP / EPS_LMAX_LP
    assert c == approx((A280 - a_lp_280) / EPS280_MAB_MASS)


# ---- Ellman -----------------------------------------------------------------
def test_ellman():
    assert core.ellman_free_thiols(0.5, 0.6, 203000.0) == approx(11.955241460541815)


# ---- HIC --------------------------------------------------------------------
def test_dar_hic():
    assert core.dar_hic({0: 5, 1: 15, 2: 45, 3: 25, 4: 10}) == approx(2.2)


def test_dar_species_fractions_sum_to_100():
    frac = core.dar_species_fractions({0: 5, 1: 15, 2: 45, 3: 25, 4: 10})
    assert sum(frac.values()) == approx(100.0)
    assert frac[2] == approx(45.0)


# ---- LC-MS ------------------------------------------------------------------
def test_dar_lcms_reduced():
    dar = core.dar_lcms_reduced({0: 10, 1: 90}, {0: 5, 1: 20, 2: 75})
    assert dar == approx(5.2)


def test_dar_lcms_intact():
    # non-reducing/native intensity-weighted mean
    assert core.dar_lcms_intact({0: 10, 2: 30, 4: 60}) == approx((0*10+2*30+4*60)/100)


def test_dar_lcms_reduced_breakdown():
    lc = {0: 10, 1: 90}
    hc = {0: 5, 1: 20, 2: 75}
    bd = core.dar_lcms_reduced_breakdown(lc, hc)
    # light chain intermediates
    assert bd["light"]["sum_i"] == approx(100)
    assert bd["light"]["sum_ki"] == approx(90)
    assert bd["light"]["mean"] == approx(0.9)
    assert bd["light"]["contribution"] == approx(1.8)
    # heavy chain intermediates
    assert bd["heavy"]["sum_i"] == approx(100)
    assert bd["heavy"]["sum_ki"] == approx(170)
    assert bd["heavy"]["mean"] == approx(1.7)
    assert bd["heavy"]["contribution"] == approx(3.4)
    # total equals the scalar function and the sum of contributions
    assert bd["dar"] == approx(5.2)
    assert bd["dar"] == approx(core.dar_lcms_reduced(lc, hc))
    assert bd["dar"] == approx(bd["light"]["contribution"] + bd["heavy"]["contribution"])


def test_dar_lcms_reduced_breakdown_real_case():
    # The worked example from the discrepancy investigation.
    lc = {0: 19186.7, 1: 6824.9, 2: 494.7}
    hc = {0: 7732.4, 1: 5823.1, 2: 6875.4, 3: 3895.1, 4: 901.7, 5: 85.7}
    bd = core.dar_lcms_reduced_breakdown(lc, hc)
    lc_mean = (0*19186.7 + 1*6824.9 + 2*494.7) / (19186.7 + 6824.9 + 494.7)
    hc_mean = (0*7732.4 + 1*5823.1 + 2*6875.4 + 3*3895.1 + 4*901.7 + 5*85.7) \
              / (7732.4 + 5823.1 + 6875.4 + 3895.1 + 901.7 + 85.7)
    assert bd["light"]["mean"] == approx(lc_mean)
    assert bd["heavy"]["mean"] == approx(hc_mean)
    assert bd["dar"] == approx(2 * lc_mean + 2 * hc_mean)


# ---- Conjugation designer ---------------------------------------------------
def test_moles_and_volume():
    # 10 mg IgG1 -> 10/1000/145000 mol
    n = core.moles_from_mass(10.0, MW_MAB)
    assert n == approx(10.0/1000/145000)
    # 6 equiv TCEP from a 5 mM stock
    moles = 6 * n
    v = core.reagent_volume_uL(moles, 5.0)
    assert v == approx(moles / (5.0e-3) * 1e6)


def test_plan_conjugation_thiomab():
    plan = core.plan_conjugation(
        "thiomab", mass_mab_mg=10.0, mw_mab=MW_MAB, mab_conc_g_per_L=5.0,
        reagents=[
            {"name": "TCEP", "equiv": 10.0, "basis": "mAb", "stock_mM": 5.0},
            {"name": "maleimide-LP", "equiv": 1.5, "basis": "thiol", "stock_mM": 10.0},
        ],
    )
    # Thiomab full reduction -> 8 free thiols per mAb
    assert plan.free_thiols_per_mab == 8
    n_mab = 10.0/1000/145000
    tcep = next(a for a in plan.additions if a.name == "TCEP")
    assert tcep.moles == approx(10.0 * n_mab)
    mal = next(a for a in plan.additions if a.name == "maleimide-LP")
    assert mal.moles == approx(1.5 * n_mab * 8)      # per thiol
    # base volume: 10 mg / 5 g/L = 2 mL = 2000 uL
    assert plan.reaction_volume_uL >= 2000.0


# ---- Yield & formulation ----------------------------------------------------
def test_yield_and_formulation():
    y = core.yield_and_formulation(conc_mgml=5.0, volume_mL=2.0, mw_adc=148765.0,
                                   starting_mass_mg=15.0, target_conc_mgml=10.0)
    assert y.recovered_mass_mg == approx(10.0)
    assert y.yield_pct == approx(10.0/15.0*100)
    assert y.molar_amount_nmol == approx(10.0/148765.0*1e6)
    assert y.v_final_mL == approx(1.0)          # concentrate 2 mL -> 1 mL
    assert y.v_change_mL == approx(-1.0)


# ---- Dye library ------------------------------------------------------------
def test_dye_library_integrity():
    assert len(core.DYE_LIBRARY) == 22
    for d in core.DYE_LIBRARY:
        # e280 must equal e_lmax * cf280 (spreadsheet definition)
        assert d["e280"] == approx(d["e_lmax"] * d["cf280"])
    af647 = core.get_dye("af647")
    assert af647 is not None and af647["lmax"] == 647


# ---- Registry ---------------------------------------------------------------
def test_make_registry_record_normalizes():
    r = core.make_registry_record(
        {"id": "ADC-001", "assay": "uv", "dar": "2.414", "conc_mgml": "", "notes": "  lot A "})
    assert r["id"] == "ADC-001"
    assert r["dar"] == approx(2.414)      # coerced from string
    assert r["conc_mgml"] is None          # blank -> None
    assert r["notes"] == "lot A"           # stripped
    # every canonical field present
    assert set(r.keys()) == set(core.REGISTRY_FIELDS)


def test_make_registry_record_requires_id():
    with pytest.raises(ValueError):
        core.make_registry_record({"assay": "uv", "dar": 2.0})


def test_registry_to_csv_escapes_and_orders():
    recs = [core.make_registry_record({"id": "ADC-001", "assay": "uv", "dar": 2.2,
                                       "notes": "AF647, lot A"})]
    csv = core.registry_to_csv(recs)
    lines = csv.split("\n")
    assert lines[0] == ",".join(core.REGISTRY_FIELDS)
    # comma inside notes forces quoting; integer-valued DAR prints as "2.2"
    assert '"AF647, lot A"' in lines[1]
    assert lines[1].startswith("ADC-001,,,uv,,2.2,")


def test_summarize_registry_stats():
    recs = [core.make_registry_record(r) for r in [
        {"id": "ADC-001", "assay": "uv", "dar": 2.414},
        {"id": "ADC-001", "assay": "hic", "dar": 2.2},
        {"id": "ADC-002", "assay": "lcms", "dar": 5.2},
        {"id": "ADC-003", "assay": "conjugation"},   # no DAR
    ]]
    s = core.summarize_registry(recs)
    assert s["n_total"] == 4
    assert s["n_ids"] == 3
    assert s["dar_n"] == 3
    assert s["dar_mean"] == approx(3.271333333333333)
    assert s["dar_sd"] == approx(1.6736981010126448)
    assert s["dar_min"] == approx(2.2) and s["dar_max"] == approx(5.2)
    assert s["by_assay"]["uv"] == 1 and s["by_id"]["ADC-001"] == 2


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
