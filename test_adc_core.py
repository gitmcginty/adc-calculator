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


def test_ellman_blank_subtraction_lowers_ratio():
    raw = core.ellman_free_thiols(0.5, 0.6, 203000.0)
    corrected = core.ellman_free_thiols(0.5, 0.6, 203000.0, 14150.0, 0.02, 0.01)
    assert corrected == approx(11.671557764867941)
    # blank background biases both channels; net effect lowers the ratio here
    assert corrected < raw


def test_ellman_zero_blank_reduces_to_raw():
    assert core.ellman_free_thiols(0.5, 0.6, 203000.0, 14150.0, 0.0, 0.0) == approx(
        core.ellman_free_thiols(0.5, 0.6, 203000.0)
    )


def test_ellman_rejects_nonpositive_net_a280():
    with pytest.raises(ValueError):
        core.ellman_free_thiols(0.5, 0.6, 203000.0, 14150.0, 0.0, 0.6)


# ---- Scattering / turbidity correction --------------------------------------
def test_scatter_absorbance_rayleigh():
    assert core.scatter_absorbance(0.05, 320.0, 280.0, 4.0) == approx(0.08529779258642231)
    # shorter target wavelength scatters more strongly than a longer one
    assert core.scatter_absorbance(0.05, 320.0, 280.0) > core.scatter_absorbance(0.05, 320.0, 495.0)


def test_scatter_corrected_absorbance():
    r = core.scatter_corrected_absorbance(1.0, 0.5, 0.05, 320.0, 495.0, 4.0)
    assert r["a280_corrected"] == approx(0.9147022074135777)
    assert r["a_lmax_corrected"] == approx(0.49126728831613614)
    assert r["scatter_280"] == approx(0.08529779258642231)
    assert r["scatter_lmax"] == approx(0.008732711683863857)


def test_scatter_correction_zero_reference_is_identity():
    r = core.scatter_corrected_absorbance(1.0, 0.5, 0.0)
    assert r["a280_corrected"] == approx(1.0)
    assert r["a_lmax_corrected"] == approx(0.5)


def test_scatter_corrected_absorbance_floors_at_zero():
    # a huge scatter reference drives the correction negative -> floored to 0
    r = core.scatter_corrected_absorbance(0.1, 0.05, 1.0, 320.0, 495.0, 4.0)
    assert r["a280_corrected"] == 0.0
    assert r["a_lmax_corrected"] == 0.0


def test_scatter_absorbance_rejects_nonpositive_wavelength():
    with pytest.raises(ValueError):
        core.scatter_absorbance(0.05, 0.0, 280.0)
    with pytest.raises(ValueError):
        core.scatter_absorbance(0.05, 320.0, -280.0)


# ---- HIC --------------------------------------------------------------------
def test_dar_hic():
    assert core.dar_hic({0: 5, 1: 15, 2: 45, 3: 25, 4: 10}) == approx(2.2)


def test_dar_species_fractions_sum_to_100():
    frac = core.dar_species_fractions({0: 5, 1: 15, 2: 45, 3: 25, 4: 10})
    assert sum(frac.values()) == approx(100.0)
    assert frac[2] == approx(45.0)


def test_dar_hic_corrected_lowers_dar():
    # Payload absorbs at 280 nm, so uncorrected area% over-weights high-DAR
    # peaks. Correcting by per-species eps280 must lower the reported DAR.
    auc = {0: 10.0, 2: 50.0, 4: 30.0, 6: 10.0}
    uncorr = core.dar_hic(auc)
    corr = core.dar_hic_corrected(auc, 203000.0, 5000.0)
    assert uncorr == approx(2.8)
    assert corr == approx(2.741596270777426)
    assert corr < uncorr


def test_dar_hic_corrected_reduces_to_uncorrected_when_payload_transparent():
    # eps280_payload == 0 -> every species shares eps280_mab -> identical to
    # the plain area-weighted mean.
    auc = {0: 5.0, 1: 15.0, 2: 45.0, 3: 25.0, 4: 10.0}
    assert core.dar_hic_corrected(auc, 203000.0, 0.0) == approx(core.dar_hic(auc))


def test_dar_species_fractions_corrected_sum_to_100():
    auc = {0: 10.0, 2: 50.0, 4: 30.0, 6: 10.0}
    frac = core.dar_species_fractions_corrected(auc, 203000.0, 5000.0)
    assert sum(frac.values()) == approx(100.0)
    assert frac[2] == approx(50.870418158189466)


def test_dar_hic_corrected_rejects_nonpositive_eps():
    with pytest.raises(ValueError):
        core.dar_hic_corrected({0: 1.0, 2: 1.0}, 0.0, 0.0)


# ---- LC-MS ------------------------------------------------------------------
def test_dar_lcms_reduced():
    dar = core.dar_lcms_reduced({0: 10, 1: 90}, {0: 5, 1: 20, 2: 75})
    assert dar == approx(5.2)


def test_dar_lcms_intact():
    # non-reducing/native intensity-weighted mean
    assert core.dar_lcms_intact({0: 10, 2: 30, 4: 60}) == approx((0*10+2*30+4*60)/100)


def test_dar_lcms_intact_corrected_reduces_to_raw():
    intensity = {0: 10, 2: 30, 4: 60}
    assert core.dar_lcms_intact_corrected(intensity) == approx(core.dar_lcms_intact(intensity))
    assert core.dar_lcms_intact_corrected(intensity, None) == approx(3.0)


def test_dar_lcms_intact_corrected_raises_dar_when_high_load_underionizes():
    intensity = {0: 10, 2: 30, 4: 60}
    response = {0: 1.0, 2: 0.8, 4: 0.6}  # higher drug load ionizes worse
    corrected = core.dar_lcms_intact_corrected(intensity, response)
    assert corrected == approx(3.2203389830508473)
    # under-counted high-load species are up-weighted -> DAR rises
    assert corrected > core.dar_lcms_intact(intensity)


def test_dar_lcms_reduced_corrected_reduces_to_raw():
    lc, hc = {0: 10, 1: 90}, {0: 5, 1: 20, 2: 75}
    assert core.dar_lcms_reduced_corrected(lc, hc) == approx(core.dar_lcms_reduced(lc, hc))


def test_dar_lcms_reduced_corrected_with_response():
    lc, hc = {0: 10, 1: 90}, {0: 5, 1: 20, 2: 75}
    lc_r, hc_r = {0: 1.0, 1: 0.9}, {0: 1.0, 1: 0.9, 2: 0.8}
    assert core.dar_lcms_reduced_corrected(lc, hc, lc_r, hc_r) == approx(5.285460807848867)


def test_dar_lcms_response_rejects_nonpositive_factor():
    with pytest.raises(ValueError):
        core.dar_lcms_intact_corrected({0: 10, 2: 30}, {0: 1.0, 2: 0.0})


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


# ---- In vitro dosing / serial dilution --------------------------------------
def test_plate_well_volume_lookup():
    assert core.plate_well_volume("96") == approx(100.0)
    assert core.plate_well_volume("384") == approx(50.0)
    assert core.plate_well_volume(6) == approx(2000.0)   # int accepted
    with pytest.raises(ValueError):
        core.plate_well_volume("100")


def test_plan_serial_dilution_spike_golden():
    """Spec §8b golden fixture: 1000 µM stock, top 10 µM, 4-fold, 4 pts, 3 reps,
    96-well (100 µL), 10x spike, no overage."""
    p = core.plan_serial_dilution(
        stock_uM=1000, top_uM=10, fold=4, n_points=4, replicates=3,
        well_volume_uL=100.0, dose_mode="spike", dose_factor=10.0, overage=1.0)
    assert p.final_concs_uM == [approx(10.0), approx(2.5), approx(0.625), approx(0.15625)]
    assert p.working_concs_uM[0] == approx(100.0)         # 10x the final
    assert p.dose_factor == approx(10.0)
    assert p.v_add_per_well_uL == approx(100.0 / 9.0)     # well / (10-1)
    assert p.final_well_volume_uL == approx(100.0 + 100.0 / 9.0)
    assert p.v_use_per_tube_uL == approx(100.0 / 3.0)     # add * reps * overage
    assert p.v_transfer_uL == approx(100.0 / 9.0)         # v_use / (fold-1)
    assert p.v_tube_total_uL == approx(100.0 / 3.0 + 100.0 / 9.0)
    assert p.total_stock_uL == approx(4.444444444444444)
    assert p.total_wells == 12
    # per-tube self-consistency: transfer + diluent == prep total for every tube
    for t in p.tubes:
        assert t.total_uL == approx(p.v_tube_total_uL)
        assert t.transfer_in_uL + t.diluent_uL == approx(t.total_uL)
    # downstream tubes each receive the serial transfer and top up with medium
    assert p.tubes[1].transfer_in_uL == approx(100.0 / 9.0)
    assert p.tubes[1].diluent_uL == approx(100.0 / 3.0)
    assert p.warnings == []


def test_plan_serial_dilution_spike_final_conc_back_check():
    """After spiking v_add of a df-fold solution into the well, the in-well
    concentration must equal the intended final value."""
    p = core.plan_serial_dilution(
        stock_uM=1000, top_uM=10, fold=3, n_points=8, replicates=3,
        well_volume_uL=100.0, dose_mode="spike", dose_factor=10.0)
    va = p.v_add_per_well_uL
    for t in p.tubes:
        in_well = t.working_conc_uM * va / (100.0 + va)
        assert in_well == approx(t.final_conc_uM)


def test_plan_serial_dilution_replace_mode():
    p = core.plan_serial_dilution(
        stock_uM=500, top_uM=50, fold=2, n_points=6, replicates=4,
        well_volume_uL=200.0, dose_mode="replace", overage=1.0)
    # replace -> dosing solution already at final conc, full-volume exchange
    assert p.dose_factor == approx(1.0)
    assert p.working_concs_uM == p.final_concs_uM
    assert p.v_add_per_well_uL == approx(200.0)
    assert p.final_well_volume_uL == approx(200.0)
    assert p.total_wells == 24


def test_plan_serial_dilution_warns_dilute_stock():
    p = core.plan_serial_dilution(
        stock_uM=50, top_uM=10, fold=3, n_points=6, replicates=3,
        well_volume_uL=100.0, dose_mode="spike", dose_factor=10.0)
    # top working conc = 100 µM > 50 µM stock
    assert any("too dilute" in w for w in p.warnings)


def test_plan_serial_dilution_rejects_bad_input():
    with pytest.raises(ValueError):
        core.plan_serial_dilution(stock_uM=0, top_uM=10, fold=3, n_points=4,
                                  replicates=3, well_volume_uL=100.0)
    with pytest.raises(ValueError):
        core.plan_serial_dilution(stock_uM=1000, top_uM=10, fold=1, n_points=4,
                                  replicates=3, well_volume_uL=100.0)  # fold must be > 1
    with pytest.raises(ValueError):
        core.plan_serial_dilution(stock_uM=1000, top_uM=10, fold=3, n_points=4,
                                  replicates=3, well_volume_uL=100.0,
                                  dose_mode="spike", dose_factor=1.0)  # bad spike factor


# ---- Plate map --------------------------------------------------------------
def test_plate_map_by_column_layout():
    concs = [10.0, 3.333, 1.111, 0.37]
    pm = core.plate_map(concs, 3, "96", "by_column")
    assert pm.orientation == "by_column"
    assert (pm.rows, pm.cols) == (8, 12)
    assert pm.used == 12 and pm.capacity == 96
    assert pm.warnings == []
    # dose increases across columns; replicates fill rows within a column
    assert (pm.wells[0].row_label, pm.wells[0].col_label) == ("A", 1)
    assert pm.wells[0].point == 0 and pm.wells[0].replicate == 1
    assert (pm.wells[1].row_label, pm.wells[1].col_label) == ("B", 1)  # replicate 2, same dose
    assert pm.wells[1].final_conc_uM == approx(10.0)
    assert (pm.wells[3].row_label, pm.wells[3].col_label) == ("A", 2)  # next dose, column 2
    assert pm.wells[3].point == 1 and pm.wells[3].final_conc_uM == approx(3.333)


def test_plate_map_by_row_layout():
    pm = core.plate_map([10.0, 5.0, 2.5], 4, "96", "by_row")
    assert pm.orientation == "by_row"
    # each point occupies a row; replicates fill columns
    assert (pm.wells[0].row_label, pm.wells[0].col_label) == ("A", 1)
    assert (pm.wells[1].row_label, pm.wells[1].col_label) == ("A", 2)  # replicate 2
    assert (pm.wells[4].row_label, pm.wells[4].col_label) == ("B", 1)  # next dose row


def test_plate_map_row_labels_wrap():
    assert core._row_label(0) == "A"
    assert core._row_label(25) == "Z"
    assert core._row_label(26) == "AA"
    assert core._row_label(31) == "AF"  # 1536-well last row


def test_plate_map_sequential_fallback_and_warning():
    # 14 points > 12 columns -> cannot use by_column, falls back to sequential
    pm = core.plate_map(list(range(1, 15)), 3, "96", "by_column")
    assert pm.orientation == "sequential"
    assert any("do not fit" in w for w in pm.warnings)
    assert pm.used == 42


def test_plate_map_overflow_warning():
    # 96-well, 40 points x 3 = 120 > 96 capacity
    pm = core.plate_map(list(range(40)), 3, "96", "by_column")
    assert pm.used == 120
    assert any("exceeds" in w for w in pm.warnings)
    # off-plate wells are flagged, not silently dropped
    assert any(not w.in_bounds for w in pm.wells)


def test_plate_map_csv_roundtrip():
    pm = core.plate_map([10.0, 5.0], 2, "96", "by_column")
    csv = core.plate_map_to_csv(pm)
    lines = csv.splitlines()
    assert lines[0] == "well,row,column,point,replicate,final_conc_uM"
    assert lines[1].startswith("A1,A,1,1,1,")
    assert len(lines) == 1 + pm.used


def test_plate_map_rejects_bad_input():
    with pytest.raises(ValueError):
        core.plate_map([10.0], 3, "100")           # unknown plate
    with pytest.raises(ValueError):
        core.plate_map([], 3, "96")                # no points
    with pytest.raises(ValueError):
        core.plate_map([10.0], 0, "96")            # no replicates


# ---- Concentration units ----------------------------------------------------
def test_convert_concentration():
    assert core.convert_concentration(5.0, "uM", "nM") == approx(5000.0)
    assert core.convert_concentration(250.0, "nM", "uM") == approx(0.25)
    assert core.convert_concentration(1.0, "uM", "pM") == approx(1e6)
    assert core.convert_concentration(3.0, "nM", "nM") == approx(3.0)
    with pytest.raises(ValueError):
        core.convert_concentration(1.0, "mM", "nM")


# ---- Selection-driven dosing ------------------------------------------------
def test_series_shape_from_selection():
    assert core.series_shape_from_selection(3, 4, "by_column") == (4, 3)
    assert core.series_shape_from_selection(3, 4, "by_row") == (3, 4)
    with pytest.raises(ValueError):
        core.series_shape_from_selection(0, 4)


def test_assign_selection_by_row_decreasing_down():
    # dose decreases top->bottom; replicates across columns
    cells = [(r, c) for r in range(2, 6) for c in range(5, 8)]  # offset 4x3 block
    concs = [10.0, 2.5, 0.625, 0.15625]
    asg = core.assign_selection(cells, concs, "by_row")
    assert asg.rectangular and (asg.rows, asg.cols) == (4, 3)
    assert asg.origin == (2, 5) and asg.warnings == []
    top = next(w for w in asg.wells if (w.row, w.col) == (2, 5))
    bot = next(w for w in asg.wells if (w.row, w.col) == (5, 5))
    assert top.point == 0 and top.final_conc_uM == approx(10.0) and top.replicate == 1
    assert bot.point == 3 and bot.final_conc_uM == approx(0.15625)
    rep2 = next(w for w in asg.wells if (w.row, w.col) == (2, 6))
    assert rep2.point == 0 and rep2.replicate == 2


def test_assign_selection_by_column():
    cells = [(r, c) for r in range(3) for c in range(4)]
    concs = [10.0, 2.5, 0.625, 0.15625]
    asg = core.assign_selection(cells, concs, "by_column")
    a1 = next(w for w in asg.wells if (w.row, w.col) == (0, 0))
    a2 = next(w for w in asg.wells if (w.row, w.col) == (0, 1))
    assert a1.point == 0 and a1.replicate == 1
    assert a2.point == 1 and a2.final_conc_uM == approx(2.5)


def test_assign_selection_non_rectangular_and_overflow():
    cells = [(0, 0), (0, 1), (0, 2), (1, 0)]  # L-shape, not a full box
    asg = core.assign_selection(cells, [10.0, 5.0], "by_column")
    assert asg.rectangular is False
    assert any("along the dose axis" in w for w in asg.warnings)  # 2 pts, span 3
    # well past the series length gets no concentration
    off = next(w for w in asg.wells if (w.row, w.col) == (0, 2))
    assert off.point == 2 and off.final_conc_uM is None


def test_assign_selection_rejects_empty():
    with pytest.raises(ValueError):
        core.assign_selection([], [10.0])


def test_aggregate_dosing_plans():
    p1 = core.plan_serial_dilution(stock_uM=1000, top_uM=10, fold=4, n_points=4,
                                   replicates=3, well_volume_uL=100)
    p2 = core.plan_serial_dilution(stock_uM=1000, top_uM=5, fold=2, n_points=3,
                                   replicates=2, well_volume_uL=100)
    agg = core.aggregate_dosing_plans([p1, p2])
    assert agg["n_groups"] == 2
    assert agg["total_wells"] == p1.total_wells + p2.total_wells
    assert agg["total_stock_uL"] == approx(p1.total_stock_uL + p2.total_stock_uL)


# ---- Dye library ------------------------------------------------------------
def test_dye_library_integrity():
    assert len(core.DYE_LIBRARY) == 22
    for d in core.DYE_LIBRARY:
        # e280 must equal e_lmax * cf280 (spreadsheet definition)
        assert d["e280"] == approx(d["e_lmax"] * d["cf280"])
    af647 = core.get_dye("af647")
    assert af647 is not None and af647["lmax"] == 647


# ---- Payload–linker library -------------------------------------------------
def test_payload_library_integrity():
    assert len(core.PAYLOAD_LIBRARY) == 9
    sourced = [p for p in core.PAYLOAD_LIBRARY if p["sourced"]]
    assert len(sourced) == 3
    # every entry has the full schema
    fields = {"name", "cls", "mw_free", "mw_conj", "lambda_max",
              "eps_lmax_free", "eps280_free", "eps_lmax_conj", "eps280_conj",
              "sourced", "source"}
    for p in core.PAYLOAD_LIBRARY:
        assert fields <= set(p.keys())
        # a reference-only entry must not carry any ε value (no fabricated ε)
        if not p["sourced"]:
            assert p["eps_lmax_free"] is None and p["eps280_free"] is None
            assert p["eps_lmax_conj"] is None and p["eps280_conj"] is None


def test_get_payload_vcmmae():
    p = core.get_payload("vc-MMAE")
    assert p is not None
    assert p["eps_lmax_conj"] == approx(15900.0)
    assert p["eps280_conj"] == approx(1500.0)
    assert p["lambda_max"] == 248
    assert p["mw_conj"] == approx(1316.63)
    assert p["cls"] == "auristatin"
    # case-insensitive
    assert core.get_payload("VC-MMAE") is p


def test_get_payload_maytansinoids():
    dm1 = core.get_payload("DM1")
    assert dm1["eps_lmax_free"] == approx(26790.0)
    assert dm1["eps280_free"] == approx(0.0)
    assert dm1["lambda_max"] == 252
    dm4 = core.get_payload("DM4")
    assert dm4["eps_lmax_free"] == approx(28044.0)
    assert dm4["eps280_free"] == approx(5700.0)


def test_get_payload_reference_only_and_missing():
    mcmmaf = core.get_payload("mc-MMAF")
    assert mcmmaf is not None and mcmmaf["sourced"] is False
    assert mcmmaf["eps_lmax_conj"] is None
    assert core.get_payload("nonexistent") is None


# ---- Product-quality aggregates: SEC purity & free-drug (§9c) ---------------
def test_sec_purity_area_normalised():
    r = core.sec_purity(1900.0, 62.0, 38.0)
    assert r["pct_monomer"] == pytest.approx(95.0, rel=REL)
    assert r["pct_hmw"] == pytest.approx(3.1, rel=REL)
    assert r["pct_lmw"] == pytest.approx(1.9, rel=REL)
    assert r["total"] == pytest.approx(2000.0, rel=REL)
    # percentages sum to 100
    assert r["pct_monomer"] + r["pct_hmw"] + r["pct_lmw"] == pytest.approx(100.0, rel=REL)


def test_sec_purity_two_peak_default_lmw():
    r = core.sec_purity(4900.0, 100.0)
    assert r["pct_monomer"] == pytest.approx(98.0, rel=REL)
    assert r["pct_hmw"] == pytest.approx(2.0, rel=REL)
    assert r["pct_lmw"] == 0.0


def test_sec_purity_rejects_negative_and_zero():
    with pytest.raises(ValueError):
        core.sec_purity(-1.0, 1.0)
    with pytest.raises(ValueError):
        core.sec_purity(1.0, -1.0)
    with pytest.raises(ValueError):
        core.sec_purity(0.0, 0.0)


def test_free_drug_percent():
    r = core.free_drug_percent(1.5, 598.5)
    assert r["pct_free"] == pytest.approx(0.25, rel=REL)
    assert r["pct_conjugated"] == pytest.approx(99.75, rel=REL)
    assert r["total"] == pytest.approx(600.0, rel=REL)
    assert r["pct_free"] + r["pct_conjugated"] == pytest.approx(100.0, rel=REL)


def test_free_drug_percent_rejects_negative_and_zero():
    with pytest.raises(ValueError):
        core.free_drug_percent(-0.1, 1.0)
    with pytest.raises(ValueError):
        core.free_drug_percent(1.0, -0.1)
    with pytest.raises(ValueError):
        core.free_drug_percent(0.0, 0.0)


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


# ---- Extinction-coefficient determination (regression) ----------------------
CONC_M = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6]
AB280 = [0.01 + 50000.0 * c for c in CONC_M]     # eps=50000, intercept=0.01
AB_LMAX = [120000.0 * c for c in CONC_M]          # eps=120000, intercept=0


def test_linear_regression_exact_line():
    fit = core.linear_regression(CONC_M, AB280)
    assert fit["slope"] == approx(50000.0)
    assert fit["intercept"] == approx(0.01)
    assert fit["r_squared"] == approx(1.0)
    assert fit["n"] == 5


def test_extinction_coefficient_path_length_1():
    ec = core.extinction_coefficient(CONC_M, AB280, 1.0)
    assert ec["eps"] == approx(50000.0)
    assert ec["r_squared"] == approx(1.0)


def test_extinction_coefficient_path_length_half():
    ec = core.extinction_coefficient(CONC_M, AB280, 0.5)
    assert ec["eps"] == approx(100000.0)


def test_extinction_coefficient_lmax_series():
    ec = core.extinction_coefficient(CONC_M, AB_LMAX, 1.0)
    assert ec["eps"] == approx(120000.0)
    assert ec["intercept"] == pytest.approx(0.0, abs=1e-9)


def test_linear_regression_rejects_degenerate_input():
    with pytest.raises(ValueError):
        core.linear_regression([1.0], [1.0])          # < 2 points
    with pytest.raises(ValueError):
        core.linear_regression([2.0, 2.0], [1.0, 3.0])  # no distinct x
    with pytest.raises(ValueError):
        core.extinction_coefficient(CONC_M, AB280, 0.0)  # bad path length


# ---- Uncertainty / error propagation ----------------------------------------
def test_dar_uv_uncertainty_matches_point_estimate():
    u = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                                sigma_a280=0.01, sigma_a_lmax=0.01)
    assert u["dar"] == approx(2.4140809554780795)
    assert u["sigma_dar"] == approx(0.03370113601940267)
    assert u["ci95_low"] == approx(2.4140809554780795 - 1.96 * 0.03370113601940267)
    assert u["ci95_high"] == approx(2.4140809554780795 + 1.96 * 0.03370113601940267)


def test_dar_uv_uncertainty_zero_sigma_is_degenerate():
    u = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP)
    assert u["sigma_dar"] == 0.0
    assert u["ci95_low"] == approx(u["dar"])
    assert u["ci95_high"] == approx(u["dar"])


def test_dar_uv_uncertainty_eps_terms_backward_compatible():
    # eps sigmas default to 0 -> identical to the read-only propagation
    u = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                                sigma_a280=0.01, sigma_a_lmax=0.01)
    assert u["sigma_dar"] == approx(0.03370113601940267)


def test_dar_uv_uncertainty_with_eps_terms():
    u = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                                sigma_a280=0.01, sigma_a_lmax=0.01,
                                sigma_eps280_mab=2030.0, sigma_eps280_lp=364.35,
                                sigma_eps_lmax_lp=481.2)
    assert u["dar"] == approx(2.4140809554780795)
    assert u["sigma_dar"] == approx(0.13795625467477654)
    assert u["ci95_low"] == approx(2.4140809554780795 - 1.96 * 0.13795625467477654)
    assert u["ci95_high"] == approx(2.4140809554780795 + 1.96 * 0.13795625467477654)


def test_dar_uv_uncertainty_eps_terms_increase_sigma():
    base = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                                   sigma_a280=0.01, sigma_a_lmax=0.01)
    more = core.dar_uv_uncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                                   sigma_a280=0.01, sigma_a_lmax=0.01,
                                   sigma_eps280_lp=364.35)
    assert more["sigma_dar"] > base["sigma_dar"]


def test_distribution_dispersion_hic():
    d = core.distribution_dispersion({0: 5, 1: 15, 2: 45, 3: 25, 4: 10})
    assert d["mean"] == approx(2.2)
    assert d["variance"] == approx(0.96)
    assert d["sd"] == approx(0.9797958971132712)


def test_dar_lcms_reduced_uncertainty():
    u = core.dar_lcms_reduced_uncertainty({0: 10, 1: 90}, {0: 5, 1: 20, 2: 75})
    assert u["dar"] == approx(5.2)
    assert u["light_sd"] == approx(0.3)
    assert u["heavy_sd"] == approx(0.5567764362830022)
    assert u["sigma_dar"] == approx(0.8944271909999159)


def test_distribution_dispersion_rejects_empty():
    with pytest.raises(ValueError):
        core.distribution_dispersion({})


# ---- DAR distribution predictor (binomial site-occupancy) -------------------
def test_predict_dar_distribution_cysteine_even_ladder():
    # Cysteine/TCEP: 4 interchain disulfides, each carries 2 drugs when reduced.
    # feed 2.5 x eff 0.8 = 2.0 drugs offered -> p_site = 2.0 / (4*2) = 0.25.
    r = core.predict_dar_distribution(4, feed_ratio=2.5, efficiency=0.8)
    assert r["p_site"] == approx(0.25)
    assert r["drugs_per_site"] == 2
    assert r["max_dar"] == 8
    assert r["mean_dar"] == approx(2.0)          # mean = feed*eff
    assert r["variance"] == approx(3.0)          # d^2 * n * p * (1-p) = 4*4*.25*.75
    assert r["sd"] == approx(3.0 ** 0.5)
    assert sum(r["distribution"].values()) == approx(1.0)
    # DAR ladder is even only — no odd species exist
    assert set(r["distribution"].keys()) == {0, 2, 4, 6, 8}
    assert 1 not in r["distribution"]
    assert 3 not in r["distribution"]
    assert r["distribution"][2] == approx(0.421875)
    assert r["distribution"][0] == approx(0.31640625)


def test_predict_dar_distribution_cysteine_full_ladder():
    # n=4, p_site=0.5 directly -> canonical DAR-4 ladder with even support.
    r = core.predict_dar_distribution(4, p_site=0.5)
    assert r["mean_dar"] == approx(4.0)          # d*n*p = 2*4*0.5
    assert r["variance"] == approx(4.0)          # d^2*n*p*(1-p) = 4*4*.25
    assert set(r["distribution"].keys()) == {0, 2, 4, 6, 8}
    assert r["distribution"][4] == approx(0.375)
    assert r["distribution"][0] == approx(0.0625)
    assert r["distribution"][8] == approx(0.0625)


def test_predict_dar_distribution_lysine_binomial():
    # Stochastic amine: single drug per site, smooth binomial over 0..n.
    r = core.predict_dar_distribution(8, p_site=0.5, drugs_per_site=1)
    assert r["drugs_per_site"] == 1
    assert r["mean_dar"] == approx(4.0)          # n*p
    assert r["variance"] == approx(2.0)          # n*p*(1-p)
    assert set(r["distribution"].keys()) == set(range(9))  # 0..8, all integers
    assert r["distribution"][3] == approx(0.21875)


def test_predict_dar_distribution_efficiency_scales_p():
    full = core.predict_dar_distribution(4, feed_ratio=4.0, efficiency=1.0)
    half = core.predict_dar_distribution(4, feed_ratio=4.0, efficiency=0.5)
    assert full["p_site"] == approx(0.5)         # 4 / (4*2)
    assert half["p_site"] == approx(0.25)
    assert full["mean_dar"] == approx(4.0)
    assert half["mean_dar"] == approx(2.0)


def test_predict_dar_distribution_rejects_bad_input():
    with pytest.raises(ValueError):
        core.predict_dar_distribution(0, p_site=0.5)          # n_sites < 1
    with pytest.raises(ValueError):
        core.predict_dar_distribution(4)                      # no p or feed
    with pytest.raises(ValueError):
        core.predict_dar_distribution(4, feed_ratio=9.0)      # p_site > 1 (>capacity 8)
    with pytest.raises(ValueError):
        core.predict_dar_distribution(4, p_site=0.5, drugs_per_site=0)  # bad d


def test_measured_dar_distribution_skewed_profile():
    r = core.measured_dar_distribution({0: 5, 2: 35, 4: 40, 6: 15, 8: 5})
    assert r["mean_dar"] == pytest.approx(3.6, rel=REL)
    assert r["sd"] == pytest.approx(1.854723699099141, rel=REL)
    assert r["variance"] == pytest.approx(3.44, rel=REL)
    assert r["total"] == pytest.approx(100.0, rel=REL)
    assert r["fractions"][2] == pytest.approx(0.35, rel=REL)
    assert r["fractions"][4] == pytest.approx(0.40, rel=REL)
    assert sum(r["fractions"].values()) == pytest.approx(1.0, rel=REL)


def test_measured_dar_distribution_matches_dispersion_and_normalises():
    # abundances need not sum to 100; mean is unchanged by a common scale factor
    a = core.measured_dar_distribution({0: 10, 2: 70, 4: 80, 6: 30, 8: 10})
    b = core.measured_dar_distribution({0: 5, 2: 35, 4: 40, 6: 15, 8: 5})
    assert a["mean_dar"] == pytest.approx(b["mean_dar"], rel=REL)
    assert a["sd"] == pytest.approx(b["sd"], rel=REL)
    assert a["fractions"][2] == pytest.approx(0.35, rel=REL)


def test_measured_dar_distribution_rejects_bad_input():
    with pytest.raises(ValueError):
        core.measured_dar_distribution({0: -1, 2: 5})
    with pytest.raises(ValueError):
        core.measured_dar_distribution({0: 0, 2: 0})


# ---------------------------------------------------------------------------
# Physical-bounds guards (spec §10b)
# ---------------------------------------------------------------------------
def test_check_physical_bounds_clean_case_no_warnings():
    assert core.check_physical_bounds(
        dar=2.5, r=0.1, eps280=203000, eps_lmax=210000,
        concentrations={"conc": 5.0},
    ) == []


def test_check_physical_bounds_returns_empty_when_nothing_passed():
    assert core.check_physical_bounds() == []


def test_check_physical_bounds_negative_dar():
    w = core.check_physical_bounds(dar=-1.2)
    assert [x["code"] for x in w] == ["dar_negative"]


def test_check_physical_bounds_high_dar():
    w = core.check_physical_bounds(dar=20.0)
    assert [x["code"] for x in w] == ["dar_high"]
    # boundary is inclusive-below: exactly DAR_MAX_PLAUSIBLE is allowed
    assert core.check_physical_bounds(dar=core.DAR_MAX_PLAUSIBLE) == []


def test_check_physical_bounds_negative_r():
    assert [x["code"] for x in core.check_physical_bounds(r=-0.5)] == ["r_negative"]


def test_check_physical_bounds_high_r():
    assert [x["code"] for x in core.check_physical_bounds(r=60.0)] == ["r_high"]
    assert core.check_physical_bounds(r=core.R_MAX_PLAUSIBLE) == []


def test_check_physical_bounds_eps_inconsistent():
    w = core.check_physical_bounds(eps280=203000, eps_lmax=5000)
    assert [x["code"] for x in w] == ["eps_inconsistent"]
    # equal is allowed (not strictly greater)
    assert core.check_physical_bounds(eps280=5000, eps_lmax=5000) == []


def test_check_physical_bounds_negative_concentration():
    w = core.check_physical_bounds(concentrations={"mAb conc (mg/mL)": -3.0})
    assert [x["code"] for x in w] == ["conc_negative"]
    assert "mAb conc (mg/mL)" in w[0]["message"]


def test_check_physical_bounds_collects_multiple():
    w = core.check_physical_bounds(
        dar=-1.0, r=-0.5, eps280=1000, eps_lmax=100,
        concentrations={"c1": -2.0, "c2": 5.0},
    )
    codes = sorted(x["code"] for x in w)
    assert codes == ["conc_negative", "dar_negative", "eps_inconsistent", "r_negative"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
