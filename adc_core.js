/*
 * adc_core.js — pure calculation core for the ADC Calculator (ES module).
 *
 * Direct 1:1 port of adc_core.py. Zero DOM, zero I/O. Every function maps to a
 * formula in adc_spec.md and is validated against the same golden values as the
 * Python test suite (see selfCheck() at the bottom / node adc_core.js).
 *
 * Units: molar epsilon in M^-1 cm^-1; mass epsilon in (mg/mL)^-1 cm^-1.
 */

export const ELLMAN_E412 = 14150.0;          // DTNB product at 412 nm (PDF p.10)
export const IGG1_INTERCHAIN_DISULFIDES = 4; // -> 8 free thiols on full reduction

// ── 1. Extinction-coefficient conversions ─────────────────────────────────
export function epsMassToMolar(epsMass, mw) { return epsMass * mw; }
export function epsMolarToMass(epsMolar, mw) {
  if (mw === 0) throw new Error("MW must be non-zero");
  return epsMolar / mw;
}

// ── 2. DAR by UV/Vis (spreadsheet I16) ─────────────────────────────────────
export function darUV(a280, aLmax, eps280Mab, epsLmaxMab, eps280Lp, epsLmaxLp) {
  if (a280 === 0) throw new Error("A280 must be non-zero");
  const r = aLmax / a280;
  const denom = r * eps280Lp - epsLmaxLp;
  if (denom === 0) throw new Error("Degenerate optics: R*eps280_LP == eps_lmax_LP");
  return (epsLmaxMab - r * eps280Mab) / denom;
}

// ── 3. ADC molecular properties ────────────────────────────────────────────
export function adcProperties(dar, mwMab, mwLp, eps280Mab, eps280Lp,
                              a280 = null, dilutionFactor = 1.0) {
  const mwAdc = mwMab + mwLp * dar;
  const epsMolar = eps280Mab + dar * eps280Lp;
  const epsMass = epsMolar / mwAdc;
  let concMgml = null, concUM = null;
  if (a280 !== null) {
    concMgml = (a280 * dilutionFactor) / epsMass;
    concUM = (concMgml / mwAdc) * 1e6;
  }
  return { dar, mwAdc, eps280AdcMolar: epsMolar, eps280AdcMass: epsMass, concMgml, concUM };
}

export function concUnknownPayloadEps(a280, aLmax, eps280Lp, epsLmaxLp,
                                      eps280MabMass, dilutionFactor = 1.0) {
  const aLp280 = (aLmax * eps280Lp) / epsLmaxLp;
  return ((a280 - aLp280) * dilutionFactor) / eps280MabMass;
}

// ── 4. Ellman's free-thiol assay (PDF p.10) ────────────────────────────────
export function ellmanFreeThiols(a412, a280, eps280Mab, eps412 = ELLMAN_E412) {
  if (a280 === 0) throw new Error("A280 must be non-zero");
  return (a412 / eps412) / (a280 / eps280Mab);
}

// ── 5. DAR by analytical HIC (PDF p.14-15) ─────────────────────────────────
function _total(obj) { return Object.values(obj).reduce((s, v) => s + v, 0); }

export function darHIC(aucBySpecies) {
  const total = _total(aucBySpecies);
  if (total === 0) throw new Error("Total AUC must be non-zero");
  let num = 0;
  for (const [k, a] of Object.entries(aucBySpecies)) num += Number(k) * a;
  return num / total;
}

export function darSpeciesFractions(aucBySpecies) {
  const total = _total(aucBySpecies);
  if (total === 0) throw new Error("Total AUC must be non-zero");
  const out = {};
  for (const [k, a] of Object.entries(aucBySpecies)) out[k] = (a / total) * 100.0;
  return out;
}

// ── 6. DAR by LC-MS (PDF p.17-18) ──────────────────────────────────────────
function _weightedMean(intensityByK) {
  const total = _total(intensityByK);
  if (total === 0) throw new Error("Total intensity must be non-zero");
  let num = 0;
  for (const [k, i] of Object.entries(intensityByK)) num += Number(k) * i;
  return num / total;
}

export function darLcmsReduced(lightChain, heavyChain) {
  return 2 * _weightedMean(lightChain) + 2 * _weightedMean(heavyChain);
}

export function darLcmsIntact(intensityByK) { return _weightedMean(intensityByK); }

// Per-chain intermediates for the reduced LC-MS DAR (sumI, sumKI, mean, 2*mean).
function _chainStats(intensityByK) {
  const sumI = _total(intensityByK);
  if (sumI === 0) throw new Error("Total intensity must be non-zero");
  let sumKI = 0;
  for (const [k, i] of Object.entries(intensityByK)) sumKI += Number(k) * i;
  const mean = sumKI / sumI;
  return { sum_i: sumI, sum_ki: sumKI, mean: mean, contribution: 2 * mean };
}

export function darLcmsReducedBreakdown(lightChain, heavyChain) {
  const light = _chainStats(lightChain);
  const heavy = _chainStats(heavyChain);
  return { light: light, heavy: heavy,
           dar: light.contribution + heavy.contribution };
}

// ── 7. Conjugation designer (PDF p.3-6) ────────────────────────────────────
export function molesFromMass(massMg, mw) { return (massMg / 1000.0) / mw; }

export function reagentVolumeUL(moles, stockConcMM) {
  if (stockConcMM === 0) throw new Error("Stock concentration must be non-zero");
  const liters = moles / (stockConcMM * 1e-3);
  return liters * 1e6;
}

const ROUTE_NOTES = {
  lysine: [
    "Activated ester 5-20 equiv (optimize by screening).",
    "mAb at 5 g/L; add 10% v/v 1.0 M NaHCO3; final 10% v/v DMSO.",
    "2 h at 25 C; quench with L-lysine 100 equiv, 20 min.",
  ],
  cys_native: [
    "1 equiv TCEP per disulfide; 6.0 equiv for full IgG1 reduction.",
    "Maleimide 1.5 equiv per free thiol; +5 mM EDTA; 10% v/v DMSO.",
    "Maleimide 1 h / 25 C; quench L-cysteine 100 equiv, 10 min.",
  ],
  cys_rebridge: [
    "Full reduction: 6.0 equiv TCEP, 1 h / 37 C.",
    "Rebridging reagent equivalents by screening; 16 h / 25 C.",
    "Purify by preparative HIC for defined DAR.",
  ],
  thiomab: [
    "Full reduction 10 equiv TCEP -> desalt to remove capping cysteines.",
    "Reoxidize with dhAA 30 equiv, 4 h / 25 C.",
    "Maleimide 1.5 equiv per free thiol, 1 h / 25 C.",
  ],
  mtgase: [
    "PNGase F 10 uL/mg to deglycosylate; tag 80 equiv (e.g. azido-PEG3-amine).",
    "MTGase 5.5 U/mg, 16 h / 25 C -> site-specific DAR 2.0.",
    "Click reagent (e.g. DBCO-LP) 3.0 equiv, 6 h / 25 C.",
  ],
};

export function planConjugation({
  route, massMabMg, mwMab = 145000.0, mabConcGPerL = 5.0,
  reagents = [], freeThiolsPerMab = null, dmsoFraction = 0.10,
}) {
  const nMab = molesFromMass(massMabMg, mwMab);
  if (freeThiolsPerMab === null) {
    freeThiolsPerMab = ["cys_native", "cys_rebridge", "thiomab"].includes(route)
      ? IGG1_INTERCHAIN_DISULFIDES * 2 : 0.0;
  }
  const nThiol = nMab * freeThiolsPerMab;
  const baseVolUL = (massMabMg / mabConcGPerL) * 1000.0;

  const additions = reagents.map((rg) => {
    const basis = rg.basis || "mAb";
    const nRef = basis === "thiol" ? nThiol : nMab;
    const moles = rg.equiv * nRef;
    return {
      name: rg.name, equivalents: rg.equiv,
      basis: basis === "thiol" ? "per free thiol" : "per mAb",
      moles, stockConcMM: rg.stock_mM, volumeUL: reagentVolumeUL(moles, rg.stock_mM),
    };
  });

  const reagentVol = additions.reduce((s, a) => s + a.volumeUL, 0);
  const reactionVolume = Math.max(baseVolUL / (1 - dmsoFraction), baseVolUL + reagentVol);
  return {
    route, massMabMg, nMabNmol: nMab * 1e9, freeThiolsPerMab, additions,
    reactionVolumeUL: reactionVolume, dmsoVolumeUL: reactionVolume * dmsoFraction,
    notes: ROUTE_NOTES[route] || [],
  };
}

// ── 8. Yield & formulation ─────────────────────────────────────────────────
export function yieldAndFormulation(concMgml, volumeML, mwAdc,
                                    startingMassMg = null, targetConcMgml = null) {
  const recovered = concMgml * volumeML;
  const yieldPct = startingMassMg ? (recovered / startingMassMg) * 100.0 : null;
  const molarNmol = (recovered / mwAdc) * 1e6;
  let vFinal = null, vChange = null;
  if (targetConcMgml) { vFinal = recovered / targetConcMgml; vChange = vFinal - volumeML; }
  return { recoveredMassMg: recovered, yieldPct, molarAmountNmol: molarNmol,
           vFinalML: vFinal, vChangeML: vChange };
}

// ── 10. Embedded dye / payload library (spreadsheet cols L-R) ───────────────
export const DYE_LIBRARY = [
  { name: "AF350", mw: 410, e_lmax: 19000, cf280: 0.19, e280: 3610, lmax: 350, comment: null },
  { name: "AF405", mw: 1028, e_lmax: 34000, cf280: 0.70, e280: 23800, lmax: 405, comment: null },
  { name: "AF430", mw: 702, e_lmax: 16000, cf280: 0.28, e280: 4480, lmax: 430, comment: null },
  { name: "AF488", mw: 643, e_lmax: 71000, cf280: 0.11, e280: 7810, lmax: 488, comment: null },
  { name: "AF500", mw: 700, e_lmax: 71000, cf280: 0.18, e280: 12780, lmax: 500, comment: null },
  { name: "AF514", mw: 714, e_lmax: 80000, cf280: 0.18, e280: 14400, lmax: 514, comment: null },
  { name: "AF532", mw: 721, e_lmax: 81000, cf280: 0.09, e280: 7290, lmax: 532, comment: null },
  { name: "AF546", mw: 1260, e_lmax: 112000, cf280: 0.12, e280: 13440, lmax: 546, comment: null },
  { name: "AF555", mw: 1250, e_lmax: 150000, cf280: 0.08, e280: 12000, lmax: 555, comment: null },
  { name: "AF568", mw: 792, e_lmax: 91300, cf280: 0.46, e280: 41998, lmax: 568, comment: null },
  { name: "AF594", mw: 820, e_lmax: 90000, cf280: 0.56, e280: 50400, lmax: 594, comment: null },
  { name: "AF610X", mw: 1285, e_lmax: 132000, cf280: 0.44, e280: 58080, lmax: 610, comment: null },
  { name: "AF633", mw: 1200, e_lmax: 100000, cf280: 0.55, e280: 55000, lmax: 633, comment: null },
  { name: "AF647", mw: 1300, e_lmax: 239000, cf280: 0.03, e280: 7170, lmax: 647, comment: null },
  { name: "AF660", mw: 1100, e_lmax: 132000, cf280: 0.10, e280: 13200, lmax: 660, comment: null },
  { name: "AF680", mw: 1150, e_lmax: 184000, cf280: 0.05, e280: 9200, lmax: 680, comment: null },
  { name: "AF700", mw: 1400, e_lmax: 192000, cf280: 0.07, e280: 13440, lmax: 700, comment: null },
  { name: "AF750", mw: 1300, e_lmax: 240000, cf280: 0.04, e280: 9600, lmax: 750, comment: null },
  { name: "AF790", mw: 1750, e_lmax: 260000, cf280: 0.08, e280: 20800, lmax: 790, comment: null },
  { name: "pHRodo deep red", mw: 1300, e_lmax: 140000, cf280: 0.33, e280: 46200, lmax: null, comment: "Dilute 1:1 with H2SO4 1M" },
  { name: "pHRodo iFL Red", mw: 1000, e_lmax: 65000, cf280: 0.12, e280: 7800, lmax: 560, comment: "Dilute 1:1 with H2SO4 1M" },
  { name: "pHRodo iFL Green", mw: 1000, e_lmax: 74500, cf280: 0.20, e280: 14900, lmax: null, comment: "Dilute 1:1 with H2SO4 1M" },
];

export function getDye(name) {
  const key = String(name).trim().toLowerCase();
  return DYE_LIBRARY.find((d) => d.name.toLowerCase() === key) || null;
}

// ============================================================================
// ADC registry: log each preparation / assay and aggregate across batches.
// Pure functions mirroring adc_core.py (make_registry_record / registry_to_csv
// / summarizeRegistry). Persistence (localStorage/files) lives in the UI.
// ============================================================================
export const REGISTRY_FIELDS = [
  "id", "date", "operator", "assay", "route",
  "dar", "conc_mgml", "conc_uM", "free_thiols", "notes",
];
const REGISTRY_NUMERIC = new Set(["dar", "conc_mgml", "conc_uM", "free_thiols"]);

function coerceNum(value) {
  if (value === null || value === undefined || typeof value === "boolean") return null;
  if (typeof value === "number") return Number.isFinite(value) ? value : null;
  const s = String(value).trim();
  if (s === "") return null;
  const n = Number(s);
  return Number.isFinite(n) ? n : null;
}

export function makeRegistryRecord(fields) {
  const f = fields || {};
  const rec = {};
  for (const key of REGISTRY_FIELDS) {
    const val = f[key];
    if (REGISTRY_NUMERIC.has(key)) rec[key] = coerceNum(val);
    else rec[key] = (val === null || val === undefined) ? "" : String(val).trim();
  }
  if (rec.id === "") throw new Error("registry record requires a non-empty 'id'");
  return rec;
}

function fmtCell(value) {
  if (value === null || value === undefined) return "";
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : String(value);
  return String(value);
}
function csvEscape(cell) {
  return /[",\n\r]/.test(cell) ? '"' + cell.replace(/"/g, '""') + '"' : cell;
}
export function registryToCsv(records) {
  const lines = [REGISTRY_FIELDS.join(",")];
  for (const r of records) {
    lines.push(REGISTRY_FIELDS.map((k) => csvEscape(fmtCell(r[k]))).join(","));
  }
  return lines.join("\n");
}

export function summarizeRegistry(records) {
  const byAssay = {}, byId = {}, dars = [];
  for (const r of records) {
    const a = r.assay || "(none)"; byAssay[a] = (byAssay[a] || 0) + 1;
    const i = r.id || "(none)"; byId[i] = (byId[i] || 0) + 1;
    const d = r.dar;
    if (typeof d === "number" && Number.isFinite(d)) dars.push(d);
  }
  const out = {
    n_total: records.length, n_ids: Object.keys(byId).length,
    by_assay: byAssay, by_id: byId, dar_n: dars.length,
    dar_mean: null, dar_sd: null, dar_min: null, dar_max: null,
  };
  if (dars.length) {
    const mean = dars.reduce((a, b) => a + b, 0) / dars.length;
    out.dar_mean = mean; out.dar_min = Math.min(...dars); out.dar_max = Math.max(...dars);
    out.dar_sd = dars.length > 1
      ? Math.sqrt(dars.reduce((a, x) => a + (x - mean) ** 2, 0) / (dars.length - 1))
      : 0.0;
  }
  return out;
}

// ── Self-check: same golden values as the Python test suite ────────────────
export function selfCheck() {
  const REL = 1e-6;
  const ok = (got, want, label) => {
    const rel = want === 0 ? Math.abs(got) : Math.abs(got - want) / Math.abs(want);
    if (rel > REL) throw new Error(`FAIL ${label}: got ${got}, want ${want}`);
    return `PASS ${label}`;
  };
  const MW_MAB = 145000.0, EPS280_MAB = 203000.0, A280 = 7.4315, A_LMAX = 0.7827;
  const MW_LP = 1559.62, EPS280_LP = 7287.0, EPS_LMAX_LP = 9624.0;
  const R = A_LMAX / A280;
  const results = [];
  results.push(ok(epsMassToMolar(1.4, MW_MAB), 203000.0, "eps conv"));
  results.push(ok(darUV(A280, A_LMAX, EPS280_MAB, R * EPS280_MAB, EPS280_LP, EPS_LMAX_LP), 0.0, "DAR case A"));
  const darB = darUV(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP);
  results.push(ok(darB, 2.4140809554780795, "DAR case B"));
  const pA = adcProperties(0.0, MW_MAB, MW_LP, EPS280_MAB, EPS280_LP, A280, 1.0);
  results.push(ok(pA.concMgml, 5.3082142857142856, "conc mg/mL case A"));
  results.push(ok(pA.concUM, 36.608374384236456, "conc uM case A"));
  const pB = adcProperties(darB, MW_MAB, MW_LP, EPS280_MAB, EPS280_LP, A280, 1.0);
  results.push(ok(pB.mwAdc, 148765.04893978272, "MW_ADC case B"));
  results.push(ok(pB.concMgml, 5.011743075614535, "conc mg/mL case B"));
  results.push(ok(ellmanFreeThiols(0.5, 0.6, 203000.0), 11.955241460541815, "Ellman"));
  results.push(ok(darHIC({ 0: 5, 1: 15, 2: 45, 3: 25, 4: 10 }), 2.2, "HIC"));
  results.push(ok(darLcmsReduced({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 }), 5.2, "LC-MS reduced"));
  const bd = darLcmsReducedBreakdown({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 });
  results.push(ok(bd.light.contribution, 1.8, "LC-MS breakdown LC 2*mean"));
  results.push(ok(bd.heavy.contribution, 3.4, "LC-MS breakdown HC 2*mean"));
  results.push(ok(bd.dar, 5.2, "LC-MS breakdown DAR"));
  // registry
  const rr = [
    makeRegistryRecord({ id: "ADC-001", assay: "uv", dar: "2.414" }),
    makeRegistryRecord({ id: "ADC-002", assay: "hic", dar: 2.2 }),
    makeRegistryRecord({ id: "ADC-002", assay: "lcms", dar: 5.2 }),
  ];
  const rs = summarizeRegistry(rr);
  ok(rs.n_total, 3, "registry n"); ok(rs.n_ids, 2, "registry ids");
  ok(rs.dar_mean, 3.271333333333333, "registry DAR mean");
  ok(rs.dar_sd, 1.6736981010126448, "registry DAR sd");
  results.push("PASS registry (3 records)");
  if (DYE_LIBRARY.length !== 22) throw new Error("dye library length != 22");
  results.push(`PASS dye library (${DYE_LIBRARY.length} entries)`);
  return results;
}

// Run self-check when executed directly under Node (node adc_core.js)
if (typeof process !== "undefined" && process.argv && process.argv[1] &&
    process.argv[1].endsWith("adc_core.js")) {
  for (const line of selfCheck()) console.log(line);
  console.log("All self-checks passed.");
}
