/*
 * adc_core.js — pure calculation core for the ADC Calculator (ES module).
 *
 * Direct 1:1 port of adc_core.py. Zero DOM, zero I/O. Every function maps to a
 * formula in adc_spec.md and is pinned to the same golden values as the Python
 * test suite (see selfCheck() at the bottom / node adc_core.js), so the two
 * cores and the embedded HTML copy cannot drift apart.
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

export function scatterAbsorbance(aRef, refWavelength, targetWavelength, exponent = 4.0) {
  if (refWavelength <= 0 || targetWavelength <= 0) throw new Error("wavelengths must be positive");
  return aRef * Math.pow(refWavelength / targetWavelength, exponent);
}

export function scatterCorrectedAbsorbance(a280, aLmax, aRef,
                                           refWavelength = 320.0,
                                           lmaxWavelength = 495.0,
                                           exponent = 4.0) {
  const s280 = scatterAbsorbance(aRef, refWavelength, 280.0, exponent);
  const sLmax = scatterAbsorbance(aRef, refWavelength, lmaxWavelength, exponent);
  return {
    a280Corrected: Math.max(a280 - s280, 0.0),
    aLmaxCorrected: Math.max(aLmax - sLmax, 0.0),
    scatter280: s280,
    scatterLmax: sLmax,
  };
}

// ── 4. Ellman's free-thiol assay (PDF p.10) ────────────────────────────────
export function ellmanFreeThiols(a412, a280, eps280Mab, eps412 = ELLMAN_E412, a412Blank = 0.0, a280Blank = 0.0) {
  // DTNB/TNB and residual reductant absorb at 412 and 280 nm; a matched
  // reagent blank measures that background. Blanks default to 0 so the call
  // reduces exactly to the raw ratio when no blank is supplied.
  const netA280 = a280 - a280Blank;
  if (netA280 <= 0) throw new Error("net A280 (A280 - blank) must be positive");
  return ((a412 - a412Blank) / eps412) / (netA280 / eps280Mab);
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

function _speciesMolarAbundance(aucBySpecies, eps280Mab, eps280Payload) {
  const mol = {};
  for (const [k, a] of Object.entries(aucBySpecies)) {
    const epsK = eps280Mab + Number(k) * eps280Payload;
    if (epsK <= 0) throw new Error("per-species eps280 must be > 0");
    mol[k] = a / epsK;
  }
  return mol;
}

export function darHICCorrected(aucBySpecies, eps280Mab, eps280Payload) {
  const mol = _speciesMolarAbundance(aucBySpecies, eps280Mab, eps280Payload);
  const total = _total(mol);
  if (total === 0) throw new Error("Total corrected abundance must be non-zero");
  let num = 0;
  for (const [k, m] of Object.entries(mol)) num += Number(k) * m;
  return num / total;
}

export function darSpeciesFractionsCorrected(aucBySpecies, eps280Mab, eps280Payload) {
  const mol = _speciesMolarAbundance(aucBySpecies, eps280Mab, eps280Payload);
  const total = _total(mol);
  if (total === 0) throw new Error("Total corrected abundance must be non-zero");
  const out = {};
  for (const [k, m] of Object.entries(mol)) out[k] = (m / total) * 100.0;
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

function _responseCorrectedMean(intensityByK, responseByK = null) {
  // Optional ionization-response correction: ESI efficiency changes with
  // drug load, so molar abundance is recovered as I_k / r_k before the
  // weighted mean. responseByK=null (or a missing key -> r=1) reproduces the
  // uncorrected intensity-weighted mean exactly.
  const molar = {};
  for (const [k, i] of Object.entries(intensityByK)) {
    if (responseByK == null) { molar[k] = i; continue; }
    const r = k in responseByK ? responseByK[k] : 1.0;
    if (r <= 0) throw new Error("response factor must be positive");
    molar[k] = i / r;
  }
  let total = 0, num = 0;
  for (const [k, m] of Object.entries(molar)) { total += m; num += Number(k) * m; }
  if (total === 0) throw new Error("Total corrected abundance must be non-zero");
  return num / total;
}

export function darLcmsReduced(lightChain, heavyChain) {
  return 2 * _weightedMean(lightChain) + 2 * _weightedMean(heavyChain);
}

export function darLcmsIntact(intensityByK) { return _weightedMean(intensityByK); }

export function darLcmsIntactCorrected(intensityByK, responseByK = null) {
  return _responseCorrectedMean(intensityByK, responseByK);
}

export function darLcmsReducedCorrected(lightChain, heavyChain, lightResponse = null, heavyResponse = null) {
  return 2 * _responseCorrectedMean(lightChain, lightResponse) + 2 * _responseCorrectedMean(heavyChain, heavyResponse);
}

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

// ── 8b. In vitro dosing / serial dilution ───────────────────────────────────
// Default per-well working (assay) volumes in uL, keyed by plate format.
export const PLATE_WELL_VOLUME_UL = {
  "6": 2000.0, "12": 1000.0, "24": 500.0, "48": 250.0,
  "96": 100.0, "384": 50.0, "1536": 10.0,
};

export function plateWellVolume(plateType) {
  const key = String(plateType).trim();
  if (!(key in PLATE_WELL_VOLUME_UL)) throw new Error(`unknown plate type ${plateType}`);
  return PLATE_WELL_VOLUME_UL[key];
}

export function planSerialDilution({
  stockUM, topUM, fold, nPoints, replicates, wellVolumeUL,
  doseMode = "spike", doseFactor = 10.0, doseVolumeUL = null,
  overage = 1.1, extraDeadUL = 0.0,
}) {
  if (stockUM <= 0) throw new Error("stock concentration must be > 0");
  if (topUM <= 0) throw new Error("top concentration must be > 0");
  if (fold <= 1) throw new Error("fold-dilution must be > 1");
  if (nPoints < 1) throw new Error("need at least 1 concentration point");
  if (replicates < 1) throw new Error("need at least 1 replicate");
  if (wellVolumeUL <= 0) throw new Error("well volume must be > 0");
  if (doseMode !== "spike" && doseMode !== "replace") throw new Error("doseMode must be 'spike' or 'replace'");
  if (doseMode === "spike" && doseFactor <= 1) throw new Error("doseFactor must be > 1 for spike mode");

  const n = Math.trunc(nPoints);
  const finalConcs = [];
  for (let i = 0; i < n; i++) finalConcs.push(topUM / fold ** i);
  const df = doseMode === "replace" ? 1.0 : Number(doseFactor);
  const workingConcs = finalConcs.map((c) => c * df);

  let vAdd, finalWellVol;
  if (doseMode === "replace") { vAdd = wellVolumeUL; finalWellVol = wellVolumeUL; }
  else {
    vAdd = doseVolumeUL ? doseVolumeUL : wellVolumeUL / (df - 1.0);
    finalWellVol = wellVolumeUL + vAdd;
  }

  const vUse = vAdd * replicates * overage + extraDeadUL;
  const vXfer = vUse / (fold - 1.0);
  const vTubeTotal = vUse + vXfer;

  const stockVolTube1 = vTubeTotal * workingConcs[0] / stockUM;
  const tubes = [];
  for (let i = 0; i < n; i++) {
    let source, transferIn, diluent;
    if (i === 0) { source = "ADC stock"; transferIn = stockVolTube1; diluent = vTubeTotal - stockVolTube1; }
    else { source = `tube ${i} (transfer)`; transferIn = vXfer; diluent = vUse; }
    tubes.push({
      index: i + 1, finalConcUM: finalConcs[i], workingConcUM: workingConcs[i],
      source, transferInUL: transferIn, diluentUL: diluent,
      totalUL: transferIn + diluent, dispenseToWellsUL: vUse,
    });
  }

  const totalDiluent = tubes.reduce((s, t) => s + t.diluentUL, 0);
  const warnings = [];
  if (workingConcs[0] > stockUM) {
    warnings.push(`Top working conc ${workingConcs[0].toPrecision(4)} uM exceeds stock `
      + `${stockUM.toPrecision(4)} uM - stock too dilute for this design.`);
  }
  if (stockVolTube1 < 1.0) {
    warnings.push(`Tube-1 stock volume ${stockVolTube1.toPrecision(3)} uL < 1 uL - hard to `
      + `pipette accurately; use an intermediate dilution or a larger prep volume.`);
  }

  return {
    finalConcsUM: finalConcs, workingConcsUM: workingConcs, doseMode, doseFactor: df,
    vAddPerWellUL: vAdd, finalWellVolumeUL: finalWellVol,
    vUsePerTubeUL: vUse, vTransferUL: vXfer, vTubeTotalUL: vTubeTotal,
    tubes, totalStockUL: stockVolTube1, totalDiluentUL: totalDiluent,
    totalWells: n * replicates, warnings,
  };
}

// ── 8c. Plate map ────────────────────────────────────────────────────────────
export const PLATE_GEOMETRY = {
  "6": [2, 3], "12": [3, 4], "24": [4, 6], "48": [6, 8],
  "96": [8, 12], "384": [16, 24], "1536": [32, 48],
};

export function rowLabel(i) {
  let s = ""; i += 1;
  while (i > 0) { const r = (i - 1) % 26; s = String.fromCharCode(65 + r) + s; i = Math.floor((i - 1) / 26); }
  return s;
}

export function plateMap(finalConcsUM, replicates, plateType, orientation = "by_column") {
  const key = String(plateType).trim();
  if (!(key in PLATE_GEOMETRY)) throw new Error(`unknown plate type ${plateType}`);
  if (replicates < 1) throw new Error("need at least 1 replicate");
  const concs = Array.from(finalConcsUM);
  const n = concs.length;
  if (n < 1) throw new Error("need at least 1 concentration point");

  const [rows, cols] = PLATE_GEOMETRY[key];
  const capacity = rows * cols;
  const used = n * replicates;
  const warnings = [];

  let orient;
  if (orientation === "by_column" && n <= cols && replicates <= rows) orient = "by_column";
  else if (orientation === "by_row" && n <= rows && replicates <= cols) orient = "by_row";
  else {
    orient = "sequential";
    if (orientation === "by_column" || orientation === "by_row") {
      warnings.push(`${n} points x ${replicates} replicates do not fit the ${rows}x${cols} `
        + `grid in ${orientation} layout; using sequential fill.`);
    }
  }

  const placed = [];
  if (orient === "by_column") {
    for (let i = 0; i < n; i++) for (let r = 0; r < replicates; r++) placed.push([r, i, i, r]);
  } else if (orient === "by_row") {
    for (let i = 0; i < n; i++) for (let r = 0; r < replicates; r++) placed.push([i, r, i, r]);
  } else {
    let pos = 0;
    for (let i = 0; i < n; i++) for (let r = 0; r < replicates; r++) { placed.push([Math.floor(pos / cols), pos % cols, i, r]); pos++; }
  }

  if (used > capacity) warnings.push(`${used} wells needed exceeds ${key}-well capacity (${capacity}).`);

  const wells = placed.map(([rr, cc, pi, rep0]) => {
    const inside = rr < rows && cc < cols;
    return {
      row: rr, col: cc, rowLabel: inside ? rowLabel(rr) : null, colLabel: cc + 1,
      point: pi, replicate: rep0 + 1, finalConcUM: concs[pi], inBounds: inside,
    };
  });

  return { plateType: key, rows, cols, orientation: orient, capacity, used, wells, warnings };
}

export function plateMapToCsv(pmap) {
  const lines = ["well,row,column,point,replicate,final_conc_uM"];
  for (const w of pmap.wells) {
    const well = w.rowLabel ? `${w.rowLabel}${w.colLabel}` : "(off-plate)";
    lines.push(`${well},${w.rowLabel || ""},${w.colLabel},${w.point + 1},`
      + `${w.replicate},${Number(w.finalConcUM).toPrecision(6)}`);
  }
  return lines.join("\n");
}

// ── 8e. Concentration units + selection-driven dosing ───────────────────────
export const CONC_UNIT_TO_UM = { uM: 1.0, nM: 1e-3, pM: 1e-6 };

export function convertConcentration(value, fromUnit, toUnit) {
  if (!(fromUnit in CONC_UNIT_TO_UM)) throw new Error(`unknown unit ${fromUnit}`);
  if (!(toUnit in CONC_UNIT_TO_UM)) throw new Error(`unknown unit ${toUnit}`);
  return value * CONC_UNIT_TO_UM[fromUnit] / CONC_UNIT_TO_UM[toUnit];
}

export function seriesShapeFromSelection(nRows, nCols, orientation = "by_column") {
  if (nRows < 1 || nCols < 1) throw new Error("selection must be at least 1x1");
  if (orientation === "by_row") return [nRows, nCols];
  return [nCols, nRows];
}

export function assignSelection(cells, finalConcsUM, orientation = "by_column") {
  const seen = new Set();
  const cleaned = [];
  for (const [r, c] of cells) {
    const key = `${r},${c}`;
    if (!seen.has(key)) { seen.add(key); cleaned.push([Math.trunc(r), Math.trunc(c)]); }
  }
  cleaned.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  if (!cleaned.length) throw new Error("no cells selected");
  const concs = Array.from(finalConcsUM);
  const n = concs.length;

  const rows = cleaned.map((x) => x[0]);
  const cols = cleaned.map((x) => x[1]);
  const r0 = Math.min(...rows), c0 = Math.min(...cols);
  const h = Math.max(...rows) - r0 + 1;
  const w = Math.max(...cols) - c0 + 1;
  const rectangular = cleaned.length === h * w;

  const doseSpan = orientation === "by_row" ? h : w;
  const warnings = [];
  if (n !== doseSpan) {
    warnings.push(`${n} concentration points but selection spans ${doseSpan} well(s) along the dose axis.`);
  }

  const wells = cleaned.map(([r, c]) => {
    let pt, rep;
    if (orientation === "by_row") { pt = r - r0; rep = c - c0 + 1; }
    else { pt = c - c0; rep = r - r0 + 1; }
    const conc = (pt >= 0 && pt < n) ? concs[pt] : null;
    return { row: r, col: c, point: pt, replicate: rep, finalConcUM: conc };
  });

  return { wells, rows: h, cols: w, origin: [r0, c0], rectangular, orientation, warnings };
}

export function aggregateDosingPlans(plans) {
  const arr = Array.from(plans);
  return {
    nGroups: arr.length,
    totalStockUL: arr.reduce((s, p) => s + p.totalStockUL, 0),
    totalDiluentUL: arr.reduce((s, p) => s + p.totalDiluentUL, 0),
    totalWells: arr.reduce((s, p) => s + p.totalWells, 0),
  };
}

// ── 9b. Extinction-coefficient determination ────────────────────────────────
// Ordinary least-squares fit of a Beer-Lambert dilution series. For known
// concentrations c (mol/L) and absorbances A at fixed wavelength and path
// length L (cm), Beer-Lambert gives A = eps*L*c, so a line A vs c has
// slope = eps*L, hence eps = slope / L.
export function linearRegression(xs, ys) {
  xs = xs.map(Number);
  ys = ys.map(Number);
  const n = xs.length;
  if (n !== ys.length) throw new Error("xs and ys must have equal length");
  if (n < 2) throw new Error("need at least 2 points");
  const mx = xs.reduce((a, b) => a + b, 0) / n;
  const my = ys.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sxx += (xs[i] - mx) ** 2; sxy += (xs[i] - mx) * (ys[i] - my); }
  if (sxx === 0) throw new Error("need at least two distinct x values");
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  let ssTot = 0;
  for (let i = 0; i < n; i++) ssTot += (ys[i] - my) ** 2;
  let rSquared;
  if (ssTot === 0) { rSquared = 1.0; }
  else {
    let ssRes = 0;
    for (let i = 0; i < n; i++) ssRes += (ys[i] - (slope * xs[i] + intercept)) ** 2;
    rSquared = 1.0 - ssRes / ssTot;
  }
  return { slope, intercept, rSquared, n };
}

export function extinctionCoefficient(concentrationsM, absorbances, pathLengthCm = 1.0) {
  if (pathLengthCm <= 0) throw new Error("path length must be positive");
  const fit = linearRegression(concentrationsM, absorbances);
  return {
    eps: fit.slope / pathLengthCm,
    slope: fit.slope,
    intercept: fit.intercept,
    rSquared: fit.rSquared,
    n: fit.n,
  };
}

// ── 9b. Uncertainty / error propagation ─────────────────────────────────────
// UV/Vis DAR error from the two absorbance reads. DAR depends on them only via
// R = aLmax/a280, so sigma_R = |R|*sqrt((sA280/a280)^2+(sALmax/aLmax)^2) and
// sigma_DAR = |dDAR/dR|*sigma_R with the analytic derivative. 95% CI = 1.96 sd.
export function darUVUncertainty(a280, aLmax, eps280Mab, epsLmaxMab, eps280Lp, epsLmaxLp,
                                 sigmaA280 = 0.0, sigmaALmax = 0.0,
                                 sigmaEps280Mab = 0.0, sigmaEpsLmaxMab = 0.0,
                                 sigmaEps280Lp = 0.0, sigmaEpsLmaxLp = 0.0) {
  if (a280 === 0) throw new Error("A280 must be non-zero");
  if (aLmax === 0 && sigmaALmax !== 0.0) throw new Error("A_lmax must be non-zero to propagate its error");
  const r = aLmax / a280;
  const denom = r * eps280Lp - epsLmaxLp;
  if (denom === 0) throw new Error("Degenerate optics: R*eps280_LP == eps_lmax_LP");
  const dar = (epsLmaxMab - r * eps280Mab) / denom;
  const numer = epsLmaxMab - r * eps280Mab;
  const dDarDr = (-eps280Mab * denom - numer * eps280Lp) / (denom * denom);
  let relVar = 0.0;
  if (sigmaA280) relVar += (sigmaA280 / a280) ** 2;
  if (sigmaALmax) relVar += (sigmaALmax / aLmax) ** 2;
  const sigmaR = Math.abs(r) * Math.sqrt(relVar);
  let variance = (dDarDr * sigmaR) ** 2;
  if (sigmaEpsLmaxMab) variance += ((1.0 / denom) * sigmaEpsLmaxMab) ** 2;
  if (sigmaEps280Mab) variance += ((-r / denom) * sigmaEps280Mab) ** 2;
  if (sigmaEps280Lp) variance += ((-numer * r / (denom * denom)) * sigmaEps280Lp) ** 2;
  if (sigmaEpsLmaxLp) variance += ((numer / (denom * denom)) * sigmaEpsLmaxLp) ** 2;
  const sigmaDar = Math.sqrt(variance);
  return {
    dar, sigmaDar,
    ci95Low: dar - 1.96 * sigmaDar,
    ci95High: dar + 1.96 * sigmaDar,
    sigmaR,
    relSigma: dar !== 0 ? sigmaDar / dar : null,
  };
}

// Population mean/variance/SD of a weighted drug-load distribution (k -> weight).
// The SD is the heterogeneity (spread of drug load), not a measurement error.
export function distributionDispersion(weights) {
  const ks = Object.keys(weights);
  let total = 0;
  for (const k of ks) total += weights[k];
  if (total === 0) throw new Error("Total weight must be non-zero");
  let mean = 0;
  for (const k of ks) mean += Number(k) * weights[k];
  mean /= total;
  let variance = 0;
  for (const k of ks) variance += weights[k] * (Number(k) - mean) ** 2;
  variance /= total;
  return { mean, variance, sd: Math.sqrt(variance) };
}

// Exact binomial coefficient C(n,k) via a multiplicative loop (no factorial
// overflow for the small site counts used in ADC conjugation).
function _binomCoeff(n, k) {
  if (k < 0 || k > n) return 0;
  if (k === 0 || k === n) return 1;
  k = Math.min(k, n - k);
  let c = 1;
  for (let i = 0; i < k; i++) c = (c * (n - i)) / (i + 1);
  return Math.round(c);
}

// Predict the full DAR distribution from a binomial site-occupancy model:
// n_sites independent sites each occupied with probability p_site, so drug
// count ~ Binomial(n_sites, p_site). Set p_site directly, or derive it from
// feed_ratio*efficiency/n_sites. Moments via distributionDispersion.
export function predictDarDistribution(nSites, { pSite = null, feedRatio = null, efficiency = 1.0, drugsPerSite = 2 } = {}) {
  // Site-occupancy DAR model keyed to conjugation chemistry.
  // drugsPerSite=2 (default): cysteine/TCEP interchain conjugation — each of
  //   nSites=4 interchain disulfides carries a PAIR of drugs, giving the even
  //   DAR ladder 0,2,4,6,8 seen for real thiol-linked ADCs.
  // drugsPerSite=1: stochastic amine/lysine — single drug per site, classic
  //   smooth binomial over 0,1,2,...
  if (nSites < 1) throw new Error("n_sites must be >= 1");
  if (drugsPerSite < 1) throw new Error("drugs_per_site must be >= 1");
  const maxDar = nSites * drugsPerSite;
  let p = pSite;
  if (p === null || p === undefined) {
    if (feedRatio === null || feedRatio === undefined) throw new Error("provide either pSite or feedRatio");
    if (feedRatio < 0 || efficiency < 0) throw new Error("feedRatio and efficiency must be >= 0");
    p = (feedRatio * efficiency) / maxDar;
  }
  if (!(p >= 0.0 && p <= 1.0)) {
    throw new Error("implied per-site probability outside [0,1]; feedRatio*efficiency exceeds the drug-carrying capacity (n_sites * drugs_per_site)");
  }
  const dist = {};
  for (let j = 0; j <= nSites; j++) {
    dist[drugsPerSite * j] = _binomCoeff(nSites, j) * p ** j * (1.0 - p) ** (nSites - j);
  }
  const disp = distributionDispersion(dist);
  return { pSite: p, distribution: dist, meanDar: disp.mean, variance: disp.variance, sd: disp.sd, drugsPerSite, maxDar };
}

// Reduced LC-MS DAR SD from per-chain load heterogeneity. DAR = 2mean(LC)+2mean(HC);
// treating the 2 light + 2 heavy chains as independent, var = 2var(LC)+2var(HC).
export function darLcmsReducedUncertainty(lightChain, heavyChain) {
  const lc = distributionDispersion(lightChain);
  const hc = distributionDispersion(heavyChain);
  const dar = 2 * lc.mean + 2 * hc.mean;
  const variance = 2 * lc.variance + 2 * hc.variance;
  const sigmaDar = Math.sqrt(variance);
  return {
    dar, sigmaDar,
    lightSd: lc.sd, heavySd: hc.sd,
    ci95Low: dar - 1.96 * sigmaDar,
    ci95High: dar + 1.96 * sigmaDar,
  };
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
  results.push(ok(ellmanFreeThiols(0.5, 0.6, 203000.0, 14150.0, 0.02, 0.01), 11.671557764867941, "Ellman blank-subtracted"));
  results.push(ok(ellmanFreeThiols(0.5, 0.6, 203000.0, 14150.0, 0.0, 0.0), 11.955241460541815, "Ellman zero blank == raw"));
  results.push(ok(scatterAbsorbance(0.05, 320.0, 280.0, 4.0), 0.08529779258642231, "scatter@280"));
  const scc = scatterCorrectedAbsorbance(1.0, 0.5, 0.05, 320.0, 495.0, 4.0);
  results.push(ok(scc.a280Corrected, 0.9147022074135777, "scatter-corrected A280"));
  results.push(ok(scc.aLmaxCorrected, 0.49126728831613614, "scatter-corrected A_lmax"));
  results.push(ok(scc.scatterLmax, 0.008732711683863857, "scatter@lmax"));
  results.push(ok(darHIC({ 0: 5, 1: 15, 2: 45, 3: 25, 4: 10 }), 2.2, "HIC"));
  results.push(ok(darHICCorrected({ 0: 10, 2: 50, 4: 30, 6: 10 }, 203000.0, 5000.0), 2.741596270777426, "HIC corrected"));
  results.push(ok(darHICCorrected({ 0: 5, 1: 15, 2: 45, 3: 25, 4: 10 }, 203000.0, 0.0), 2.2, "HIC corrected reduces to uncorrected"));
  const hicFrac = darSpeciesFractionsCorrected({ 0: 10, 2: 50, 4: 30, 6: 10 }, 203000.0, 5000.0);
  results.push(ok(hicFrac[2], 50.870418158189466, "HIC corrected fraction DAR2"));
  results.push(ok(darLcmsReduced({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 }), 5.2, "LC-MS reduced"));
  results.push(ok(darLcmsIntactCorrected({ 0: 10, 2: 30, 4: 60 }), 3.0, "LC-MS intact uncorrected"));
  results.push(ok(darLcmsIntactCorrected({ 0: 10, 2: 30, 4: 60 }, { 0: 1.0, 2: 0.8, 4: 0.6 }), 3.2203389830508473, "LC-MS intact response-corrected"));
  results.push(ok(darLcmsReducedCorrected({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 }), 5.2, "LC-MS reduced corrected == raw when no response"));
  results.push(ok(darLcmsReducedCorrected({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 }, { 0: 1.0, 1: 0.9 }, { 0: 1.0, 1: 0.9, 2: 0.8 }), 5.285460807848867, "LC-MS reduced response-corrected"));
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
  // extinction-coefficient regression
  const CONC = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6];
  const AB280 = CONC.map((c) => 0.01 + 50000.0 * c);
  const fit = linearRegression(CONC, AB280);
  ok(fit.slope, 50000.0, "linreg slope");
  ok(fit.intercept, 0.01, "linreg intercept");
  ok(fit.rSquared, 1.0, "linreg r2");
  const ec = extinctionCoefficient(CONC, AB280, 1.0);
  ok(ec.eps, 50000.0, "eps L=1");
  ok(extinctionCoefficient(CONC, AB280, 0.5).eps, 100000.0, "eps L=0.5");
  results.push("PASS extinction coefficient (regression)");
  // uncertainty / error propagation
  const uu = darUVUncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP, 0.01, 0.01);
  ok(uu.dar, 2.4140809554780795, "DAR unc point");
  ok(uu.sigmaDar, 0.03370113601940267, "DAR unc sigma");
  const uue = darUVUncertainty(A280, A_LMAX, EPS280_MAB, 0.0, EPS280_LP, EPS_LMAX_LP,
                               0.01, 0.01, 2030.0, 0.0, 364.35, 481.2);
  ok(uue.sigmaDar, 0.13795625467477654, "DAR unc sigma (with eps terms)");
  const dd = distributionDispersion({ 0: 5, 1: 15, 2: 45, 3: 25, 4: 10 });
  ok(dd.sd, 0.9797958971132712, "HIC dispersion sd");
  const pd = predictDarDistribution(4, { feedRatio: 2.5, efficiency: 0.8 });
  results.push(ok(pd.pSite, 0.25, "DAR dist p_site (cys)"));
  results.push(ok(pd.meanDar, 2.0, "DAR dist mean (cys)"));
  results.push(ok(pd.sd, 1.7320508075688772, "DAR dist sd (cys)"));
  results.push(ok(pd.distribution[2], 0.421875, "DAR dist P(DAR=2) (cys)"));
  results.push(ok(pd.distribution[1] === undefined ? 1 : 0, 1, "DAR dist has no odd species (cys)"));
  const pd2 = predictDarDistribution(4, { pSite: 0.5 });
  results.push(ok(pd2.meanDar, 4.0, "DAR dist via p_site (cys)"));
  const pdLys = predictDarDistribution(8, { pSite: 0.5, drugsPerSite: 1 });
  results.push(ok(pdLys.meanDar, 4.0, "DAR dist lysine mean"));
  results.push(ok(pdLys.distribution[3], 0.21875, "DAR dist lysine P(DAR=3)"));
  const lu = darLcmsReducedUncertainty({ 0: 10, 1: 90 }, { 0: 5, 1: 20, 2: 75 });
  ok(lu.sigmaDar, 0.8944271909999159, "LC-MS reduced sigma");
  results.push("PASS uncertainty (UV / HIC / LC-MS)");
  // in vitro dosing / serial dilution (spec §8b golden fixture)
  const dp = planSerialDilution({ stockUM: 1000, topUM: 10, fold: 4, nPoints: 4,
    replicates: 3, wellVolumeUL: 100.0, doseMode: "spike", doseFactor: 10.0, overage: 1.0 });
  ok(dp.finalConcsUM[3], 0.15625, "dosing final conc[3]");
  ok(dp.workingConcsUM[0], 100.0, "dosing working conc[0]");
  ok(dp.vAddPerWellUL, 100.0 / 9.0, "dosing spike volume");
  ok(dp.vUsePerTubeUL, 100.0 / 3.0, "dosing dispense per tube");
  ok(dp.vTransferUL, 100.0 / 9.0, "dosing transfer volume");
  ok(dp.totalStockUL, 4.444444444444444, "dosing tube-1 stock");
  ok(dp.tubes[1].diluentUL, 100.0 / 3.0, "dosing tube-2 diluent");
  ok(dp.totalWells, 12, "dosing total wells");
  ok(plateWellVolume("96"), 100.0, "plate 96 well volume");
  const dpr = planSerialDilution({ stockUM: 500, topUM: 50, fold: 2, nPoints: 6,
    replicates: 4, wellVolumeUL: 200.0, doseMode: "replace", overage: 1.0 });
  ok(dpr.workingConcsUM[0], 50.0, "replace working==final");
  ok(dpr.vAddPerWellUL, 200.0, "replace add == well volume");
  results.push("PASS in vitro dosing (serial dilution)");
  const pm = plateMap([10, 3.333, 1.111, 0.37], 3, "96", "by_column");
  if (pm.orientation !== "by_column") throw new Error("plate map orientation");
  if (pm.used !== 12) throw new Error("plate map used count");
  if (pm.wells[0].rowLabel + pm.wells[0].colLabel !== "A1") throw new Error("plate map A1");
  if (pm.wells[1].rowLabel + pm.wells[1].colLabel !== "B1") throw new Error("plate map B1 replicate");
  if (pm.wells[3].rowLabel + pm.wells[3].colLabel !== "A2") throw new Error("plate map A2 next dose");
  if (rowLabel(26) !== "AA") throw new Error("row label AA");
  const pmSeq = plateMap([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 3, "96", "by_column");
  if (pmSeq.orientation !== "sequential" || pmSeq.warnings.length === 0) throw new Error("plate map fallback");
  results.push("PASS plate map (layout + labels)");
  if (convertConcentration(5, "uM", "nM") !== 5000) throw new Error("unit uM->nM");
  if (convertConcentration(250, "nM", "uM") !== 0.25) throw new Error("unit nM->uM");
  const [sp, sr] = seriesShapeFromSelection(3, 4, "by_row");
  if (sp !== 3 || sr !== 4) throw new Error("series shape by_row");
  const [sp2, sr2] = seriesShapeFromSelection(3, 4, "by_column");
  if (sp2 !== 4 || sr2 !== 3) throw new Error("series shape by_column");
  // dose decreasing down rows (top->bottom), reps across columns
  const cells = [];
  for (let r = 2; r < 6; r++) for (let c = 5; c < 8; c++) cells.push([r, c]);
  const asg = assignSelection(cells, [10, 2.5, 0.625, 0.15625], "by_row");
  if (!asg.rectangular || asg.rows !== 4 || asg.cols !== 3) throw new Error("selection bbox");
  const wTop = asg.wells.find((w) => w.row === 2 && w.col === 5);
  const wBot = asg.wells.find((w) => w.row === 5 && w.col === 5);
  if (wTop.point !== 0 || wTop.finalConcUM !== 10) throw new Error("selection top dose");
  if (wBot.point !== 3 || wBot.finalConcUM !== 0.15625) throw new Error("selection bottom dose");
  if (wTop.replicate !== 1) throw new Error("selection replicate");
  results.push("PASS units + selection assignment");
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
