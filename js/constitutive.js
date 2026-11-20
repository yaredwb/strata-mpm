/*
 * STRATA — Computational Geomechanics Constitutive Models & Continuum Analytics
 * Exact 2x2 analytic SVD, Hencky logarithmic finite-strain hyperelasticity,
 * non-associated plane-strain Mohr-Coulomb return mapping, Nordic strain-softening
 * sensitive quick-clay Tresca return mapping, effective-stress liquefaction coupling,
 * and analytical Prandtl / Terzaghi / Rankine / Mohr-circle tools.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory();
  } else {
    root.Strata = root.Strata || {};
    root.Strata.Constitutive = factory();
  }
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  const PI = Math.PI;
  const DEG2RAD = PI / 180;

  // Material model enum
  const MAT = {
    SAND_MC: 0,       // Mohr-Coulomb frictional-cohesive granular soil
    QUICK_CLAY: 1,    // Strain-softening sensitive marine quick clay (Tresca s_u(kappa))
    LIQ_SAND: 2,      // Effective-stress liquefiable saturated loose sand
    CRUST_CLAY: 3,    // Stiff dry weathered clay crust (higher c', mild softening)
    RIGID_STRUCT: 4,  // Stiff elastic concrete/steel structure (footing / sheet-pile wall)
    BEDROCK: 5        // Immobile stiff moraine / bedrock base
  };

  function clamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
  }

  // Exact analytic 2x2 SVD: F = U * diag(s1, s2) * V^T, det(U)=det(V)=+1
  // F = [[f00, f01], [f10, f11]]
  function svd2x2(f00, f01, f10, f11, out) {
    const E = 0.5 * (f00 + f11);
    const F = 0.5 * (f00 - f11);
    const G = 0.5 * (f01 + f10);
    const H = 0.5 * (f10 - f01);

    const Q = Math.hypot(E, H);
    const R = Math.hypot(F, G);
    let s1 = Q + R;
    let s2 = Q - R;

    const a1 = Math.atan2(G, F);
    const a2 = Math.atan2(H, E);
    const thU = 0.5 * (a1 + a2);
    const thV = 0.5 * (a1 - a2);

    const cu = Math.cos(thU);
    const su = Math.sin(thU);
    const cv = Math.cos(thV);
    const sv = Math.sin(thV);

    // Guard against inverted elements in extreme impact
    if (s1 < 0.2) s1 = 0.2;
    if (s2 < 0.2) s2 = 0.2;
    if (s1 > 4.0) s1 = 4.0;
    if (s2 > 4.0) s2 = 4.0;

    out.u00 = cu;  out.u01 = -su;
    out.u10 = su;  out.u11 = cu;
    out.v00 = cv;  out.v01 = -sv;
    out.v10 = sv;  out.v11 = cv;
    out.s1 = s1;   out.s2 = s2;
    return out;
  }

  // Reconstruct 2x2 matrix A = U * diag(d1, d2) * V^T
  function mulUDVT(u00, u01, u10, u11, d1, d2, v00, v01, v10, v11, out) {
    const a00 = u00 * d1;
    const a01 = u01 * d2;
    const a10 = u10 * d1;
    const a11 = u11 * d2;
    out.m00 = a00 * v00 + a01 * v01;
    out.m01 = a00 * v10 + a01 * v11;
    out.m10 = a10 * v00 + a11 * v01;
    out.m11 = a10 * v10 + a11 * v11;
    return out;
  }

  // Reconstruct symmetric Cauchy/Kirchhoff tensor T = U * diag(t1, t2) * U^T
  function symUDUT(u00, u01, u10, u11, t1, t2, out) {
    out.t00 = u00 * u00 * t1 + u01 * u01 * t2;
    out.t11 = u10 * u10 * t1 + u11 * u11 * t2;
    out.t01 = u00 * u10 * t1 + u01 * u11 * t2;
    return out;
  }

  // Sensitive quick-clay exponential strain-softening law:
  // s_u(kappa) = s_ur + (s_u0 - s_ur) * exp(-kappa / kappa95)
  function quickClayStrength(kappa, su0, sensitivity, kappa95) {
    const st = Math.max(1.0, sensitivity);
    const sur = su0 / st;
    const k95 = Math.max(1e-4, kappa95);
    return sur + (su0 - sur) * Math.exp(-Math.max(0, kappa) / k95);
  }

  // Plane-strain Lamé parameters from Young's modulus E and Poisson's ratio nu
  function lameParams(E, nu) {
    const mu = E / (2 * (1 + nu));
    const lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    const k2d = lambda + mu; // 2D plane-strain bulk modulus: p_2D = -K_2D * (e1 + e2)
    return { mu, lambda, k2d };
  }

  // Exact plane-strain non-associated Mohr-Coulomb return mapping in Hencky log-strain space.
  // Supports effective stress via excess pore pressure uPore >= 0 (liquefaction coupling).
  // Inputs:
  //   eps1, eps2: trial principal Hencky strains ln(s1), ln(s2)
  //   mu, k2d: Lamé shear modulus & 2D plane-strain bulk modulus (Pa)
  //   c: cohesion (Pa), phiDeg: friction angle (deg), psiDeg: dilation angle (deg)
  //   uPore: excess pore-water pressure (Pa) reducing effective mean stress p' = p_tot - uPore
  // Outputs written into `res`:
  //   eps1New, eps2New, tau1Tot, tau2Tot, pEff, qShear, dLambda, yielded
  function returnMapMohrCoulomb(eps1, eps2, mu, k2d, c, phiDeg, psiDeg, uPore, res) {
    const sinPhi = Math.sin(clamp(phiDeg, 0, 55) * DEG2RAD);
    const cosPhi = Math.cos(clamp(phiDeg, 0, 55) * DEG2RAD);
    const sinPsi = Math.sin(clamp(psiDeg, 0, phiDeg) * DEG2RAD);

    const evTr = eps1 + eps2;
    const edTr = eps1 - eps2;
    const pTotTr = -k2d * evTr;        // Total in-plane mean compressive stress (Pa, + in compression)
    const qTr = mu * Math.abs(edTr);   // In-plane max shear stress (tau1 - tau2)/2 >= 0
    const sgn = edTr >= 0 ? 1 : -1;

    // Effective confining stress (cannot go below tensile apex)
    const pEffTr = pTotTr - uPore;
    const cCos = c * cosPhi;
    const fTr = qTr - pEffTr * sinPhi - cCos;

    if (fTr <= 0 && (sinPhi <= 1e-6 || pEffTr >= -(cCos / Math.max(sinPhi, 1e-6)))) {
      // Purely elastic step
      res.eps1New = eps1;
      res.eps2New = eps2;
      res.tau1Tot = -pTotTr + sgn * qTr;
      res.tau2Tot = -pTotTr - sgn * qTr;
      res.pEff = pEffTr;
      res.qShear = qTr;
      res.dLambda = 0;
      res.yielded = false;
      return res;
    }

    // Tensile apex check (cohesionless sand apex at p' = 0; cohesive soil at p' = -c*cot(phi))
    const pApex = sinPhi > 1e-5 ? -(cCos / sinPhi) : -Math.max(c * 0.6, 500);

    // Smooth cone return mapping
    const denom = mu + k2d * sinPhi * sinPsi;
    let dLambda = fTr / denom;
    let qNew = qTr - mu * dLambda;
    let pEffNew = pEffTr + k2d * dLambda * sinPsi;
    let pTotNew = pTotTr + k2d * dLambda * sinPsi;

    if (qNew < 0 || pEffNew < pApex) {
      // Apex projection (tensile detachment / zero effective confinement liquefaction)
      qNew = 0;
      pEffNew = pApex;
      // Preserve existing compressive total stress without artificially injecting upward free-surface splash pressure
      pTotNew = uPore > 0 ? Math.max(0, pTotTr) : pApex;
      dLambda = Math.max(0, qTr / mu);
    }

    const evNew = -pTotNew / k2d;
    const edNew = sgn * (qNew / mu);
    res.eps1New = 0.5 * (evNew + edNew);
    res.eps2New = 0.5 * (evNew - edNew);
    res.tau1Tot = -pTotNew + sgn * qNew;
    res.tau2Tot = -pTotNew - sgn * qNew;
    res.pEff = pEffNew;
    res.qShear = qNew;
    res.dLambda = Math.max(0, dLambda);
    res.yielded = true;
    return res;
  }

  // Exact strain-softening Tresca (undrained sensitive quick clay) return mapping.
  // Solves q_tr - mu * dLambda - s_u(kappa + dLambda) = 0 via Newton-Raphson.
  function returnMapQuickClay(eps1, eps2, mu, k2d, kappa, su0, sensitivity, kappa95, res) {
    const evTr = eps1 + eps2;
    const edTr = eps1 - eps2;
    const pTotTr = -k2d * evTr;
    const qTr = mu * Math.abs(edTr);
    const sgn = edTr >= 0 ? 1 : -1;

    const suCur = quickClayStrength(kappa, su0, sensitivity, kappa95);
    const tensionCap = -0.35 * suCur; // Small tensile strength allows crisp backscarp tension cracks

    let pNew = pTotTr;
    let tensionCut = false;
    if (pNew < tensionCap) {
      pNew = tensionCap;
      tensionCut = true;
    }

    const fTr = qTr - suCur;
    if (fTr <= 0) {
      const evNew = -pNew / k2d;
      res.eps1New = 0.5 * (evNew + edTr);
      res.eps2New = 0.5 * (evNew - edTr);
      res.tau1Tot = -pNew + sgn * qTr;
      res.tau2Tot = -pNew - sgn * qTr;
      res.pEff = pNew;
      res.qShear = qTr;
      res.dLambda = 0;
      res.suNew = suCur;
      res.yielded = tensionCut;
      return res;
    }

    // 4-iteration exact Newton-Raphson for exponential softening
    const sur = su0 / Math.max(1.0, sensitivity);
    const k95 = Math.max(1e-4, kappa95);
    let dl = fTr / mu;
    for (let it = 0; it < 4; it++) {
      const eTerm = (su0 - sur) * Math.exp(-(kappa + dl) / k95);
      const suAt = sur + eTerm;
      const gVal = qTr - mu * dl - suAt;
      const gPrime = -mu + eTerm / k95; // Always < 0 since mu >> su0 / k95
      const step = gVal / Math.min(gPrime, -mu * 0.1);
      dl = Math.max(0, dl - step);
    }

    const suNew = quickClayStrength(kappa + dl, su0, sensitivity, kappa95);
    const qNew = Math.max(0, suNew);
    const evNew = -pNew / k2d;
    const edNew = sgn * (qNew / mu);

    res.eps1New = 0.5 * (evNew + edNew);
    res.eps2New = 0.5 * (evNew - edNew);
    res.tau1Tot = -pNew + sgn * qNew;
    res.tau2Tot = -pNew - sgn * qNew;
    res.pEff = pNew;
    res.qShear = qNew;
    res.dLambda = dl;
    res.suNew = suNew;
    res.yielded = true;
    return res;
  }

  // Exact analytical Prandtl / Vesic bearing capacity factors (Nc, Nq, Ngamma)
  function bearingCapacityFactors(phiDeg) {
    if (phiDeg < 1e-4) {
      return { Nc: 2 + PI, Nq: 1.0, Ngamma: 0.0 };
    }
    const phi = phiDeg * DEG2RAD;
    const tanPhi = Math.tan(phi);
    const Nq = Math.exp(PI * tanPhi) * Math.pow(Math.tan(PI / 4 + 0.5 * phi), 2);
    const Nc = (Nq - 1) / tanPhi;
    const Ngamma = 2 * (Nq + 1) * tanPhi; // Vesic (1973)
    return { Nc, Nq, Ngamma };
  }

  // Exact analytical ultimate strip footing bearing capacity q_u (Pa)
  function prandtlBearingCapacity(c, phiDeg, gamma, B, qSurcharge) {
    const f = bearingCapacityFactors(phiDeg);
    return c * f.Nc + (qSurcharge || 0) * f.Nq + 0.5 * gamma * B * f.Ngamma;
  }

  // Compute Mohr's circle parameters from 2D effective stress tensor (sxx, syy, txy; compression positive!)
  function mohrCircleFromStress(sxxComp, syyComp, txy) {
    const center = 0.5 * (sxxComp + syyComp);
    const radius = Math.hypot(0.5 * (sxxComp - syyComp), txy);
    const sigma1 = center + radius;
    const sigma3 = center - radius;
    const angleRad = 0.5 * Math.atan2(2 * txy, sxxComp - syyComp);
    const mobSinPhi = center > 100 ? clamp(radius / center, 0, 1) : 0;
    const mobPhiDeg = Math.asin(mobSinPhi) / DEG2RAD;
    return { center, radius, sigma1, sigma3, angleRad, mobPhiDeg };
  }

  // Compute initial elastic Hencky deformation gradient F0 for K0 geostatic equilibrium:
  // sigma_yy = -sigmaV (vertical compression), sigma_xx = -K0 * sigmaV (horizontal compression)
  function geostaticInitialF(sigmaV, K0, mu, k2d) {
    const sxx = -K0 * sigmaV;
    const syy = -sigmaV;
    const p = -0.5 * (sxx + syy);
    const q = 0.5 * (sxx - syy); // positive when K0 < 1
    const ev = -p / k2d;
    const ed = q / mu;
    const e1 = 0.5 * (ev + ed);
    const e2 = 0.5 * (ev - ed);
    return { f00: Math.exp(e1), f11: Math.exp(e2) };
  }

  return {
    MAT,
    svd2x2,
    mulUDVT,
    symUDUT,
    quickClayStrength,
    lameParams,
    returnMapMohrCoulomb,
    returnMapQuickClay,
    bearingCapacityFactors,
    prandtlBearingCapacity,
    mohrCircleFromStress,
    geostaticInitialF
  };
});
