// Node.js quantitative verification suite for STRATA constitutive models & MLS-MPM solver
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const C = require('../js/constitutive.js');
const MPMSolver = require('../js/mpm-solver.js');

function approx(got, exp, tol, label) {
  const err = Math.abs(got - exp);
  assert.ok(err <= tol, `${label}: got ${got}, expected ${exp} (err=${err.toExponential(3)} > ${tol})`);
  console.log(`[PASS] ${label.padEnd(46)} got=${got.toFixed(6)} err=${err.toExponential(2)}`);
}

console.log('=== STRATA JS Geomechanics Verification Suite ===');

// 1. Exact analytic 2x2 SVD reconstruction & orthogonality
const svd = C.svd2x2(1.22, -0.31, 0.19, 0.86, {});
const rec = C.mulUDVT(svd.u00, svd.u01, svd.u10, svd.u11, svd.s1, svd.s2, svd.v00, svd.v01, svd.v10, svd.v11, {});
approx(Math.max(Math.abs(rec.m00 - 1.22), Math.abs(rec.m01 + 0.31), Math.abs(rec.m10 - 0.19), Math.abs(rec.m11 - 0.86)), 0, 1e-13, '2x2 SVD reconstruction ||U*S*V^T - F||');
approx(Math.abs(svd.u00 * svd.u01 + svd.u10 * svd.u11), 0, 1e-14, '2x2 SVD U orthogonality');

// 2. Non-associated Mohr-Coulomb return mapping consistency
const { mu, k2d } = C.lameParams(10e6, 0.3);
const retMC = C.returnMapMohrCoulomb(0.008, -0.024, mu, k2d, 10000, 30, 5, 0, {});
assert.equal(retMC.yielded, true);
const fMC = retMC.qShear - retMC.pEff * Math.sin((30 * Math.PI) / 180) - 10000 * Math.cos((30 * Math.PI) / 180);
approx(fMC, 0, 1e-7, 'Mohr-Coulomb return mapping |f(sig)| (Pa)');

// 3. Liquefaction effective-stress strength collapse when uPore >= pTot (115.4 kPa -> r_u = 1.0)
const retLiq = C.returnMapMohrCoulomb(-0.004, -0.008, mu, k2d, 0, 32, 0, 120000, {});
approx(retLiq.qShear, 0, 1e-6, 'Liquefied sand shear capacity collapse (q->0)');

// 4. Strain-softening Nordic quick clay exact Newton-Raphson consistency
const retQC = C.returnMapQuickClay(0.014, -0.024, mu, k2d, 0.04, 26000, 32, 0.16, {});
const suExp = C.quickClayStrength(0.04 + retQC.dLambda, 26000, 32, 0.16);
approx(retQC.qShear - suExp, 0, 1e-6, 'Quick-clay softening |q - s_u(k)| (Pa)');

// 5. Analytical Prandtl / Vesic bearing capacity factors
const fUndrained = C.bearingCapacityFactors(0);
approx(fUndrained.Nc, 2 + Math.PI, 1e-12, 'Prandtl undrained Nc = 2 + pi (5.14159)');
const fPhi30 = C.bearingCapacityFactors(30);
approx(fPhi30.Nq, 18.401122, 1e-4, 'Prandtl Nq (phi=30 deg)');
approx(fPhi30.Nc, 30.139625, 1e-4, 'Prandtl Nc (phi=30 deg)');

// 6. K0 geostatic initialization & solver stability test
const solver = new MPMSolver(12000);
solver.loadScenario('quickclay');
assert.ok(solver.count > 4000, `Expected >4000 material points, got ${solver.count}`);
for (let s = 0; s < 120; s++) solver.step();
const intactMetrics = solver.computeMetrics();
assert.ok(Number.isFinite(intactMetrics.ke) && intactMetrics.maxSpeed < 0.45, `Intact slope unstable: maxSpeed=${intactMetrics.maxSpeed}`);
console.log(`[PASS] ${'JS MLS-MPM intact K0 slope equilibrium'.padEnd(46)} MPs=${solver.count} max|v|=${intactMetrics.maxSpeed.toFixed(3)} m/s`);

// 7. Toe excavation triggers retrogressive quick-clay softening & mobilization
solver.triggerExcavationPreset();
for (let s = 0; s < 650; s++) solver.step();
const postExcav = solver.computeMetrics();
assert.ok(postExcav.maxPlastic > 0.25, `Expected plastic shear band after toe cut, got maxPlastic=${postExcav.maxPlastic}`);
console.log(`[PASS] ${'Toe cut triggers quick-clay localization'.padEnd(46)} maxEpsP=${postExcav.maxPlastic.toFixed(3)} remolded=${postExcav.remoldedPct.toFixed(1)}%`);

// 8. Mohr circle probe works cleanly inside soil continuum
const probe = solver.probeAt(10.0, 7.0);
assert.ok(probe && probe.mohr.sigma1 >= probe.mohr.sigma3, 'Mohr circle probe returned invalid principal stresses');
console.log(`[PASS] ${'Live Mohr circle continuum probe'.padEnd(46)} sig1=${(probe.mohr.sigma1 / 1e3).toFixed(1)}kPa sig3=${(probe.mohr.sigma3 / 1e3).toFixed(1)}kPa`);

console.log('>>> ALL STRATA JS GEOMECHANICS TESTS PASSED! <<<');
