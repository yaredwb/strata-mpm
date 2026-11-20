/*
 * STRATA — Large-Deformation Elastoplastic MLS-MPM + APIC Geomechanics Solver
 * Features:
 *   - Quadratic B-Spline MLS-MPM with APIC affine transfer (zero cell-crossing noise)
 *   - Structure-of-Arrays Float64/Float32 buffers for 60 FPS in-browser simulation
 *   - Coupled excess pore-water pressure generation & 2D Eulerian Darcy consolidation
 *   - Exact K0 geostatic stress initialization
 *   - 6 Computational Geomechanics Scenarios + Live Interactive Excavation/Pour/Footing/Seismic tools
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('./constitutive.js'));
  } else {
    root.Strata = root.Strata || {};
    root.Strata.MPMSolver = factory(root.Strata.Constitutive);
  }
})(typeof globalThis !== 'undefined' ? globalThis : this, function (C) {
  'use strict';

  const MAT = C.MAT;

  const SCENARIOS = {
    quickclay: {
      id: 'quickclay',
      name: 'Nordic Quick-Clay Retrogressive Landslide',
      domainW: 44.0,
      domainH: 18.0,
      nx: 110,
      ny: 45,
      gravity: 9.81,
      su0: 42000,         // 42 kPa intact undrained shear strength (FS ~ 1.6 intact!)
      sensitivity: 40.0,  // St = 40 -> remolded s_ur = 1.05 kPa!
      kappa95: 0.12,
      phiDeg: 32.0,
      psiDeg: 2.0,
      cohesion: 24000,    // Weathered dry clay crust cohesion
      youngs: 8.5e6,      // 8.5 MPa
      poisson: 0.34,
      liqBeta: 0.0,
      darcyCv: 0.0,
      defaultField: 0,
      autoTriggerHint: 'Click "Excavate Toe" or drag with the Excavator on the slope toe to trigger a retrogressive quick-clay failure cascade.'
    },
    liquefaction: {
      id: 'liquefaction',
      name: 'Seismic Liquefaction & Structure Settlement',
      domainW: 32.0,
      domainH: 16.0,
      nx: 96,
      ny: 48,
      gravity: 9.81,
      su0: 25000,
      sensitivity: 10.0,
      kappa95: 0.20,
      phiDeg: 30.0,
      psiDeg: 0.0,
      cohesion: 200,
      youngs: 7.2e6,
      poisson: 0.30,
      liqBeta: 34.0,      // Excess pore-pressure generation coupling
      darcyCv: 0.035,     // Consolidation diffusivity (m^2/s)
      defaultField: 4,    // Pore pressure ratio r_u
      autoTriggerHint: 'Press "Trigger Earthquake" to apply cyclic laminar-box base shaking — watch r_u climb to 1.0 and the heavy footing tilt & settle into liquefied sand.'
    },
    footing: {
      id: 'footing',
      name: 'Strip Footing Bearing Capacity (Prandtl–Terzaghi)',
      domainW: 24.0,
      domainH: 12.0,
      nx: 96,
      ny: 48,
      gravity: 9.81,
      su0: 22000,
      sensitivity: 1.0,
      kappa95: 0.25,
      phiDeg: 20.0,
      psiDeg: 4.0,
      cohesion: 12000,    // 12 kPa
      youngs: 9.5e6,
      poisson: 0.30,
      liqBeta: 0.0,
      darcyCv: 0.0,
      footingWidth: 3.6,
      defaultField: 1,    // Plastic deviatoric strain (highlights Prandtl slip lines!)
      autoTriggerHint: 'Constant-velocity rigid footing indentation — watch Prandtl logarithmic-spiral slip wedges nucleate and compare load q vs analytical q_u.'
    },
    column: {
      id: 'column',
      name: 'Granular Column Collapse (Lajeunesse Benchmark)',
      domainW: 28.0,
      domainH: 12.0,
      nx: 112,
      ny: 48,
      gravity: 9.81,
      su0: 25000,
      sensitivity: 1.0,
      kappa95: 0.25,
      phiDeg: 32.0,
      psiDeg: 2.0,
      cohesion: 0.0,      // Pure dry cohesionless sand
      youngs: 8.5e6,
      poisson: 0.28,
      liqBeta: 0.0,
      darcyCv: 0.0,
      colR0: 3.6,
      colH0: 6.2,
      defaultField: 1,
      autoTriggerHint: 'Instantaneous dam-break release of cohesionless sand (a = H0/R0 = 1.72) — observe static dead-zone wedge and experimental runout scaling.'
    },
    retaining: {
      id: 'retaining',
      name: 'Cantilever Sheet-Pile Retaining Wall',
      domainW: 32.0,
      domainH: 15.0,
      nx: 96,
      ny: 45,
      gravity: 9.81,
      su0: 25000,
      sensitivity: 1.0,
      kappa95: 0.25,
      phiDeg: 30.0,
      psiDeg: 3.0,
      cohesion: 1500,
      youngs: 9.0e6,
      poisson: 0.30,
      liqBeta: 0.0,
      darcyCv: 0.0,
      defaultField: 1,
      autoTriggerHint: 'Use the Excavator to dig in front of the sheet-pile wall (or click "Excavate Front") and watch the Rankine active wedge mobilize at 45° + phi/2.'
    },
    biaxial: {
      id: 'biaxial',
      name: 'Biaxial Compression Shear-Band Localization',
      domainW: 12.0,
      domainH: 12.0,
      nx: 72,
      ny: 72,
      gravity: 0.0,       // Weightless element test with confining pressure
      su0: 35000,
      sensitivity: 2.2,
      kappa95: 0.14,
      phiDeg: 28.0,
      psiDeg: 5.0,
      cohesion: 16000,
      youngs: 11.0e6,
      poisson: 0.28,
      liqBeta: 0.0,
      darcyCv: 0.0,
      defaultField: 1,
      autoTriggerHint: 'Displacement-controlled biaxial compression test — watch conjugate X-shaped shear bands spontaneously localize at 45° + phi/2.'
    }
  };

  class MPMSolver {
    constructor(maxParticles) {
      this.maxParticles = maxParticles || 16000;
      this.count = 0;

      const N = this.maxParticles;
      this.x = new Float64Array(N);
      this.y = new Float64Array(N);
      this.vx = new Float64Array(N);
      this.vy = new Float64Array(N);
      this.mass = new Float64Array(N);
      this.vol0 = new Float64Array(N);

      // Elastic deformation gradient F (2x2)
      this.f00 = new Float64Array(N);
      this.f01 = new Float64Array(N);
      this.f10 = new Float64Array(N);
      this.f11 = new Float64Array(N);

      // APIC affine velocity gradient C (2x2)
      this.c00 = new Float64Array(N);
      this.c01 = new Float64Array(N);
      this.c10 = new Float64Array(N);
      this.c11 = new Float64Array(N);

      // Stress tensor (compression positive for reporting: sxx, syy, txy) & state variables
      this.sxx = new Float32Array(N);
      this.syy = new Float32Array(N);
      this.txy = new Float32Array(N);
      this.pEff = new Float32Array(N);
      this.qShear = new Float32Array(N);
      this.kappa = new Float32Array(N);     // Equivalent plastic deviatoric strain
      this.suCur = new Float32Array(N);     // Current undrained strength (Pa)
      this.uPore = new Float32Array(N);     // Excess pore-water pressure (Pa)
      this.sigV0 = new Float32Array(N);     // Initial vertical effective overburden (Pa)
      this.ru = new Float32Array(N);        // Excess pore-pressure ratio u / sigV0
      this.mat = new Uint8Array(N);         // Material ID
      this.strataBand = new Uint8Array(N);  // Visual geological layer band index
      this.active = new Uint8Array(N);      // 1 = live, 0 = excavated

      // Reusable scratch objects (zero GC allocations inside substep loop!)
      this._svd = { u00: 1, u01: 0, u10: 0, u11: 1, v00: 1, v01: 0, v10: 0, v11: 1, s1: 1, s2: 1 };
      this._ret = { eps1New: 0, eps2New: 0, tau1Tot: 0, tau2Tot: 0, pEff: 0, qShear: 0, dLambda: 0, suNew: 0, yielded: false };
      this._matTmp = { m00: 1, m01: 0, m10: 0, m11: 1 };
      this._symTmp = { t00: 0, t11: 0, t01: 0 };
      this._wx = new Float64Array(3);
      this._wy = new Float64Array(3);

      // Simulation & seismic parameters
      this.params = {
        phiDeg: 32,
        psiDeg: 2,
        cohesion: 12000,
        su0: 26000,
        sensitivity: 32,
        kappa95: 0.16,
        youngs: 8.0e6,
        poisson: 0.32,
        gravity: 9.81,
        liqBeta: 14.5,
        darcyCv: 0.18,
        seismicActive: false,
        seismicPgaG: 0.35,
        seismicFreqHz: 2.2,
        seismicDuration: 5.5,
        seismicStartTime: 0
      };

      // Domain & Eulerian background grid
      this.scenarioId = 'quickclay';
      this.domainW = 44.0;
      this.domainH = 18.0;
      this.nx = 110;
      this.ny = 45;
      this.dx = 0.4;
      this.invDx = 2.5;
      this.dt = 0.0012;
      this.time = 0;
      this.stepCount = 0;

      // Driven footing / platen state
      this.platenY = 0;
      this.platenY0 = 0;
      this.platenXMin = 0;
      this.platenXMax = 0;
      this.platenVel = 0;
      this.platenForceY = 0;
      this.platenArea = 1;
      this.biaxialH0 = 6.0;
      this.biaxialConfine = 35000;

      // Telemetry history ring buffer for live engineering charts
      this.history = [];
      this.maxHistory = 260;
      this.lastHistoryTime = -1;

      this.initGrid(110, 45);
      this.loadScenario('quickclay');
    }

    initGrid(nx, ny) {
      this.nx = nx;
      this.ny = ny;
      this.dx = this.domainW / nx;
      this.invDx = 1 / this.dx;
      const gn = (nx + 1) * (ny + 1);
      this.gMass = new Float64Array(gn);
      this.gVx = new Float64Array(gn);
      this.gVy = new Float64Array(gn);
      this.gFixed = new Uint8Array(gn);
      this.gPoreNum = new Float64Array(gn);
      this.gPore = new Float64Array(gn);
      this.gPoreNext = new Float64Array(gn);
    }

    computeStableDt() {
      const E = this.params.youngs * 2.2; // Account for stiffer structural/crust elements
      const nu = Math.min(0.45, this.params.poisson);
      const { mu, lambda } = C.lameParams(E, nu);
      const cp = Math.sqrt((lambda + 2 * mu) / 1850);
      return Math.min(0.0018, (0.34 * this.dx) / cp);
    }

    clearParticles() {
      this.count = 0;
      this.time = 0;
      this.stepCount = 0;
      this.history = [];
      this.lastHistoryTime = -1;
      this.params.seismicActive = false;
    }

    addParticle(px, py, matId, rho, vol0, sigmaV, K0, strataBand) {
      if (this.count >= this.maxParticles) return -1;
      const i = this.count++;
      this.x[i] = px;
      this.y[i] = py;
      this.vx[i] = 0;
      this.vy[i] = 0;
      this.mass[i] = rho * vol0;
      this.vol0[i] = vol0;
      this.mat[i] = matId;
      this.strataBand[i] = strataBand || 0;
      this.active[i] = 1;

      this.c00[i] = 0; this.c01[i] = 0;
      this.c10[i] = 0; this.c11[i] = 0;
      this.kappa[i] = 0;
      this.uPore[i] = 0;
      this.ru[i] = 0;
      this.sigV0[i] = Math.max(sigmaV, 3500);
      this.suCur[i] = matId === MAT.QUICK_CLAY ? this.params.su0 : this.params.cohesion;

      // Exact K0 geostatic elastic initialization so slopes start in calm equilibrium!
      if (sigmaV > 1 && matId !== MAT.BEDROCK) {
        const E = this.materialYoungs(matId);
        const nu = this.materialPoisson(matId);
        const { mu, k2d } = C.lameParams(E, nu);
        const f0 = C.geostaticInitialF(sigmaV, K0, mu, k2d);
        this.f00[i] = f0.f00; this.f01[i] = 0;
        this.f10[i] = 0;      this.f11[i] = f0.f11;
        this.sxx[i] = K0 * sigmaV;
        this.syy[i] = sigmaV;
        this.txy[i] = 0;
        this.pEff[i] = 0.5 * (1 + K0) * sigmaV;
        this.qShear[i] = 0.5 * Math.abs(1 - K0) * sigmaV;
      } else {
        this.f00[i] = 1; this.f01[i] = 0;
        this.f10[i] = 0; this.f11[i] = 1;
        this.sxx[i] = 0; this.syy[i] = 0; this.txy[i] = 0;
        this.pEff[i] = 0; this.qShear[i] = 0;
      }
      return i;
    }

    materialYoungs(matId) {
      if (matId === MAT.RIGID_STRUCT) return this.params.youngs * 4.5;
      if (matId === MAT.CRUST_CLAY) return this.params.youngs * 1.4;
      if (matId === MAT.BEDROCK) return this.params.youngs * 3.0;
      return this.params.youngs;
    }

    materialPoisson(matId) {
      if (matId === MAT.RIGID_STRUCT) return 0.24;
      if (matId === MAT.QUICK_CLAY) return 0.38; // Near-undrained Poisson's ratio
      return this.params.poisson;
    }

    loadScenario(id) {
      const cfg = SCENARIOS[id] || SCENARIOS.quickclay;
      this.scenarioId = cfg.id;
      this.domainW = cfg.domainW;
      this.domainH = cfg.domainH;
      this.initGrid(cfg.nx, cfg.ny);
      this.clearParticles();

      this.params.gravity = cfg.gravity;
      this.params.phiDeg = cfg.phiDeg;
      this.params.psiDeg = cfg.psiDeg;
      this.params.cohesion = cfg.cohesion;
      this.params.su0 = cfg.su0;
      this.params.sensitivity = cfg.sensitivity;
      this.params.kappa95 = cfg.kappa95;
      this.params.youngs = cfg.youngs;
      this.params.poisson = cfg.poisson;
      this.params.liqBeta = cfg.liqBeta;
      this.params.darcyCv = cfg.darcyCv;
      this.dt = this.computeStableDt();

      const dx = this.dx;
      const ppc = 2; // 2x2 = 4 material points per grid cell
      const sub = dx / ppc;
      const vol0 = sub * sub;

      this.platenVel = 0;
      this.platenForceY = 0;

      if (cfg.id === 'quickclay') {
        // Nordic marine quick-clay river slope:
        // Left plateau height y=12.4m (x in [1.2, 16.5]), sloping down to valley floor y=4.6m (x in [24.0, 42.5])
        // Dry crust top 1.5m over sensitive quick clay over stiff moraine base (y < 2.2m)
        const rho = 1900;
        const K0 = 0.75;
        for (let gx = 3; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 3; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                // Smooth realistic Norwegian valley slope profile
                let ySurf;
                if (px <= 15.0) ySurf = 11.8;
                else if (px >= 26.0) ySurf = 4.4;
                else {
                  const t = (px - 15.0) / 11.0;
                  ySurf = 11.8 - 7.4 * (0.5 - 0.5 * Math.cos(t * Math.PI));
                }
                if (py > ySurf) continue;
                const depth = Math.max(0, ySurf - py);
                let matId = MAT.QUICK_CLAY;
                let band = Math.floor(py / 1.35) % 4;
                if (py < 2.1) {
                  matId = MAT.BEDROCK;
                  band = 5;
                } else if (depth < 1.25 && px < 23.0) {
                  matId = MAT.CRUST_CLAY;
                  band = 4;
                }
                const sigmaV = rho * this.params.gravity * depth;
                this.addParticle(px, py, matId, rho, vol0, sigmaV, K0, band);
              }
            }
          }
        }
      } else if (cfg.id === 'liquefaction') {
        // Saturated loose contractive sand layer supporting a heavy rigid concrete building block
        const rhoSand = 1920;
        const rhoStruct = 2750;
        const yGround = 9.2;
        const bMin = 13.0;
        const bMax = 19.0;
        const bTop = 12.2;
        const K0 = 0.52;
        for (let gx = 2; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 3; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                if (py <= yGround) {
                  const depth = yGround - py;
                  const matId = py < 1.8 ? MAT.BEDROCK : MAT.LIQ_SAND;
                  const band = py < 1.8 ? 5 : Math.floor(py / 1.2) % 4;
                  const gammaEff = (rhoSand - 1000) * this.params.gravity;
                  let sigV = gammaEff * depth;
                  if (px >= bMin - 0.8 && px <= bMax + 0.8) {
                    sigV += rhoStruct * this.params.gravity * (bTop - yGround) * Math.exp(-depth * 0.14);
                  }
                  this.addParticle(px, py, matId, rhoSand, vol0, sigV, K0, band);
                } else if (py <= bTop && px >= bMin && px <= bMax) {
                  const extraRho = px < 15.4 ? rhoStruct * 1.28 : rhoStruct * 0.86;
                  const depth = bTop - py;
                  this.addParticle(px, py, MAT.RIGID_STRUCT, extraRho, vol0, extraRho * this.params.gravity * depth, 0.3, 6);
                }
              }
            }
          }
        }
      } else if (cfg.id === 'footing') {
        // Classic Prandtl-Terzaghi strip footing bearing capacity test
        const rho = 1900;
        const yGround = 8.4;
        const B = cfg.footingWidth;
        const cx = this.domainW * 0.5;
        this.platenXMin = cx - 0.5 * B;
        this.platenXMax = cx + 0.5 * B;
        this.platenY = yGround;
        this.platenY0 = yGround;
        this.platenVel = -0.34; // Controlled downward indentation speed (m/s)
        this.platenArea = B;
        const K0 = 0.62;

        for (let gx = 2; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 3; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                if (py > yGround) continue;
                const depth = yGround - py;
                const band = Math.floor(py / 0.85) % 4;
                const sigV = rho * this.params.gravity * depth;
                this.addParticle(px, py, MAT.SAND_MC, rho, vol0, sigV, K0, band);
              }
            }
          }
        }
      } else if (cfg.id === 'column') {
        // Granular column collapse (Lajeunesse / Staron-Hinch MPM benchmark)
        const rho = 1850;
        const xMin = 2.0 * dx;
        const xMax = xMin + cfg.colR0;
        const yMin = 2.0 * dx;
        const yMax = yMin + cfg.colH0;
        const K0 = 0.48;
        for (let gx = 2; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 2; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                if (px < xMin || px > xMax || py < yMin || py > yMax) continue;
                const depth = yMax - py;
                const band = Math.floor((py - yMin) / 0.72) % 4;
                const sigV = rho * this.params.gravity * depth;
                this.addParticle(px, py, MAT.SAND_MC, rho, vol0, sigV, K0, band);
              }
            }
          }
        }
      } else if (cfg.id === 'retaining') {
        // Cantilevered stiff sheet-pile wall retaining cohesionless-cohesive backfill
        const rhoSoil = 1900;
        const rhoWall = 2600;
        const wallXMin = 15.6;
        const wallXMax = 16.5;
        const yBack = 11.4;
        const yFront = 6.2;
        const K0 = 0.50;

        for (let gx = 2; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 3; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                const inWall = px >= wallXMin && px <= wallXMax && py >= 2.2 && py <= yBack + 0.6;
                if (inWall) {
                  this.addParticle(px, py, MAT.RIGID_STRUCT, rhoWall, vol0, rhoWall * this.params.gravity * (yBack + 0.6 - py), 0.3, 6);
                } else if (px < wallXMin && py <= yBack) {
                  const depth = yBack - py;
                  const matId = py < 2.2 ? MAT.BEDROCK : MAT.SAND_MC;
                  const band = py < 2.2 ? 5 : Math.floor(py / 1.1) % 4;
                  this.addParticle(px, py, matId, rhoSoil, vol0, rhoSoil * this.params.gravity * depth, K0, band);
                } else if (px > wallXMax && py <= yFront) {
                  const depth = yFront - py;
                  const matId = py < 2.2 ? MAT.BEDROCK : MAT.SAND_MC;
                  const band = py < 2.2 ? 5 : Math.floor(py / 1.1) % 4;
                  this.addParticle(px, py, matId, rhoSoil, vol0, rhoSoil * this.params.gravity * depth, K0, band);
                }
              }
            }
          }
        }
      } else if (cfg.id === 'biaxial') {
        // Weightless biaxial compression sample (4.4m x 7.4m) with tiny center heterogeneity
        // so conjugate X-shaped shear bands nucleate cleanly!
        const rho = 1900;
        const cx = this.domainW * 0.5;
        const halfW = 2.2;
        const yBot = 2.0;
        const yTop = 9.4;
        this.biaxialH0 = yTop - yBot;
        this.platenXMin = cx - halfW - 0.6;
        this.platenXMax = cx + halfW + 0.6;
        this.platenY = yTop;
        this.platenY0 = yTop;
        this.platenVel = -0.22;
        this.platenArea = 2 * halfW;

        for (let gx = 2; gx < this.nx - 2; gx++) {
          for (let gy = 2; gy < this.ny - 2; gy++) {
            for (let sy = 0; sy < ppc; sy++) {
              for (let sx = 0; sx < ppc; sx++) {
                const px = (gx + (sx + 0.5) / ppc) * dx;
                const py = (gy + (sy + 0.5) / ppc) * dx;
                if (px < cx - halfW || px > cx + halfW || py < yBot || py > yTop) continue;
                const band = Math.floor((py - yBot) / 0.74) % 4;
                const idx = this.addParticle(px, py, MAT.SAND_MC, rho, vol0, this.biaxialConfine, 1.0, band);
                // Seed a 1.5% micro-imperfection at the specimen center to break exact symmetry
                if (idx >= 0 && Math.hypot(px - cx, py - 0.5 * (yBot + yTop)) < 0.45) {
                  this.kappa[idx] = 0.025;
                }
              }
            }
          }
        }
      }

      // Quasi-static kinetic-damping gravity-turn-on equilibration (FLAC / Anura3D style)
      if (cfg.id === 'quickclay' || cfg.id === 'liquefaction' || cfg.id === 'retaining') {
        this.geostaticRelax(240);
      }
      this.recordHistory(true);
    }

    // Quasi-static Cundall kinetic-damping gravity-turn-on equilibration (Anura3D / FLAC style)
    // Reaches true (sxx, syy, txy) slope equilibrium in <250 steps with minimal residual wave.
    geostaticRelax(nSteps) {
      const savedSens = this.params.sensitivity;
      const savedBeta = this.params.liqBeta;
      this.params.sensitivity = 1.0;
      this.params.liqBeta = 0.0;

      for (let s = 0; s < nSteps; s++) {
        this.step(this.dt);
        const resetVel = (s + 1) % 45 === 0;
        const dampV = resetVel ? 0.0 : 0.965;
        for (let i = 0; i < this.count; i++) {
          if (!this.active[i] || this.mat[i] === MAT.BEDROCK) continue;
          this.vx[i] *= dampV; this.vy[i] *= dampV;
          this.c00[i] *= dampV; this.c01[i] *= dampV;
          this.c10[i] *= dampV; this.c11[i] *= dampV;
          this.kappa[i] = 0;
          this.uPore[i] = 0;
          this.ru[i] = 0;
        }
      }

      this.params.sensitivity = savedSens;
      this.params.liqBeta = savedBeta;
      for (let i = 0; i < this.count; i++) {
        this.vx[i] = 0; this.vy[i] = 0;
        this.c00[i] = 0; this.c01[i] = 0; this.c10[i] = 0; this.c11[i] = 0;
        this.kappa[i] = 0;
        this.uPore[i] = 0;
        this.ru[i] = 0;
        if (this.mat[i] === MAT.QUICK_CLAY) this.suCur[i] = this.params.su0;
      }
      this.time = 0;
      this.stepCount = 0;
      this.history = [];
      this.lastHistoryTime = -1;
    }

    // Trigger scripted toe excavation for Quick-Clay or Retaining Wall presets
    triggerExcavationPreset() {
      if (this.scenarioId === 'quickclay') {
        // Excavate river-bank toe notch -> removes passive toe support & initiates retrogressive quick-clay landslide!
        this.excavateCircle(24.2, 4.4, 3.4);
        this.excavateCircle(21.0, 5.2, 3.5);
        this.excavateCircle(18.6, 6.8, 3.1);
      } else if (this.scenarioId === 'retaining') {
        // Excavate passive toe in front of sheet-pile wall
        this.excavateRect(16.6, 3.3, 27.5, 7.2);
      } else {
        this.excavateCircle(this.domainW * 0.5, this.domainH * 0.4, 2.0);
      }
    }

    triggerEarthquake() {
      this.params.seismicActive = true;
      this.params.seismicStartTime = this.time;
    }

    excavateCircle(cx, cy, radius) {
      const r2 = radius * radius;
      let removed = 0;
      for (let i = 0; i < this.count; i++) {
        if (!this.active[i] || this.mat[i] === MAT.BEDROCK) continue;
        const dx = this.x[i] - cx;
        const dy = this.y[i] - cy;
        if (dx * dx + dy * dy <= r2) {
          this.active[i] = 0;
          removed++;
        }
      }
      return removed;
    }

    excavateRect(xMin, yMin, xMax, yMax) {
      let removed = 0;
      for (let i = 0; i < this.count; i++) {
        if (!this.active[i] || this.mat[i] === MAT.BEDROCK || this.mat[i] === MAT.RIGID_STRUCT) continue;
        if (this.x[i] >= xMin && this.x[i] <= xMax && this.y[i] >= yMin && this.y[i] <= yMax) {
          this.active[i] = 0;
          removed++;
        }
      }
      return removed;
    }

    pourParticles(cx, cy, radius, matId) {
      const sub = this.dx * 0.5;
      const vol0 = sub * sub;
      const rho = matId === MAT.RIGID_STRUCT ? 2550 : 1900;
      let added = 0;
      for (let ox = -radius; ox <= radius; ox += sub) {
        for (let oy = -radius; oy <= radius; oy += sub) {
          if (ox * ox + oy * oy > radius * radius) continue;
          const px = cx + ox;
          const py = cy + oy;
          if (px < 2 * this.dx || px > this.domainW - 2 * this.dx || py < 2 * this.dx || py > this.domainH - 2 * this.dx) continue;
          const idx = this.addParticle(px, py, matId, rho, vol0, 2000, 0.55, Math.floor(py / 1.1) % 4);
          if (idx >= 0) {
            this.vy[idx] = -0.6;
            added++;
          }
        }
      }
      return added;
    }

    placeSurchargeBlock(cx, cy, width, height) {
      const sub = this.dx * 0.5;
      const vol0 = sub * sub;
      const rho = 2750;
      const x0 = cx - 0.5 * width;
      const y0 = cy - 0.5 * height;
      for (let px = x0; px <= x0 + width; px += sub) {
        for (let py = y0; py <= y0 + height; py += sub) {
          if (px < 2 * this.dx || px > this.domainW - 2 * this.dx || py < 2 * this.dx || py > this.domainH - 2 * this.dx) continue;
          const idx = this.addParticle(px, py, MAT.RIGID_STRUCT, rho, vol0, rho * this.params.gravity * (y0 + height - py), 0.3, 6);
          if (idx >= 0) this.vy[idx] = -0.25;
        }
      }
    }

    currentSeismicAccelX() {
      if (!this.params.seismicActive) return 0;
      const elapsed = this.time - this.params.seismicStartTime;
      if (elapsed < 0 || elapsed > this.params.seismicDuration) {
        this.params.seismicActive = false;
        return 0;
      }
      // Realistic Husid/Saragoni-Hart tapered earthquake envelope + 2-harmonic motion
      const frac = elapsed / this.params.seismicDuration;
      const env = Math.sin(frac * Math.PI) * Math.pow(Math.sin(Math.min(1, frac * 1.35) * Math.PI), 0.5);
      const w = 2 * Math.PI * this.params.seismicFreqHz;
      const wave = 0.72 * Math.sin(w * elapsed) + 0.28 * Math.sin(1.63 * w * elapsed + 0.8);
      return this.params.seismicPgaG * 9.81 * env * wave;
    }

    // Single explicit Quadratic B-Spline MLS-MPM + APIC substep
    step(dtOverride) {
      const dt = dtOverride || this.dt;
      const nx = this.nx;
      const ny = this.ny;
      const stride = nx + 1;
      const dx = this.dx;
      const invDx = this.invDx;
      const stressScale = (-4.0 * dt) * invDx * invDx;

      const gMass = this.gMass;
      const gVx = this.gVx;
      const gVy = this.gVy;
      const gFixed = this.gFixed;
      const gPoreNum = this.gPoreNum;
      const gPore = this.gPore;
      const gPoreNext = this.gPoreNext;

      gMass.fill(0);
      gVx.fill(0);
      gVy.fill(0);
      gFixed.fill(0);
      const doPoreDiffusion = this.params.darcyCv > 1e-4 && this.scenarioId === 'liquefaction';
      if (doPoreDiffusion) gPoreNum.fill(0);

      const wx = this._wx;
      const wy = this._wy;
      const svd = this._svd;
      const ret = this._ret;
      const mTmp = this._matTmp;
      const sTmp = this._symTmp;

      const P = this.params;
      const count = this.count;

      // 1. P2G: Scatter mass, APIC momentum, and Kirchhoff elastoplastic stress to Eulerian grid
      for (let p = 0; p < count; p++) {
        if (!this.active[p]) continue;

        const px = this.x[p];
        const py = this.y[p];
        const gx = px * invDx;
        const gy = py * invDx;
        const bx = Math.floor(gx - 0.5);
        const by = Math.floor(gy - 0.5);
        if (bx < 0 || by < 0 || bx > nx - 2 || by > ny - 2) continue;

        const fx = gx - bx;
        const fy = gy - by;
        wx[0] = 0.5 * (1.5 - fx) * (1.5 - fx);
        wx[1] = 0.75 - (fx - 1.0) * (fx - 1.0);
        wx[2] = 0.5 * (fx - 0.5) * (fx - 0.5);
        wy[0] = 0.5 * (1.5 - fy) * (1.5 - fy);
        wy[1] = 0.75 - (fy - 1.0) * (fy - 1.0);
        wy[2] = 0.5 * (fy - 0.5) * (fy - 0.5);

        const mp = this.mass[p];
        const vp0 = this.vol0[p];
        const matId = this.mat[p];

        // Rigid bedrock particles mark Eulerian nodes as no-slip Dirichlet base
        if (matId === MAT.BEDROCK) {
          for (let oy = 0; oy < 3; oy++) {
            const rowOff = (by + oy) * stride + bx;
            for (let ox = 0; ox < 3; ox++) {
              gFixed[rowOff + ox] = 1;
            }
          }
          continue;
        }

        // Compute Kirchhoff stress from current F (tension positive!)
        const E = this.materialYoungs(matId);
        const nu = this.materialPoisson(matId);
        const { mu, k2d } = C.lameParams(E, nu);

        C.svd2x2(this.f00[p], this.f01[p], this.f10[p], this.f11[p], svd);
        const e1 = Math.log(svd.s1);
        const e2 = Math.log(svd.s2);

        let tau1, tau2;
        if (matId === MAT.RIGID_STRUCT || matId === MAT.BEDROCK) {
          const ev = e1 + e2;
          const pTot = -k2d * ev;
          const q = mu * (e1 - e2);
          tau1 = -pTot + q;
          tau2 = -pTot - q;
          this.pEff[p] = pTot;
          this.qShear[p] = Math.abs(q);
        } else {
          // Already returned onto yield surface during G2P; evaluate Hencky Kirchhoff stress
          const ev = e1 + e2;
          const pTot = -k2d * ev;
          const q = mu * (e1 - e2);
          tau1 = -pTot + q;
          tau2 = -pTot - q;
        }

        C.symUDUT(svd.u00, svd.u01, svd.u10, svd.u11, tau1, tau2, sTmp);
        // Save compression-positive Cauchy stress for live Mohr circle inspection
        const invJ = 1 / Math.max(0.25, svd.s1 * svd.s2);
        this.sxx[p] = -sTmp.t00 * invJ;
        this.syy[p] = -sTmp.t11 * invJ;
        this.txy[p] = -sTmp.t01 * invJ;

        // Combine APIC affine matrix m_p*C_p and MLS-MPM stress force (-4*dt/dx^2 * V0 * tau)
        const sFac = stressScale * vp0;
        const q00 = sFac * sTmp.t00 + mp * this.c00[p];
        const q01 = sFac * sTmp.t01 + mp * this.c01[p];
        const q10 = sFac * sTmp.t01 + mp * this.c10[p];
        const q11 = sFac * sTmp.t11 + mp * this.c11[p];

        const pvx = this.vx[p];
        const pvy = this.vy[p];
        const uP = this.uPore[p];

        for (let oy = 0; oy < 3; oy++) {
          const wyv = wy[oy];
          const dyNode = (by + oy) * dx - py;
          const rowOff = (by + oy) * stride + bx;
          for (let ox = 0; ox < 3; ox++) {
            const w = wx[ox] * wyv;
            const dxNode = (bx + ox) * dx - px;
            const idx = rowOff + ox;
            gMass[idx] += w * mp;
            gVx[idx] += w * (mp * pvx + q00 * dxNode + q01 * dyNode);
            gVy[idx] += w * (mp * pvy + q10 * dxNode + q11 * dyNode);
            if (doPoreDiffusion) gPoreNum[idx] += w * mp * uP;
          }
        }
      }

      // 2. Eulerian Grid Update: Gravity, Seismic Shaking, Rigid Platens & Boundary Conditions
      const axSeis = this.currentSeismicAccelX();
      const grav = P.gravity;
      const damp = 1.0 - 0.14 * dt; // Subtle numerical high-frequency damping
      const isFooting = this.scenarioId === 'footing';
      const isBiaxial = this.scenarioId === 'biaxial';
      const isLiq = this.scenarioId === 'liquefaction';
      let platenReactionImpulse = 0;

      if (isFooting && this.platenY > this.platenY0 - 2.3) {
        this.platenY += this.platenVel * dt;
      } else if (isBiaxial && this.platenY > this.platenY0 - 1.35) {
        this.platenY += this.platenVel * dt;
      }

      for (let gy = 0; gy <= ny; gy++) {
        const yNode = gy * dx;
        const rowOff = gy * stride;
        for (let gx = 0; gx <= nx; gx++) {
          const idx = rowOff + gx;
          const m = gMass[idx];
          if (m <= 1e-12) continue;

          let vx = (gVx[idx] / m + axSeis * dt) * damp;
          let vy = (gVy[idx] / m - grav * dt) * damp;
          const xNode = gx * dx;

          // Domain & rigid bedrock boundary conditions
          if (gFixed[idx]) {
            vx = 0;
            vy = 0;
          }
          if (!isLiq && (gx <= 2 || gx >= nx - 2)) {
            vx = 0; // Roller vertical walls
          }
          if (isLiq && (gx <= 1 || gx >= nx - 1)) {
            vx *= 0.85;
          }
          if (gy <= 2) {
            // Rough no-slip bedrock base
            vx = 0;
            vy = Math.max(0, vy);
          }
          if (gy >= ny - 2) {
            vy = Math.min(0, vy);
          }

          // Rigid strip footing or biaxial compression top/bottom platens
          if (isFooting && xNode >= this.platenXMin && xNode <= this.platenXMax && yNode >= this.platenY - 0.55 * dx) {
            const vBefore = vy;
            vy = Math.min(vy, this.platenVel);
            vx *= 0.88; // Semi-rough concrete footing base
            platenReactionImpulse += m * (vBefore - vy);
          } else if (isBiaxial) {
            if (yNode <= 2.0 + 0.5 * dx) {
              vy = Math.max(0, vy);
            }
            if (yNode >= this.platenY - 0.5 * dx) {
              const vBefore = vy;
              vy = Math.min(vy, this.platenVel);
              platenReactionImpulse += m * (vBefore - vy);
            }
            // Lateral confining stress approximation on left/right edges of biaxial sample
            if (yNode > 2.15 && yNode < this.platenY - 0.15) {
              const cx = this.domainW * 0.5;
              if (Math.abs(xNode - cx) > 1.95 && gMass[idx + (xNode > cx ? 1 : -1)] < 1e-6) {
                const dirIn = xNode > cx ? -1 : 1;
                vx += dirIn * (this.biaxialConfine * dx * dt) / Math.max(m, 40);
              }
            }
          }

          gVx[idx] = vx;
          gVy[idx] = vy;
          if (doPoreDiffusion) gPore[idx] = gPoreNum[idx] / m;
        }
      }

      if (isLiq) {
        // Tied-DOF laminar shear-box lateral boundary (standard seismic liquefaction BC)
        for (let gy = 3; gy < ny - 2; gy++) {
          const idxL = gy * stride + 2;
          const idxR = gy * stride + (nx - 2);
          if (gMass[idxL] > 1e-12 && gMass[idxR] > 1e-12) {
            const vAvg = 0.5 * (gVx[idxL] + gVx[idxR]);
            gVx[idxL] = vAvg;
            gVx[idxR] = vAvg;
          }
        }
      }

      if (isFooting || isBiaxial) {
        const instForce = Math.abs(platenReactionImpulse) / dt;
        // Exponential moving average smooths micro-step contact jitter for crisp load-settlement plot
        this.platenForceY = this.platenForceY * 0.92 + instForce * 0.08;
      }

      // 2b. Optional 2D Eulerian Darcy Consolidation Diffusion (FLIP delta_u transfer -> zero boundary leakage!)
      if (doPoreDiffusion) {
        const alphaD = Math.min(0.18, P.darcyCv * dt * invDx * invDx);
        gPoreNext.fill(0);
        for (let gy = 2; gy < ny - 2; gy++) {
          const rowOff = gy * stride;
          for (let gx = 2; gx < nx - 2; gx++) {
            const idx = rowOff + gx;
            if (gMass[idx] <= 1e-12) continue;
            const uC = gPore[idx];
            if (gMass[idx + stride] <= 1e-12) {
              gPoreNext[idx] = -2.2 * alphaD * uC;
              continue;
            }
            const uL = gMass[idx - 1] > 1e-12 ? gPore[idx - 1] : uC;
            const uR = gMass[idx + 1] > 1e-12 ? gPore[idx + 1] : uC;
            const uD = gMass[idx - stride] > 1e-12 ? gPore[idx - stride] : uC;
            const uU = gPore[idx + stride];
            gPoreNext[idx] = alphaD * (uL + uR + uD + uU - 4 * uC);
          }
        }
      }

      // 3. G2P: Gather velocity + APIC C_p, advect, and perform exact elastoplastic return mapping
      const xLimMin = 2.05 * dx;
      const xLimMax = this.domainW - 2.05 * dx;
      const yLimMin = 2.05 * dx;
      const yLimMax = this.domainH - 2.05 * dx;

      for (let p = 0; p < count; p++) {
        if (!this.active[p]) continue;
        const matId = this.mat[p];
        if (matId === MAT.BEDROCK) {
          this.vx[p] = 0;
          this.vy[p] = 0;
          continue;
        }

        const px = this.x[p];
        const py = this.y[p];
        const gx = px * invDx;
        const gy = py * invDx;
        const bx = Math.floor(gx - 0.5);
        const by = Math.floor(gy - 0.5);
        if (bx < 0 || by < 0 || bx > nx - 2 || by > ny - 2) continue;

        const fx = gx - bx;
        const fy = gy - by;
        wx[0] = 0.5 * (1.5 - fx) * (1.5 - fx);
        wx[1] = 0.75 - (fx - 1.0) * (fx - 1.0);
        wx[2] = 0.5 * (fx - 0.5) * (fx - 0.5);
        wy[0] = 0.5 * (1.5 - fy) * (1.5 - fy);
        wy[1] = 0.75 - (fy - 1.0) * (fy - 1.0);
        wy[2] = 0.5 * (fy - 0.5) * (fy - 0.5);

        let nvx = 0, nvy = 0;
        let c00 = 0, c01 = 0, c10 = 0, c11 = 0;
        let uGrid = 0;
        const cFac = 4.0 * invDx * invDx;

        for (let oy = 0; oy < 3; oy++) {
          const wyv = wy[oy];
          const dyNode = (by + oy) * dx - py;
          const rowOff = (by + oy) * stride + bx;
          for (let ox = 0; ox < 3; ox++) {
            const w = wx[ox] * wyv;
            const dxNode = (bx + ox) * dx - px;
            const idx = rowOff + ox;
            const gvx = gVx[idx];
            const gvy = gVy[idx];
            nvx += w * gvx;
            nvy += w * gvy;
            c00 += w * gvx * dxNode;
            c01 += w * gvx * dyNode;
            c10 += w * gvy * dxNode;
            c11 += w * gvy * dyNode;
            if (doPoreDiffusion) uGrid += w * gPoreNext[idx];
          }
        }

        c00 *= cFac; c01 *= cFac;
        c10 *= cFac; c11 *= cFac;

        this.vx[p] = nvx;
        this.vy[p] = nvy;
        this.c00[p] = c00; this.c01[p] = c01;
        this.c10[p] = c10; this.c11[p] = c11;

        let nxPos = px + nvx * dt;
        let nyPos = py + nvy * dt;
        if (nxPos < xLimMin) { nxPos = xLimMin; this.vx[p] = 0; }
        if (nxPos > xLimMax) { nxPos = xLimMax; this.vx[p] = 0; }
        if (nyPos < yLimMin) { nyPos = yLimMin; this.vy[p] = 0; }
        if (nyPos > yLimMax) { nyPos = yLimMax; this.vy[p] = 0; }
        this.x[p] = nxPos;
        this.y[p] = nyPos;

        // Trial elastic deformation gradient F_tr = (I + dt*C) * F
        const a00 = 1 + dt * c00;
        const a01 = dt * c01;
        const a10 = dt * c10;
        const a11 = 1 + dt * c11;
        const of00 = this.f00[p], of01 = this.f01[p];
        const of10 = this.f10[p], of11 = this.f11[p];
        const nf00 = a00 * of00 + a01 * of10;
        const nf01 = a00 * of01 + a01 * of11;
        const nf10 = a10 * of00 + a11 * of10;
        const nf11 = a10 * of01 + a11 * of11;

        if (matId === MAT.RIGID_STRUCT) {
          // Stiff elastic concrete/steel structure (no plastic yielding)
          this.f00[p] = nf00; this.f01[p] = nf01;
          this.f10[p] = nf10; this.f11[p] = nf11;
          continue;
        }

        // Exact 2x2 SVD & Hencky return-mapping plasticity
        const E = this.materialYoungs(matId);
        const nu = this.materialPoisson(matId);
        const { mu, k2d } = C.lameParams(E, nu);

        C.svd2x2(nf00, nf01, nf10, nf11, svd);
        const e1 = Math.log(svd.s1);
        const e2 = Math.log(svd.s2);

        if (matId === MAT.QUICK_CLAY) {
          C.returnMapQuickClay(e1, e2, mu, k2d, this.kappa[p], P.su0, P.sensitivity, P.kappa95, ret);
          this.kappa[p] += ret.dLambda;
          this.suCur[p] = ret.suNew;
          this.pEff[p] = ret.pEff;
          this.qShear[p] = ret.qShear;
        } else if (matId === MAT.CRUST_CLAY) {
          // Weathered dry clay crust: moderate cohesion + gentle brittle fracture softening
          const cEff = C.quickClayStrength(this.kappa[p], P.cohesion, 4.0, 0.24);
          C.returnMapMohrCoulomb(e1, e2, mu, k2d, cEff, 24.0, 1.0, 0, ret);
          this.kappa[p] += ret.dLambda;
          this.suCur[p] = cEff;
          this.pEff[p] = ret.pEff;
          this.qShear[p] = ret.qShear;
        } else if (matId === MAT.LIQ_SAND) {
          // Saturated loose sand with FLIP Darcy consolidation + cyclic excess pore pressure buildup
          let uP = Math.max(0, this.uPore[p] + (doPoreDiffusion ? uGrid : 0));
          C.returnMapMohrCoulomb(e1, e2, mu, k2d, P.cohesion, P.phiDeg, P.psiDeg, uP, ret);
          const dLam = ret.dLambda;
          if (dLam > 0 && P.liqBeta > 0) {
            const pRef = Math.max(this.sigV0[p], 4000);
            const ruCur = Math.min(0.96, uP / pRef);
            const du = P.liqBeta * Math.max(0, ret.pEff) * dLam * Math.max(0, 1 - ruCur / 0.96);
            uP = Math.min(0.96 * pRef, uP + du);
          }
          if (P.seismicActive && P.liqBeta > 0) {
            const shearRate = Math.abs(c01 + c10) + Math.abs(c00 - c11);
            const pRef = Math.max(this.sigV0[p], 4000);
            const ruCur = Math.min(0.96, uP / pRef);
            const duCyc = P.liqBeta * 1.45 * pRef * shearRate * dt * Math.max(0, 1 - ruCur / 0.96);
            uP = Math.min(0.96 * pRef, uP + duCyc);
          }
          this.uPore[p] = uP;
          this.ru[p] = Math.min(1.0, uP / Math.max(this.sigV0[p], 3500));
          this.kappa[p] += dLam;
          this.pEff[p] = ret.pEff;
          this.qShear[p] = ret.qShear;
        } else {
          // Standard Mohr-Coulomb frictional-cohesive soil (with mild softening in biaxial mode)
          const cEff = this.scenarioId === 'biaxial'
            ? C.quickClayStrength(this.kappa[p], P.cohesion, P.sensitivity, P.kappa95)
            : P.cohesion;
          C.returnMapMohrCoulomb(e1, e2, mu, k2d, cEff, P.phiDeg, P.psiDeg, 0, ret);
          this.kappa[p] += ret.dLambda;
          this.pEff[p] = ret.pEff;
          this.qShear[p] = ret.qShear;
        }

        // Reconstruct updated elastic F = U * diag(exp(e1New), exp(e2New)) * V^T
        C.mulUDVT(
          svd.u00, svd.u01, svd.u10, svd.u11,
          Math.exp(ret.eps1New), Math.exp(ret.eps2New),
          svd.v00, svd.v01, svd.v10, svd.v11,
          mTmp
        );
        this.f00[p] = mTmp.m00; this.f01[p] = mTmp.m01;
        this.f10[p] = mTmp.m10; this.f11[p] = mTmp.m11;
      }

      this.time += dt;
      this.stepCount++;
      if (this.time - this.lastHistoryTime >= 0.04) {
        this.recordHistory(false);
      }
    }

    computeMetrics() {
      let ke = 0;
      let maxSpeed = 0;
      let maxPlastic = 0;
      let maxRu = 0;
      let activeCount = 0;
      let quickCount = 0;
      let remoldedCount = 0;
      let maxRunoutX = 0;
      let structMinY = 1e9;
      let structTiltRad = 0;
      let sxSum = 0, sySum = 0, sCount = 0;

      for (let i = 0; i < this.count; i++) {
        if (!this.active[i]) continue;
        activeCount++;
        const v2 = this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i];
        ke += 0.5 * this.mass[i] * v2;
        const spd = Math.sqrt(v2);
        if (spd > maxSpeed) maxSpeed = spd;
        if (this.kappa[i] > maxPlastic) maxPlastic = this.kappa[i];
        if (this.ru[i] > maxRu) maxRu = this.ru[i];
        if (this.x[i] > maxRunoutX && this.mat[i] !== MAT.BEDROCK) maxRunoutX = this.x[i];
        if (this.mat[i] === MAT.QUICK_CLAY) {
          quickCount++;
          if (this.suCur[i] < this.params.su0 * 0.25) remoldedCount++;
        }
        if (this.mat[i] === MAT.RIGID_STRUCT) {
          sxSum += this.x[i];
          sySum += this.y[i];
          sCount++;
          if (this.y[i] < structMinY) structMinY = this.y[i];
        }
      }

      // Compute rigid building tilt angle via principal axis covariance
      if (sCount > 10) {
        const mx = sxSum / sCount;
        const my = sySum / sCount;
        let cxx = 0, cxy = 0, cyy = 0;
        for (let i = 0; i < this.count; i++) {
          if (!this.active[i] || this.mat[i] !== MAT.RIGID_STRUCT) continue;
          const dx = this.x[i] - mx;
          const dy = this.y[i] - my;
          cxx += dx * dx; cxy += dx * dy; cyy += dy * dy;
        }
        structTiltRad = 0.5 * Math.atan2(2 * cxy, cxx - cyy);
      }

      const footingSettlement = Math.max(0, this.platenY0 - this.platenY);
      const footingQ = this.platenArea > 0 ? this.platenForceY / this.platenArea : 0;
      const prandtlQu = this.scenarioId === 'footing'
        ? C.prandtlBearingCapacity(this.params.cohesion, this.params.phiDeg, 1900 * this.params.gravity, SCENARIOS.footing.footingWidth, 0)
        : 0;

      const remoldedPct = quickCount > 0 ? (100 * remoldedCount) / quickCount : 0;
      const runoutRatio = this.scenarioId === 'column'
        ? Math.max(1.0, (maxRunoutX - 2 * this.dx) / SCENARIOS.column.colR0)
        : 0;
      const axialStrainPct = this.scenarioId === 'biaxial'
        ? (100 * footingSettlement) / this.biaxialH0
        : 0;

      return {
        time: this.time,
        activeCount,
        ke,
        maxSpeed,
        maxPlastic,
        maxRu,
        remoldedPct,
        footingSettlement,
        footingQ,
        prandtlQu,
        runoutRatio,
        axialStrainPct,
        structTiltDeg: (structTiltRad * 180) / Math.PI
      };
    }

    recordHistory(force) {
      if (!force && this.time - this.lastHistoryTime < 0.04) return;
      this.lastHistoryTime = this.time;
      const m = this.computeMetrics();
      this.history.push(m);
      if (this.history.length > this.maxHistory) {
        this.history.shift();
      }
    }

    // Spatial stress/state probe around world point (wx, wy) using Gaussian kernel
    probeAt(wx, wy, radius) {
      const r = radius || this.dx * 1.6;
      const r2 = r * r;
      let wSum = 0;
      let sxx = 0, syy = 0, txy = 0;
      let pEff = 0, qShear = 0, kappa = 0, uPore = 0, ru = 0, suCur = 0, vx = 0, vy = 0;
      let dominantMat = -1;
      let bestDist2 = r2;

      for (let i = 0; i < this.count; i++) {
        if (!this.active[i]) continue;
        const dx = this.x[i] - wx;
        const dy = this.y[i] - wy;
        const d2 = dx * dx + dy * dy;
        if (d2 > r2) continue;
        const w = Math.exp(-2.5 * (d2 / r2));
        wSum += w;
        sxx += w * this.sxx[i];
        syy += w * this.syy[i];
        txy += w * this.txy[i];
        pEff += w * this.pEff[i];
        qShear += w * this.qShear[i];
        kappa += w * this.kappa[i];
        uPore += w * this.uPore[i];
        ru += w * this.ru[i];
        suCur += w * this.suCur[i];
        vx += w * this.vx[i];
        vy += w * this.vy[i];
        if (d2 < bestDist2) {
          bestDist2 = d2;
          dominantMat = this.mat[i];
        }
      }

      if (wSum < 1e-6) return null;
      const inv = 1 / wSum;
      sxx *= inv; syy *= inv; txy *= inv;
      uPore *= inv;
      // Convert total stresses to effective stresses for Mohr's circle display
      const sxxEff = sxx - uPore;
      const syyEff = syy - uPore;
      const mohr = C.mohrCircleFromStress(sxxEff, syyEff, txy);

      return {
        x: wx,
        y: wy,
        mat: dominantMat,
        sxxEff,
        syyEff,
        txy,
        pEff: pEff * inv,
        qShear: qShear * inv,
        kappa: kappa * inv,
        uPore,
        ru: ru * inv,
        suCur: suCur * inv,
        speed: Math.hypot(vx * inv, vy * inv),
        mohr
      };
    }
  }

  MPMSolver.SCENARIOS = SCENARIOS;
  return MPMSolver;
});
