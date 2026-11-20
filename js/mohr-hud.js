/*
 * STRATA — Live Engineering Instruments HUD
 *   1. Real-Time Mohr's Circle + Mohr-Coulomb / Strain-Softening Tresca Failure Envelope
 *      (shows Total vs Effective Mohr circle shifting left under pore-pressure buildup!)
 *   2. Adaptive Real-Time Engineering Strip Chart (Load-Settlement vs Analytical Prandtl q_u,
 *      Lajeunesse Runout vs Time, Liquefaction r_u & Tilt, Remolded Quick-Clay %, Biaxial q-eps)
 */
(function (root) {
  'use strict';

  root.Strata = root.Strata || {};
  const C = root.Strata.Constitutive;
  const MAT = C.MAT;

  class EngineeringHud {
    constructor(mohrCanvas, chartCanvas) {
      this.mohrCanvas = mohrCanvas;
      this.mohrCtx = mohrCanvas.getContext('2d');
      this.chartCanvas = chartCanvas;
      this.chartCtx = chartCanvas.getContext('2d');
    }

    syncCanvas(canvas) {
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const w = Math.max(220, canvas.clientWidth || 300);
      const h = Math.max(140, canvas.clientHeight || 170);
      const pw = Math.round(w * dpr);
      const ph = Math.round(h * dpr);
      if (canvas.width !== pw || canvas.height !== ph) {
        canvas.width = pw;
        canvas.height = ph;
      }
      return { w, h, dpr };
    }

    drawMohrCircle(solver, probe) {
      const { w, h, dpr } = this.syncCanvas(this.mohrCanvas);
      const ctx = this.mohrCtx;
      ctx.save();
      ctx.scale(dpr, dpr);

      ctx.fillStyle = '#070b14';
      ctx.fillRect(0, 0, w, h);

      const padL = 34;
      const padR = 14;
      const padT = 20;
      const padB = 24;
      const plotW = w - padL - padR;
      const plotH = h - padT - padB;

      // Dynamic stress scale (kPa)
      let maxSig = 160;
      if (probe) {
        maxSig = Math.max(120, Math.ceil((Math.max(probe.mohr.sigma1 + probe.uPore, 80000) / 1000) / 40) * 40 + 40);
      }
      const maxTau = maxSig * (plotH / plotW);
      const ox = padL + plotW * 0.10; // Leave small tension side
      const oy = h - padB;
      const kpaScale = (plotW * 0.88) / maxSig;

      const toX = (kpa) => ox + kpa * kpaScale;
      const toY = (kpa) => oy - kpa * kpaScale;

      // Subtle axes & grid
      ctx.strokeStyle = 'rgba(130, 165, 220, 0.10)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      const gridStep = maxSig > 200 ? 60 : 40;
      for (let s = 0; s <= maxSig; s += gridStep) {
        const px = toX(s);
        ctx.moveTo(px, padT);
        ctx.lineTo(px, oy);
      }
      for (let t = gridStep; t <= maxTau; t += gridStep) {
        const py = toY(t);
        ctx.moveTo(padL, py);
        ctx.lineTo(w - padR, py);
      }
      ctx.stroke();

      // Main σ' and τ axes
      ctx.strokeStyle = 'rgba(175, 200, 235, 0.35)';
      ctx.beginPath();
      ctx.moveTo(padL, oy);
      ctx.lineTo(w - padR, oy);
      ctx.moveTo(ox, padT);
      ctx.lineTo(ox, oy);
      ctx.stroke();

      // Failure Envelope for probed material (or active scenario default)
      const matId = probe ? probe.mat : (solver.scenarioId === 'quickclay' ? MAT.QUICK_CLAY : MAT.SAND_MC);
      ctx.save();
      if (matId === MAT.QUICK_CLAY) {
        // Intact Tresca envelope s_u0 (dashed)
        const su0Kpa = solver.params.su0 * 1e-3;
        ctx.setLineDash([4, 4]);
        ctx.strokeStyle = 'rgba(255, 165, 86, 0.45)';
        ctx.lineWidth = 1.1;
        ctx.beginPath();
        ctx.moveTo(padL, toY(su0Kpa));
        ctx.lineTo(w - padR, toY(su0Kpa));
        ctx.stroke();

        // Current softened Tresca envelope s_u(kappa)
        const suCurKpa = probe ? probe.suCur * 1e-3 : su0Kpa;
        ctx.setLineDash([]);
        ctx.strokeStyle = '#ff546c';
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        ctx.moveTo(padL, toY(suCurKpa));
        ctx.lineTo(w - padR, toY(suCurKpa));
        ctx.stroke();
      } else {
        // Mohr-Coulomb Failure Envelope: tau_f = c' + sigma' * tan(phi')
        const cKpa = (matId === MAT.CRUST_CLAY ? solver.params.cohesion : solver.params.cohesion) * 1e-3;
        const phiRad = ((matId === MAT.CRUST_CLAY ? 24 : solver.params.phiDeg) * Math.PI) / 180;
        const tanPhi = Math.tan(phiRad);
        const sStart = Math.max(-maxSig * 0.1, tanPhi > 1e-4 ? -cKpa / tanPhi : 0);
        ctx.strokeStyle = '#ff546c';
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        ctx.moveTo(toX(sStart), toY(Math.max(0, cKpa + sStart * tanPhi)));
        ctx.lineTo(toX(maxSig), toY(cKpa + maxSig * tanPhi));
        ctx.stroke();
      }
      ctx.restore();

      if (probe) {
        const cEff = probe.mohr.center * 1e-3;
        const rad = probe.mohr.radius * 1e-3;
        const uKpa = probe.uPore * 1e-3;

        // Total stress Mohr circle (dashed slate when excess pore pressure u > 1 kPa)
        if (uKpa > 1.0) {
          const cTot = cEff + uKpa;
          ctx.save();
          ctx.setLineDash([4, 3]);
          ctx.strokeStyle = 'rgba(160, 185, 220, 0.50)';
          ctx.lineWidth = 1.3;
          ctx.beginPath();
          ctx.arc(toX(cTot), oy, Math.max(1.5, rad * kpaScale), Math.PI, 0, false);
          ctx.stroke();

          // Effective-stress shift arrow (u_e pushes circle LEFT toward failure envelope!)
          ctx.setLineDash([]);
          ctx.strokeStyle = '#ff6b7d';
          ctx.lineWidth = 1.4;
          const ay = oy - Math.max(10, rad * kpaScale * 0.55);
          ctx.beginPath();
          ctx.moveTo(toX(cTot), ay);
          ctx.lineTo(toX(cEff) + 3, ay);
          ctx.lineTo(toX(cEff) + 8, ay - 3.5);
          ctx.moveTo(toX(cEff) + 3, ay);
          ctx.lineTo(toX(cEff) + 8, ay + 3.5);
          ctx.stroke();
          ctx.restore();
        }

        // Effective stress Mohr circle (solid glowing cyan/amber)
        ctx.save();
        const yielded = probe.kappa > 0.015;
        ctx.fillStyle = yielded ? 'rgba(255, 152, 56, 0.16)' : 'rgba(56, 214, 255, 0.15)';
        ctx.strokeStyle = yielded ? '#ff9e3d' : '#38d6ff';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.arc(toX(cEff), oy, Math.max(2.0, rad * kpaScale), Math.PI, 0, false);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();

        // Principal stress points sigma'_3 and sigma'_1
        ctx.fillStyle = '#ffffff';
        ctx.beginPath();
        ctx.arc(toX(probe.mohr.sigma1 * 1e-3), oy, 2.8, 0, Math.PI * 2);
        ctx.arc(toX(probe.mohr.sigma3 * 1e-3), oy, 2.8, 0, Math.PI * 2);
        ctx.fill();
        ctx.restore();
      }

      // Labels
      ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
      ctx.fillStyle = 'rgba(212, 226, 248, 0.85)';
      ctx.fillText("MOHR'S CIRCLE (σ', τ) [kPa]", 8, 13);
      ctx.fillStyle = 'rgba(255, 92, 114, 0.9)';
      ctx.fillText(matId === MAT.QUICK_CLAY ? 'Tresca s_u(κ)' : 'MC Envelope', w - 82, 13);
      ctx.fillStyle = 'rgba(155, 178, 210, 0.72)';
      ctx.fillText("σ' (kPa)", w - 50, h - 6);
      ctx.fillText('0', ox - 3, h - 6);
      ctx.fillText(String(maxSig), toX(maxSig) - 18, h - 6);

      ctx.restore();
    }

    drawStripChart(solver) {
      const { w, h, dpr } = this.syncCanvas(this.chartCanvas);
      const ctx = this.chartCtx;
      ctx.save();
      ctx.scale(dpr, dpr);

      ctx.fillStyle = '#070b14';
      ctx.fillRect(0, 0, w, h);

      const padL = 38;
      const padR = 12;
      const padT = 22;
      const padB = 22;
      const plotW = w - padL - padR;
      const plotH = h - padT - padB;
      const hist = solver.history;

      // Subtle frame & grid
      ctx.strokeStyle = 'rgba(130, 165, 220, 0.12)';
      ctx.lineWidth = 1;
      ctx.strokeRect(padL, padT, plotW, plotH);
      for (let i = 1; i <= 3; i++) {
        const py = padT + (plotH * i) / 4;
        ctx.beginPath();
        ctx.moveTo(padL, py);
        ctx.lineTo(padL + plotW, py);
        ctx.stroke();
      }

      const id = solver.scenarioId;
      let title = 'KINETIC ENERGY & MOBILIZATION';
      let xMax = Math.max(4.0, solver.time);
      let yMax = 100;
      let refLine = null;
      let refLabel = '';
      let series1 = [];
      let series2 = [];
      let col1 = '#38d6ff';
      let col2 = '#ff9c3b';
      let xUnit = 't (s)';
      let yUnit = '';

      if (id === 'footing') {
        title = 'BEARING PRESSURE q vs s/B';
        const B = root.Strata.MPMSolver.SCENARIOS.footing.footingWidth;
        const quKpa = hist.length ? hist[hist.length - 1].prandtlQu * 1e-3 : 240;
        xMax = 0.65;
        yMax = Math.max(quKpa * 1.35, 160);
        refLine = quKpa;
        refLabel = `Prandtl q_u = ${quKpa.toFixed(0)} kPa`;
        xUnit = 's/B';
        yUnit = 'kPa';
        for (const m of hist) {
          series1.push([m.footingSettlement / B, m.footingQ * 1e-3]);
        }
      } else if (id === 'quickclay') {
        title = 'REMOLDED QUICK CLAY (%) & SPEED';
        yMax = 100;
        yUnit = '%';
        for (const m of hist) {
          series1.push([m.time, m.remoldedPct]);
          series2.push([m.time, Math.min(100, m.maxSpeed * 12)]);
        }
      } else if (id === 'liquefaction') {
        title = 'PEAK PORE RATIO r_u & TILT (°)';
        yMax = 1.1;
        refLine = 1.0;
        refLabel = 'Liquefaction r_u = 1.0';
        yUnit = 'r_u';
        for (const m of hist) {
          series1.push([m.time, m.maxRu]);
          series2.push([m.time, Math.min(1.1, Math.abs(m.structTiltDeg) / 25)]);
        }
      } else if (id === 'column') {
        title = 'NORMALIZED RUNOUT R(t)/R0';
        yMax = 4.2;
        refLine = 1 + 1.24 * (6.2 / 3.6);
        refLabel = `Lajeunesse = ${refLine.toFixed(2)}`;
        yUnit = 'R/R0';
        for (const m of hist) {
          series1.push([m.time, m.runoutRatio]);
        }
      } else if (id === 'biaxial') {
        title = 'AXIAL STRESS q (kPa) vs STRAIN ε_a';
        xMax = 22;
        yMax = 140;
        xUnit = 'ε_a (%)';
        yUnit = 'kPa';
        for (const m of hist) {
          series1.push([m.axialStrainPct, m.footingQ * 1e-3]);
        }
      } else {
        title = 'MAX PLASTIC STRAIN εp & SPEED';
        yMax = 2.0;
        for (const m of hist) {
          series1.push([m.time, Math.min(yMax, m.maxPlastic)]);
          series2.push([m.time, Math.min(yMax, m.maxSpeed * 0.4)]);
        }
      }

      const toX = (xv) => padL + (Math.min(xMax, Math.max(0, xv)) / Math.max(xMax, 1e-4)) * plotW;
      const toY = (yv) => padT + plotH - (Math.min(yMax, Math.max(0, yv)) / Math.max(yMax, 1e-4)) * plotH;

      // Analytical reference benchmark line
      if (refLine !== null && refLine <= yMax) {
        const ry = toY(refLine);
        ctx.save();
        ctx.setLineDash([4, 4]);
        ctx.strokeStyle = 'rgba(255, 198, 72, 0.85)';
        ctx.lineWidth = 1.3;
        ctx.beginPath();
        ctx.moveTo(padL, ry);
        ctx.lineTo(padL + plotW, ry);
        ctx.stroke();
        ctx.fillStyle = 'rgba(255, 206, 98, 0.92)';
        ctx.font = '9.5px ui-monospace, SFMono-Regular, Menlo, monospace';
        ctx.fillText(refLabel, padL + 6, Math.max(padT + 10, ry - 4));
        ctx.restore();
      }

      function plotLine(pts, color, fillAlpha) {
        if (pts.length < 2) return;
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        for (let i = 0; i < pts.length; i++) {
          const px = toX(pts[i][0]);
          const py = toY(pts[i][1]);
          if (i === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        }
        ctx.stroke();
        if (fillAlpha) {
          ctx.lineTo(toX(pts[pts.length - 1][0]), padT + plotH);
          ctx.lineTo(toX(pts[0][0]), padT + plotH);
          ctx.closePath();
          ctx.fillStyle = fillAlpha;
          ctx.fill();
        }
        ctx.restore();
      }

      plotLine(series2, col2, null);
      plotLine(series1, col1, 'rgba(56, 214, 255, 0.12)');

      ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
      ctx.fillStyle = 'rgba(212, 226, 248, 0.85)';
      ctx.fillText(title, 8, 14);
      ctx.fillStyle = 'rgba(155, 178, 210, 0.72)';
      ctx.fillText(yMax.toFixed(yMax < 5 ? 1 : 0), 4, padT + 8);
      ctx.fillText('0', 22, padT + plotH);
      ctx.fillText(xUnit, padL + plotW - 36, h - 6);
      if (yUnit) ctx.fillText(yUnit, 4, padT + plotH * 0.55);

      ctx.restore();
    }
  }

  root.Strata.EngineeringHud = EngineeringHud;
})(typeof globalThis !== 'undefined' ? globalThis : this);
