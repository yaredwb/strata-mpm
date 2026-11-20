/*
 * STRATA — High-DPI 60 FPS Computational Geomechanics Continuum Renderer
 * Renders MLS-MPM material points across 7 scientific colormaps, plus optional
 * oriented Principal Stress Crosses (sigma1/sigma3), Velocity Vectors, Eulerian MPM Mesh,
 * Rigid Footing/Biaxial Platens, Calibrated Colorbar, and Interactive Tool Cursors.
 */
(function (root) {
  'use strict';

  root.Strata = root.Strata || {};
  const C = root.Strata.Constitutive;
  const MAT = C.MAT;

  const FIELD_META = [
    { id: 0, name: 'Stratigraphy & Quick-Clay Phase', unit: 'phase', max: 1.0, cmap: 'phase' },
    { id: 1, name: 'Plastic Deviatoric Strain εp_q', unit: '—', max: 0.65, cmap: 'inferno' },
    { id: 2, name: 'Deviatoric Shear Stress q', unit: 'kPa', max: 55.0, cmap: 'inferno' },
    { id: 3, name: 'Mean Effective Stress p\'', unit: 'kPa', max: 140.0, cmap: 'viridis' },
    { id: 4, name: 'Pore-Pressure Ratio r_u = Δu/σ\'v0', unit: 'r_u', max: 1.0, cmap: 'liquefaction' },
    { id: 5, name: 'Softening Ratio s_u / s_u0', unit: 'ratio', max: 1.0, cmap: 'softening' },
    { id: 6, name: 'Velocity Magnitude |v|', unit: 'm/s', max: 4.5, cmap: 'turbo' }
  ];

  // Precomputed 256-entry RGB LUTs for ultra-fast 60 FPS particle rendering
  function buildLUT(stops) {
    const lut = new Array(256);
    const n = stops.length - 1;
    for (let i = 0; i < 256; i++) {
      const t = (i / 255) * n;
      const idx = Math.min(n - 1, Math.floor(t));
      const f = t - idx;
      const a = stops[idx];
      const b = stops[idx + 1];
      const r = Math.round(a[0] + (b[0] - a[0]) * f);
      const g = Math.round(a[1] + (b[1] - a[1]) * f);
      const bl = Math.round(a[2] + (b[2] - a[2]) * f);
      lut[i] = `rgb(${r},${g},${bl})`;
    }
    return lut;
  }

  const LUTS = {
    inferno: buildLUT([
      [14, 16, 30], [72, 16, 108], [156, 40, 102],
      [226, 86, 54], [252, 168, 16], [252, 252, 168]
    ]),
    viridis: buildLUT([
      [32, 18, 72], [54, 78, 138], [34, 140, 141],
      [96, 200, 96], [250, 232, 42]
    ]),
    liquefaction: buildLUT([
      [26, 72, 118], [36, 144, 158], [226, 188, 46],
      [246, 92, 30], [255, 32, 68]
    ]),
    softening: buildLUT([
      [255, 42, 134], [244, 132, 46], [190, 214, 88],
      [62, 178, 162], [46, 104, 152]
    ]),
    turbo: buildLUT([
      [34, 44, 92], [32, 146, 218], [62, 224, 150],
      [234, 212, 48], [248, 86, 36]
    ])
  };

  const STRATA_SAND = ['#be965d', '#a57c4b', '#cca66c', '#946e42'];
  const STRATA_CLAY = [
    [76, 112, 144], [92, 130, 158], [66, 99, 128], [102, 140, 168]
  ];

  function clamp(v, lo, hi) {
    return v < lo ? lo : (v > hi ? hi : v);
  }

  class ContinuumRenderer {
    constructor(canvas) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d', { alpha: false });
      this.viewTransform = { ox: 0, oy: 0, scale: 1, w: 800, h: 500 };
    }

    syncSize() {
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const w = Math.max(320, this.canvas.clientWidth || 800);
      const h = Math.max(240, this.canvas.clientHeight || 500);
      const pw = Math.round(w * dpr);
      const ph = Math.round(h * dpr);
      if (this.canvas.width !== pw || this.canvas.height !== ph) {
        this.canvas.width = pw;
        this.canvas.height = ph;
      }
      return { w, h, dpr };
    }

    computeTransform(solver, w, h) {
      const padL = 44;
      const padR = 24;
      const padT = 26;
      const padB = 38;
      const scale = Math.min((w - padL - padR) / solver.domainW, (h - padT - padB) / solver.domainH);
      const ox = padL + 0.5 * ((w - padL - padR) - solver.domainW * scale);
      const oy = h - padB - 0.5 * ((h - padT - padB) - solver.domainH * scale);
      this.viewTransform = { ox, oy, scale, w, h };
      return this.viewTransform;
    }

    screenToWorld(clientX, clientY) {
      const rect = this.canvas.getBoundingClientRect();
      const sx = clientX - rect.left;
      const sy = clientY - rect.top;
      const t = this.viewTransform;
      return {
        x: (sx - t.ox) / Math.max(t.scale, 1e-4),
        y: (t.oy - sy) / Math.max(t.scale, 1e-4)
      };
    }

    render(solver, uiState, probeData) {
      const { w, h, dpr } = this.syncSize();
      const ctx = this.ctx;
      const t = this.computeTransform(solver, w, h);
      const { ox, oy, scale } = t;

      ctx.save();
      ctx.scale(dpr, dpr);

      // Background studio gradient
      const bg = ctx.createLinearGradient(0, 0, 0, h);
      bg.addColorStop(0, '#090d16');
      bg.addColorStop(1, '#05070d');
      ctx.fillStyle = bg;
      ctx.fillRect(0, 0, w, h);

      // Domain boundary frame
      const domX = ox;
      const domY = oy - solver.domainH * scale;
      const domW = solver.domainW * scale;
      const domH = solver.domainH * scale;
      ctx.fillStyle = 'rgba(13, 19, 32, 0.55)';
      ctx.fillRect(domX, domY, domW, domH);

      // Optional Eulerian MPM computational background mesh
      if (uiState.showMesh) {
        ctx.strokeStyle = 'rgba(84, 152, 224, 0.14)';
        ctx.lineWidth = 0.75;
        ctx.beginPath();
        for (let gx = 0; gx <= solver.nx; gx += 2) {
          const x = ox + gx * solver.dx * scale;
          ctx.moveTo(x, domY);
          ctx.lineTo(x, oy);
        }
        for (let gy = 0; gy <= solver.ny; gy += 2) {
          const y = oy - gy * solver.dx * scale;
          ctx.moveTo(domX, y);
          ctx.lineTo(domX + domW, y);
        }
        ctx.stroke();
      }

      // Metric coordinate grid & axis labels
      const stepM = solver.domainW > 30 ? 4 : 2;
      ctx.strokeStyle = 'rgba(135, 175, 230, 0.08)';
      ctx.lineWidth = 1;
      ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
      ctx.fillStyle = 'rgba(148, 172, 208, 0.68)';
      ctx.beginPath();
      for (let xm = 0; xm <= solver.domainW + 1e-4; xm += stepM) {
        const px = ox + xm * scale;
        ctx.moveTo(px, domY);
        ctx.lineTo(px, oy);
        ctx.fillText(xm.toFixed(0) + 'm', px - 8, oy + 16);
      }
      for (let ym = 0; ym <= solver.domainH + 1e-4; ym += stepM) {
        const py = oy - ym * scale;
        ctx.moveTo(domX, py);
        ctx.lineTo(domX + domW, py);
        ctx.fillText(ym.toFixed(0) + 'm', domX - 28, py + 3);
      }
      ctx.stroke();

      // Subtle water table indicator for liquefaction scenario
      if (solver.scenarioId === 'liquefaction') {
        const yWater = oy - 9.2 * scale;
        ctx.save();
        ctx.setLineDash([6, 4]);
        ctx.strokeStyle = 'rgba(56, 200, 255, 0.52)';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(domX, yWater);
        ctx.lineTo(domX + domW, yWater);
        ctx.stroke();
        ctx.fillStyle = 'rgba(84, 214, 255, 0.78)';
        ctx.fillText('GWT (Saturated)', domX + 8, yWater - 6);
        ctx.restore();
      }

      // Draw Material Points
      const fieldId = uiState.fieldId;
      const meta = FIELD_META[fieldId] || FIELD_META[0];
      const lut = LUTS[meta.cmap] || LUTS.inferno;
      const pSize = Math.max(2.3, solver.dx * scale * 0.62);
      const halfP = 0.5 * pSize;
      const su0 = Math.max(1, solver.params.su0);

      for (let i = 0; i < solver.count; i++) {
        if (!solver.active[i]) continue;
        const sx = ox + solver.x[i] * scale;
        const sy = oy - solver.y[i] * scale;
        const matId = solver.mat[i];

        if (matId === MAT.BEDROCK) {
          ctx.fillStyle = '#242a37';
        } else if (matId === MAT.RIGID_STRUCT) {
          ctx.fillStyle = '#a4b2c7';
        } else if (fieldId === 0) {
          if (matId === MAT.CRUST_CLAY) {
            const k = clamp(solver.kappa[i] / 0.45, 0, 1);
            ctx.fillStyle = k > 0.3 ? '#d87a38' : '#946b44';
          } else if (matId === MAT.QUICK_CLAY) {
            const remold = clamp(1 - solver.suCur[i] / su0, 0, 1);
            const base = STRATA_CLAY[solver.strataBand[i] & 3];
            const remPow = Math.pow(remold, 0.72);
            let r = base[0] + (36 - base[0]) * remPow;
            let g = base[1] + (214 - base[1]) * remPow;
            let b = base[2] + (248 - base[2]) * remPow;
            const sb = clamp(solver.kappa[i] / 0.55, 0, 1);
            if (sb > 0.15) {
              const hs = 0.55 * sb * sb;
              r += (255 - r) * hs;
              g += (124 - g) * hs;
              b += (34 - b) * hs;
            }
            ctx.fillStyle = `rgb(${r | 0},${g | 0},${b | 0})`;
          } else {
            ctx.fillStyle = STRATA_SAND[solver.strataBand[i] & 3];
          }
        } else {
          let norm = 0;
          if (fieldId === 1) norm = Math.pow(clamp(solver.kappa[i] / meta.max, 0, 1), 0.62);
          else if (fieldId === 2) norm = Math.pow(clamp((solver.qShear[i] * 1e-3) / meta.max, 0, 1), 0.72);
          else if (fieldId === 3) norm = clamp((solver.pEff[i] * 1e-3) / meta.max, 0, 1);
          else if (fieldId === 4) norm = clamp(solver.ru[i], 0, 1);
          else if (fieldId === 5) {
            const ratio = matId === MAT.QUICK_CLAY ? clamp(solver.suCur[i] / su0, 0, 1) : 1.0;
            norm = ratio;
          } else {
            const spd = Math.hypot(solver.vx[i], solver.vy[i]);
            norm = clamp(spd / meta.max, 0, 1);
          }
          ctx.fillStyle = lut[Math.min(255, Math.max(0, (norm * 255) | 0))];
        }

        ctx.fillRect(sx - halfP, sy - halfP, pSize, pSize);
      }

      // Draw Rigid Indentation Footing or Biaxial Compression Platens
      if (solver.scenarioId === 'footing' || solver.scenarioId === 'biaxial') {
        const fx0 = ox + solver.platenXMin * scale;
        const fx1 = ox + solver.platenXMax * scale;
        const fyBot = oy - solver.platenY * scale;
        const fh = Math.max(14, 1.15 * scale);
        ctx.fillStyle = '#4e5d75';
        ctx.strokeStyle = '#c4d5ee';
        ctx.lineWidth = 1.6;
        ctx.fillRect(fx0, fyBot - fh, fx1 - fx0, fh);
        ctx.strokeRect(fx0, fyBot - fh, fx1 - fx0, fh);

        // Downward load force arrow
        const fcx = 0.5 * (fx0 + fx1);
        ctx.strokeStyle = '#ff9d3b';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(fcx, fyBot - fh - 22);
        ctx.lineTo(fcx, fyBot - fh - 4);
        ctx.lineTo(fcx - 5, fyBot - fh - 10);
        ctx.moveTo(fcx, fyBot - fh - 4);
        ctx.lineTo(fcx + 5, fyBot - fh - 10);
        ctx.stroke();

        if (solver.scenarioId === 'biaxial') {
          const fyLow = oy - 2.0 * scale;
          ctx.fillStyle = '#3c485c';
          ctx.fillRect(fx0, fyLow, fx1 - fx0, 12);
          ctx.strokeRect(fx0, fyLow, fx1 - fx0, 12);
        }
      }

      // Optional Principal Stress Tensor Crosses (sigma1 major, sigma3 minor)
      if (uiState.showStressTensors) {
        this.drawPrincipalStressCrosses(ctx, solver, ox, oy, scale);
      }

      // Optional Velocity Vector Field
      if (uiState.showVelocities) {
        this.drawVelocityVectors(ctx, solver, ox, oy, scale);
      }

      // Domain crisp border
      ctx.strokeStyle = 'rgba(145, 180, 235, 0.25)';
      ctx.lineWidth = 1.2;
      ctx.strokeRect(domX, domY, domW, domH);

      // Seismic base excitation indicator
      if (solver.params.seismicActive) {
        const ax = solver.currentSeismicAccelX() / 9.81;
        ctx.fillStyle = 'rgba(255, 72, 98, 0.18)';
        ctx.strokeStyle = 'rgba(255, 92, 114, 0.72)';
        ctx.lineWidth = 1.2;
        ctx.fillRect(domX + 10, domY + 10, 210, 26);
        ctx.strokeRect(domX + 10, domY + 10, 210, 26);
        ctx.fillStyle = '#ff98a8';
        ctx.font = 'bold 11px ui-monospace, SFMono-Regular, Menlo, monospace';
        ctx.fillText(`SEISMIC EXCITATION  a_x = ${(ax >= 0 ? '+' : '') + ax.toFixed(2)} g`, domX + 18, domY + 27);
      }

      // Active Stress Probe Reticle
      if (probeData) {
        const px = ox + probeData.x * scale;
        const py = oy - probeData.y * scale;
        ctx.save();
        ctx.strokeStyle = '#38e8ff';
        ctx.lineWidth = 1.6;
        ctx.beginPath();
        ctx.arc(px, py, 8, 0, Math.PI * 2);
        ctx.moveTo(px - 14, py); ctx.lineTo(px - 5, py);
        ctx.moveTo(px + 5, py);  ctx.lineTo(px + 14, py);
        ctx.moveTo(px, py - 14); ctx.lineTo(px, py - 5);
        ctx.moveTo(px, py + 5);  ctx.lineTo(px, py + 14);
        ctx.stroke();
        ctx.restore();
      }

      // Interactive Tool Brush Preview
      if (uiState.pointerInCanvas && uiState.tool !== 'probe') {
        const px = ox + uiState.pointerWorldX * scale;
        const py = oy - uiState.pointerWorldY * scale;
        ctx.save();
        if (uiState.tool === 'excavate') {
          ctx.strokeStyle = 'rgba(255, 92, 112, 0.88)';
          ctx.setLineDash([5, 4]);
          ctx.lineWidth = 1.6;
          ctx.beginPath();
          ctx.arc(px, py, uiState.brushRadius * scale, 0, Math.PI * 2);
          ctx.stroke();
        } else if (uiState.tool === 'pour_sand' || uiState.tool === 'pour_clay') {
          ctx.strokeStyle = uiState.tool === 'pour_clay' ? 'rgba(56, 214, 255, 0.85)' : 'rgba(255, 185, 82, 0.85)';
          ctx.setLineDash([4, 4]);
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          ctx.arc(px, py, uiState.brushRadius * scale, 0, Math.PI * 2);
          ctx.stroke();
        } else if (uiState.tool === 'footing') {
          const bw = 3.2 * scale;
          const bh = 2.0 * scale;
          ctx.strokeStyle = 'rgba(200, 222, 255, 0.85)';
          ctx.setLineDash([5, 3]);
          ctx.strokeRect(px - 0.5 * bw, py - 0.5 * bh, bw, bh);
        }
        ctx.restore();
      }

      // Calibrated Field Legend / Colorbar (top-right inside domain)
      this.drawColorbar(ctx, domX + domW - 216, domY + 12, 204, meta);

      ctx.restore();
    }

    drawPrincipalStressCrosses(ctx, solver, ox, oy, scale) {
      const stepX = Math.max(1.2, solver.domainW / 24);
      const stepY = Math.max(0.9, solver.domainH / 12);
      const maxLen = Math.min(stepX, stepY) * scale * 0.44;
      const normRef = 140000; // 140 kPa reference full cross length

      ctx.save();
      ctx.lineWidth = 1.3;
      for (let wx = stepX; wx < solver.domainW - stepX * 0.5; wx += stepX) {
        for (let wy = stepY; wy < solver.domainH - stepY * 0.5; wy += stepY) {
          const pr = solver.probeAt(wx, wy, solver.dx * 1.35);
          if (!pr || pr.mat === MAT.BEDROCK) continue;
          const m = pr.mohr;
          if (m.sigma1 < 1500) continue;

          const cx = ox + wx * scale;
          const cy = oy - wy * scale;
          const ang = -m.angleRad; // Screen Y is inverted relative to world Y
          const l1 = clamp(m.sigma1 / normRef, 0.15, 1.0) * maxLen;
          const l3 = clamp(Math.max(0, m.sigma3) / normRef, 0.08, 1.0) * maxLen;

          // Major principal stress sigma'_1 (warm amber-white)
          ctx.strokeStyle = 'rgba(255, 205, 115, 0.85)';
          ctx.beginPath();
          ctx.moveTo(cx - Math.cos(ang) * l1, cy - Math.sin(ang) * l1);
          ctx.lineTo(cx + Math.cos(ang) * l1, cy + Math.sin(ang) * l1);
          ctx.stroke();

          // Minor principal stress sigma'_3 (cool cyan)
          const ang3 = ang + Math.PI * 0.5;
          ctx.strokeStyle = 'rgba(78, 224, 255, 0.78)';
          ctx.beginPath();
          ctx.moveTo(cx - Math.cos(ang3) * l3, cy - Math.sin(ang3) * l3);
          ctx.lineTo(cx + Math.cos(ang3) * l3, cy + Math.sin(ang3) * l3);
          ctx.stroke();
        }
      }
      ctx.restore();
    }

    drawVelocityVectors(ctx, solver, ox, oy, scale) {
      const stride = solver.nx + 1;
      const skip = 3;
      ctx.save();
      ctx.strokeStyle = 'rgba(110, 245, 255, 0.78)';
      ctx.lineWidth = 1.25;
      ctx.beginPath();
      for (let gy = 3; gy < solver.ny - 2; gy += skip) {
        for (let gx = 3; gx < solver.nx - 2; gx += skip) {
          const idx = gy * stride + gx;
          if (solver.gMass[idx] < 1e-6) continue;
          const vx = solver.gVx[idx];
          const vy = solver.gVy[idx];
          const spd = Math.hypot(vx, vy);
          if (spd < 0.08) continue;
          const x0 = ox + gx * solver.dx * scale;
          const y0 = oy - gy * solver.dx * scale;
          const len = Math.min(24, spd * 6.5);
          const dx = (vx / spd) * len;
          const dy = -(vy / spd) * len;
          ctx.moveTo(x0, y0);
          ctx.lineTo(x0 + dx, y0 + dy);
        }
      }
      ctx.stroke();
      ctx.restore();
    }

    drawColorbar(ctx, x, y, width, meta) {
      ctx.save();
      ctx.fillStyle = 'rgba(8, 12, 22, 0.78)';
      ctx.strokeStyle = 'rgba(145, 180, 235, 0.22)';
      ctx.lineWidth = 1;
      ctx.fillRect(x, y, width, 44);
      ctx.strokeRect(x, y, width, 44);

      ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
      ctx.fillStyle = '#dce8fb';
      ctx.fillText(meta.name, x + 8, y + 14);

      if (meta.cmap === 'phase') {
        const items = [
          ['#946b44', 'Crust'], ['#4d7492', 'Intact Clay'],
          ['#22d8f8', 'Remolded'], ['#ff8224', 'Shear']
        ];
        let cx = x + 8;
        for (const [col, lbl] of items) {
          ctx.fillStyle = col;
          ctx.fillRect(cx, y + 23, 8, 8);
          ctx.fillStyle = '#a8bad6';
          ctx.fillText(lbl, cx + 11, y + 31);
          cx += lbl.length * 5.8 + 16;
        }
      } else {
        const lut = LUTS[meta.cmap] || LUTS.inferno;
        const barW = width - 16;
        for (let i = 0; i < barW; i++) {
          const idx = Math.min(255, Math.floor((i / barW) * 255));
          ctx.fillStyle = lut[idx];
          ctx.fillRect(x + 8 + i, y + 20, 1.5, 9);
        }
        ctx.fillStyle = '#9bb0d0';
        ctx.fillText('0', x + 8, y + 40);
        const mid = (meta.max * 0.5).toFixed(meta.max <= 1.5 ? 2 : 0);
        const hi = meta.max.toFixed(meta.max <= 1.5 ? 2 : 0) + (meta.unit !== '—' ? ' ' + meta.unit : '');
        ctx.fillText(mid, x + 8 + barW * 0.45, y + 40);
        ctx.fillText(hi, x + width - 8 - hi.length * 6, y + 40);
      }
      ctx.restore();
    }
  }

  ContinuumRenderer.FIELD_META = FIELD_META;
  root.Strata.ContinuumRenderer = ContinuumRenderer;
})(typeof globalThis !== 'undefined' ? globalThis : this);
