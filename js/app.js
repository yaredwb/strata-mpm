/*
 * STRATA — Interactive MLS-MPM Geomechanics Laboratory Controller
 */
(function (root) {
  'use strict';

  const C = root.Strata.Constitutive;
  const MPMSolver = root.Strata.MPMSolver;
  const ContinuumRenderer = root.Strata.ContinuumRenderer;
  const EngineeringHud = root.Strata.EngineeringHud;

  function $(id) {
    return document.getElementById(id);
  }

  const MAT_NAMES = ['Mohr-Coulomb Sand', 'Sensitive Quick Clay', 'Liquefiable Sand', 'Weathered Crust', 'Concrete Structure', 'Bedrock'];

  function init() {
    const simCanvas = $('sim-canvas');
    const mohrCanvas = $('mohr-canvas');
    const chartCanvas = $('chart-canvas');

    const solver = new MPMSolver(16000);
    const renderer = new ContinuumRenderer(simCanvas);
    const hud = new EngineeringHud(mohrCanvas, chartCanvas);

    const uiState = {
      theme: 'dark',
      paused: false,
      fieldId: 0,
      tool: 'probe', // 'probe' | 'excavate' | 'pour_sand' | 'pour_clay' | 'footing'
      brushRadius: 1.9,
      showStressTensors: false,
      showVelocities: false,
      showMesh: false,
      substepsPerFrame: 12,
      pointerInCanvas: false,
      pointerWorldX: 12.0,
      pointerWorldY: 7.0,
      probePinned: false,
      probeWorldX: 14.0,
      probeWorldY: 7.5
    };

    function applyTheme(nextTheme) {
      uiState.theme = nextTheme === 'light' ? 'light' : 'dark';
      document.documentElement.dataset.theme = uiState.theme;
      $('btn-theme').textContent = uiState.theme === 'light' ? 'Dark Mode' : 'Light Mode';
    }

    let activeProbe = null;
    let pointerDown = false;

    function syncScenarioUi(id) {
      const cfg = MPMSolver.SCENARIOS[id];
      if (!cfg) return;
      uiState.fieldId = cfg.defaultField;
      uiState.probePinned = false;
      uiState.probeWorldX = solver.domainW * 0.36;
      uiState.probeWorldY = solver.domainH * 0.40;

      for (const btn of $('scenario-strip').querySelectorAll('.pill')) {
        btn.classList.toggle('active', btn.dataset.scenario === id);
      }
      for (const btn of $('field-grid').querySelectorAll('.mode-btn')) {
        btn.classList.toggle('active', Number(btn.dataset.field) === uiState.fieldId);
      }

      // Sync parameter sliders
      $('ctrl-su0').value = Math.round(solver.params.su0 * 1e-3);
      $('out-su0').textContent = Math.round(solver.params.su0 * 1e-3) + ' kPa';
      $('ctrl-sens').value = solver.params.sensitivity;
      $('out-sens').textContent = 'St = ' + solver.params.sensitivity.toFixed(0);
      $('ctrl-phi').value = solver.params.phiDeg;
      $('out-phi').textContent = solver.params.phiDeg.toFixed(0) + '°';
      $('ctrl-cohesion').value = Math.round(solver.params.cohesion * 1e-3);
      $('out-cohesion').textContent = Math.round(solver.params.cohesion * 1e-3) + ' kPa';
      $('ctrl-youngs').value = (solver.params.youngs * 1e-6).toFixed(1);
      $('out-youngs').textContent = (solver.params.youngs * 1e-6).toFixed(1) + ' MPa';
      $('ctrl-pga').value = solver.params.seismicPgaG;
      $('out-pga').textContent = solver.params.seismicPgaG.toFixed(2) + ' g';

      $('scenario-hint').textContent = cfg.autoTriggerHint;
      $('btn-excavate-preset').classList.toggle('hidden', id !== 'quickclay' && id !== 'retaining');
      $('btn-excavate-preset').textContent = id === 'retaining' ? 'Excavate Passive Toe' : 'Excavate Slope Toe';
    }

    function selectTool(tool) {
      uiState.tool = tool;
      for (const btn of $('tool-strip').querySelectorAll('.tool-pill')) {
        btn.classList.toggle('active', btn.dataset.tool === tool);
      }
    }

    function applyActiveTool(wx, wy, isSingleClick) {
      if (uiState.tool === 'excavate') {
        solver.excavateCircle(wx, wy, uiState.brushRadius);
      } else if (uiState.tool === 'pour_sand') {
        solver.pourParticles(wx, wy, uiState.brushRadius * 0.72, C.MAT.SAND_MC);
      } else if (uiState.tool === 'pour_clay') {
        solver.pourParticles(wx, wy, uiState.brushRadius * 0.72, C.MAT.QUICK_CLAY);
      } else if (uiState.tool === 'footing' && isSingleClick) {
        solver.placeSurchargeBlock(wx, wy, 3.2, 2.0);
      }
    }

    // Bind Scenario Presets
    $('scenario-strip').addEventListener('click', (e) => {
      const btn = e.target.closest('.pill');
      if (!btn || !btn.dataset.scenario) return;
      solver.loadScenario(btn.dataset.scenario);
      syncScenarioUi(btn.dataset.scenario);
    });

    // Bind Field Visualization Modes
    $('field-grid').addEventListener('click', (e) => {
      const btn = e.target.closest('.mode-btn');
      if (!btn) return;
      uiState.fieldId = Number(btn.dataset.field);
      for (const b of $('field-grid').querySelectorAll('.mode-btn')) {
        b.classList.toggle('active', b === btn);
      }
    });

    // Bind Interactive Tools
    $('tool-strip').addEventListener('click', (e) => {
      const btn = e.target.closest('.tool-pill');
      if (!btn || !btn.dataset.tool) return;
      selectTool(btn.dataset.tool);
    });

    // Bind Action Buttons
    $('btn-excavate-preset').addEventListener('click', () => {
      solver.triggerExcavationPreset();
    });

    $('btn-earthquake').addEventListener('click', () => {
      solver.triggerEarthquake();
    });

    $('btn-reset').addEventListener('click', () => {
      solver.loadScenario(solver.scenarioId);
      syncScenarioUi(solver.scenarioId);
    });

    $('btn-theme').addEventListener('click', () => {
      applyTheme(uiState.theme === 'dark' ? 'light' : 'dark');
    });

    $('btn-pause').addEventListener('click', () => {
      uiState.paused = !uiState.paused;
      $('btn-pause').textContent = uiState.paused ? 'Resume' : 'Pause';
      $('btn-pause').classList.toggle('active', uiState.paused);
    });

    $('btn-tensors').addEventListener('click', () => {
      uiState.showStressTensors = !uiState.showStressTensors;
      $('btn-tensors').classList.toggle('active', uiState.showStressTensors);
    });

    $('btn-vectors').addEventListener('click', () => {
      uiState.showVelocities = !uiState.showVelocities;
      $('btn-vectors').classList.toggle('active', uiState.showVelocities);
    });

    $('btn-mesh').addEventListener('click', () => {
      uiState.showMesh = !uiState.showMesh;
      $('btn-mesh').classList.toggle('active', uiState.showMesh);
    });

    $('btn-snap').addEventListener('click', () => {
      const url = simCanvas.toDataURL('image/png');
      const a = document.createElement('a');
      a.href = url;
      a.download = 'strata-' + solver.scenarioId + '-t' + solver.time.toFixed(1) + 's.png';
      document.body.appendChild(a);
      a.click();
      a.remove();
    });

    $('btn-csv').addEventListener('click', () => {
      const rows = ['time_s,active_mps,kinetic_energy_J,max_speed_ms,max_plastic_strain,max_ru,remolded_pct,footing_settlement_m,footing_q_kPa,runout_ratio'];
      for (const m of solver.history) {
        rows.push([
          m.time.toFixed(4), m.activeCount, m.ke.toFixed(2), m.maxSpeed.toFixed(4),
          m.maxPlastic.toFixed(4), m.maxRu.toFixed(4), m.remoldedPct.toFixed(2),
          m.footingSettlement.toFixed(4), (m.footingQ * 1e-3).toFixed(2), m.runoutRatio.toFixed(4)
        ].join(','));
      }
      const blob = new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'strata-' + solver.scenarioId + '-telemetry.csv';
      document.body.appendChild(a);
      a.click();
      a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });

    // Bind Sliders
    $('ctrl-su0').addEventListener('input', (e) => {
      const kpa = Number(e.target.value);
      solver.params.su0 = kpa * 1e3;
      $('out-su0').textContent = kpa + ' kPa';
    });
    $('ctrl-sens').addEventListener('input', (e) => {
      const st = Number(e.target.value);
      solver.params.sensitivity = st;
      $('out-sens').textContent = 'St = ' + st.toFixed(0);
    });
    $('ctrl-phi').addEventListener('input', (e) => {
      const v = Number(e.target.value);
      solver.params.phiDeg = v;
      $('out-phi').textContent = v.toFixed(0) + '°';
    });
    $('ctrl-cohesion').addEventListener('input', (e) => {
      const kpa = Number(e.target.value);
      solver.params.cohesion = kpa * 1e3;
      $('out-cohesion').textContent = kpa + ' kPa';
    });
    $('ctrl-youngs').addEventListener('input', (e) => {
      const mpa = Number(e.target.value);
      solver.params.youngs = mpa * 1e6;
      solver.dt = solver.computeStableDt();
      $('out-youngs').textContent = mpa.toFixed(1) + ' MPa';
    });
    $('ctrl-pga').addEventListener('input', (e) => {
      const g = Number(e.target.value);
      solver.params.seismicPgaG = g;
      $('out-pga').textContent = g.toFixed(2) + ' g';
    });
    $('ctrl-brush').addEventListener('input', (e) => {
      const r = Number(e.target.value);
      uiState.brushRadius = r;
      $('out-brush').textContent = r.toFixed(1) + ' m';
    });

    // Canvas Pointer Interaction
    simCanvas.addEventListener('pointerdown', (e) => {
      pointerDown = true;
      simCanvas.setPointerCapture(e.pointerId);
      const wpt = renderer.screenToWorld(e.clientX, e.clientY);
      uiState.pointerWorldX = wpt.x;
      uiState.pointerWorldY = wpt.y;
      if (uiState.tool === 'probe') {
        uiState.probePinned = !uiState.probePinned;
        uiState.probeWorldX = wpt.x;
        uiState.probeWorldY = wpt.y;
      } else {
        applyActiveTool(wpt.x, wpt.y, true);
      }
    });

    simCanvas.addEventListener('pointermove', (e) => {
      uiState.pointerInCanvas = true;
      const wpt = renderer.screenToWorld(e.clientX, e.clientY);
      uiState.pointerWorldX = wpt.x;
      uiState.pointerWorldY = wpt.y;
      if (uiState.tool === 'probe' && !uiState.probePinned) {
        uiState.probeWorldX = wpt.x;
        uiState.probeWorldY = wpt.y;
      } else if (pointerDown && uiState.tool !== 'footing') {
        applyActiveTool(wpt.x, wpt.y, false);
      }
    });

    simCanvas.addEventListener('pointerup', (e) => {
      pointerDown = false;
      if (simCanvas.hasPointerCapture(e.pointerId)) {
        simCanvas.releasePointerCapture(e.pointerId);
      }
    });

    simCanvas.addEventListener('pointerleave', () => {
      uiState.pointerInCanvas = false;
      pointerDown = false;
    });

    window.addEventListener('keydown', (e) => {
      if (e.target && e.target.tagName === 'INPUT') return;
      if (e.code === 'Space') {
        e.preventDefault();
        $('btn-pause').click();
      } else if (e.key === '1') selectTool('probe');
      else if (e.key === '2') selectTool('excavate');
      else if (e.key === '3') selectTool('pour_sand');
      else if (e.key === '4') selectTool('pour_clay');
      else if (e.key === '5') selectTool('footing');
      else if (e.key === 'r' || e.key === 'R') $('btn-reset').click();
      else if (e.key === 't' || e.key === 'T') $('btn-theme').click();
    });

    function updateTelemetry(metrics, probe) {
      $('val-time').textContent = metrics.time.toFixed(2) + ' s';
      $('val-mps').textContent = String(metrics.activeCount);
      $('val-dt').textContent = (solver.dt * 1e3).toFixed(2) + ' ms';
      $('val-maxeps').textContent = metrics.maxPlastic.toFixed(3);
      $('val-maxru').textContent = metrics.maxRu.toFixed(2);
      $('val-remold').textContent = metrics.remoldedPct.toFixed(1) + '%';

      if (probe) {
        $('probe-mat').textContent = MAT_NAMES[probe.mat] || 'Continuum';
        $('val-sig1').textContent = (probe.mohr.sigma1 * 1e-3).toFixed(1) + ' kPa';
        $('val-sig3').textContent = (probe.mohr.sigma3 * 1e-3).toFixed(1) + ' kPa';
        $('val-peff').textContent = (probe.pEff * 1e-3).toFixed(1) + ' kPa';
        $('val-qshear').textContent = (probe.qShear * 1e-3).toFixed(1) + ' kPa';
        $('val-upore').textContent = (probe.uPore * 1e-3).toFixed(1) + ' kPa (' + probe.ru.toFixed(2) + ')';
        $('val-mobphi').textContent = probe.mat === C.MAT.QUICK_CLAY
          ? 's_u ' + (probe.suCur * 1e-3).toFixed(1) + ' kPa'
          : probe.mohr.mobPhiDeg.toFixed(1) + '°';
      } else {
        $('probe-mat').textContent = 'Hover Soil';
        $('val-sig1').textContent = '—';
        $('val-sig3').textContent = '—';
        $('val-peff').textContent = '—';
        $('val-qshear').textContent = '—';
        $('val-upore').textContent = '—';
        $('val-mobphi').textContent = '—';
      }
    }

    applyTheme('dark');
    syncScenarioUi('quickclay');

    let lastMs = performance.now();
    let fpsAccum = 0;
    let fpsCount = 0;

    function tick(now) {
      const dtReal = Math.min(0.05, (now - lastMs) * 0.001);
      lastMs = now;

      if (!uiState.paused) {
        for (let s = 0; s < uiState.substepsPerFrame; s++) {
          solver.step();
        }
      }

      activeProbe = solver.probeAt(uiState.probeWorldX, uiState.probeWorldY);
      renderer.render(solver, uiState, activeProbe);
      hud.drawMohrCircle(solver, activeProbe, uiState.theme);
      hud.drawStripChart(solver, uiState.theme);

      const metrics = solver.computeMetrics();
      updateTelemetry(metrics, activeProbe);

      fpsAccum += dtReal;
      fpsCount++;
      if (fpsAccum >= 0.45) {
        $('val-fps').textContent = Math.round(fpsCount / fpsAccum) + ' FPS';
        fpsAccum = 0;
        fpsCount = 0;
      }

      requestAnimationFrame(tick);
    }

    requestAnimationFrame(tick);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})(typeof globalThis !== 'undefined' ? globalThis : this);
