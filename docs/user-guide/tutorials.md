# Tutorial library

This page provides step-by-step tutorials for common IMFPY workflows.  Each tutorial includes goals, prerequisites, and verification checkpoints so you can confirm that you are on track.

---

## Tutorial 1 – GUI deep dive

**Goal:** Understand every control in the desktop application and run a custom simulation.

### Prerequisites
- IMFPY installed (see [installation](installation.md)).
- Optional: a game-ready mouse for 3D navigation.

### Steps

1. **Launch the GUI:** `python -m imfpy.gui.main`.
2. **Explore the sidebar:** Hover over each field to read the built-in tooltips.  Open **Help → Quick Start Guide** if you need a refresher on any control.
3. **Adjust global parameters:**
   - Set `Particles` to 8.
   - Change `Steps` to 8000 and `Δt` to 20 seconds.
   - Enter `1.5e20` for `GM` to approximate a heavier central body.
4. **Edit the ensemble:**
   - Click **Populate Circular Orbit** to refresh the defaults.
   - Manually tweak one particle’s velocity to create a retrograde orbit.
5. **Apply fields:**
   - Set an electric field of `(2e-4, 0, 0)` V/m.
   - Set a magnetic field of `(0, 0, 5e-9)` T.
6. **Run & review:**
   - Press **Run Simulation**.
   - When finished, rotate the 3D view and inspect the legend to identify your retrograde particle.
   - Check the status bar for the backend that executed the run.
7. **Save a snapshot:** Use **File → Export Result…** to create `elliptical_run.npz`.

### Verification
- The 3D viewer should show diverging trajectories (one retrograde).
- The status bar message ends with “Simulation finished (backend=<name>)”.

---

## Tutorial 2 – Physics validation workflow

**Goal:** Compare the Fortran, Python, and GPU backends for numerical agreement.

### Prerequisites
- All backends installed (`python -m imfpy.simulation.runner --list-backends`).
- CLI familiarity.

### Steps

1. Generate a baseline config:
   ```bash
   python -m imfpy.simulation.runner \
       --backend python \
       --particles 4 \
       --steps 4000 \
       --dt 60 \
       --output baseline_python.npz
   ```
2. Re-run with the Fortran backend:
   ```bash
   python -m imfpy.simulation.runner \
       --backend fortran \
       --particles 4 \
       --steps 4000 \
       --dt 60 \
       --output baseline_fortran.npz
   ```
3. Re-run with the GPU backend (if available) and save as `baseline_gpu.npz`.
4. Compare results in a notebook:
   ```python
   import numpy as np

   def load(path):
       data = np.load(path)
       return data["positions"], data["velocities"]

   p_pos, p_vel = load("baseline_python.npz")
   f_pos, f_vel = load("baseline_fortran.npz")

   print("Max |Δpos|:", np.max(np.abs(p_pos - f_pos)))
   print("Max |Δvel|:", np.max(np.abs(p_vel - f_vel)))
   ```
5. Investigate discrepancies by reducing the time step or adjusting OpenMP thread counts.

### Verification
- Maximum position and velocity differences should be within numerical tolerances (< 1e-6 for identical initial conditions).
- If GPU results diverge, consult [troubleshooting](troubleshooting.md#gpu-specific-issues).

---

## Tutorial 3 – Automating with the CLI

**Goal:** Run a parameter sweep headlessly and aggregate the results.

### Prerequisites
- Familiarity with shell scripting.
- Optional: GNU Parallel or similar tooling for concurrency.

### Steps

1. Create a configuration YAML (`configs/sweep.yaml`):
   ```yaml
   runs:
     - name: beta-low
       particles: 8
       steps: 6000
       dt: 30
       beta: 0.1
     - name: beta-mid
       particles: 8
       steps: 6000
       dt: 30
       beta: 0.5
     - name: beta-high
       particles: 8
       steps: 6000
       dt: 30
       beta: 0.9
   ```
2. Write a helper script (`scripts/run_sweep.sh`):
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail

   BACKEND=${1:-fortran}
   CONFIG=configs/sweep.yaml

   python - <<'PY'
   import yaml
   from pathlib import Path
   from subprocess import run

   cfg = yaml.safe_load(Path("configs/sweep.yaml").read_text())
   for run_cfg in cfg["runs"]:
       name = run_cfg["name"]
       cmd = [
           "python", "-m", "imfpy.simulation.runner",
           "--backend", "${BACKEND}",
           "--particles", str(run_cfg["particles"]),
           "--steps", str(run_cfg["steps"]),
           "--dt", str(run_cfg["dt"]),
           "--beta", str(run_cfg["beta"]),
           "--output", f"outputs/{name}.npz",
       ]
       print("Running", " ".join(cmd))
       run(cmd, check=True)
   PY
   ```
3. Run the sweep and monitor outputs in `outputs/`.
4. Post-process results in a notebook or script to build plots comparing β values.

### Verification
- Each run produces an `.npz` file with the expected naming scheme.
- Scripts exit with code 0 (no exceptions printed).

---

## Tutorial 4 – Embedding IMFPY in notebooks

**Goal:** Use the Python API directly for custom analysis.

```python
from imfpy.simulation import ParticleEnsemble, SimulationConfig, SimulationRunner
import numpy as np

positions = np.zeros((4, 3))
velocities = np.zeros((4, 3))
charges = np.full(4, 1.6e-19)
masses = np.full(4, 1e-18)
betas = np.full(4, 0.5)

ensemble = ParticleEnsemble(
    positions=positions,
    velocities=velocities,
    charges=charges,
    masses=masses,
    betas=betas,
)

config = SimulationConfig(
    particles=ensemble,
    num_steps=2000,
    time_step=60.0,
    gm_sun=1.32712440018e20,
    electric_field=(0.0, 0.0, 0.0),
    magnetic_field=(0.0, 0.0, 5e-9),
    backend="python",
)

runner = SimulationRunner()
result = runner.run(config)

print("Trajectory array shape:", result.positions().shape)
```

See the `notebooks/` directory for complete Jupyter examples, including visualisation recipes with Plotly and animation exports.

---

## What next?

- Run through the [troubleshooting matrix](troubleshooting.md) if you encounter build/runtime issues.
- Consult the [FAQ](faq.md) for numerical stability tips and backend selection advice.
- Share your own tutorials!  Pull requests with walkthroughs, videos, or screenshots are very welcome.
