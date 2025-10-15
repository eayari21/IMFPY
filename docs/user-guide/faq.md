# Frequently asked questions

## General

### What problems does IMFPY solve?
IMFPY simulates charged dust grain trajectories in the heliosphere or interstellar medium by integrating the coupled gravitational, radiation pressure, and Lorentz forces.  It is useful for mission design studies, astrophysics research, and classroom demonstrations.

### Do I need the Fortran backend?
No.  The Python backend works everywhere.  The Fortran backend provides a significant speed boost for large ensembles but requires a compatible compiler.

### Does IMFPY ship with sample data?
Yes.  The GUI boots with a circular orbit ensemble and the tutorials provide additional configuration snippets.  For more complex scenarios, use the CLI to load YAML/JSON files.

## Installation & configuration

### Which Python versions are supported?
Python 3.10 through 3.12.  Earlier versions are untested and may fail to import PyQt6.

### How do I install GPU support?
Install a CuPy wheel that matches your CUDA toolkit (e.g. `cupy-cuda12x`).  Then open **Help → Diagnostics** to confirm that the GPU backend is recognised.

### Can I run IMFPY inside Docker or on a cluster?
Yes.  Use the CLI for headless environments.  Provide a `--backend python` flag if CUDA/Fortran toolchains are unavailable.  When running the GUI on remote clusters, use X forwarding or virtual desktops.

## Usage

### How do I reset the GUI to defaults?
Click **Reset** in the control panel.  This repopulates the particle table with a four-particle circular orbit and restores default global parameters.

### Where are simulation results stored?
CLI runs write `.npz` archives wherever you point `--output`.  In the GUI, exported runs default to the last-used directory.  The in-memory result is accessible until you start a new simulation.

### Can I run multiple simulations at once?
The GUI executes one simulation at a time to keep the interface responsive.  Use the CLI in combination with GNU Parallel or your scheduler to run multiple jobs concurrently.

### How do I change the integrator?
IMFPY currently ships with a Runge-Kutta 4 solver.  Advanced users can extend `imfpy/fortran` or implement new Python backends.  See `imfpy/simulation/runner.py` for the dispatch logic.

## Troubleshooting

### The GUI froze during a run.
Check the terminal window for stack traces.  If the worker crashed, the status bar will show “Error”.  Press **Help → Quick Start Guide** to review correct input ranges and rerun with a smaller time step.

### The Fortran build keeps recompiling.
Delete the `imfpy/fortran/_build` directory to force a clean build.  Ensure the Python process has permission to write to this directory.

### GPU runs are slower than CPU runs.
For small particle counts (< 1000), CPU backends often outperform GPUs due to launch overhead.  Increase the ensemble size or stick with the Fortran backend.

## Project & contributing

### Where can I request features?
Open an issue on GitHub with the `enhancement` label.  Outline the use case, proposed behaviour, and any relevant references.

### How do I contribute documentation?
Edit the files in `docs/user-guide/` and submit a pull request.  Screenshots should be added to `media/figures/` and referenced with relative paths.

### Is there a roadmap?
Yes!  Check the GitHub Projects board linked from the repository README.  Contributions and feedback help shape future releases.
