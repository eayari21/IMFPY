# IMFPY Dust Trajectory Toolkit

This repository now ships a cross-platform PyQt6 application that drives a Fortran integrator to study interstellar dust grain trajectories. The GUI makes it easy to configure ensembles of particles, choose between CPU and GPU backends, and quickly visualise the resulting trajectories.

## Features

- **High-performance Fortran backend** compiled via `numpy.f2py` with OpenMP support for multi-core CPUs.
- **CuPy GPU backend** (optional) for running the same equations of motion on supported NVIDIA hardware.
- **NumPy/Numba reference backend** for validation, testing, and fallback when accelerators are unavailable.
- **Interactive GUI** written with PyQt6 that lets you:
  - Configure particle counts, initial state vectors, masses, charges, and radiation pressure coefficients (β).
  - Adjust integration length, time step, gravitational parameter (GM), and interplanetary electric/magnetic field vectors.
  - Launch simulations in a background worker while keeping the interface responsive.
  - Visualise multi-particle trajectories in 3D with Matplotlib.

## Getting Started

1. **Install dependencies**

   ```bash
   python -m pip install -r requirements.txt  # create this file if you package the project
   ```

   Key Python dependencies include `numpy`, `scipy`, `matplotlib`, `PyQt6`, and optionally `cupy` for GPU support. A Fortran compiler such as `gfortran` is required to build the backend.

2. **Launch the GUI**

   ```bash
   python -m imfpy.gui.main
   ```

   On first run the Fortran source will be compiled automatically. Subsequent launches reuse the cached shared library under `imfpy/fortran/_build/`.

3. **Command line runner**

   If you prefer automated workflows, you can run the same engine via CLI:

   ```bash
   python -m imfpy.simulation.runner --backend fortran --particles 32 --steps 5000 --dt 5.0
   ```

## Repository Layout

- `imfpy/fortran/` – Fortran RK4 integrator and build helpers.
- `imfpy/simulation/` – Python orchestration layer, configuration dataclasses, and backends.
- `imfpy/gui/` – PyQt6 graphical application.
- `tests/` – Space for automated tests (add your own as the project grows).

## Development Tips

- The Fortran module is automatically rebuilt if you delete the contents of `imfpy/fortran/_build/` or pass `force=True` to `imfpy.fortran.build.build_module`.
- Use the Python backend during algorithm prototyping; once satisfied switch to the Fortran backend for production-scale runs.
- The GPU backend requires CuPy compiled against the CUDA toolkit version available on your system. If CuPy is missing the GUI will fall back to the CPU options.

## License

MIT (or specify the actual licence applicable to your project).
