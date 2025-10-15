# IMFPY Dust Trajectory Toolkit

IMFPY is an end-to-end toolkit for exploring interstellar and heliospheric dust grain dynamics.  
It combines a high-performance Fortran integrator with pure-Python and GPU backends and ships a full-featured PyQt6 desktop application for crafting ensembles, running simulations, and visualising the resulting trajectories.

---

## Table of contents

1. [Features at a glance](#features-at-a-glance)
2. [System requirements](#system-requirements)
3. [Installation](#installation)
   - [Python virtual environment](#1-create-a-python-virtual-environment)
   - [Install IMFPY](#2-install-imfpy)
   - [Optional GPU support](#optional-gpu-support)
4. [Quick start](#quick-start)
   - [Launch the desktop app](#launch-the-desktop-app)
   - [Command line runner](#command-line-runner)
5. [Walkthrough tutorials](#walkthrough-tutorials)
6. [Project layout](#project-layout)
7. [Troubleshooting & diagnostics](#troubleshooting--diagnostics)
8. [Additional documentation](#additional-documentation)
9. [Contributing](#contributing)
10. [License](#license)

---

## Features at a glance

- **High-performance Fortran backend** compiled through `numpy.f2py` with OpenMP support for multi-core CPUs.
- **CuPy GPU backend** (optional) that mirrors the CPU equations of motion on compatible NVIDIA hardware.
- **NumPy/Numba reference backend** for validation, testing, and execution on systems without a Fortran compiler.
- **Interactive PyQt6 GUI** with:
  - Particle ensemble editor (positions, velocities, masses, charges, radiation-pressure β coefficients).
  - Global integration controls (time step, step count, gravitational parameter, electric & magnetic field vectors).
  - Background execution worker to keep the interface responsive.
  - Built-in help viewer, quick-start overlays, and direct links to the full documentation.
  - 3D Matplotlib trajectory viewer with mouse navigation tools.
- **Automation-friendly CLI** for batch experiments, CI pipelines, or headless environments.
- **Extensive documentation** (installation, tutorials, troubleshooting, FAQs) in [`docs/user-guide`](docs/user-guide/index.md) and mirrored inside the application itself.

---

## System requirements

| Component | Recommendation |
| --- | --- |
| Operating system | Linux, macOS 12+, or Windows 10/11 |
| Python | 3.10 – 3.12 |
| CPU backend | Any modern x86-64/ARM64 CPU with OpenMP support |
| GPU backend (optional) | NVIDIA GPU with CUDA 11.8+ and a CuPy build that matches your CUDA toolkit |
| Memory | ≥ 8 GB RAM recommended for large ensembles |
| Additional tooling | `gfortran` (or compatible) for building the native backend |

> **Tip:** If you only need the pure Python backend you can skip the Fortran toolchain and still run the GUI.

---

## Installation

### 1. Create a Python virtual environment

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
python -m pip install --upgrade pip
```

### 2. Install IMFPY

Install the project dependencies and build the Fortran extension on first launch.

```bash
python -m pip install -r requirements.txt
```

On Linux you may need the system packages `gfortran`, `libopenblas-dev`, and `python3-dev`.  macOS users can install `gcc` via Homebrew (`brew install gcc`).  Windows users can install the [MSYS2](https://www.msys2.org) toolchain.

After installing the requirements you can launch the GUI or CLI directly (see [Quick start](#quick-start)).  The Fortran backend is auto-built the first time you invoke it and the compiled module is cached under `imfpy/fortran/_build/`.

### Optional GPU support

To enable the GPU backend install a CuPy build that matches your CUDA toolkit.  Example for CUDA 11.8:

```bash
python -m pip install cupy-cuda11x
```

You can confirm GPU availability from the GUI status bar (`Help → Diagnostics`) or via the CLI:

```bash
python -m imfpy.simulation.runner --list-backends
```

---

## Quick start

### Launch the desktop app

```bash
python -m imfpy.gui.main
```

What happens on first launch:

1. IMFPY compiles the Fortran integrator (if available) and validates the Python/GPU backends.
2. The main window opens with a pre-populated four-particle circular orbit ensemble.
3. Press **Run Simulation** to execute with the default parameters.  Results appear in the 3D viewer and the status bar will indicate completion.
4. Open the **Help → Quick Start Guide** menu (or press `F1`) for a step-by-step tour of every control.

### Command line runner

The CLI exposes the same solver stack for scripted runs:

```bash
python -m imfpy.simulation.runner \
    --backend fortran \
    --particles 32 \
    --steps 5000 \
    --dt 5.0 \
    --output trajectories.npz
```

Use `--help` for the full set of options and see [`docs/user-guide/tutorials.md`](docs/user-guide/tutorials.md#automating-with-the-cli) for complete automation workflows.

---

## Walkthrough tutorials

- **Guided tour of the GUI** – learn how to configure particle ensembles, tweak global fields, and interpret plots.
- **Physics validation workflow** – cross-check Fortran, Python, and GPU backends using the same initial conditions.
- **Batch experiments** – create YAML/JSON configs, run them from the CLI, and post-process the results in notebooks.

Each tutorial includes screenshots, step-by-step instructions, and troubleshooting checkpoints.  Start with [`docs/user-guide/tutorials.md`](docs/user-guide/tutorials.md).

---

## Project layout

```text
imfpy/
├── fortran/        # Native integrator sources and build helpers
├── gui/            # PyQt6 application (widgets, windows, help system)
└── simulation/     # Backend orchestration, data classes, CLI runner
```

Additional directories:

- `docs/user-guide/` – Markdown documentation surfaced on GitHub and referenced by the in-app help viewer.
- `media/` – Figures used in talks, documentation, and presentations.
- `notebooks/` – Example scientific workflows (exploratory analysis, plotting).
- `tests/` – Space for automated unit tests.

---

## Troubleshooting & diagnostics

1. Launch the app and open **Help → Diagnostics** for a live system report (Python version, backend availability, CUDA status).
2. Review the troubleshooting matrix in [`docs/user-guide/troubleshooting.md`](docs/user-guide/troubleshooting.md) for fixes to common compiler, GUI, and runtime issues.
3. Run the CLI self-test to validate all backends:

   ```bash
   python -m imfpy.simulation.runner --self-test
   ```

4. Still stuck?  File an issue with logs, OS details, and the contents of `imfpy/fortran/_build/build.log`.

---

## Additional documentation

- [Full user guide](docs/user-guide/index.md)
- [Installation deep dive](docs/user-guide/installation.md)
- [Tutorial library](docs/user-guide/tutorials.md)
- [Troubleshooting & diagnostics](docs/user-guide/troubleshooting.md)
- [Frequently asked questions](docs/user-guide/faq.md)

These files are packaged with the application and can be opened from inside the GUI via **Help → Open Documentation Folder**.

---

## Contributing

1. Fork and clone the repository.
2. Create a new branch for your feature or fix.
3. Install development dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

4. Run the linters/tests (add your own in `tests/`).
5. Submit a pull request summarising the change and linking relevant issues.

---

## License

MIT (or specify the actual licence applicable to your project).
