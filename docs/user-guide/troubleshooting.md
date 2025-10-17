# Troubleshooting & diagnostics

This page lists common issues and actionable fixes.  When reporting problems, include your OS, Python version, and the log snippets referenced below.

---

## Quick diagnostic commands

```bash
python -m imfpy.simulation.runner --self-test   # Exercises every available backend
python -m imfpy.simulation.runner --list-backends
python -c "import imfpy, sys; print(sys.version); print(imfpy.__version__)"
```

Inside the GUI, open **Help → Diagnostics** to see the same information in a dialog.

---

## Installation issues

| Symptom | Likely cause | Resolution |
| --- | --- | --- |
| `ModuleNotFoundError: No module named 'PyQt6'` | Dependencies missing | Re-run `python -m pip install -r requirements.txt` in the active virtual environment. |
| `error: f951: Fatal Error: no input files` | Fortran compiler not found | Install `gfortran` (see [installation guide](installation.md#step-by-step-installation)); restart the shell so your `PATH` updates. |
| `numpy.f2py` build fails with OpenMP errors | Compiler lacks OpenMP support | Install a newer GCC/Clang build. On macOS use `brew install gcc` and set `export FC=gfortran`. |
| CuPy wheel fails to install | CUDA mismatch | Choose the wheel that matches your CUDA toolkit (e.g. `cupy-cuda12x` for CUDA 12). |

---

## GUI/runtime issues

| Symptom | Likely cause | Resolution |
| --- | --- | --- |
| GUI does not launch; errors mention `xcb` (Linux) | Missing Qt platform plugins | Install system package `libxcb-xinerama0` (Debian/Ubuntu) or `qt6-base` (Arch). Alternatively, set `QT_DEBUG_PLUGINS=1` and inspect output. |
| Simulation never finishes; status bar stuck on *Running…* | Backend exception in worker thread | Check the terminal for stack traces.  Open **Help → Diagnostics** and run `python -m imfpy.simulation.runner --self-test` to isolate the failing backend. |
| 3D viewer blank | Matplotlib backend not initialised | Ensure `matplotlib` is installed (comes via requirements).  Try running `python -c "import matplotlib"` to confirm. |
| Particle table does not update when changing counts | Pending validation error | Ensure all numeric fields contain valid numbers (no empty cells).  Press **Reset** to repopulate defaults. |
| `TypeError: scaledToHeight(... AspectRatioMode)` when opening HDF plotter | Qt 6 tightened the `QPixmap.scaledToHeight` signature; passing an aspect ratio enum now raises | Replace the call with `pixmap.scaledToHeight(height, Qt.TransformationMode.SmoothTransformation)` or use `pixmap.scaled(width, height, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)` so the last argument is a `Qt.TransformationMode`. |

---

## CLI issues

| Symptom | Likely cause | Resolution |
| --- | --- | --- |
| `ValueError: Unknown backend 'gpu'` | CuPy not installed or incompatible | Install the correct CuPy wheel and re-run `--list-backends`. |
| `FileNotFoundError` when using `--output` | Directory missing | Create the output directory first (`mkdir -p outputs`). |
| `yaml` module missing in automation tutorial | PyYAML not installed | `python -m pip install pyyaml` (optional dependency for batch scripts). |

---

## GPU-specific issues

| Symptom | Likely cause | Resolution |
| --- | --- | --- |
| `cupy.cuda.runtime.CUDARuntimeError: cudaErrorInsufficientDriver` | Driver/toolkit mismatch | Update NVIDIA drivers and ensure CUDA runtime >= toolkit required by CuPy. |
| GPU backend slower than CPU | Small particle counts | GPUs shine for ≥ 1e3 particles.  For small ensembles stick with CPU/Fortran. |
| GPU backend missing in diagnostics | CuPy import failed | Run `python -c "import cupy; cupy.show_config()"` for detailed errors. |

---

## Logging & support

- GUI logs: check the terminal that launched the app.  Worker exceptions appear there with full tracebacks.
- Fortran build logs: `imfpy/fortran/_build/build.log`.
- CLI logs: run with `--verbose` for extra diagnostics.

When opening a GitHub issue include:

1. Operating system & version.
2. Python version (`python --version`).
3. IMFPY commit hash or release tag.
4. Backend used (`fortran`, `python`, or `gpu`).
5. Steps to reproduce and relevant logs.

We are happy to help – the more details you provide, the faster we can reproduce the bug.
