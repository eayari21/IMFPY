# Installation

Whether you are a researcher exploring interstellar dust dynamics or a student completing a class project, IMFPY is designed to install cleanly on Windows, macOS, and Linux.  This page walks you through prerequisites, step-by-step installation, and optional GPU acceleration.

---

## Prerequisites

| Requirement | Notes |
| --- | --- |
| Python 3.10 – 3.12 | Install from [python.org](https://www.python.org/downloads/) or your package manager. |
| Compiler toolchain | `gfortran` for Linux/macOS, MSYS2 `mingw-w64-gfortran` for Windows. |
| Build tools | Ensure `pip`, `setuptools`, and `wheel` are up to date: `python -m pip install --upgrade pip setuptools wheel`. |
| Optional CUDA | CUDA 11.8+ with matching CuPy wheels for GPU acceleration. |

> **No compiler?** The Python backend still works.  Install the requirements and choose the **python** backend inside the GUI or CLI.

---

## Step-by-step installation

### 1. Clone the repository (developers)

```bash
git clone https://github.com/<your-org>/IMFPY.git
cd IMFPY
```

> If you are installing from a release archive, unzip it and `cd` into the extracted folder instead.

### 2. Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
```

### 3. Install dependencies

```bash
python -m pip install -r requirements.txt
```

This pulls in PyQt6, NumPy/SciPy, Matplotlib, and backend-specific dependencies.  To avoid binary compatibility issues ensure that your Python installation and `pip` are the latest versions for your platform.

### 4. Verify your compiler (optional)

```bash
gfortran --version
```

If the command is not found, install the toolchain:

- **Ubuntu/Debian:** `sudo apt install gfortran libopenblas-dev python3-dev`
- **Fedora:** `sudo dnf install gcc-gfortran python3-devel openblas-devel`
- **macOS:** `brew install gcc`
- **Windows:** install [MSYS2](https://www.msys2.org) and then `pacman -S mingw-w64-x86_64-gcc` from the MSYS2 MinGW shell.

### 5. Launch IMFPY

```bash
python -m imfpy.gui.main
```

The first launch builds the Fortran shared library (if a compiler is present).  Progress and errors appear both in the terminal and inside the GUI status bar.

---

## Optional GPU support

IMFPY can offload integrations to NVIDIA GPUs via CuPy.  Install the wheel that matches your CUDA toolkit.

| CUDA version | Install command |
| --- | --- |
| CUDA 11.x | `python -m pip install cupy-cuda11x` |
| CUDA 12.x | `python -m pip install cupy-cuda12x` |
| ROCm (experimental) | Follow the [CuPy ROCm docs](https://docs.cupy.dev/en/stable/install.html#install-cupy-for-rocm-platforms). |

After installation run:

```bash
python -m imfpy.simulation.runner --list-backends
```

You should see `gpu` listed.  Inside the GUI open **Help → Diagnostics** to confirm GPU availability.

---

## Keeping IMFPY up to date

1. Pull the latest changes: `git pull` (or download the newest release archive).
2. Activate your virtual environment.
3. Re-run `python -m pip install -r requirements.txt` to capture dependency updates.
4. Clear the Fortran build cache if you suspect stale binaries: delete `imfpy/fortran/_build/`.

---

## Uninstalling

- Deactivate and delete the virtual environment.
- Remove the project folder.
- Optional: delete cached builds under your system’s temporary directory (e.g., `%TEMP%` on Windows).

That’s it!  You are ready to dive into the tutorials and start exploring dust trajectories.
