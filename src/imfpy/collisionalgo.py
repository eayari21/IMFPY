from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

"""
Python script to process FORTRAN integrator data simulating charged particles
in the Interplanetary Magnetic Field (IMF).
__author__      = Ethan Ayari,
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust
"""


# %%CALL FORTRAN 95 SCRIPT
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
FORTRAN_DIR = REPO_ROOT / "src" / "fortran"
FORTRAN_SOURCE = FORTRAN_DIR / "earthcopy.f95"
FORTRAN_BINARY = FORTRAN_DIR / "a.out"
FORTRAN_OUTPUT = FORTRAN_DIR / "subtest.8"


def _prepare_fortran_source(q: str, m: str, pitch: str, x: str, y: str, z: str, KE: str) -> str:
    """Return the modified FORTRAN source ready for compilation."""

    if not FORTRAN_SOURCE.exists():
        raise FileNotFoundError(f"Could not locate FORTRAN source at {FORTRAN_SOURCE}")

    scriptlines = FORTRAN_SOURCE.read_text().splitlines(keepends=True)

    scriptlines[10] = "       l_p=" + pitch + "\n"
    scriptlines[53] = "      l_p=" + pitch + "\n"
    scriptlines[11] = "       pm=" + m + "\n"
    scriptlines[9] = "       q=" + q + "\n"
    scriptlines[40] = "      y=" + y + "\n"
    scriptlines[41] = "      z=" + z + "\n"
    scriptlines[39] = "      x=" + x + "\n"
    scriptlines[43] = "      xe=" + KE + "\n"

    return "".join(scriptlines)


def _compile_fortran(source_code: str) -> None:
    """Compile the modified FORTRAN source in-place using gfortran."""

    FORTRAN_DIR.mkdir(parents=True, exist_ok=True)
    original_source = FORTRAN_SOURCE.read_text()

    try:
        FORTRAN_SOURCE.write_text(source_code)

        gfortran = shutil.which("gfortran")
        if gfortran is None:
            raise RuntimeError("gfortran compiler not found on PATH")

        subprocess.run([gfortran, FORTRAN_SOURCE.name], cwd=FORTRAN_DIR, check=True)
    finally:
        FORTRAN_SOURCE.write_text(original_source)


def _run_fortran_binary() -> None:
    """Execute the compiled FORTRAN binary within the FORTRAN directory."""

    if not FORTRAN_BINARY.exists():
        raise RuntimeError("FORTRAN binary not found; ensure compilation succeeded")

    subprocess.run([str(FORTRAN_BINARY)], cwd=FORTRAN_DIR, check=True)


def _cleanup_fortran_artifacts() -> None:
    """Remove build artefacts produced by the FORTRAN workflow."""

    if FORTRAN_BINARY.exists():
        FORTRAN_BINARY.unlink()
    if FORTRAN_OUTPUT.exists():
        FORTRAN_OUTPUT.unlink()


def fly_particle(q: str, m: str, pitch: str, x: str, y: str, z: str, KE: str) -> pd.DataFrame:
    """Run the FORTRAN integrator and return the trajectory as a DataFrame."""

    source_code = _prepare_fortran_source(q, m, pitch, x, y, z, KE)
    _compile_fortran(source_code)
    try:
        _run_fortran_binary()
        data = pd.read_csv(FORTRAN_OUTPUT, header=None, delim_whitespace=True)

        raw_output = REPO_ROOT / "data" / "raw" / FORTRAN_OUTPUT.name
        raw_output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(FORTRAN_OUTPUT, raw_output)
    finally:
        _cleanup_fortran_artifacts()

    return data.transpose()


# %%PROCESS OUTPUT FROM FORTRAN 95 SCRIPT
def process_data(data, col, lab):
    data = data  # .transpose()
    # t = data.iloc[0]
    x = data.iloc[1]
    y = data.iloc[2]
    z = data.iloc[3]
    # vx = data.iloc[4]
    # vy = data.iloc[5]
    # vz = data.iloc[6]

    u = np.linspace(0, np.pi, 30)
    v = np.linspace(0, 2 * np.pi, 30)
    xsphere = 2.5*np.outer(np.sin(u), np.sin(v))
    ysphere = 2.5*np.outer(np.sin(u), np.cos(v))
    zsphere = 2.5*np.outer(np.cos(u), np.ones_like(v))

    # magv = np.sqrt(vx**2 + vy**2 + vz**2)
    """
    plt.plot((x+y)**2,z)
    plt.xlabel("Radius")
    plt.ylabel("Z")
    plt.show()
    """
    # mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = Axes3D(fig)
    # ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(xsphere, ysphere, zsphere, color='gray',
                      label='Aluminum (Al) Sphere')
    ax.plot(x, y, z, c='r')
    ax.set_xlabel('X', fontsize=20, fontweight='bold')
    ax.set_ylabel('Y', fontsize=20, fontweight='bold')
    ax.set_zlabel('Z', fontsize=20, fontweight='bold')
    ax.set_title(lab, fontsize=20, fontweight='bold')
    ax.legend()
    plt.show()


"""
    if os.path.exists("F:\\Terella\\fortranintegrator\\a.exe"):
        os.remove("F:\\Terella\\fortranintegrator\\a.exe")
    if os.path.exists("F:\\Terella\\fortranintegrator\\subtest.8"):
        os.remove("F:\\Terella\\fortranintegrator\\subtest.8")

"""


# %%CALL SIMULATION
def main() -> None:
    """Execute the FORTRAN integrator workflow and create an animation."""

    df1 = fly_particle("1", "1", "pi/2.", "6", "6", "1", "1")

    process_data(df1, 'm', 'Proton w/ angle of pi/2')

    t = df1.iloc[0]
    x = df1.iloc[1]
    y = df1.iloc[2]
    z = df1.iloc[3]
    num_data_points = len(t)

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    def animate_func(num):
        ax.clear()  # Clears the figure to update the line, point, title, and axes

        # Updating Trajectory Line (num+1 due to Python indexing)
        ax.plot3D(x[:num + 1], y[:num + 1], z[:num + 1], c='red')

        # Updating Point Location
        ax.scatter(x[num], y[num], z[num], c='purple', marker='o', label="Electron")

        # Adding Constant Origin
        ax.plot3D(x[0], y[0], z[0], c='black', marker='o')

        # Setting Axes Limits
        ax.set_xlim3d([-8, 8])
        ax.set_ylim3d([-8, 8])
        ax.set_zlim3d([-3, 3])

        # Adding Figure Labels
        ax.set_title(
            'Trajectory \nTime = ' + str(np.round(t[num], decimals=7)) + ' sec',
            fontsize=20,
            fontweight='bold',
        )
        ax.set_xlabel('X', fontsize=20, fontweight='bold')
        ax.set_ylabel('Y', fontsize=20, fontweight='bold')
        ax.set_zlabel('Z', fontsize=20, fontweight='bold')

        ax.legend()

    line_ani = animation.FuncAnimation(fig, animate_func, interval=10, frames=num_data_points)
    plt.show()

    video_path = REPO_ROOT / "media" / "videos" / "imfanimation.mp4"
    video_path.parent.mkdir(parents=True, exist_ok=True)

    writergif = animation.FFMpegWriter(fps=num_data_points / 60)
    line_ani.save(video_path, writer=writergif)


if __name__ == "__main__":
    main()
