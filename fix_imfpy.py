from pathlib import Path
import re
import subprocess, sys

ROOT = Path(__file__).resolve().parent
PKG = (ROOT / "src" / "imfpy") if (ROOT / "src" / "imfpy").exists() else (ROOT / "imfpy")

def patch_numba_norm():
    """Replace np.linalg.norm(positions, axis=1) with numba-friendly form."""
    target = PKG / "simulation" / "backends" / "python_backend.py"
    txt = target.read_text(encoding="utf-8")

    # Replace r = np.linalg.norm(positions, axis=1) with sqrt((positions*positions).sum(axis=1))
    pattern = r"np\.linalg\.norm\s*\(\s*positions\s*,\s*axis\s*=\s*1\s*\)"
    repl = "np.sqrt((positions * positions).sum(axis=1))"
    new_txt, n = re.subn(pattern, repl, txt)
    if n:
        target.write_text(new_txt, encoding="utf-8")
        print(f"[patch] Replaced {n} occurrence(s) of np.linalg.norm(..., axis=1) in {target}")
    else:
        print(f"[patch] No np.linalg.norm(..., axis=1) occurrence found in {target} (already patched?)")

def reinstall_editable():
    print("[pip] reinstall editable imfpy …")
    subprocess.run([sys.executable, "-m", "pip", "install", "-e", "."], check=True)

def run_tests():
    print("[pytest] running tests/test_simulation_runner.py …")
    subprocess.run([sys.executable, "-m", "pytest", "tests/test_simulation_runner.py"], check=False)

if __name__ == "__main__":
    patch_numba_norm()
    reinstall_editable()
    run_tests()
    print("\n[ok] Patch complete. In the GUI, pick backend = 'python' to run without Fortran.\n")


