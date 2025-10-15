#!/usr/bin/env bash
set -euo pipefail

blue(){ printf "\033[1;34m==>\033[0m %s\n" "$*"; }
warn(){ printf "\033[1;33m[WARN]\033[0m %s\n" "$*"; }
err(){  printf "\033[1;31m[ERROR]\033[0m %s\n" "$*" >&2; }

# --- sanity: repo root ---
[[ -f requirements.txt ]] || { err "Run from repo root (where requirements.txt lives)."; exit 1; }

PROJECT_ROOT="$PWD"
VENV_PATH="${PROJECT_ROOT}/.venv"

blue "Starting IMFPY bootstrap in: ${PROJECT_ROOT}"

# --- CLT ---
if ! xcode-select -p >/dev/null 2>&1; then
  warn "Command Line Tools missing; opening installer..."
  xcode-select --install || true
  err "Install CLT, then re-run this script."
  exit 1
fi

# --- Homebrew ---
if command -v brew >/dev/null 2>&1; then
  BREW_BIN="$(command -v brew)"
  blue "Homebrew: $("$BREW_BIN" --version | head -n1)"
else
  blue "Installing Homebrew (non-interactive)..."
  NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  if [[ -x "/opt/homebrew/bin/brew" ]]; then BREW_BIN="/opt/homebrew/bin/brew"
  elif [[ -x "/usr/local/bin/brew" ]]; then BREW_BIN="/usr/local/bin/brew"
  else err "brew not on PATH after install"; exit 1; fi
fi
eval "$("$BREW_BIN" shellenv)" || true
grep -q 'brew shellenv' "${HOME}/.zprofile" 2>/dev/null || echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> "${HOME}/.zprofile"

# --- Python + gcc ---
for pkg in python@3.11 gcc; do
  "$BREW_BIN" list --versions "$pkg" >/dev/null 2>&1 || { blue "Installing $pkg ..."; "$BREW_BIN" install "$pkg"; }
  blue "$pkg OK"
done
PY311="$("$BREW_BIN" --prefix python@3.11)/bin/python3.11"
[[ -x "$PY311" ]] || { err "python3.11 not found"; exit 1; }

# Ensure a "gfortran" shim (Homebrew often installs gfortran-14, etc.)
if ! command -v gfortran >/dev/null 2>&1; then
  CAND="$(ls /opt/homebrew/bin/gfortran-* 2>/dev/null | head -n1 || true)"
  if [[ -n "${CAND}" ]]; then
    blue "Linking ${CAND} -> /opt/homebrew/bin/gfortran (may prompt for sudo)"
    sudo ln -sf "${CAND}" /opt/homebrew/bin/gfortran
  fi
fi
command -v gfortran >/dev/null 2>&1 || { err "gfortran still not found; run: brew install gcc"; exit 1; }
export FC="$(command -v gfortran)"; export F77="$FC"; export F90="$FC"
blue "Using gfortran: $FC"

# --- venv ---
if [[ ! -d "$VENV_PATH" ]]; then
  blue "Creating venv: $VENV_PATH"
  "$PY311" -m venv "$VENV_PATH"
else
  blue "venv exists: $VENV_PATH"
fi
# shellcheck disable=SC1090
source "$VENV_PATH/bin/activate"
blue "Python: $(python -c 'import sys; print(sys.executable)')"

python -m pip install --upgrade pip setuptools wheel

# --- pyproject (leave existing alone) ---
USE_SRC_LAYOUT=0
if [[ -f "src/imfpy/__init__.py" ]]; then USE_SRC_LAYOUT=1
elif [[ -f "imfpy/__init__.py" ]]; then USE_SRC_LAYOUT=0
else err "imfpy package not found (looked for src/imfpy/__init__.py or imfpy/__init__.py)"; exit 1
fi
if [[ ! -f "pyproject.toml" ]]; then
  blue "Creating pyproject.toml"
  if [[ $USE_SRC_LAYOUT -eq 1 ]]; then
    cat > pyproject.toml <<'EOF'
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "imfpy"
version = "0.0.0"
description = "Field & particle simulation tools"
readme = "README.md"
requires-python = ">=3.11"

[tool.setuptools]
package-dir = {"" = "src"}
packages = ["imfpy"]
EOF
  else
    cat > pyproject.toml <<'EOF'
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "imfpy"
version = "0.0.0"
description = "Field & particle simulation tools"
readme = "README.md"
requires-python = ">=3.11"

[tool.setuptools]
packages = ["imfpy"]
EOF
  fi
else
  blue "pyproject.toml exists; keeping it."
fi

# --- setup.cfg (optional, harmless if exists) ---
if [[ ! -f "setup.cfg" ]]; then
  blue "Creating setup.cfg"
  cat > setup.cfg <<'EOF'
[metadata]
license = MIT
author = Ethan Ayari

[options]
include_package_data = True
install_requires =
    numpy
    scipy
    matplotlib
    PyQt6
    numba
EOF
else
  blue "setup.cfg exists; keeping it."
fi

# --- deps ---
blue "Installing requirements.txt"
python -m pip install -r requirements.txt

# --- Patch imfpy/fortran/build.py so imports are clean (no f2py at import) ---
if [[ -f "imfpy/fortran/build.py" ]]; then
  blue "Patching imfpy/fortran/build.py"
  python - <<'PY'
from pathlib import Path
import re
p = Path("imfpy/fortran/build.py")
s = p.read_text(encoding="utf-8")

# A) Ensure `from __future__ import annotations` is first (after docstring)
FUT = "from __future__ import annotations"
lines = [ln for ln in s.splitlines() if ln.strip() != FUT]
def doc_end(ls):
    i=0
    while i<len(ls) and (ls[i].strip()=="" or ls[i].lstrip().startswith("#")): i+=1
    if i<len(ls) and ls[i].lstrip().startswith(('"""',"'''")):
        q=ls[i].lstrip()[:3]
        if ls[i].strip().endswith(q) and len(ls[i].strip())>3: return i+1
        j=i+1
        while j<len(ls):
            if ls[j].strip().endswith(q): return j+1
            j+=1
    return 0
end = doc_end(lines)
body = lines[end:]
body = [ln for ln in body if ln.strip() != "import sys"]
new = lines[:end] + [FUT, "import sys"] + body
s = "\n".join(new) + "\n"

# B) Remove dangerous f2py imports
s = re.sub(r"^\s*from\s+numpy\.f2py\.__main__\s+import\s+main\s+as\s+f2py_main\s*$","",s,flags=re.MULTILINE)
s = re.sub(r"^\s*import\s+numpy\.f2py\.__main__(\s+as\s+f2py_main)?\s*$","",s,flags=re.MULTILINE)

# C) Ensure subprocess import
if "import subprocess" not in s:
    s = s.replace("import sys", "import sys\nimport subprocess", 1)

# D) Safe wrapper that only runs when explicitly called
if "def _f2py_run(" not in s:
    s += (
        "\n\ndef _f2py_run(args):\n"
        "    cmd = [sys.executable, '-m', 'numpy.f2py', *args]\n"
        "    subprocess.run(cmd, check=True)\n"
    )

# E) Replace any f2py_main(...) uses
s = re.sub(r"\bf2py_main\s*\(", "_f2py_run(", s)

# F) Remove top-level side-effect calls (autobuilds)
def comment_top_calls(text, names=("build_module","_f2py_run")):
    out=[]
    for ln in text.splitlines():
        if any(ln.startswith(n + "(") for n in names):
            out.append("# DISABLED_AUTOBUILD: " + ln)
        else:
            out.append(ln)
    return "\n".join(out) + "\n"
s = comment_top_calls(s)

# G) Neutralize any __main__ block
s = re.sub(r"if\s+__name__\s*==\s*['\"]__main__['\"]\s*:\s*[\s\S]*$", "if __name__ == '__main__':\n    pass\n", s)

p.write_text(s, encoding="utf-8")
print("Patched", p)
PY
else
  warn "imfpy/fortran/build.py not found; skipping patch."
fi

# --- install editable ---
blue "Installing imfpy (editable)"
python -m pip install -e . >/dev/null

# --- smoke import ---
blue "Smoke import"
python - <<'PY'
import imfpy
print("imfpy:", imfpy.__file__)
import imfpy.fortran.build as b
print("build.py imports cleanly")
PY

# --- tests (no -q, to avoid argv shenanigans if any remain) ---
if [[ -f "tests/test_simulation_runner.py" ]]; then
  blue "Running pytest tests/test_simulation_runner.py"
  python -m pytest tests/test_simulation_runner.py
else
  warn "tests/test_simulation_runner.py not found; skipping tests."
fi

blue "Done!"
echo
echo "Next:"
echo "  source \"$VENV_PATH/bin/activate\""
echo "  python -m imfpy.gui.main"
