#!/usr/bin/env bash

set -euo pipefail

log() {
    printf "\033[1;34m==>\033[0m %s\n" "$*"
}

warn() {
    printf "\033[1;33m[WARN]\033[0m %s\n" "$*"
}

error() {
    printf "\033[1;31m[ERROR]\033[0m %s\n" "$*" >&2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
VENV_PATH="${PROJECT_ROOT}/.venv"

log "IMFPY macOS setup starting..."

if [[ "$(uname -s)" != "Darwin" ]]; then
    error "This setup script is intended to run on macOS."
    exit 1
fi

if ! xcode-select -p >/dev/null 2>&1; then
    warn "Xcode Command Line Tools are not installed. A dialog will prompt you to install them."
    xcode-select --install || true
    error "Command Line Tools must be installed before continuing. Re-run this script after installation completes."
    exit 1
fi

if command -v brew >/dev/null 2>&1; then
    BREW_BIN="$(command -v brew)"
    log "Homebrew detected: $(${BREW_BIN} --version | head -n1)"
else
    log "Homebrew not found. Installing (non-interactively)..."
    NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    if [[ -x "/opt/homebrew/bin/brew" ]]; then
        BREW_BIN="/opt/homebrew/bin/brew"
    elif [[ -x "/usr/local/bin/brew" ]]; then
        BREW_BIN="/usr/local/bin/brew"
    else
        error "Homebrew installation did not place brew on PATH."
        exit 1
    fi
fi

eval "$(${BREW_BIN} shellenv)"

BREW_PACKAGES=(python@3.11 gcc)
for package in "${BREW_PACKAGES[@]}"; do
    if ! ${BREW_BIN} list --versions "${package}" >/dev/null 2>&1; then
        log "Installing ${package} via Homebrew..."
        ${BREW_BIN} install "${package}"
    else
        log "${package} already installed via Homebrew."
    fi
done

PYTHON_BIN="$(${BREW_BIN} --prefix python@3.11)/bin/python3.11"
if [[ ! -x "${PYTHON_BIN}" ]]; then
    error "Unable to locate python3.11 from Homebrew at ${PYTHON_BIN}."
    exit 1
fi

log "Using Python interpreter: ${PYTHON_BIN}"

if [[ ! -d "${VENV_PATH}" ]]; then
    log "Creating virtual environment at ${VENV_PATH}"
    "${PYTHON_BIN}" -m venv "${VENV_PATH}"
else
    log "Virtual environment already exists at ${VENV_PATH}"
fi

# shellcheck disable=SC1091
source "${VENV_PATH}/bin/activate"

log "Upgrading pip, setuptools, and wheel"
pip install --upgrade pip setuptools wheel

log "Installing Python requirements"
pip install -r "${PROJECT_ROOT}/requirements.txt"

log "Running test suite smoke check"
pytest tests/test_simulation_runner.py

log "Attempting to compile the Fortran backend (dust_integrator)"
python - <<'PY'
from imfpy.fortran import build

try:
    module = build.build_module(force=True)
except build.FortranBuildError as exc:
    raise SystemExit(f"Fortran build failed: {exc}")
else:
    print(f"Fortran backend ready: {module.__file__}")
PY

log "Setup complete! Activate your environment with:"
log "    source ${VENV_PATH}/bin/activate"
log "Then launch the GUI with:"
log "    python -m imfpy.gui.main"
