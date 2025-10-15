"""Utilities to build and load the Fortran dust integrator module."""
from __future__ import annotations

import contextlib
import importlib.util
import io
import os
from pathlib import Path
from types import ModuleType
from typing import Optional

from numpy.f2py.__main__ import main as f2py_main

MODULE_NAME = "dust_integrator"
_BUILD_DIR = Path(__file__).with_name("_build")
_SOURCE_FILE = Path(__file__).with_name("dust_integrator.f90")


class FortranBuildError(RuntimeError):
    """Raised when the Fortran backend cannot be compiled."""


def _existing_module_path() -> Optional[Path]:
    if not _BUILD_DIR.exists():
        return None
    candidates = sorted(_BUILD_DIR.glob(f"{MODULE_NAME}*.so"))
    if not candidates:
        candidates = sorted(_BUILD_DIR.glob(f"{MODULE_NAME}*.pyd"))
    return candidates[-1] if candidates else None


def _import_module_from_path(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(MODULE_NAME, path)
    if spec is None or spec.loader is None:
        raise FortranBuildError(f"Cannot load module spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_module(force: bool = False) -> ModuleType:
    """Build and import the Fortran module."""
    _BUILD_DIR.mkdir(exist_ok=True)

    if not force:
        existing = _existing_module_path()
        if existing is not None:
            return _import_module_from_path(existing)

    if not _SOURCE_FILE.exists():
        raise FortranBuildError(f"Missing source file: {_SOURCE_FILE}")

    args = [
        "-c",
        str(_SOURCE_FILE),
        "-m",
        MODULE_NAME,
        "--quiet",
        "--backend",
        "distutils",
        "--opt=-O3",
        "--opt=-fopenmp",
        "--f90flags=-fopenmp",
        "--link-args=-fopenmp",
    ]

    env = os.environ.copy()
    env.setdefault("NPY_DISTUTILS_APPEND_FLAGS", "1")

    cwd = os.getcwd()
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    try:
        os.chdir(_BUILD_DIR)
        try:
            with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(
                stderr_buffer
            ):
                f2py_main(args)
        except SystemExit as exc:  # pragma: no cover - depends on numpy behaviour
            code = exc.code if isinstance(exc.code, int) else 1
            if code not in (0, None):
                output = stdout_buffer.getvalue().strip()
                errors = stderr_buffer.getvalue().strip()
                combined = "\n".join(part for part in (output, errors) if part)
                message = "F2PY failed to compile the Fortran backend"
                if combined:
                    message = f"{message}:\n{combined}"
                raise FortranBuildError(message) from exc
        finally:
            stdout_buffer.close()
            stderr_buffer.close()
    finally:
        os.chdir(cwd)

    module_path = _existing_module_path()
    if module_path is None:
        raise FortranBuildError("Compilation succeeded but module not found.")

    return _import_module_from_path(module_path)


def load_module() -> ModuleType:
    """Load the Fortran module, building it if necessary."""
    try:
        return build_module(force=False)
    except FortranBuildError as exc:  # pragma: no cover - informative message
        raise FortranBuildError(
            "Unable to compile the Fortran dust integrator. "
            "Ensure that gfortran (or another supported Fortran compiler) "
            "is installed and available on PATH."
        ) from exc


__all__ = ["MODULE_NAME", "build_module", "load_module", "FortranBuildError"]
