"""Tests for Fortran build helpers."""
from __future__ import annotations

import textwrap

import pytest

pytest.importorskip("numpy")

from imfpy.fortran import build


def test_build_module_raises_friendly_error_on_f2py_failure(monkeypatch, tmp_path):
    source = tmp_path / "dust_integrator.f90"
    source.write_text(textwrap.dedent(
        """
        subroutine dummy()
        return
        end
        """
    ))

    build_dir = tmp_path / "build"

    monkeypatch.setattr(build, "_SOURCE_FILE", source)
    monkeypatch.setattr(build, "_BUILD_DIR", build_dir)

    def fake_f2py_main(args):
        print("usage: f2py ...")
        raise SystemExit(2)

    monkeypatch.setattr(build, "f2py_main", fake_f2py_main)

    with pytest.raises(build.FortranBuildError) as excinfo:
        build.build_module(force=True)

    assert "F2PY failed" in str(excinfo.value)
    assert "usage: f2py" in str(excinfo.value)
