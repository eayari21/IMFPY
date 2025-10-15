"""In-application help content for the IMFPY GUI."""
from __future__ import annotations

import html
import platform
import sys
import webbrowser
from pathlib import Path
from textwrap import dedent
from typing import Final, Iterable

from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QMessageBox,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from ..fortran import FortranBuildError
from ..simulation import SimulationRunner
from ..simulation.backends import GPUNotAvailableError

REPO_ROOT_FALLBACK: Final[Path] = Path(__file__).resolve().parents[2]
DOCS_FOLDER: Final[Path] = REPO_ROOT_FALLBACK / "docs" / "user-guide"


HELP_HTML: Final[str] = dedent(
    """
    <h1>IMFPY Quick Start Guide</h1>
    <p>
        Welcome to the IMFPY Dust Trajectory Simulator! This guide walks through the most
        common actions directly from inside the application. For detailed documentation,
        open the full manual via the <strong>Help → Open Documentation Folder</strong> menu.
    </p>
    <h2>1. Interface overview</h2>
    <ul>
        <li><strong>Controls sidebar</strong> – configure particle counts, step counts, time step, and GM.</li>
        <li><strong>Field panels</strong> – set the electric (V/m) and magnetic (T) field vectors.</li>
        <li><strong>Particle table</strong> – edit positions, velocities, charges, masses, and β values.</li>
        <li><strong>Toolbar & canvas</strong> – interact with the 3D Matplotlib viewer once a run finishes.</li>
        <li><strong>Status bar</strong> – live feedback about builds, running simulations, and errors.</li>
    </ul>
    <h2>2. Run your first simulation</h2>
    <ol>
        <li>Review the defaults (four particles, circular orbit).</li>
        <li>Press <em>Run Simulation</em>. The status bar updates to show progress.</li>
        <li>When complete, rotate the 3D plot with the left mouse button and zoom with the scroll wheel.</li>
        <li>Use the legend to distinguish individual particles.</li>
    </ol>
    <h2>3. Customise scenarios</h2>
    <ul>
        <li>Change <em>Particles</em> to expand or shrink the table dynamically.</li>
        <li>Use <em>Populate Circular Orbit</em> to reset the ensemble, then tweak individual rows for custom trajectories.</li>
        <li>Tweak <em>Δt</em>, <em>Steps</em>, and <em>GM</em> for different central bodies.</li>
        <li>Experiment with electric/magnetic fields to reproduce solar wind conditions.</li>
        <li>Use <em>Reset</em> to restore defaults if anything goes wrong.</li>
    </ul>
    <h2>4. Choose backends</h2>
    <p>
        The backend dropdown selects between the optimised Fortran engine, the portable Python implementation,
        and the optional GPU mode. The list automatically hides unavailable backends; install CuPy to enable GPU support.
    </p>
    <h2>5. Export & share results</h2>
    <ul>
        <li>After a run, open <strong>File → Export Result…</strong> to save an <code>.npz</code> archive.</li>
        <li>The archive contains positions, velocities, metadata, and units.</li>
        <li>Load the file in Jupyter with <code>numpy.load</code> for post-processing.</li>
    </ul>
    <h2>6. Diagnostics & troubleshooting</h2>
    <ul>
        <li>Open <strong>Help → Diagnostics</strong> to see Python, backend, and CUDA availability.</li>
        <li>Check the terminal for detailed tracebacks when a simulation fails.</li>
        <li>Consult the troubleshooting matrix in the bundled documentation if builds fail or plots are empty.</li>
    </ul>
    <h2>7. Command-line workflows</h2>
    <p>
        Prefer scripting? Run <code>python -m imfpy.simulation.runner --help</code> to discover the CLI options.
        The <em>Tutorial library</em> inside the documentation folder includes automation examples and batch scripts.
    </p>
    <h2>8. Need more help?</h2>
    <p>
        The documentation lives in <code>docs/user-guide</code> next to this application. Share feedback or
        questions via the GitHub issue tracker – we love hearing about your projects!
    </p>
    """
).strip()


class HelpDialog(QDialog):
    """Dialog window that renders the bundled quick-start HTML."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("IMFPY Quick Start Guide")
        self.setModal(True)
        self.resize(820, 640)

        layout = QVBoxLayout(self)
        self.browser = QTextBrowser(self)
        self.browser.setOpenExternalLinks(True)
        self.browser.setHtml(HELP_HTML)
        layout.addWidget(self.browser)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def showEvent(self, event) -> None:  # type: ignore[override]
        """Ensure the cursor starts at the top when the dialog is shown."""
        super().showEvent(event)
        self.browser.moveCursor(self.browser.textCursor().Start)


def open_docs_folder(parent: QWidget | None = None) -> None:
    """Open the user-guide directory in the system file browser or warn if missing."""

    docs_path = DOCS_FOLDER
    if not docs_path.exists():
        QMessageBox.warning(
            parent,
            "Documentation not found",
            (
                "Could not locate the bundled documentation folder at\n"
                f"{docs_path}. If you installed IMFPY via a package manager, check the"
                " online README for the latest guide."
            ),
        )
        return

    webbrowser.open(docs_path.as_uri())


def _format_items(pairs: Iterable[tuple[str, str]]) -> str:
    items = [f"<li><strong>{html.escape(label)}:</strong> {html.escape(value)}" for label, value in pairs]
    return "\n".join(items)


def _backend_status_html() -> str:
    runner = SimulationRunner()
    entries: list[str] = []
    for name in ("fortran", "python", "gpu"):
        try:
            runner.registry.get(name)  # may raise if backend unavailable
        except GPUNotAvailableError as exc:
            entries.append(
                f"<li><strong>{name}</strong>: ❌ {html.escape(str(exc))}</li>"
            )
        except FortranBuildError as exc:
            entries.append(
                f"<li><strong>{name}</strong>: ❌ Fortran build failed – {html.escape(str(exc))}</li>"
            )
        except Exception as exc:  # pragma: no cover - defensive
            entries.append(
                f"<li><strong>{name}</strong>: ⚠️ {html.escape(type(exc).__name__ + ': ' + str(exc))}</li>"
            )
        else:
            entries.append(f"<li><strong>{name}</strong>: ✅ Available</li>")
    return "\n".join(entries)


def _build_diagnostics_html() -> str:
    python_line = sys.version.splitlines()[0]
    try:
        import numpy as np
    except Exception:  # pragma: no cover - optional
        numpy_version = "Not installed"
    else:
        numpy_version = np.__version__

    try:  # pragma: no cover - optional dependency
        import cupy as cp
    except Exception:
        cupy_version = "Not installed"
    else:
        cupy_version = cp.__version__

    env_pairs = [
        ("Python", python_line),
        ("Executable", sys.executable),
        ("Platform", platform.platform()),
        ("NumPy", numpy_version),
        ("CuPy", cupy_version),
    ]

    html_chunks = [
        "<h1>IMFPY Diagnostics</h1>",
        "<p>Environment summary and backend availability.</p>",
        "<h2>Environment</h2>",
        f"<ul>{_format_items(env_pairs)}</ul>",
        "<h2>Backends</h2>",
        f"<ul>{_backend_status_html()}</ul>",
        "<p>Use the command line self-test (<code>python -m imfpy.simulation.runner --self-test</code>) for deeper checks.</p>",
    ]
    return "\n".join(html_chunks)


class DiagnosticsDialog(QDialog):
    """Dialog that displays runtime environment and backend availability."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("IMFPY Diagnostics")
        self.resize(720, 520)

        layout = QVBoxLayout(self)
        self.browser = QTextBrowser(self)
        self.browser.setOpenExternalLinks(True)
        layout.addWidget(self.browser)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.refresh()

    def refresh(self) -> None:
        """Refresh the diagnostics content."""

        self.browser.setHtml(_build_diagnostics_html())
        self.browser.moveCursor(self.browser.textCursor().Start)


__all__ = ["DiagnosticsDialog", "HelpDialog", "open_docs_folder"]
