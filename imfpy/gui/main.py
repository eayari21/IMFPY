"""PyQt6 GUI for running dust trajectory simulations."""
from __future__ import annotations

import sys
import traceback

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal
from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QFileDialog,
    QPushButton,
    QSpinBox,
    QSplitter,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from imfpy.simulation import ParticleEnsemble, SimulationConfig, SimulationRunner
from imfpy.simulation.results import SimulationResult

from .help import DiagnosticsDialog, HelpDialog, open_docs_folder
from .widgets import ParticleTableWidget, VectorInput


class WorkerSignals(QObject):
    finished = pyqtSignal(SimulationResult)
    error = pyqtSignal(str)


class SimulationWorker(QRunnable):
    def __init__(self, config: SimulationConfig) -> None:
        super().__init__()
        self.config = config
        self.signals = WorkerSignals()

    def run(self) -> None:  # pragma: no cover - executed in thread
        try:
            runner = SimulationRunner()
            result = runner.run(self.config)
        except Exception as exc:  # noqa: BLE001
            tb = traceback.format_exc()
            self.signals.error.emit(f"{exc}\n{tb}")
        else:
            self.signals.finished.emit(result)


class SimulationCanvas(FigureCanvasQTAgg):
    def __init__(self, parent: QWidget | None = None) -> None:
        self.figure = Figure(figsize=(6, 5))
        super().__init__(self.figure)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.axes = self.figure.add_subplot(111, projection="3d")
        self.figure.tight_layout()

    def plot_result(self, result: SimulationResult) -> None:
        self.axes.clear()
        positions = result.positions()
        num_particles = result.num_particles
        for idx in range(num_particles):
            traj = positions[:, idx, :]
            self.axes.plot(traj[:, 0], traj[:, 1], traj[:, 2], label=f"Particle {idx+1}")
            self.axes.scatter(traj[0, 0], traj[0, 1], traj[0, 2], marker="o")
        self.axes.scatter([0], [0], [0], color="gold", marker="*", s=200, label="Star")
        self.axes.set_xlabel("x (m)")
        self.axes.set_ylabel("y (m)")
        self.axes.set_zlabel("z (m)")
        self.axes.legend(loc="upper right")
        self.axes.set_title(f"Trajectories ({result.backend})")
        self.figure.canvas.draw_idle()


class SimulationWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("IMFPY Dust Trajectory Simulator")
        self.thread_pool = QThreadPool.globalInstance()
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)
        self._help_dialog: HelpDialog | None = None
        self._diagnostics_dialog: DiagnosticsDialog | None = None
        self._last_result: SimulationResult | None = None
        status = QStatusBar()
        self.status_label = QLabel("Welcome")
        status.addWidget(self.status_label)
        self.setStatusBar(status)
        self._build_menus()
        self._build_pages()

    def _build_pages(self) -> None:
        self.launch_widget = self._create_launch_screen()
        self.main_widget = self._create_main_screen()
        self.stack.addWidget(self.launch_widget)
        self.stack.addWidget(self.main_widget)
        self.stack.setCurrentWidget(self.launch_widget)

    def _create_launch_screen(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(48, 48, 48, 48)
        layout.setSpacing(24)
        layout.addStretch(1)

        title = QLabel("Welcome to the IMFPY Dust Trajectory Simulator")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setWordWrap(True)
        title.setStyleSheet("font-size: 24pt; font-weight: 600;")
        layout.addWidget(title)

        subtitle = QLabel(
            "Craft ensembles of particles, run the integrator, and inspect the trajectories."
        )
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        button_row = QHBoxLayout()
        button_row.setSpacing(16)
        button_row.addStretch(1)

        start_button = QPushButton("Start a Simulation")
        start_button.setMinimumWidth(220)
        start_button.setDefault(True)
        start_button.clicked.connect(self._show_main_view)
        button_row.addWidget(start_button)

        quickstart_button = QPushButton("Quick Start Guide")
        quickstart_button.clicked.connect(self._show_help)
        button_row.addWidget(quickstart_button)

        button_row.addStretch(1)
        layout.addLayout(button_row)
        layout.addStretch(2)

        return widget

    def _create_main_screen(self) -> QWidget:
        main = QWidget()
        main_layout = QVBoxLayout(main)
        main_layout.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)

        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        control_layout.setContentsMargins(12, 12, 12, 12)
        control_layout.setSpacing(12)

        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        self.particle_spin = QSpinBox()
        self.particle_spin.setRange(1, 10_000)
        self.particle_spin.setValue(4)
        self.particle_spin.setToolTip("Number of dust particles to integrate.")
        form.addRow("Particles", self.particle_spin)

        self.step_spin = QSpinBox()
        self.step_spin.setRange(10, 200_000)
        self.step_spin.setValue(2000)
        self.step_spin.setToolTip("Total integration steps to perform.")
        form.addRow("Steps", self.step_spin)

        self.dt_spin = QDoubleSpinBox()
        self.dt_spin.setDecimals(6)
        self.dt_spin.setRange(1e-6, 1e9)
        self.dt_spin.setValue(60.0)
        self.dt_spin.setToolTip("Time step (seconds) between integration points.")
        form.addRow("Δt (s)", self.dt_spin)

        self.gm_spin = QDoubleSpinBox()
        self.gm_spin.setDecimals(6)
        self.gm_spin.setRange(1e-6, 1e30)
        self.gm_spin.setValue(1.32712440018e20)
        self.gm_spin.setToolTip("Gravitational parameter GM of the central body (m^3/s^2).")
        form.addRow("GM (m^3/s^2)", self.gm_spin)

        self.backend_combo = QComboBox()
        self.backend_combo.addItems(["fortran", "python", "gpu"])
        self.backend_combo.setToolTip(
            "Select the simulation backend (availability depends on installed dependencies)."
        )
        form.addRow("Backend", self.backend_combo)

        control_layout.addLayout(form)

        self.e_field_input = VectorInput("Electric field (V/m)")
        self.b_field_input = VectorInput("Magnetic field (T)")
        self.e_field_input.setToolTip("Uniform electric field applied to all particles (V/m).")
        self.b_field_input.setToolTip("Uniform magnetic field applied to all particles (tesla).")
        control_layout.addWidget(self.e_field_input)
        control_layout.addWidget(self.b_field_input)

        self.table = ParticleTableWidget()
        self.table.setToolTip("Edit particle initial conditions, mass, charge, and radiation-pressure β.")
        control_layout.addWidget(self.table)
        self.table.set_particle_count(self.particle_spin.value())
        self.table.populate_circular_orbit()

        button_layout = QHBoxLayout()
        button_layout.setSpacing(12)
        self.run_button = QPushButton("Run Simulation")
        self.reset_button = QPushButton("Reset")
        self.run_button.setToolTip("Start the simulation with the current settings.")
        self.reset_button.setToolTip("Restore default parameters and circular orbit ensemble.")
        self.run_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.reset_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.reset_button)
        control_layout.addLayout(button_layout)
        control_layout.addStretch(1)

        self.run_button.clicked.connect(self._start_simulation)
        self.reset_button.clicked.connect(self._reset)
        self.particle_spin.valueChanged.connect(self.table.set_particle_count)

        control_container = QScrollArea()
        control_container.setWidgetResizable(True)
        control_container.setWidget(control_widget)
        control_container.setMinimumWidth(360)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(12, 12, 12, 12)
        plot_layout.setSpacing(12)
        self.canvas = SimulationCanvas()
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        splitter.addWidget(control_container)
        splitter.addWidget(plot_widget)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)

        main_layout.addWidget(splitter)
        return main

    def _show_main_view(self) -> None:
        self.stack.setCurrentWidget(self.main_widget)
        self.status_label.setText("Ready")
        self.run_button.setFocus()

    def _build_menus(self) -> None:
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("&File")
        self.export_action = QAction("Export Result…", self)
        self.export_action.setEnabled(False)
        self.export_action.triggered.connect(self._export_result)
        file_menu.addAction(self.export_action)

        help_menu = menu_bar.addMenu("&Help")
        self.quick_start_action = QAction("Quick Start Guide", self)
        self.quick_start_action.setShortcut(QKeySequence.StandardKey.HelpContents)
        self.quick_start_action.triggered.connect(self._show_help)
        help_menu.addAction(self.quick_start_action)

        diagnostics_action = QAction("Diagnostics", self)
        diagnostics_action.triggered.connect(self._show_diagnostics)
        help_menu.addAction(diagnostics_action)

        docs_action = QAction("Open Documentation Folder", self)
        docs_action.triggered.connect(lambda: open_docs_folder(self))
        help_menu.addAction(docs_action)

    def _reset(self) -> None:
        self.table.populate_circular_orbit()
        self.status_label.setText("Defaults restored")

    def _start_simulation(self) -> None:
        try:
            config = self._build_config()
        except ValueError as exc:
            QMessageBox.critical(self, "Invalid configuration", str(exc))
            return
        self.status_label.setText("Running simulation...")
        self.run_button.setEnabled(False)
        worker = SimulationWorker(config)
        worker.signals.finished.connect(self._simulation_finished)
        worker.signals.error.connect(self._simulation_failed)
        self.thread_pool.start(worker)

    def _build_config(self) -> SimulationConfig:
        table_data = self.table.data()
        ensemble = ParticleEnsemble(
            positions=table_data.positions,
            velocities=table_data.velocities,
            charges=table_data.charges,
            masses=table_data.masses,
            betas=table_data.betas,
        )
        config = SimulationConfig(
            particles=ensemble,
            num_steps=self.step_spin.value(),
            time_step=self.dt_spin.value(),
            gm_sun=self.gm_spin.value(),
            electric_field=self.e_field_input.vector(),
            magnetic_field=self.b_field_input.vector(),
            backend=self.backend_combo.currentText(),
        )
        return config

    def _simulation_finished(self, result: SimulationResult) -> None:
        self.canvas.plot_result(result)
        self.status_label.setText("Simulation finished")
        self.run_button.setEnabled(True)
        self._last_result = result
        self.export_action.setEnabled(True)

    def _simulation_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Simulation failed", message)
        self.status_label.setText("Error")
        self.run_button.setEnabled(True)

    def _show_help(self) -> None:
        if self._help_dialog is None:
            self._help_dialog = HelpDialog(self)
        self._help_dialog.show()
        self._help_dialog.raise_()
        self._help_dialog.activateWindow()

    def _show_diagnostics(self) -> None:
        if self._diagnostics_dialog is None:
            self._diagnostics_dialog = DiagnosticsDialog(self)
        self._diagnostics_dialog.refresh()
        self._diagnostics_dialog.show()
        self._diagnostics_dialog.raise_()
        self._diagnostics_dialog.activateWindow()

    def _export_result(self) -> None:
        if self._last_result is None:
            QMessageBox.information(self, "No result", "Run a simulation before exporting.")
            return

        default_name = "trajectory_run.npz"
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Simulation Result",
            default_name,
            "NumPy archive (*.npz)",
        )
        if not path:
            return

        np.savez(
            path,
            times=self._last_result.times,
            state=self._last_result.state,
            backend=self._last_result.backend,
        )
        self.status_label.setText(f"Exported result to {path}")


def run() -> None:
    app = QApplication(sys.argv)
    window = SimulationWindow()
    screen = app.primaryScreen()
    if screen is not None:
        available = screen.availableGeometry()
        width = int(available.width() * 0.8)
        height = int(available.height() * 0.8)
        window.resize(width, height)
    else:
        window.resize(1280, 720)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":  # pragma: no cover
    run()
