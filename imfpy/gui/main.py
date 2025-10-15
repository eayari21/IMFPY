"""PyQt6 GUI for running dust trajectory simulations."""
from __future__ import annotations

import sys
import traceback
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from ..simulation import ParticleEnsemble, SimulationConfig, SimulationRunner
from ..simulation.results import SimulationResult
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
        self._build_ui()

    def _build_ui(self) -> None:
        central = QSplitter(Qt.Orientation.Horizontal)
        self.setCentralWidget(central)

        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        control_layout.setContentsMargins(8, 8, 8, 8)

        form = QFormLayout()
        self.particle_spin = QSpinBox()
        self.particle_spin.setRange(1, 10_000)
        self.particle_spin.setValue(4)
        form.addRow("Particles", self.particle_spin)

        self.step_spin = QSpinBox()
        self.step_spin.setRange(10, 200_000)
        self.step_spin.setValue(2000)
        form.addRow("Steps", self.step_spin)

        self.dt_spin = QDoubleSpinBox()
        self.dt_spin.setDecimals(6)
        self.dt_spin.setRange(1e-6, 1e9)
        self.dt_spin.setValue(60.0)
        form.addRow("Δt (s)", self.dt_spin)

        self.gm_spin = QDoubleSpinBox()
        self.gm_spin.setDecimals(6)
        self.gm_spin.setRange(1e-6, 1e30)
        self.gm_spin.setValue(1.32712440018e20)
        form.addRow("GM (m^3/s^2)", self.gm_spin)

        self.backend_combo = QComboBox()
        self.backend_combo.addItems(["fortran", "python", "gpu"])
        form.addRow("Backend", self.backend_combo)

        control_layout.addLayout(form)

        self.e_field_input = VectorInput("Electric field (V/m)")
        self.b_field_input = VectorInput("Magnetic field (T)")
        control_layout.addWidget(self.e_field_input)
        control_layout.addWidget(self.b_field_input)

        self.table = ParticleTableWidget()
        control_layout.addWidget(self.table)
        self.table.set_particle_count(self.particle_spin.value())
        self.table.populate_circular_orbit()

        button_layout = QHBoxLayout()
        self.run_button = QPushButton("Run Simulation")
        self.reset_button = QPushButton("Reset")
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.reset_button)
        control_layout.addLayout(button_layout)
        control_layout.addStretch(1)

        self.run_button.clicked.connect(self._start_simulation)
        self.reset_button.clicked.connect(self._reset)
        self.particle_spin.valueChanged.connect(self.table.set_particle_count)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.canvas = SimulationCanvas()
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        central.addWidget(control_widget)
        central.addWidget(plot_widget)
        central.setSizes([400, 800])

        status = QStatusBar()
        self.status_label = QLabel("Ready")
        status.addWidget(self.status_label)
        self.setStatusBar(status)

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

    def _simulation_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Simulation failed", message)
        self.status_label.setText("Error")
        self.run_button.setEnabled(True)


def run() -> None:
    app = QApplication(sys.argv)
    window = SimulationWindow()
    window.resize(1280, 720)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":  # pragma: no cover
    run()
