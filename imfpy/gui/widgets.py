"""Reusable widgets for the IMFPY GUI."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
from PyQt6.QtCore import QRegularExpression
from PyQt6.QtGui import QRegularExpressionValidator
from PyQt6.QtWidgets import (
    QGridLayout,
    QGroupBox,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)


@dataclass
class ParticleTableData:
    positions: np.ndarray
    velocities: np.ndarray
    masses: np.ndarray
    charges: np.ndarray
    betas: np.ndarray


class ParticleTableWidget(QGroupBox):
    """Table widget allowing per-particle configuration."""

    headers = ["x", "y", "z", "vx", "vy", "vz", "mass", "charge", "beta"]

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__("Particle Parameters", parent)
        self.table = QTableWidget(0, len(self.headers))
        self.table.setHorizontalHeaderLabels(self.headers)
        header = self.table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table)

        button_row = QWidget()
        button_layout = QGridLayout(button_row)
        self.populate_button = QPushButton("Populate Circular Orbit")
        self.clear_button = QPushButton("Clear")
        button_layout.addWidget(QLabel("Presets:"), 0, 0)
        button_layout.addWidget(self.populate_button, 0, 1)
        button_layout.addWidget(self.clear_button, 0, 2)
        layout.addWidget(button_row)

        self.populate_button.clicked.connect(self.populate_circular_orbit)
        self.clear_button.clicked.connect(self.clear)

    def set_particle_count(self, count: int) -> None:
        current = self.table.rowCount()
        if current < count:
            for _ in range(count - current):
                self.table.insertRow(self.table.rowCount())
        elif current > count:
            for _ in range(current - count):
                self.table.removeRow(self.table.rowCount() - 1)

        validator = QRegularExpressionValidator(QRegularExpression(r"^[-+]?\d*(?:\.\d*)?(?:[eE][-+]?\d+)?$"))
        for row in range(self.table.rowCount()):
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                if item is None:
                    item = QTableWidgetItem("0.0")
                    self.table.setItem(row, col, item)
                line_edit = self.table.cellWidget(row, col)
                if line_edit is None:
                    line_edit = QLineEdit(self.table)
                    line_edit.setValidator(validator)
                    line_edit.setText(item.text())
                    self.table.setCellWidget(row, col, line_edit)
                else:
                    line_edit.setValidator(validator)
                    if not line_edit.text():
                        line_edit.setText("0.0")

    def populate_circular_orbit(self) -> None:
        rows = self.table.rowCount()
        if rows == 0:
            return
        radius = 1.5e11
        speed = 26_000.0
        mass = 1e-16
        charge = 1.6e-19
        beta = 0.1
        for i in range(rows):
            angle = 2 * np.pi * i / max(rows, 1)
            pos = radius * np.array([np.cos(angle), np.sin(angle), 0.0])
            vel = speed * np.array([-np.sin(angle), np.cos(angle), 0.0])
            values = [
                *pos.tolist(),
                *vel.tolist(),
                mass,
                charge,
                beta,
            ]
            for col, value in enumerate(values):
                widget = self.table.cellWidget(i, col)
                if widget is not None:
                    widget.setText(f"{value:.6e}")

    def clear(self) -> None:
        for row in range(self.table.rowCount()):
            for col in range(self.table.columnCount()):
                widget = self.table.cellWidget(row, col)
                if widget is not None:
                    widget.setText("0.0")

    def data(self) -> ParticleTableData:
        rows = self.table.rowCount()
        values = np.zeros((rows, len(self.headers)), dtype=float)
        for row in range(rows):
            for col in range(len(self.headers)):
                widget = self.table.cellWidget(row, col)
                text = widget.text() if widget is not None else "0.0"
                try:
                    values[row, col] = float(text)
                except ValueError:
                    raise ValueError(f"Invalid numeric value at row {row + 1}, column {self.headers[col]}")
        positions = values[:, 0:3]
        velocities = values[:, 3:6]
        masses = values[:, 6]
        charges = values[:, 7]
        betas = values[:, 8]
        return ParticleTableData(positions, velocities, masses, charges, betas)


class VectorInput(QWidget):
    """Simple widget for editing 3-component vectors."""

    def __init__(self, label: str, default: List[float] | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QGridLayout(self)
        layout.addWidget(QLabel(label), 0, 0)
        self.fields: List[QLineEdit] = []
        values = default or [0.0, 0.0, 0.0]
        validator = QRegularExpressionValidator(QRegularExpression(r"^[-+]?\d*(?:\.\d*)?(?:[eE][-+]?\d+)?$"))
        for i, axis in enumerate("xyz"):
            layout.addWidget(QLabel(axis), 0, i + 1)
            field = QLineEdit(self)
            field.setValidator(validator)
            field.setText(f"{values[i]:.6g}")
            layout.addWidget(field, 1, i + 1)
            self.fields.append(field)
        layout.setColumnStretch(0, 2)
        for col in range(1, 4):
            layout.setColumnStretch(col, 1)

    def vector(self) -> np.ndarray:
        values = []
        for field in self.fields:
            text = field.text() or "0.0"
            values.append(float(text))
        return np.array(values, dtype=float)


__all__ = ["ParticleTableWidget", "ParticleTableData", "VectorInput"]
