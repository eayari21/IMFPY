"""Simulation configuration models."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Literal, Sequence

import numpy as np

BackendName = Literal["fortran", "python", "gpu"]


def _to_array(values: Sequence[float], *, length: int, dtype=float) -> np.ndarray:
    array = np.asarray(values, dtype=dtype)
    if array.shape != (length,):
        raise ValueError(f"Expected sequence of length {length}, got shape {array.shape}")
    return array


@dataclass
class ParticleEnsemble:
    """Container describing the initial state for a set of particles."""

    positions: np.ndarray
    velocities: np.ndarray
    charges: np.ndarray
    masses: np.ndarray
    betas: np.ndarray

    def __post_init__(self) -> None:
        for name, value in (
            ("positions", self.positions),
            ("velocities", self.velocities),
            ("charges", self.charges),
            ("masses", self.masses),
            ("betas", self.betas),
        ):
            if not isinstance(value, np.ndarray):
                raise TypeError(f"{name} must be a numpy.ndarray")
        if self.positions.shape != self.velocities.shape:
            raise ValueError("positions and velocities must share the same shape")
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError("positions must be shaped (n_particles, 3)")
        n_particles = self.positions.shape[0]
        for name, vector in (
            ("velocities", self.velocities),
            ("charges", self.charges),
            ("masses", self.masses),
            ("betas", self.betas),
        ):
            if vector.shape[0] != n_particles:
                raise ValueError(f"{name} length {vector.shape[0]} does not match {n_particles}")
        if np.any(self.masses <= 0):
            raise ValueError("Mass values must be strictly positive")

    @property
    def count(self) -> int:
        return self.positions.shape[0]

    @classmethod
    def from_sequences(
        cls,
        positions: Iterable[Sequence[float]],
        velocities: Iterable[Sequence[float]],
        charges: Sequence[float],
        masses: Sequence[float],
        betas: Sequence[float],
    ) -> "ParticleEnsemble":
        pos_array = np.array(list(positions), dtype=float)
        vel_array = np.array(list(velocities), dtype=float)
        charge_array = np.asarray(charges, dtype=float)
        mass_array = np.asarray(masses, dtype=float)
        beta_array = np.asarray(betas, dtype=float)
        return cls(pos_array, vel_array, charge_array, mass_array, beta_array)

    def as_initial_state(self) -> np.ndarray:
        """Return concatenated position/velocity matrix shaped (n_particles, 6)."""
        return np.hstack((self.positions, self.velocities))


@dataclass
class SimulationConfig:
    """Full configuration for a dust ensemble integration."""

    particles: ParticleEnsemble
    num_steps: int
    time_step: float
    gm_sun: float = 1.32712440018e20  # m^3 s^-2
    electric_field: np.ndarray = field(default_factory=lambda: np.zeros(3))
    magnetic_field: np.ndarray = field(default_factory=lambda: np.zeros(3))
    backend: BackendName = "fortran"
    gpu_device: int | None = None

    def __post_init__(self) -> None:
        if self.num_steps <= 1:
            raise ValueError("num_steps must be > 1")
        if self.time_step <= 0:
            raise ValueError("time_step must be positive")
        self.electric_field = _to_array(self.electric_field, length=3)
        self.magnetic_field = _to_array(self.magnetic_field, length=3)
        if self.backend not in ("fortran", "python", "gpu"):
            raise ValueError(f"Unsupported backend: {self.backend}")

    @property
    def time_span(self) -> float:
        return self.time_step * float(self.num_steps - 1)

    @property
    def times(self) -> np.ndarray:
        return np.linspace(0.0, self.time_span, self.num_steps, dtype=float)

    def q_over_m(self) -> np.ndarray:
        return self.particles.charges / self.particles.masses


__all__ = ["BackendName", "ParticleEnsemble", "SimulationConfig"]
