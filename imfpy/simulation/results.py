"""Result containers for simulations."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SimulationResult:
    times: np.ndarray
    state: np.ndarray  # shape (num_steps, num_particles, 6)
    backend: str

    def positions(self) -> np.ndarray:
        return self.state[:, :, :3]

    def velocities(self) -> np.ndarray:
        return self.state[:, :, 3:]

    @property
    def num_particles(self) -> int:
        return self.state.shape[1]

    @property
    def num_steps(self) -> int:
        return self.state.shape[0]


__all__ = ["SimulationResult"]
