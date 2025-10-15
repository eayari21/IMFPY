"""NumPy/Numba reference backend."""
from __future__ import annotations

import logging
from typing import Callable

import numpy as np

from ..config import SimulationConfig
from ..results import SimulationResult

LOGGER = logging.getLogger(__name__)

try:  # pragma: no cover - optional dependency
    from numba import njit
except Exception:  # pragma: no cover - optional dependency
    njit = None


def _derivatives_numpy(
    state: np.ndarray,
    gm_sun: float,
    q_over_m: np.ndarray,
    beta: np.ndarray,
    electric_field: np.ndarray,
    magnetic_field: np.ndarray,
) -> np.ndarray:
    positions = state[:, :3]
    velocities = state[:, 3:]
    r = np.linalg.norm(positions, axis=1)
    r3 = np.where(r > 1e-12, r**3, np.inf)
    grav_coeff = -gm_sun * (1.0 - beta) / r3
    grav = positions * grav_coeff[:, None]
    cross = np.cross(velocities, magnetic_field)
    lorentz = q_over_m[:, None] * (electric_field + cross)
    acc = grav + lorentz
    return np.hstack((velocities, acc))


if njit is not None:  # pragma: no cover - compiled path
    _derivatives_numba = njit(_derivatives_numpy, fastmath=True, cache=True)
else:
    _derivatives_numba = None


def _rk4_step(state: np.ndarray, dt: float, deriv: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    k1 = deriv(state)
    k2 = deriv(state + 0.5 * dt * k1)
    k3 = deriv(state + 0.5 * dt * k2)
    k4 = deriv(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


class PythonBackend:
    """Pure Python/NumPy backend for validation and GPU fallbacks."""

    def run(self, config: SimulationConfig) -> SimulationResult:
        particles = config.particles
        dt = config.time_step
        num_steps = config.num_steps
        state = np.zeros((num_steps, particles.count, 6), dtype=float)
        state[0] = particles.as_initial_state()

        q_over_m = config.q_over_m()
        beta = particles.betas
        electric_field = config.electric_field
        magnetic_field = config.magnetic_field

        if _derivatives_numba is not None:
            deriv = lambda data: _derivatives_numba(data, config.gm_sun, q_over_m, beta, electric_field, magnetic_field)
        else:
            deriv = lambda data: _derivatives_numpy(
                data, config.gm_sun, q_over_m, beta, electric_field, magnetic_field
            )

        for i in range(1, num_steps):
            state[i] = _rk4_step(state[i - 1], dt, deriv)

        return SimulationResult(times=config.times, state=state, backend="python")


__all__ = ["PythonBackend"]
