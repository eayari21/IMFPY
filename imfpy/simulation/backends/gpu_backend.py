"""CuPy powered backend for GPU acceleration."""
from __future__ import annotations

import logging

from ..config import SimulationConfig
from ..results import SimulationResult
from .python_backend import _rk4_step

LOGGER = logging.getLogger(__name__)

try:  # pragma: no cover - optional dependency
    import cupy as cp
except Exception:  # pragma: no cover - optional dependency
    cp = None


class GPUNotAvailableError(RuntimeError):
    """Raised when GPU execution is requested but CuPy is unavailable."""


class GPUBackend:
    """GPU accelerated backend built on top of CuPy."""

    def __init__(self) -> None:
        if cp is None:
            raise GPUNotAvailableError(
                "CuPy is not installed. Install cupy-cudaXX to enable GPU acceleration."
            )

    def run(self, config: SimulationConfig) -> SimulationResult:
        particles = config.particles
        dt = config.time_step
        num_steps = config.num_steps

        q_over_m = cp.asarray(config.q_over_m())
        beta = cp.asarray(particles.betas)
        electric_field = cp.asarray(config.electric_field)
        magnetic_field = cp.asarray(config.magnetic_field)

        state = cp.zeros((num_steps, particles.count, 6), dtype=cp.float64)
        state[0] = cp.asarray(particles.as_initial_state())

        def deriv(data: cp.ndarray) -> cp.ndarray:
            positions = data[:, :3]
            velocities = data[:, 3:]
            r = cp.linalg.norm(positions, axis=1)
            r3 = cp.where(r > 1e-12, r ** 3, cp.inf)
            grav_coeff = -config.gm_sun * (1.0 - beta) / r3
            grav = positions * grav_coeff[:, None]
            cross = cp.cross(velocities, magnetic_field)
            lorentz = q_over_m[:, None] * (electric_field + cross)
            acc = grav + lorentz
            return cp.hstack((velocities, acc))

        for i in range(1, num_steps):
            state[i] = _rk4_step(state[i - 1], dt, deriv)

        final_state = cp.asnumpy(state)
        return SimulationResult(times=config.times, state=final_state, backend="gpu")


__all__ = ["GPUBackend", "GPUNotAvailableError"]
