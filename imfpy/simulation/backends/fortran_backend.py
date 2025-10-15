"""Interface to the Fortran dust integrator."""
from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from ...fortran import FortranBuildError, load_module
from ..config import SimulationConfig
from ..results import SimulationResult

LOGGER = logging.getLogger(__name__)


class FortranBackend:
    """Run simulations using the Fortran RK4 integrator."""

    def __init__(self) -> None:
        try:
            self._module = load_module()
        except FortranBuildError:
            LOGGER.exception("Failed to load Fortran integrator")
            raise

    def run(self, config: SimulationConfig) -> SimulationResult:
        particles = config.particles
        initial_state = particles.as_initial_state().T.copy(order="F")
        q_over_m = config.q_over_m().astype(np.float64)
        beta = particles.betas.astype(np.float64)
        e_field = config.electric_field.astype(np.float64)
        b_field = config.magnetic_field.astype(np.float64)

        num_particles = particles.count
        num_steps = config.num_steps

        results = np.zeros((6, num_particles, num_steps), dtype=np.float64, order="F")

        LOGGER.debug(
            "Starting Fortran simulation: particles=%s, steps=%s, dt=%s", num_particles, num_steps, config.time_step
        )

        integrator = getattr(self._module, "integrate_particles")
        state, status = integrator(
            num_particles,
            num_steps,
            config.time_step,
            config.gm_sun,
            q_over_m,
            beta,
            e_field,
            b_field,
            initial_state,
            results,
        )

        # f2py returns (results, status)
        if status != 0:
            raise RuntimeError(f"Fortran integrator returned non-zero status {status}")

        final_state = np.transpose(state, (2, 1, 0)).copy()

        return SimulationResult(times=config.times, state=final_state, backend="fortran")


__all__ = ["FortranBackend"]
