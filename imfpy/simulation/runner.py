"""Simulation runner orchestrating available backends."""
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from typing import Any, Callable, Dict

import numpy as np

from .backends import FortranBackend, GPUBackend, GPUNotAvailableError, PythonBackend
from .config import BackendName, ParticleEnsemble, SimulationConfig
from .results import SimulationResult

LOGGER = logging.getLogger(__name__)


@dataclass
class BackendRegistry:
    factories: Dict[BackendName, Callable[[], Any]]

    def get(self, name: BackendName):
        try:
            return self.factories[name]()
        except KeyError as exc:  # pragma: no cover - defensive
            raise ValueError(f"Unknown backend {name}") from exc


DEFAULT_BACKENDS = BackendRegistry(
    factories={
        "fortran": FortranBackend,
        "python": PythonBackend,
        "gpu": GPUBackend,
    }
)


class SimulationRunner:
    """High-level orchestrator for dust simulations."""

    def __init__(self, registry: BackendRegistry = DEFAULT_BACKENDS) -> None:
        self.registry = registry

    def run(self, config: SimulationConfig) -> SimulationResult:
        LOGGER.info("Running simulation with backend=%s", config.backend)
        try:
            backend = self.registry.get(config.backend)
        except GPUNotAvailableError as exc:
            raise RuntimeError("GPU backend requested but unavailable") from exc
        result = backend.run(config)
        LOGGER.info("Simulation finished")
        return result


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a dust trajectory simulation.")
    parser.add_argument("--backend", default="fortran", choices=["fortran", "python", "gpu"], help="Backend to use")
    parser.add_argument("--particles", type=int, default=4, help="Number of particles")
    parser.add_argument("--steps", type=int, default=200, help="Number of integration steps")
    parser.add_argument("--dt", type=float, default=60.0, help="Time step size")
    parser.add_argument("--beta", type=float, default=0.1, help="Radiation pressure coefficient")
    parser.add_argument("--charge", type=float, default=1.6e-19, help="Charge per particle")
    parser.add_argument("--mass", type=float, default=1e-16, help="Mass per particle")
    return parser.parse_args(argv)


def _build_default_config(args: argparse.Namespace) -> SimulationConfig:
    radius = 1.5e11  # 1 AU
    speeds = 26_000.0

    positions = np.zeros((args.particles, 3))
    velocities = np.zeros((args.particles, 3))

    for i in range(args.particles):
        angle = 2 * np.pi * i / max(args.particles, 1)
        positions[i] = radius * np.array([np.cos(angle), np.sin(angle), 0.0])
        velocities[i] = speeds * np.array([-np.sin(angle), np.cos(angle), 0.0])

    ensemble = ParticleEnsemble(
        positions=positions,
        velocities=velocities,
        charges=np.full(args.particles, args.charge),
        masses=np.full(args.particles, args.mass),
        betas=np.full(args.particles, args.beta),
    )

    config = SimulationConfig(
        particles=ensemble,
        num_steps=args.steps,
        time_step=args.dt,
        backend=args.backend,
    )
    return config


def main(argv: list[str] | None = None) -> SimulationResult | None:
    logging.basicConfig(level=logging.INFO)
    args = _parse_args(argv)
    runner = SimulationRunner()
    config = _build_default_config(args)
    result = runner.run(config)
    LOGGER.info(
        "Result summary: backend=%s, particles=%s, steps=%s", result.backend, result.num_particles, result.num_steps
    )
    return result


if __name__ == "__main__":  # pragma: no cover
    main()
