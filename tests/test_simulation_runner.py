import numpy as np
import pytest

from imfpy.simulation import ParticleEnsemble, SimulationConfig, SimulationRunner


def test_python_backend_keeps_circular_orbit():
    radius = 1.0
    speed = 1.0
    positions = np.array([[radius, 0.0, 0.0]])
    velocities = np.array([[0.0, speed, 0.0]])
    charges = np.array([1.0])
    masses = np.array([1.0])
    betas = np.array([0.0])

    ensemble = ParticleEnsemble(
        positions=positions,
        velocities=velocities,
        charges=charges,
        masses=masses,
        betas=betas,
    )
    config = SimulationConfig(
        particles=ensemble,
        num_steps=10,
        time_step=0.01,
        gm_sun=1.0,
        backend="python",
    )
    result = SimulationRunner().run(config)
    final_position = result.positions()[-1, 0]
    assert np.linalg.norm(final_position) == pytest.approx(radius, rel=1e-2)
