import numpy as np

# =========================
# Physical constants
# =========================
AU = 1.495978707e11        # m
G = 6.67430e-11            # m^3 kg^-1 s^-2
M_sun = 1.98847e30         # kg
GM = G * M_sun

# =========================
# Acceleration model
# =========================
def acceleration(r, beta):
    """
    Compute acceleration due to gravity + radiation pressure.

    Parameters
    ----------
    r : ndarray, shape (2,)
        Position vector [m]
    beta : float
        Radiation pressure coefficient

    Returns
    -------
    a : ndarray, shape (2,)
        Acceleration vector [m/s^2]
    """
    norm_r = np.linalg.norm(r)
    if norm_r == 0:
        return np.zeros(2)

    grav = -GM * r / norm_r**3
    rad  = beta * GM * r / norm_r**3

    return grav + rad


# =========================
# RK4 integrator
# =========================
def rk4_step(state, dt, beta):
    """
    One RK4 step.

    state = [x, y, vx, vy]
    """
    def deriv(s):
        r = s[:2]
        v = s[2:]
        a = acceleration(r, beta)
        return np.hstack((v, a))

    k1 = deriv(state)
    k2 = deriv(state + 0.5 * dt * k1)
    k3 = deriv(state + 0.5 * dt * k2)
    k4 = deriv(state + dt * k3)

    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


# =========================
# Trajectory propagation
# =========================
def propagate(state0, t_end, dt, beta):
    """
    Propagate trajectory.

    Returns
    -------
    states : ndarray, shape (N, 4)
    """
    steps = int(t_end / dt)
    states = np.zeros((steps, 4))
    states[0] = state0

    for i in range(1, steps):
        states[i] = rk4_step(states[i-1], dt, beta)

    return states


# =========================
# Example usage
# =========================
if __name__ == "__main__":
    # Initial conditions at 5 AU
    r0 = np.array([5 * AU, 0.0])
    v0 = np.array([0.0, -26e3])  # interstellar inflow speed

    beta = 0.5
    dt = 3600.0                 # 1 hour
    t_end = 2.0 * 365.25 * 24 * 3600  # 2 years

    state0 = np.hstack((r0, v0))
    traj = propagate(state0, t_end, dt, beta)

    print("Final position [AU]:", traj[-1, :2] / AU)
