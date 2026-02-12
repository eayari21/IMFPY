# isd_py/rhs.py
import numpy as np
from .forces import compute_forces

def rhs(t, y, const, grain, field_model, flags):
    """
    y = [x, y, z, vx, vy, vz, q]
    returns dy/dt
    """
    x, yv, z, vx, vy, vz, q = y

    E, B, vplas, p, tilt, thetad, theta = field_model.fields(t, x, yv, z)

    st = {"x": x, "y": yv, "z": z, "vx": vx, "vy": vy, "vz": vz, "q": q}
    fx, fy, fz, fsg, igone, icrash = compute_forces(
        const, grain, st, E, B,
        flags["ifgr"], flags["ifrp"], flags["ifl"], flags["nochar"]
    )

    if igone or icrash:
        return np.zeros_like(y), igone, icrash, p, tilt, thetad, theta

    dxdt = np.array([
        vx,
        vy,
        vz,
        fx / grain.gm,
        fy / grain.gm,
        fz / grain.gm,
        0.0
    ], dtype=float)

    return dxdt, 0, 0, p, tilt, thetad, theta
