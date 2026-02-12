# isd_py/forces.py
import numpy as np
from .constants import ConstantsCGS
from .grain import Grain

def compute_forces(const: ConstantsCGS,
                   grain: Grain,
                   state,   # dict with x,y,z,vx,vy,vz,q
                   E, B,
                   ifgr: int, ifrp: int, ifl: int,
                   nochar: int):
    """
    Returns:
      fx, fy, fz, fsg, igone, icrash
    """
    x = state["x"]; y = state["y"]; z = state["z"]
    vx = state["vx"]; vy = state["vy"]; vz = state["vz"]
    q = state["q"]

    r2 = x*x + y*y + z*z
    r = float(np.sqrt(r2))
    rxy = float(np.sqrt(x*x + y*y))

    # crash if close to Sun
    if r <= 2.5 * const.rs:
        return 0.0, 0.0, 0.0, 0.0, 0, 1

    # igone: simple criterion (Fortran uses yau>10)
    yau = y / const.oneau
    igone = 1 if (yau > 10.0) else 0
    if igone:
        return 0.0, 0.0, 0.0, 0.0, 1, 0

    # Gravity + radiation pressure
    fsgr = const.gc * const.sm * grain.gm / r2
    fac = (1.0 - ifrp * grain.beta)

    fsgrx = -fsgr * fac * x / r
    fsgry = -fsgr * fac * y / r
    fsgrz = -fsgr * fac * z / r
    fsg = float(np.sqrt(fsgrx*fsgrx + fsgry*fsgry + fsgrz*fsgrz))

    fx = ifgr * fsgrx
    fy = ifgr * fsgry
    fz = ifgr * fsgrz

    # Poynting-Robertson drag (Krivov 1998 Eq.4)
    vr = (x*vx + y*vy + z*vz) / r
    fprx = -fsgr * grain.beta * (vr*x/r + vx) / const.c
    fpry = -fsgr * grain.beta * (vr*y/r + vy) / const.c
    fprz = -fsgr * grain.beta * (vr*z/r + vz) / const.c

    fx += ifrp * fprx
    fy += ifrp * fpry
    fz += ifrp * fprz

    # Lorentz force
    if nochar != 1:
        ex, ey, ez = E
        bx, by, bz = B

        fex = q * ex
        fey = q * ey
        fez = q * ez

        fbx = q * (vy*bz - vz*by) / const.c
        fby = q * (vz*bx - vx*bz) / const.c
        fbz = q * (vx*by - vy*bx) / const.c

        fx += ifl * (fex + fbx)
        fy += ifl * (fey + fby)
        fz += ifl * (fez + fbz)

    return float(fx), float(fy), float(fz), fsg, 0, 0
