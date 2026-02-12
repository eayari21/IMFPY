# isd_py/field_model.py
from dataclasses import dataclass
import numpy as np
from .constants import ConstantsCGS

@dataclass
class FieldModel:
    const: ConstantsCGS
    t0_years: float   # Fortran t0 (years offset from 1964)
    fi0: float = 0.0  # radians

    # Fortran constants inside fields()
    omegasun: float = 2.86533e-6
    b0: float = 3.5e-5

    def fields(self, t_sec: float, x: float, y: float, z: float):
        """
        Returns:
          (ex, ey, ez), (bx, by, bz), (vplasx, vplasy, vplasz), p, tilt, thetad, theta
        """
        c = self.const.c
        rs = self.const.rs
        oneau = self.const.oneau
        pi = self.const.pi

        r = np.sqrt(x*x + y*y + z*z)
        rxy = np.sqrt(x*x + y*y)

        # Avoid singularity exactly on axis
        if rxy == 0.0:
            rxy = 1e-30

        sinfi = y / rxy
        cosfi = x / rxy
        sinth = rxy / r
        costh = z / r

        # "NEW": radial-dependent solar wind velocity replaced by constant 400 km/s
        vr = 400.0e5
        vth = 0.0
        vfi = 0.0

        # polar->cartesian (solar wind velocity)
        vswx = vr*sinth*cosfi + vth*costh*cosfi - vfi*sinfi
        vswy = vr*sinth*sinfi + vth*costh*sinfi + vfi*cosfi
        vswz = vr*costh - vth*sinth
        vsw = float(np.sqrt(vswx*vswx + vswy*vswy + vswz*vswz))

        # Source surface
        rss = 2.5 * rs

        # Dust spherical angles
        cthetad = z / r
        sthetad = rxy / r
        thetad = float(np.arccos(cthetad))
        fid = float(np.arctan2(y, x))
        if fid < 0.0:
            fid += 2.0*pi

        # Solar wind travel time from rss to r
        tv = (r - rss) / vsw

        # Tilt model (piecewise) — identical logic
        rm1 = pi/8.0
        rm2 = pi/14.0
        b2  = pi/2.0 - 4.0*rm2
        b3  = pi - 11.0*rm1
        b4  = 3.0*pi/2.0 - 15.0*rm2
        b5  = 2.0*pi - 22.0*rm1
        b6  = 5.0*pi/2.0 - 26.0*rm2
        b7  = 3.0*pi - 33.0*rm1
        b8  = 7.0*pi/2.0 - 37.0*rm2
        b9  = 4.0*pi - 44.0*rm1
        b10 = 9.0*pi/2.0 - 48.0*rm2

        tt = self.t0_years + (t_sec - tv) / 365.0 / 86400.0
        tyear = tt

        # Default
        tilt = 0.0
        p0 = 1.0

        if 0.0 <= tyear < 4.0:
            tilt = rm1*tt
            p0 = -1.0
        elif 4.0 <= tyear < 11.0:
            tilt = rm2*tt + b2
            p0 = 1.0
        elif 11.0 <= tyear < 15.0:
            tilt = rm1*tt + b3
            p0 = 1.0
        elif 15.0 <= tyear < 22.0:
            tilt = rm2*tt + b4
            p0 = -1.0
        elif 22.0 <= tyear < 26.0:
            tilt = rm1*tt + b5
            p0 = -1.0
        elif 26.0 <= tyear < 33.0:
            tilt = rm2*tt + b6
            p0 = 1.0
        elif 33.0 <= tyear < 37.0:
            tilt = rm1*tt + b7
            p0 = 1.0
        elif 37.0 <= tyear < 44.0:
            tilt = rm2*tt + b8
            p0 = -1.0
        elif 44.0 <= tyear < 48.0:
            tilt = rm1*tt + b9
            p0 = -1.0
        elif 48.0 <= tyear < 55.0:
            tilt = rm2*tt + b10
            p0 = 1.0

        sintilt = abs(np.sin(tilt))
        fi00 = fid - self.fi0 - self.omegasun*t_sec + tv*self.omegasun
        sinfi00 = np.sin(fi00)

        # Pei 2012 form (correct for all tilts)
        tgtilt = abs(np.tan(tilt))
        theta = float(pi/2.0 - np.arctan(tgtilt*sinfi00))

        # Polarity p
        p = p0 if (thetad < theta) else -p0

        # Parker spiral Br, Bphi (Burger 2005)
        r0 = oneau
        br = p*self.b0*(r0/r)**2
        if r > rss:
            bfi = -p*self.b0*(r0/r)**2 * self.omegasun*(r - rss)*sthetad/vsw
        else:
            bfi = 0.0
        bth = 0.0

        # polar->cartesian for B
        bx = br*sinth*cosfi + bth*costh*cosfi - bfi*sinfi
        by = br*sinth*sinfi + bth*costh*sinfi + bfi*cosfi
        bz = br*costh - bth*sinth

        # convective E = -v x B / c
        ex = -(vswy*bz - vswz*by)/c
        ey = -(vswz*bx - vswx*bz)/c
        ez = -(vswx*by - vswy*bx)/c

        return (float(ex), float(ey), float(ez)), (float(bx), float(by), float(bz)), (float(vswx), float(vswy), float(vswz)), float(p), float(tilt), float(thetad), float(theta)
