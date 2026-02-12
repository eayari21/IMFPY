# isd_py/constants.py
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class ConstantsCGS:
    # Fundamental constants used in the Fortran startup()
    eq: float = 4.803e-10     # statcoulomb
    em: float = 9.108e-28     # g
    pm: float = 1.660e-24     # g (proton mass proxy used in code)
    gc: float = 6.668e-8      # cgs G
    c: float  = 2.998e10      # cm/s
    bk: float = 1.380e-16     # erg/K

    # Solar / heliocentric
    rs: float = 6.96e10       # cm
    sm: float = 1.98892e33    # g
    oneau: float = 149598.0e8 # cm  (same as Fortran)
    pi: float = float(np.arccos(-1.0))

    # Derived helpers (computed on demand)
    @property
    def sqpi(self) -> float:
        return float(np.sqrt(self.pi))

    @property
    def eqpbk(self) -> float:
        return self.eq / self.bk

    @property
    def prcc1(self) -> float:
        return float(np.sqrt(2.0 * self.bk / self.pm))

    @property
    def prcc2(self) -> float:
        return float(np.sqrt(2.0 * self.eq / self.pm))
