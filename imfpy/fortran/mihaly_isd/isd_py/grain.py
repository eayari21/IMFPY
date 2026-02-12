# isd_py/grain.py
from dataclasses import dataclass
import numpy as np
from .constants import ConstantsCGS

@dataclass
class Grain:
    # size (radius) in cm
    gr: float
    # grain mass in g
    gm: float
    # beta = Frad/Fgrav
    beta: float
    # surface area
    gs: float
    # constants used in charging (kept for completeness)
    clpf: float
    cphe: float
    # Fortran uses these in current model; we keep them but default to simplified usage
    delmax: float
    emaxse: float
    emaxp4k: float
    tsetk: float
    eqgsp4: float
    elcc: float
    eqpktp: float
    f1: float

def betacalc_from_file(gm: float, path="beta-2.dat"):
    """
    Fortran betacalc reads 100 lines of (tm, pbeta) and linearly interpolates.
    It also enforces beta=0.1 if gm < 1.02e-17.
    """
    if gm < 1.02e-17:
        return 0.1

    data = np.loadtxt(path, dtype=float)
    tm = data[:, 0]
    pb = data[:, 1]

    # find bracket
    idx = np.searchsorted(tm, gm) - 1
    idx = int(np.clip(idx, 0, len(tm) - 2))

    rr = (gm - tm[idx]) / (tm[idx + 1] - tm[idx])
    beta = pb[idx] + rr * (pb[idx + 1] - pb[idx])
    return float(beta)

def dust_setup(gr_cm: float,
               const: ConstantsCGS,
               dust_type: float = 2.0,
               dau: float = 1.0,
               beta_file="beta-2.dat") -> Grain:
    """
    Port of Fortran subroutine dust() with the important parts retained.

    Notes:
      - density ro is hard-set to 2.5 g/cm^3 in the Fortran.
      - Q_pr interpolation tables exist in Fortran; we preserve them for clpf/cphe.
      - beta is computed from beta-2.dat via betacalc().
    """
    pi = const.pi
    bk = const.bk
    eq = const.eq
    em = const.em
    c = const.c

    # Meyer-Vernet / photoemission params from startup
    # hi=0.1 for dielectric (type=2), else 1.0
    hi = 0.1 if dust_type == 2.0 else 1.0
    tp = 2.0 * 1.6e-12
    eqpktp = eq / tp
    f1 = 2.5e10 * hi / (dau**2)

    # Secondary emission params from startup
    delsec = 1.5
    emaxse0 = 400.0
    delmax = 2.0 * 3.7 * delsec
    emaxse = emaxse0 * 1.6e-12
    emaxp4k = emaxse / 4.0 / bk
    tsetk = 2.5 * 1.6e-12

    # Grain physical properties
    ro = 2.5
    gs = 4.0 * pi * gr_cm**2
    eqgsp4 = eq * gs / 4.0
    elcc = eq * gr_cm**2 * np.sqrt(8.0 * bk * pi / em)
    gm = (4.0 * pi / 3.0) * gr_cm**3 * ro

    # beta
    beta = betacalc_from_file(gm, path=beta_file)

    # Qpr tables from Fortran dust()
    r = np.array([0.,1.e-5,1.44e-5,1.77e-5,2.04e-5,2.98e-5,4.51e-5,6.63e-5,
                  1.02e-4,2.36e-4,5.57e-4,1.3e-3,2.94e-3,6.51e-3,1.42e-2,3.08e-2,
                  6.66e-2,1.44e-1,3.1e-1,6.68e-1,1.0], dtype=float)

    q1 = np.array([0.00,1.54,1.88,2.03,2.07,1.97,1.73,1.52,1.38,1.18,1.06,
                   *([1.0]*10)], dtype=float)
    q2 = np.array([0.00,0.31,0.59,0.76,0.88,1.05,1.18,1.10,0.96,0.74,0.63,
                   0.65,0.74,0.87, *([1.0]*7)], dtype=float)

    # Find interpolation bin i such that gr < r[i]
    i = 1
    for ii in range(1, 18):
        if gr_cm < r[ii]:
            i = ii
            break
    else:
        i = 20  # end behavior

    if i == 1:
        rr = 1.0
    else:
        rr = (gr_cm - r[i-1]) / (r[i] - r[i-1])

    if i != 1:
        qpr = q1[i-1] + (q1[i] - q1[i-1]) * rr
        if dust_type == 2.0:
            qpr = q2[i-1] + (q2[i] - q2[i-1]) * rr
    else:
        # rare small-grain branch; keep Fortran behavior
        dr = r[1] - gr_cm
        qpr = (q1[1] - 0.3) * np.exp(-dr * 1.e5) + 0.3
        if dust_type == 2.0:
            qpr = q2[1] * np.exp(-dr * 5.5e5)

    # Light pressure force constant at 1 AU
    clpf = np.pi * gr_cm**2 / c * 1.36e6 * qpr / (dau**2)
    cphe = gs / 4.0 * eq * f1

    return Grain(
        gr=gr_cm, gm=float(gm), beta=float(beta),
        gs=float(gs), clpf=float(clpf), cphe=float(cphe),
        delmax=float(delmax), emaxse=float(emaxse), emaxp4k=float(emaxp4k),
        tsetk=float(tsetk), eqgsp4=float(eqgsp4), elcc=float(elcc),
        eqpktp=float(eqpktp), f1=float(f1)
    )
