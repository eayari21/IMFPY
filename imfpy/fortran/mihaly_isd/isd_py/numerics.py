# isd_py/numerics.py
import numpy as np

def exxp(x: float) -> float:
    """Fortran clamps |x| to 70 before exp()."""
    xx = abs(x)
    arg = x
    if xx > 70.0:
        arg = (x / xx) * 70.0
    return float(np.exp(arg))

def erfc1(x: float) -> float:
    """
    Fortran approximation used (Numerical Recipes-like).
    """
    z = abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    arg = (-z*z - 1.26551223 +
           t*(1.00002368 + t*(0.37409196 +
           t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
           t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))))
    val = t * exxp(arg)
    if x < 0.0:
        val = 2.0 - val
    return float(val)

def erf1(x: float) -> float:
    return float(1.0 - erfc1(x))
