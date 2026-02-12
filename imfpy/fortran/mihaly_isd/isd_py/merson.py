# isd_py/merson.py
import numpy as np

class MersonIntegrator:
    """
    Direct port of the Fortran Merson method logic.
    You give it an RHS function that returns (dy, igone, icrash, ...).
    We only use dy for stepping; igone/icrash are handled outside.
    """
    def __init__(self, acc=1e-5, h=1.0, hmin=1e-10, jtest=0):
        self.acc = float(acc)
        self.h = float(h)
        self.hmin = float(hmin)
        self.jtest = int(jtest)
        self.rzero = 1e-13

    def integrate(self, rhs_fun, t0, t1, y0):
        """
        Integrate from t0 to t1.
        Returns y(t1) and final step size self.h.
        """
        y = y0.astype(float).copy()
        z = float(t0)
        zend = float(t1)

        bcc = self.acc
        zmin = self.hmin
        itest = self.jtest
        s = self.h
        iswh = 0

        while True:
            hsv = s
            cof = zend - z
            if abs(s) >= abs(cof):
                s = cof
                if abs(cof / hsv) < self.rzero:
                    # effectively done
                    self.h = hsv
                    return y

                iswh = 1

            yz = y.copy()
            ht = (1.0/3.0) * s

            f1 = rhs_fun(z, yz)[0]
            z1 = z + ht
            a = ht * f1
            w = yz + a

            f2 = rhs_fun(z1, w)[0]
            a = 0.5 * a
            w = yz + 0.5*ht*f2 + a

            f3 = rhs_fun(z1, w)[0]
            z2 = z + 0.5*ht
            b = 4.5 * ht * f3
            w = yz + 0.25*b + 0.75*a

            f4 = rhs_fun(z2, w)[0]
            z3 = z + 0.5*s
            a = 2.0*ht*f4 + a
            w = yz + 3.0*a - b

            f5 = rhs_fun(z3, w)[0]

            # error estimate & accept/reject
            b = -0.5*ht*f5 - b + 2.0*a
            w_new = w - b

            # check each component
            ok = True
            for k in range(len(y)):
                ak = abs(5.0*bcc*w_new[k])
                bk = abs(b[k])
                if abs(w_new[k]) <= self.rzero:
                    continue
                if bk > ak:
                    ok = False
                    break

            if ok:
                # accept
                y = w_new
                z = z + s

                if iswh == 1:
                    self.h = hsv
                    return y

                # step doubling test
                can_double = True
                for k in range(len(y)):
                    ak = abs(5.0*bcc*y[k])
                    bk = abs(b[k])
                    if bk > 0.03125 * ak:
                        can_double = False
                        break
                if can_double:
                    s = 2.0*s

                # continue until reach zend
                continue

            # reject: halve step
            cof2 = 0.5 * s
            if abs(cof2) >= zmin:
                # rollback and retry smaller step
                y = yz
                s = cof2
                iswh = 0
                continue

            # too small
            if itest == 0:
                raise RuntimeError(f"Merson step underflow: h={s}, hmin={zmin}, t={z}")

            # jtest=1: force constant step hmin
            s = zmin if hsv >= 0 else -zmin
            if iswh == 1:
                self.h = hsv
                return y
