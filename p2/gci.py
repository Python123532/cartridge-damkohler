# Grid Convergence Index (Roache). f1 = finest grid, r = refinement ratio.
import numpy as np


def gci(f1, f2, f3, r=2.0, Fs=1.25):
    e21 = f1 - f2
    e32 = f2 - f3
    monotone = e21 * e32 > 0
    p = np.nan
    if monotone and e21 != 0:
        p = np.log(abs(e32 / e21)) / np.log(r)
    asymptotic = bool(monotone and np.isfinite(p) and 0.5 <= p <= 3.0)
    if not asymptotic:
        # not converging cleanly: assume first order, bigger safety factor
        p = 1.0
        Fs = 3.0
    f_ext = f1 + e21 / (r ** p - 1)
    g = Fs * abs(e21 / f1) / (r ** p - 1)
    return dict(p=p, f_ext=f_ext, gci=g, monotone=bool(monotone), asymptotic=asymptotic)
