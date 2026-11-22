# Mixing time vs biology time. Da = t_mix / t_bio. If Da << 1 the jar is well mixed.
import csv
from pathlib import Path
import numpy as np
from .params import JAR, G, MU, BIO_TIMESCALES

RESULTS = Path(__file__).resolve().parents[1] / "results"


def V_c_joshi_sharma(jar=JAR):
    # Joshi & Sharma 1979, SI. simplified form, not checked against the paper (no full text)
    return 1.31 * (G * jar.D * jar.U_G) ** (1 / 3)


def D_L_deckwer(jar=JAR):
    # Deckwer et al. 1974 Eq. 1: D = 2.7 d^1.4 uG^0.3, all cgs (cm, cm/s, cm^2/s).
    # their data: 15 and 20 cm columns (plus other authors'), mean error 30%
    D_cm = jar.D * 100
    UG_cms = jar.U_G * 100
    D_L_cgs = 2.7 * D_cm ** 1.4 * UG_cms ** 0.3
    return D_L_cgs * 1e-4


def D_L_baird_rice(jar=JAR):
    # Baird & Rice 1975 Eq. 11, K = 0.35 dimensionless. "order-of-magnitude" only;
    # data d = 8.2-107 cm, uG = 0.3-45 cm/s; caution for small d and L/d < 5
    return 0.35 * jar.D ** (4 / 3) * (G * jar.U_G) ** (1 / 3)


def kappa(x0, L, tol=0.05):
    # one Fourier mode, 1-D. release at x0. order of magnitude only.
    amp = 2 * abs(np.cos(np.pi * x0 / L))
    return np.log(amp / tol) / np.pi ** 2


KAPPA_SCALING = 1.0
KAPPA_END = kappa(0.0, 1.0)                          # 0.374 worst case
KAPPA_RELEASE = kappa(JAR.y_outlet - 0.001, JAR.H)   # 0.346 at 57 mm
KAPPA_5PCT = KAPPA_END


def t_mix(D_L, L=JAR.H, kappa=KAPPA_SCALING):
    return kappa * L ** 2 / D_L


def damkohler_table(jar=JAR):
    # bracket over both correlations and all three kappas
    tmixes = []
    for D in [D_L_deckwer(jar), D_L_baird_rice(jar)]:
        for k in [KAPPA_RELEASE, KAPPA_END, KAPPA_SCALING]:
            tmixes.append(t_mix(D, jar.H, k))
    lo = min(tmixes)
    hi = max(tmixes)
    return (lo, hi), damkohler_bracket((lo, hi))


def damkohler_bracket(tmix_range):
    lo, hi = tmix_range
    out = {}
    for name, (tb_lo, tb_hi) in BIO_TIMESCALES.items():
        out[name] = (lo / tb_hi, hi / tb_lo)
    return out


def settling_velocity(d=1e-6, drho=100.0):
    # Stokes settling of one cell. d and drho are guesses
    return drho * G * d ** 2 / (18 * MU)


def solver_tmix_range(path=None):
    # min and max t_uniform from run_designs.py output
    path = Path(path) if path else RESULTS / "t_uniform.csv"
    if not path.exists():
        return None
    times = []
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["t_uniform"] not in ("", "nan"):
                times.append(float(row["t_uniform"]))
    return min(times), max(times)


def L_star(D_L, t_bio, kappa=KAPPA_SCALING, Da=1.0):
    # a method, not a channel prediction: uses the jar's D_L
    return np.sqrt(Da * t_bio * D_L / kappa)


def L_star_bracket(jar=JAR):
    t_fast = min(lo for lo, hi in BIO_TIMESCALES.values())
    vals = {}
    for cname, D in [("Deckwer", D_L_deckwer(jar)), ("Baird-Rice", D_L_baird_rice(jar))]:
        vals[f"{cname}, kappa=1"] = L_star(D, t_fast, KAPPA_SCALING)
        vals[f"{cname}, kappa={KAPPA_END:.3f} (end)"] = L_star(D, t_fast, KAPPA_END)
        vals[f"{cname}, kappa={KAPPA_RELEASE:.3f} (57 mm)"] = L_star(D, t_fast, KAPPA_RELEASE)
    return t_fast, vals


def summary(jar=JAR):
    dD = D_L_deckwer(jar)
    dB = D_L_baird_rice(jar)
    lines = []
    lines.append(f"U_G = {jar.U_G*1e3:.3f} mm/s,  V_c = {V_c_joshi_sharma(jar):.4f} m/s")
    lines.append(f"kappa: end {KAPPA_END:.3f}, at 57 mm {KAPPA_RELEASE:.3f}")
    for name, d in [("Deckwer   ", dD), ("Baird-Rice", dB)]:
        lines.append(f"D_L {name} = {d:.3e} m^2/s   t_mix k=1: {t_mix(d):.1f} s   "
                     f"k=end: {t_mix(d, kappa=KAPPA_END):.1f} s   k=57mm: {t_mix(d, kappa=KAPPA_RELEASE):.1f} s")
    lines.append(f"Correlation disagreement: {dB/dD:.2f}x")

    (lo, hi), da = damkohler_table(jar)
    lines.append(f"t_mix bracket (correlations, 1-D kappa) = {lo:.1f} - {hi:.1f} s")
    for name, (a, b) in da.items():
        lines.append(f"  Da[{name:13s}] = {a:.1e} - {b:.1e}")

    sr = solver_tmix_range()
    if sr is None:
        lines.append("no results/t_uniform.csv yet, run scripts/run_designs.py")
    else:
        lines.append(f"t_mix bracket (2-D solver, all designs/sweeps) = {sr[0]:.1f} - {sr[1]:.1f} s   <- USE THIS ONE")
        for name, (a, b) in damkohler_bracket(sr).items():
            lines.append(f"  Da[{name:13s}] = {a:.1e} - {b:.1e}   (worst-case separation {1/b:,.0f}x)")

    vs = settling_velocity()
    lines.append(f"Cell settling (Stokes, d=1 um, drho=100): {vs:.1e} m/s -> {jar.H/vs/86400:.0f} d to fall H; "
                 f"vs V_c ratio {vs/V_c_joshi_sharma(jar):.0e}")

    t_fast, Ls = L_star_bracket(jar)
    lines.append(f"L* METHOD ILLUSTRATION with jar D_L (Da=1 vs {t_fast/3600:.0f} h):")
    for name, v in Ls.items():
        lines.append(f"  {name:34s} L* = {v:.2f} m")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
