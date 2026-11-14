r"""Jar-scale transport time scales and the Damkohler argument.

The project's central claim:

    Da_i = t_mix / t_bio,i   for each downstream biological process i.

If Da_i << 1 for every i, transport is not rate-limiting and the geometry that
sets t_mix cannot move the outcome. The crossover length L* is where Da = 1 for
the fastest biological process

Mixing-time prefactor kappa
---------------------------
t_mix = kappa * L^2 / D_L, kappa = 1 (scaling estimate). For
a 1-D closed domain with an END release, the slowest Fourier mode is
2 cos(pi x / L) exp(-pi^2 D t / L^2); requiring it to fall below the +/-5%
uniformity tolerance gives

    kappa_5% = ln(2/0.05) / pi^2 = 0.374.

So the scaling estimate is conservative by ~2.7x. 
"""
import numpy as np
from .params import JAR, G, BIO_TIMESCALES


def U_G(jar=JAR):
    return jar.U_G


def V_c_joshi_sharma(jar=JAR):
    """Joshi & Sharma 1979 liquid circulation velocity, m/s (SI)."""
    return 1.31 * (G * jar.D * jar.U_G) ** (1 / 3)


def D_L_deckwer(jar=JAR):
    """Deckwer et al. 1974 axial dispersion; correlation is in cgs units."""
    D_cm, UG_cms = jar.D * 100, jar.U_G * 100
    return 0.678 * D_cm ** 1.4 * UG_cms ** 0.3 * 1e-4      # cm^2/s -> m^2/s


def D_L_baird_rice(jar=JAR):
    """Baird & Rice 1975 (Chem Eng J 9:171), isotropic-turbulence derivation:
    D_L = 0.35 D^{4/3} (g U_G)^{1/3}, SI. Second correlation = computational
    control #3: if the conclusion flips between correlations it is not robust."""
    return 0.35 * jar.D ** (4 / 3) * (G * jar.U_G) ** (1 / 3)


KAPPA_SCALING = 1.0
KAPPA_5PCT = np.log(2 / 0.05) / np.pi ** 2      # 0.374


def t_mix(D_L, L=JAR.H, kappa=KAPPA_SCALING):
    return kappa * L ** 2 / D_L


def damkohler_table(jar=JAR):
    """Da = t_mix / t_bio, bracketed over correlation, kappa and bio range.
    Returns dict process -> (Da_min, Da_max) plus the t_mix bracket."""
    DLs = [D_L_deckwer(jar), D_L_baird_rice(jar)]
    tmixes = [t_mix(d, jar.H, k) for d in DLs for k in (KAPPA_5PCT, KAPPA_SCALING)]
    tm_lo, tm_hi = min(tmixes), max(tmixes)
    out = {}
    for name, (tb_lo, tb_hi) in BIO_TIMESCALES.items():
        out[name] = (tm_lo / tb_hi, tm_hi / tb_lo)
    return (tm_lo, tm_hi), out


def L_star(D_L, t_bio, kappa=KAPPA_SCALING, Da=1.0):
    r"""Length at which dispersive transport time equals Da * t_bio:
    kappa L^2 / D_L = Da t_bio  ->  L* = sqrt(Da t_bio D_L / kappa).

    CAVEAT carried with every number: this uses the JAR dispersion coefficient
    as a stand-in for a raft channel. A real channel's D_L depends on its own
    cross-section and aeration; any net recirculation flow U adds an advective
    route with L*_adv = U t_bio, which is far longer. So L* is a LOWER bound on
    the scale at which geometry can start to matter."""
    return np.sqrt(Da * t_bio * D_L / kappa)


def L_star_bracket(jar=JAR):
    t_fast = min(lo for lo, _ in BIO_TIMESCALES.values())   # 1 h root capture
    vals = {}
    for cname, D in (("Deckwer", D_L_deckwer(jar)), ("Baird-Rice", D_L_baird_rice(jar))):
        for kname, k in (("kappa=1", KAPPA_SCALING), ("kappa=0.374", KAPPA_5PCT)):
            vals[f"{cname}, {kname}"] = L_star(D, t_fast, k)
    return t_fast, vals


def summary(jar=JAR):
    lines = []
    lines.append(f"U_G            = {jar.U_G*1e3:.3f} mm/s")
    lines.append(f"V_c (J&S)      = {V_c_joshi_sharma(jar):.4f} m/s")
    dD, dB = D_L_deckwer(jar), D_L_baird_rice(jar)
    lines.append(f"D_L Deckwer    = {dD:.3e} m^2/s   t_mix(k=1) = {t_mix(dD):.1f} s   t_mix(k=0.374) = {t_mix(dD,kappa=KAPPA_5PCT):.1f} s")
    lines.append(f"D_L Baird-Rice = {dB:.3e} m^2/s   t_mix(k=1) = {t_mix(dB):.1f} s   t_mix(k=0.374) = {t_mix(dB,kappa=KAPPA_5PCT):.1f} s")
    lines.append(f"Correlation disagreement: {dB/dD:.2f}x")
    (tlo, thi), da = damkohler_table(jar)
    lines.append(f"t_mix bracket  = {tlo:.1f} - {thi:.1f} s")
    for k, (a, b) in da.items():
        lines.append(f"  Da[{k:13s}] = {a:.1e} - {b:.1e}   (separation {1/b:,.0f}x - {1/a:,.0f}x)")
    t_fast, Ls = L_star_bracket(jar)
    lines.append(f"L* (Da=1 vs fastest bio process, {t_fast/3600:.0f} h):")
    for k, v in Ls.items():
        lines.append(f"  {k:24s} L* = {v:.2f} m")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
