"""Velocity fields for the two cartridge geometries.

Each design is a stream-function SHAPE with unit amplitude. fv_solver.scale_to_U0
then scales it to the target peak speed U0. psi = 0 on every wall, so nothing
flows out of the jar.

Design A (airlift): the circulation stays below the outlet holes, at
    y_out = 58 mm.
    planar : two counter-rotating slab cells, upflow along the centreline
             psi = sin(2 pi x / W) sin(pi y / y_out),  y < y_out; 0 above
    axisym : one toroidal cell, upflow on the axis
             psi = (r/R)^2 (1 - (r/R)^2) sin(pi z / y_out),  z < y_out; 0 above
    Above y_out the liquid only diffuses (optionally with D_t cut by 10x, which
    is the pessimistic stagnant-layer case).

Design B (passive, air stone off to the side at 0.35 W): one cell over the full
    depth, with the upflow core moved over to the stone. This is done with a
    smooth monotone coordinate map
        xi(s) = s + a s (1 - s),  s = x/W,  a = (0.5 - f)/(f(1 - f)),
    which puts the peak of sin(pi xi) at s = f = 0.35 and keeps
    xi' = 1 + a(1 - 2s) > 0.
    (A power-law stretch xi = s^p looks like it should do the same thing, but
    xi' blows up at the wall. Since the field gets scaled to a peak speed, that
    singular wall velocity takes magnitude away from the bulk as the grid is
    refined, and the solution never converges. The GCI study caught this, so
    I'm leaving it here as a warning.)

SOURCE LOCATIONS ARE ASSUMPTIONS, so swap them for the real CAD coordinates:
    A: outlet holes, 4 mm off the axis, 1 mm below y_out.
    B: basket slots at mid-width, 80% of the liquid height.
"""
import numpy as np
from .params import JAR
from .fv_solver import Grid, scale_to_U0, point_source


def _xi(x, W, frac):
    # normalised position, clipped so we never go outside [0, 1]
    s = x / W
    s = np.clip(s, 0, 1)
    a = (0.5 - frac) / (frac * (1 - frac))
    assert a < 1, "map not monotone for this stone offset"
    xi = s + a * s * (1 - s)
    return xi


def design_A_planar(n=56, U0=0.05, jar=JAR):
    g = Grid(jar.D, jar.H, n, n)
    X = g.Xf
    Y = g.Yf
    yo = jar.y_outlet

    # clamp Y at yo so the sine stays in range, then zero everything above it
    Y_clamped = np.minimum(Y, yo)
    inside = np.sin(2 * np.pi * X / jar.D) * np.sin(np.pi * Y_clamped / yo)
    psi = np.where(Y <= yo, inside, 0.0)
    psi = scale_to_U0(g, psi, U0)

    c0 = point_source(g, jar.D / 2 + 0.004, yo - 0.001)
    return g, psi, c0


def design_A_axisym(n=56, U0=0.05, jar=JAR):
    R = jar.D / 2
    g = Grid(R, jar.H, n // 2, n, axisym=True)
    r = g.Xf / R
    Z = g.Yf
    yo = jar.y_outlet

    Z_clamped = np.minimum(Z, yo)
    radial = r ** 2 * (1 - r ** 2)
    inside = radial * np.sin(np.pi * Z_clamped / yo)
    psi = np.where(Z <= yo, inside, 0.0)
    psi = scale_to_U0(g, psi, U0)

    c0 = point_source(g, 0.004, yo - 0.001)
    return g, psi, c0


def design_B_planar(n=56, U0=0.05, jar=JAR):
    g = Grid(jar.D, jar.H, n, n)
    X = g.Xf
    Y = g.Yf

    xi = _xi(X, jar.D, jar.x_stone_frac)
    horizontal = np.sin(np.pi * xi)
    vertical = np.sin(np.pi * Y / jar.H)
    psi = horizontal * vertical
    psi = scale_to_U0(g, psi, U0)

    c0 = point_source(g, 0.5 * jar.D, 0.8 * jar.H)
    return g, psi, c0


def stagnant_layer_D(g, Dt, jar=JAR, factor=0.1):
    """Dt field with the diffusivity cut by `factor` (10x) above the airlift outlets."""
    D = np.full((g.nx, g.ny), Dt)
    above = g.yc > jar.y_outlet
    D[:, above] = D[:, above] * factor
    return D


def randomized(g, psi, rng, U0):
    """Control #4: a random divergence-free field with the same peak speed.
    It is built from a random smooth psi (a sum of sine modes that respect the
    walls), so there is no coherent circulation tied to either geometry. If the
    A vs B differences survive this, they come from the SOURCE position and not
    from the flow structure."""
    X = g.Xf / g.W
    Y = g.Yf / g.H
    ps = np.zeros_like(psi)
    for i in range(12):
        kx, ky = rng.integers(1, 6, size=2)
        amp = rng.normal()
        mode = np.sin(kx * np.pi * X) * np.sin(ky * np.pi * Y)
        ps = ps + amp * mode
    if g.axisym:
        ps = ps * X ** 2          # keeps psi ~ r^2 near the axis
    return scale_to_U0(g, ps, U0)


DESIGNS = {
    "A_planar": design_A_planar,
    "A_axisym": design_A_axisym,
    "B_planar": design_B_planar,
}
