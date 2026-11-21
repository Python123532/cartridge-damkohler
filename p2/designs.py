# Flow fields for the two cartridges, as a stream function psi (zero on every wall).
# Vertical velocity is -dpsi/dx, so the plume is where psi drops fastest.
# A: airlift, flow only below the outlet holes. B: plume at the stone, 0.35 W.
# Source positions are guesses until the CAD numbers are in.
import numpy as np
from .params import JAR, STAGNANT_FACTOR
from .fv_solver import Grid, scale_to_U0, point_source, face_fluxes


def stretch(x, W, frac):
    # smooth map that moves x = frac*W to the middle. (x/W)**p looks simpler but
    # blows up at the wall and the grid never converges, found that with the GCI
    s = np.clip(x / W, 0, 1)
    a = (0.5 - frac) / (frac * (1 - frac))
    return s + a * s * (1 - s)


def design_A_planar(n=56, U0=0.05, jar=JAR):
    g = Grid(jar.D, jar.H, n, n)
    yo = jar.y_outlet
    below = g.Yf <= yo
    psi = np.sin(2 * np.pi * g.Xf / jar.D) * np.sin(np.pi * np.minimum(g.Yf, yo) / yo)
    psi = np.where(below, psi, 0.0)
    psi = scale_to_U0(g, psi, U0)
    c0 = point_source(g, jar.D / 2 + 0.004, yo - 0.001)
    return g, psi, c0


def design_A_axisym(n=56, U0=0.05, jar=JAR):
    R = jar.D / 2
    g = Grid(R, jar.H, n // 2, n, axisym=True)
    yo = jar.y_outlet
    r = g.Xf / R
    below = g.Yf <= yo
    # minus sign makes the flow go UP the axis
    psi = -(r ** 2) * (1 - r ** 2) * np.sin(np.pi * np.minimum(g.Yf, yo) / yo)
    psi = np.where(below, psi, 0.0)
    psi = scale_to_U0(g, psi, U0)
    c0 = point_source(g, 0.004, yo - 0.001)
    return g, psi, c0


def design_B_planar(n=56, U0=0.05, jar=JAR):
    g = Grid(jar.D, jar.H, n, n)
    xi = stretch(g.Xf, jar.D, jar.x_stone_frac)
    psi = np.sin(2 * np.pi * xi) * np.sin(np.pi * g.Yf / jar.H)
    psi = scale_to_U0(g, psi, U0)
    c0 = point_source(g, 0.5 * jar.D, 0.8 * jar.H)
    return g, psi, c0


def vertical_velocity(g, psi):
    _, Qy = face_fluxes(g, psi)
    return 0.5 * (Qy[:, 1:] / g.Ay[:, 1:] + Qy[:, :-1] / g.Ay[:, :-1])


def stagnant_layer_D(g, Dt, jar=JAR, factor=STAGNANT_FACTOR):
    D = np.full((g.nx, g.ny), Dt)
    D[:, g.yc > jar.y_outlet] *= factor
    return D


DESIGNS = {
    "A_planar": design_A_planar,
    "A_axisym": design_A_axisym,
    "B_planar": design_B_planar,
}
