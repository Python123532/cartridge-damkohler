# Checks that the solver solves the equations right. The dye test checks the equations are right.
# run: pytest -q
import numpy as np
import pytest

from p2.fv_solver import Grid, run, face_fluxes, divergence
from p2.designs import DESIGNS, vertical_velocity

NAMES = list(DESIGNS)


def test_no_flow_decay_matches_cosine_mode():
    # u = 0: amplitude should fall like exp(-pi^2 D t / L^2)
    L, D = 0.07, 1e-4
    g = Grid(L, L, 64, 4)
    c0 = 1 + 0.5 * np.cos(np.pi * g.Xc / L)
    no_flow = np.zeros((g.nx + 1, g.ny + 1))
    r = run(g, no_flow, D, c0, t_end=5.0, tol=-1, record_every=0.05)
    start = np.abs(c0 - c0.mean()).max()
    end = np.abs(r.c_final - r.c_final.mean()).max()
    rate = np.log(start / end) / r.t[-1]
    assert rate == pytest.approx(D * (np.pi / L) ** 2, rel=0.01)


@pytest.mark.parametrize("name", NAMES)
def test_mass_is_kept(name):
    g, psi, c0 = DESIGNS[name](28, 0.05)
    Qx, Qy = face_fluxes(g, psi)
    assert np.abs(divergence(g, Qx, Qy)).max() < 1e-15
    assert (c0 * g.vol).sum() == pytest.approx(1.0, rel=1e-12)
    r = run(g, psi, 3e-5, c0, t_end=20)
    assert r.mass_err < 1e-12


# the two Designs.py bugs passed the mass test with the flow going the wrong way.
# these catch them.
def test_A_axisym_goes_up_the_axis():
    g, psi, c0 = DESIGNS["A_axisym"](56, 0.05)
    u = vertical_velocity(g, psi)
    j = np.argmin(np.abs(g.yc - 0.025))
    assert u[0, j] > 0.8 * 0.05


def test_B_plume_is_at_the_stone():
    g, psi, c0 = DESIGNS["B_planar"](56, 0.05)
    u = vertical_velocity(g, psi)
    j = np.argmin(np.abs(g.yc - 0.5 * g.H))
    i = np.argmin(np.abs(g.xc - 0.35 * g.W))
    assert u[i, j] > 0.6 * 0.05
    assert abs(g.xc[np.argmax(u[:, j])] / g.W - 0.35) < 0.05
