# Finite volume: dc/dt + div(u c) = div(D grad c), closed walls.
# Fluxes come from psi at cell corners, so every cell's net flow is exactly zero.
# Planar and axisymmetric use the same code, only areas and volumes change.
import numpy as np
from dataclasses import dataclass, field

SOURCE_RADIUS = 2.5e-3   # m, size of an outlet hole. fixed so every grid gets the same release


class Grid:
    def __init__(self, W, H, nx, ny, axisym=False, depth=1.0):
        self.W, self.H, self.nx, self.ny = W, H, nx, ny
        self.axisym = axisym
        self.xf = np.linspace(0, W, nx + 1)
        self.yf = np.linspace(0, H, ny + 1)
        self.xc = 0.5 * (self.xf[1:] + self.xf[:-1])
        self.yc = 0.5 * (self.yf[1:] + self.yf[:-1])
        self.dx = W / nx
        self.dy = H / ny

        if axisym:
            self.k = 2 * np.pi
            ring = np.pi * (self.xf[1:] ** 2 - self.xf[:-1] ** 2)
            self.vol = np.outer(ring, np.full(ny, self.dy))
            self.Ax = np.outer(2 * np.pi * self.xf, np.full(ny, self.dy))
            self.Ay = np.outer(ring, np.ones(ny + 1))
        else:
            self.k = depth
            self.vol = np.full((nx, ny), self.dx * self.dy * depth)
            self.Ax = np.full((nx + 1, ny), self.dy * depth)
            self.Ay = np.full((nx, ny + 1), self.dx * depth)

        self.Xc, self.Yc = np.meshgrid(self.xc, self.yc, indexing="ij")
        self.Xf, self.Yf = np.meshgrid(self.xf, self.yf, indexing="ij")


def face_fluxes(grid, psi):
    # Qx through vertical faces, Qy through horizontal faces
    Qx = grid.k * (psi[:, 1:] - psi[:, :-1])
    Qy = -grid.k * (psi[1:, :] - psi[:-1, :])
    return Qx, Qy


def divergence(grid, Qx, Qy):
    return (Qx[1:, :] - Qx[:-1, :]) + (Qy[:, 1:] - Qy[:, :-1])


def cell_speed(grid, Qx, Qy):
    Ax = np.where(grid.Ax > 0, grid.Ax, np.inf)   # axis face has zero area
    ux = 0.5 * (Qx[1:, :] / Ax[1:, :] + Qx[:-1, :] / Ax[:-1, :])
    uy = 0.5 * (Qy[:, 1:] / grid.Ay[:, 1:] + Qy[:, :-1] / grid.Ay[:, :-1])
    return np.sqrt(ux ** 2 + uy ** 2)


def scale_to_U0(grid, psi, U0):
    Qx, Qy = face_fluxes(grid, psi)
    top = cell_speed(grid, Qx, Qy).max()
    if top == 0:
        return psi
    return psi * U0 / top


def point_source(grid, x, y, radius=SOURCE_RADIUS):
    radius = max(radius, 0.75 * max(grid.dx, grid.dy))
    inside = np.sqrt((grid.Xc - x) ** 2 + (grid.Yc - y) ** 2) <= radius
    c = np.zeros((grid.nx, grid.ny))
    c[inside] = 1.0 / grid.vol[inside].sum()
    return c


def vanleer(r):
    return (r + abs(r)) / (1 + abs(r))


def tvd_face_values(c, Q, axis):
    # second order face value with van Leer limiter. 2 ghost cells each side
    P = np.pad(c, [(2, 2), (0, 0)] if axis == 0 else [(0, 0), (2, 2)], mode="edge")
    n = c.shape[axis]
    LL = np.take(P, range(0, n + 1), axis=axis)
    L = np.take(P, range(1, n + 2), axis=axis)
    R = np.take(P, range(2, n + 3), axis=axis)
    RR = np.take(P, range(3, n + 4), axis=axis)
    jump = R - L
    safe = np.where(jump != 0, jump, 1.0)
    r_plus = np.where(jump != 0, (L - LL) / safe, 0.0)
    r_minus = np.where(jump != 0, (RR - R) / safe, 0.0)
    from_left = L + 0.5 * vanleer(r_plus) * jump
    from_right = R - 0.5 * vanleer(r_minus) * jump
    return np.where(Q >= 0, from_left, from_right)


@dataclass
class Result:
    t: np.ndarray
    dev: np.ndarray
    t_uniform: float
    mass_err: float
    c_final: np.ndarray
    snapshots: dict = field(default_factory=dict)


def run(grid, psi, Dt, c0, t_end=600.0, tol=0.05, cfl=0.4, record_every=0.5,
        snapshot_times=(), scheme="vanleer"):
    Qx, Qy = face_fluxes(grid, psi)
    Qx[0, :] = 0
    Qx[-1, :] = 0
    Qy[:, 0] = 0
    Qy[:, -1] = 0

    # diffusion conductance. harmonic mean so the stagnant layer jump is handled right
    D = np.broadcast_to(np.asarray(Dt, float), (grid.nx, grid.ny)).copy()
    Gx = np.zeros((grid.nx + 1, grid.ny))
    Gy = np.zeros((grid.nx, grid.ny + 1))
    Gx[1:-1] = 2 * D[1:] * D[:-1] / (D[1:] + D[:-1]) * grid.Ax[1:-1] / grid.dx
    Gy[:, 1:-1] = 2 * D[:, 1:] * D[:, :-1] / (D[:, 1:] + D[:, :-1]) * grid.Ay[:, 1:-1] / grid.dy

    # time step: advection and diffusion together
    out = (np.maximum(Qx[1:], 0) + np.maximum(-Qx[:-1], 0)
           + np.maximum(Qy[:, 1:], 0) + np.maximum(-Qy[:, :-1], 0)
           + Gx[1:] + Gx[:-1] + Gy[:, 1:] + Gy[:, :-1])
    dt = cfl * np.min(grid.vol / np.maximum(out, 1e-300))

    def rate(c):
        # change per second in every cell
        Px = np.pad(c, [(1, 1), (0, 0)], mode="edge")
        Py = np.pad(c, [(0, 0), (1, 1)], mode="edge")
        if scheme == "upwind":
            cx = np.where(Qx > 0, Px[:-1], Px[1:])
            cy = np.where(Qy > 0, Py[:, :-1], Py[:, 1:])
        else:
            cx = tvd_face_values(c, Qx, 0)
            cy = tvd_face_values(c, Qy, 1)
        Fx = Qx * cx - Gx * (Px[1:] - Px[:-1])
        Fy = Qy * cy - Gy * (Py[:, 1:] - Py[:, :-1])
        return -((Fx[1:] - Fx[:-1]) + (Fy[:, 1:] - Fy[:, :-1])) / grid.vol

    c = c0.astype(float).copy()
    mass0 = (c * grid.vol).sum()
    cbar = mass0 / grid.vol.sum()
    nsteps = int(np.ceil(t_end / dt))
    rec = max(1, int(round(record_every / dt)))
    snap_left = sorted(snapshot_times)
    times, devs, snaps = [], [], {}
    mass_err = 0.0
    t_uni = np.nan

    for n in range(nsteps + 1):
        t = n * dt
        if snap_left and t >= snap_left[0]:
            snaps[snap_left.pop(0)] = c.copy()
        if n % rec == 0:
            dev = np.abs(c - cbar).max() / cbar
            times.append(t)
            devs.append(dev)
            mass_err = max(mass_err, abs((c * grid.vol).sum() - mass0) / mass0)
            if dev <= tol and np.isnan(t_uni):
                t_uni = t
                if not snap_left:
                    break
        if n == nsteps:
            break
        if scheme == "upwind":
            c = c + dt * rate(c)
        else:
            # Heun (SSP-RK2)
            c1 = c + dt * rate(c)
            c = 0.5 * c + 0.5 * (c1 + dt * rate(c1))

    return Result(np.array(times), np.array(devs), t_uni, mass_err, c, snaps)
