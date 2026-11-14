"""Finite-volume scalar transport in a prescribed, divergence-free flow.

    dc/dt + div(u c) = div(D_t grad c),   closed walls (no flux through them).

A few choices I made on purpose:

1. Flux form from a stream function. The face volumetric fluxes are just
   differences of psi at the cell CORNERS. That means the net outflow of every
   cell cancels out to zero (down to roundoff), so mass conservation comes for
   free from the structure instead of being something we hope for.
   (tests/test_verification.py is supposed to check this.)

2. One code path for two geometries. Planar (x, y) and axisymmetric (r, z) only
   differ in the face areas, the cell volumes, and the factor that turns a psi
   difference into a volumetric flux:
       planar : Q = b * dpsi           (b = depth into the page)
       axisym : Q = 2*pi * dpsi        (Stokes stream function)
   Design A's centred airlift really is axisymmetric, so (r, z) is the honest
   2-D version of it. The planar model turns a ring-shaped cell into two slab
   cells. Design B (stone off to the side) has no symmetry axis so it stays
   planar. Comparing planar-A against axisym-A tells you how much the 2-D
   reduction itself costs us.

3. Advection. scheme="vanleer" (the default now) is second-order TVD with the
   van Leer limiter, stepped with SSP-RK2 (Heun). scheme="upwind" is the
   first-order scheme from the original 56x56 solver, kept around so we can
   compare.

   Upwind adds numerical diffusion, D_num ~ |u| dx / 2. At 56x56 with
   U0 = 0.05 m/s, dx = 1.33 mm, that works out to D_num ~ 3.3e-5 m^2/s, which
   is about the same size as the physical D_t we are trying to sweep. The cell
   Peclet number U dx / D_t is around 2. So the upwind answer is partly just a
   property of the grid. The GCI study shows this and the TVD scheme mostly
   gets rid of it.

   Diffusion is central, with a harmonic-mean face diffusivity (this is the
   right way to average across the jump at the stagnant layer). The explicit
   time step limit takes advection and diffusion together.
"""
import numpy as np
from dataclasses import dataclass, field


@dataclass
class Grid:
    W: float              # horizontal size: width (planar) or radius (axisym)
    H: float              # height
    nx: int
    ny: int
    axisym: bool = False
    depth: float = 1.0    # out-of-plane depth for planar (drops out of c anyway)

    def __post_init__(self):
        # cell faces and cell centres
        self.xf = np.linspace(0, self.W, self.nx + 1)
        self.yf = np.linspace(0, self.H, self.ny + 1)
        self.xc = 0.5 * (self.xf[1:] + self.xf[:-1])
        self.yc = 0.5 * (self.yf[1:] + self.yf[:-1])
        self.dx = self.W / self.nx
        self.dy = self.H / self.ny

        if self.axisym:
            self.k = 2 * np.pi
            # area of each annular ring, shape (nx,)
            ring = np.pi * (self.xf[1:] ** 2 - self.xf[:-1] ** 2)
            ones_y = np.full(self.ny, self.dy)
            self.vol = np.outer(ring, ones_y)
            # x-faces are cylinders, y-faces are rings
            self.Ax = np.outer(2 * np.pi * self.xf, np.full(self.ny, self.dy))  # (nx+1, ny)
            self.Ay = np.outer(ring, np.ones(self.ny + 1))                       # (nx, ny+1)
        else:
            self.k = self.depth
            self.vol = np.full((self.nx, self.ny), self.dx * self.dy * self.depth)
            self.Ax = np.full((self.nx + 1, self.ny), self.dy * self.depth)
            self.Ay = np.full((self.nx, self.ny + 1), self.dx * self.depth)

        self.Xc, self.Yc = np.meshgrid(self.xc, self.yc, indexing="ij")
        self.Xf, self.Yf = np.meshgrid(self.xf, self.yf, indexing="ij")


def face_fluxes(grid, psi):
    """psi lives at the corners, shape (nx+1, ny+1). Gives back Qx (nx+1, ny)
    and Qy (nx, ny+1), the volumetric fluxes in the +x and +y directions."""
    dpsi_y = psi[:, 1:] - psi[:, :-1]
    dpsi_x = psi[1:, :] - psi[:-1, :]
    Qx = grid.k * dpsi_y
    Qy = -grid.k * dpsi_x
    return Qx, Qy


def cell_speed(grid, Qx, Qy):
    # x-faces on the axis have zero area, so swap in inf there to avoid 0/0
    # (the flux is zero there anyway)
    Ax_right = np.where(grid.Ax[1:, :] > 0, grid.Ax[1:, :], np.inf)
    Ax_left = np.where(grid.Ax[:-1, :] > 0, grid.Ax[:-1, :], np.inf)

    ux_right = Qx[1:, :] / Ax_right
    ux_left = Qx[:-1, :] / Ax_left
    ux = 0.5 * (ux_right + ux_left)

    uy_top = Qy[:, 1:] / grid.Ay[:, 1:]
    uy_bottom = Qy[:, :-1] / grid.Ay[:, :-1]
    uy = 0.5 * (uy_top + uy_bottom)

    return np.sqrt(ux ** 2 + uy ** 2)


def scale_to_U0(grid, psi_shape, U0):
    """Scale a unit psi shape so the biggest cell-centre speed comes out as U0."""
    Qx, Qy = face_fluxes(grid, psi_shape)
    speed = cell_speed(grid, Qx, Qy)
    s = speed.max()
    if s > 0:
        return psi_shape * (U0 / s)
    return psi_shape


def divergence(grid, Qx, Qy):
    div_x = Qx[1:, :] - Qx[:-1, :]
    div_y = Qy[:, 1:] - Qy[:, :-1]
    return div_x + div_y


@dataclass
class Result:
    t: np.ndarray
    dev: np.ndarray          # max |c - cbar| / cbar at each recorded time
    t_uniform: float         # first time dev <= tol (nan if it never gets there)
    mass_err: float          # worst relative mass drift
    c_final: np.ndarray
    snapshots: dict = field(default_factory=dict)


def _vanleer(r):
    return (r + np.abs(r)) / (1 + np.abs(r))


def _ratio(num, den, eps):
    # gradient ratio r = num / den, but just 0 where the denominator is ~0
    ok = np.abs(den) > eps
    safe_den = np.where(ok, den, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(ok, num / safe_den, 0.0)
    return r


def _tvd_face_x(P, Q):
    """P is c padded with 2 zero-gradient ghost cells at each end in x.
    Returns the face values (nx+1 of them) using the van Leer limiter."""
    n = P.shape[0] - 4
    LL = P[0:n + 1]
    L = P[1:n + 2]
    R = P[2:n + 3]
    RR = P[3:n + 4]
    eps = 1e-300
    dRL = R - L
    rp = _ratio(L - LL, dRL, eps)
    rm = _ratio(RR - R, dRL, eps)
    c_from_left = L + 0.5 * _vanleer(rp) * dRL
    c_from_right = R - 0.5 * _vanleer(rm) * dRL
    return np.where(Q >= 0, c_from_left, c_from_right)


def _tvd_face_y(P, Q):
    """Same thing as _tvd_face_x but along y. I kept them as two copies
    instead of one function with an axis argument."""
    n = P.shape[1] - 4
    LL = P[:, 0:n + 1]
    L = P[:, 1:n + 2]
    R = P[:, 2:n + 3]
    RR = P[:, 3:n + 4]
    eps = 1e-300
    dRL = R - L
    rp = _ratio(L - LL, dRL, eps)
    rm = _ratio(RR - R, dRL, eps)
    c_from_left = L + 0.5 * _vanleer(rp) * dRL
    c_from_right = R - 0.5 * _vanleer(rm) * dRL
    return np.where(Q >= 0, c_from_left, c_from_right)


def run(grid, psi, Dt, c0, t_end=600.0, tol=0.05, cfl=0.4, record_every=0.5,
        snapshot_times=(), scheme="vanleer"):
    """March forward until the field is uniform (dev <= tol) or t_end.
    Dt can be a scalar or an (nx, ny) array."""
    nx, ny = grid.nx, grid.ny
    Qx, Qy = face_fluxes(grid, psi)
    Dc = np.broadcast_to(np.asarray(Dt, float), (nx, ny)).copy()

    # Face diffusivities by harmonic mean, then turned into diffusive
    # conductances G = D * A / distance.
    Gx = np.zeros((nx + 1, ny))
    D_hm_x = 2 * Dc[1:] * Dc[:-1] / (Dc[1:] + Dc[:-1] + 1e-300)
    Gx[1:-1] = D_hm_x * grid.Ax[1:-1] / grid.dx

    Gy = np.zeros((nx, ny + 1))
    D_hm_y = 2 * Dc[:, 1:] * Dc[:, :-1] / (Dc[:, 1:] + Dc[:, :-1] + 1e-300)
    Gy[:, 1:-1] = D_hm_y * grid.Ay[:, 1:-1] / grid.dy

    # Walls: zero flux. This also makes sure there is no advective flux
    # through them (psi = 0 on the boundary).
    Qx[0, :] = 0.0
    Qx[-1, :] = 0.0
    Qy[:, 0] = 0.0
    Qy[:, -1] = 0.0

    # Stability limit: for every cell add up the outgoing advective flux and
    # all the conductances, then take the worst cell.
    out_adv = (np.maximum(Qx[1:], 0) + np.maximum(-Qx[:-1], 0)
               + np.maximum(Qy[:, 1:], 0) + np.maximum(-Qy[:, :-1], 0))
    out_diff = Gx[1:] + Gx[:-1] + Gy[:, 1:] + Gy[:, :-1]
    out = out_adv + out_diff
    dt_cells = grid.vol / np.maximum(out, 1e-300)
    dt = cfl * np.min(dt_cells)

    # upwind splits of the face fluxes
    Qxp = np.maximum(Qx, 0)
    Qxm = np.minimum(Qx, 0)
    Qyp = np.maximum(Qy, 0)
    Qym = np.minimum(Qy, 0)

    c = c0.astype(float).copy()
    M0 = (c * grid.vol).sum()
    cbar = M0 / grid.vol.sum()
    nsteps = int(np.ceil(t_end / dt))
    rec = max(1, int(round(record_every / dt)))
    ts = []
    devs = []
    mass_err = 0.0
    snaps = {}
    snap_iter = sorted(snapshot_times)
    t_uni = np.nan

    def diff_flux(c):
        # fresh padded arrays every call, not the most efficient but it is clear
        pad_x = np.zeros((nx + 2, ny))
        pad_y = np.zeros((nx, ny + 2))
        pad_x[1:-1] = c
        pad_y[:, 1:-1] = c
        Dx = -Gx * (pad_x[1:] - pad_x[:-1])
        Dy = -Gy * (pad_y[:, 1:] - pad_y[:, :-1])
        return Dx, Dy

    def rhs_upwind(c):
        pad_x = np.zeros((nx + 2, ny))
        pad_y = np.zeros((nx, ny + 2))
        pad_x[1:-1] = c
        pad_y[:, 1:-1] = c
        adv_x = Qxp * pad_x[:-1] + Qxm * pad_x[1:]
        adv_y = Qyp * pad_y[:, :-1] + Qym * pad_y[:, 1:]
        dif_x = -Gx * (pad_x[1:] - pad_x[:-1])
        dif_y = -Gy * (pad_y[:, 1:] - pad_y[:, :-1])
        Fx = adv_x + dif_x
        Fy = adv_y + dif_y
        net = (Fx[1:] - Fx[:-1]) + (Fy[:, 1:] - Fy[:, :-1])
        return net / grid.vol

    def rhs_tvd(c):
        # zero-gradient ghost cells (2 each side) for the limiter stencil
        Px = np.pad(c, ((2, 2), (0, 0)), mode="edge")
        Py = np.pad(c, ((0, 0), (2, 2)), mode="edge")
        Dx, Dy = diff_flux(c)
        face_x = _tvd_face_x(Px, Qx)
        face_y = _tvd_face_y(Py, Qy)
        Fx = Qx * face_x + Dx
        Fy = Qy * face_y + Dy
        net = (Fx[1:] - Fx[:-1]) + (Fy[:, 1:] - Fy[:, :-1])
        return net / grid.vol

    for n in range(nsteps + 1):
        t = n * dt

        # grab a snapshot if we just passed one of the requested times
        if len(snap_iter) > 0 and t >= snap_iter[0]:
            snaps[snap_iter.pop(0)] = c.copy()

        # record diagnostics every `rec` steps
        if n % rec == 0:
            spread = np.abs(c - cbar)
            dev = spread.max() / cbar
            ts.append(t)
            devs.append(dev)
            mass_now = (c * grid.vol).sum()
            mass_err = max(mass_err, abs(mass_now - M0) / M0)
            if dev <= tol and np.isnan(t_uni):
                t_uni = t
                if len(snap_iter) == 0:
                    break

        if n == nsteps:
            break

        if scheme == "upwind":
            c = c - dt * rhs_upwind(c)
        else:
            # SSP-RK2 (Heun)
            k1 = rhs_tvd(c)
            c1 = c - dt * k1
            k2 = rhs_tvd(c1)
            c = 0.5 * c + 0.5 * (c1 - dt * k2)

    return Result(np.array(ts), np.array(devs), t_uni, mass_err, c, snaps)


SOURCE_RADIUS = 2.5e-3   # m, physical size of an outlet hole / slot cluster


def point_source(grid, x, y, radius=SOURCE_RADIUS):
    """Unit-mass release spread over the cells within a FIXED physical radius
    of (x, y). That way the starting condition is the same on every grid,
    which the grid-convergence study needs."""
    radius = max(radius, 0.75 * max(grid.dx, grid.dy))
    dist = np.sqrt((grid.Xc - x) ** 2 + (grid.Yc - y) ** 2)
    m = dist <= radius
    if not m.any():
        # radius is smaller than a cell, so just use the nearest cell
        i = np.argmin(np.abs(grid.xc - x))
        j = np.argmin(np.abs(grid.yc - y))
        m[i, j] = True
    c = np.zeros((grid.nx, grid.ny))
    c[m] = 1.0 / grid.vol[m].sum()
    return c
