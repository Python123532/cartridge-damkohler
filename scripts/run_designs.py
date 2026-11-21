# Mixing time for every design: grid check at 28/56/112 plus the U0 x Dt sweep.
# python scripts/run_designs.py   (~8 min)
# writes results/t_uniform.csv, results/grid_convergence.csv, figures/t_uniform_designs.png
import sys, csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from p2.designs import DESIGNS, stagnant_layer_D
from p2.fv_solver import run
from p2.gci import gci
from p2.params import U0_SWEEP, DT_SWEEP

RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

GRIDS = (28, 56, 112)
U0_C = 0.05
DT_C = 3e-5
COLORS = {"A_planar": "#2a78d6", "B_planar": "#eb6834", "A_axisym": "#1baf7a"}
LABELS = {"A_planar": "A airlift, planar", "A_axisym": "A airlift, axisymmetric", "B_planar": "B passive, planar"}


def mix_time(name, n, U0, Dt, stagnant=False):
    g, psi, c0 = DESIGNS[name](n, U0)
    D = Dt
    if stagnant:
        D = stagnant_layer_D(g, Dt)
    r = run(g, psi, D, c0, t_end=1500, record_every=0.05)
    assert r.mass_err < 1e-10   # mass should never leak
    return r.t_uniform


def find(rows, name, n, U0, Dt, layer):
    for r in rows:
        if r["design"] == name and r["N"] == n and r["U0"] == U0 and r["Dt"] == Dt and r["layer"] == layer:
            return r["t_uniform"]
    return None


start = time.time()
rows = []

# grid check, central case
for name in DESIGNS:
    for n in GRIDS:
        t = mix_time(name, n, U0_C, DT_C)
        rows.append({"design": name, "N": n, "U0": U0_C, "Dt": DT_C, "layer": "none", "t_uniform": t})
        print(name, n, round(t, 2))

# sweep at N = 56
for name in DESIGNS:
    for U0 in U0_SWEEP:
        for Dt in DT_SWEEP:
            if not (U0 == U0_C and Dt == DT_C):   # already done above
                t = mix_time(name, 56, U0, Dt)
                rows.append({"design": name, "N": 56, "U0": U0, "Dt": Dt, "layer": "none", "t_uniform": t})
            if name.startswith("A"):
                t = mix_time(name, 56, U0, Dt, stagnant=True)
                rows.append({"design": name, "N": 56, "U0": U0, "Dt": Dt, "layer": "stagnant", "t_uniform": t})

with open(RES / "t_uniform.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["design", "N", "U0", "Dt", "layer", "t_uniform"])
    w.writeheader()
    for r in rows:
        r2 = dict(r)
        r2["t_uniform"] = f"{r['t_uniform']:.2f}"
        w.writerow(r2)

# GCI table
gci_rows = []
for name in DESIGNS:
    t28 = find(rows, name, 28, U0_C, DT_C, "none")
    t56 = find(rows, name, 56, U0_C, DT_C, "none")
    t112 = find(rows, name, 112, U0_C, DT_C, "none")
    o = gci(t112, t56, t28)
    gci_rows.append({"design": name, "t28": t28, "t56": t56, "t112": t112, "p_obs": o["p"],
                     "t_extrap": o["f_ext"], "GCI_fine_pct": 100 * o["gci"], "asymptotic": o["asymptotic"]})

with open(RES / "grid_convergence.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(gci_rows[0].keys()))
    w.writeheader()
    for r in gci_rows:
        r2 = {}
        for k, v in r.items():
            r2[k] = f"{v:.3f}" if isinstance(v, float) else v
        w.writerow(r2)

print("\nGRID CONVERGENCE (central case)")
for r in gci_rows:
    print(f"  {r['design']:9s} {r['t28']:6.2f} {r['t56']:6.2f} {r['t112']:6.2f} s  p={r['p_obs']:.2f}  "
          f"t_ext={r['t_extrap']:.2f}  GCI={r['GCI_fine_pct']:.1f}%  asymptotic={r['asymptotic']}")

print("\nSWEEP RANGES (N=56)")
for name in DESIGNS:
    for layer in ("none", "stagnant"):
        v = [r["t_uniform"] for r in rows if r["design"] == name and r["N"] == 56 and r["layer"] == layer]
        if v:
            print(f"  {name:9s} {layer:9s} {min(v):6.1f} - {max(v):6.1f} s")

# figure: grid check on the left, sweep on the right
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "font.size": 9, "figure.dpi": 150,
                     "savefig.bbox": "tight"})
fig, (left, right) = plt.subplots(1, 2, figsize=(9, 3.4))

for r in gci_rows:
    h = [1 / n for n in GRIDS]
    left.plot(h, [r["t28"], r["t56"], r["t112"]], "o-", color=COLORS[r["design"]], lw=1.6, ms=5,
              label=f"{LABELS[r['design']]} (GCI {r['GCI_fine_pct']:.1f}%)")
    left.plot([0], [r["t_extrap"]], "*", color=COLORS[r["design"]], ms=9)
left.set_xlabel("grid spacing h = 1/N")
left.set_ylabel("time to ±5% uniformity (s)")
left.set_xlim(left=-0.002)
left.set_title("Grid convergence, central case (★ = extrapolated)", loc="left", fontsize=10)
left.legend(frameon=False, fontsize=7)

for i, name in enumerate(DESIGNS):
    plain = [r["t_uniform"] for r in rows if r["design"] == name and r["N"] == 56 and r["layer"] == "none"]
    stag = [r["t_uniform"] for r in rows if r["design"] == name and r["N"] == 56 and r["layer"] == "stagnant"]
    right.scatter([i - 0.05] * len(plain), plain, marker="o", s=22, color=COLORS[name], alpha=0.85)
    if stag:
        right.scatter([i + 0.15] * len(stag), stag, marker="^", s=22, color=COLORS[name], alpha=0.85,
                      label="stagnant layer" if i == 0 else None)
right.set_xticks(range(3))
right.set_xticklabels([LABELS[n].replace(", ", "\n") for n in DESIGNS], fontsize=8)
right.set_yscale("log")
right.set_ylabel("time to ±5% uniformity (s)")
right.set_title(r"Sweep: $U_0$ × $D_t$ (9 cases each), N = 56", loc="left", fontsize=10)
right.legend(frameon=False, fontsize=7)

fig.savefig(FIG / "t_uniform_designs.png")
plt.close(fig)
print(f"\nruntime {time.time()-start:.0f} s")
