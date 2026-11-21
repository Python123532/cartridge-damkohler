# Time-scale ladder: mixing vs biology on one log axis. The one figure to show first.
# python scripts/ladder_figure.py   (needs results/t_uniform.csv for the solver bar)
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from p2 import timescales as ts
from p2.params import BIO_TIMESCALES

bars = [("mixing: correlations", ts.damkohler_table()[0], "#2a78d6")]
solver = ts.solver_tmix_range()
if solver:
    bars.append(("mixing: 2-D solver", solver, "#2a78d6"))
for name, rng in BIO_TIMESCALES.items():
    bars.append((name, rng, "#8a8984"))

plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e5e1", "font.size": 9, "figure.dpi": 150, "savefig.bbox": "tight"})
fig, ax = plt.subplots(figsize=(7, 2.8))
for i, (label, (lo, hi), color) in enumerate(reversed(bars)):
    ax.plot([lo, hi], [i, i], color=color, lw=7, solid_capstyle="round")
    ax.text(hi * 1.6, i, label, va="center", fontsize=8)
for t, label in [(60, "1 min"), (3600, "1 h"), (86400, "1 d"), (604800, "1 wk")]:
    ax.axvline(t, color="#c3c2b7", lw=0.8, zorder=0)
    ax.text(t * 1.15, -0.55, label, fontsize=7, color="#52514e")
ax.set_xscale("log")
ax.set_yticks([])
ax.set_xlim(0.5, 2e8)
ax.set_ylim(-0.8, len(bars) - 0.5)
ax.set_xlabel("time scale (s)")
worst = min(lo for lo, hi in BIO_TIMESCALES.values()) / (solver[1] if solver else ts.damkohler_table()[0][1])
ax.set_title(f"Jar mixing vs biology: worst-case separation {worst:.0f}x", loc="left", fontsize=10)
(ROOT / "figures").mkdir(exist_ok=True)
fig.savefig(ROOT / "figures" / "fig1_timescale_ladder.png")
print("wrote figures/fig1_timescale_ladder.png")
