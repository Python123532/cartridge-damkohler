# Measured (dye) vs predicted (solver) mixing time.
# needs results/t_uniform.csv; adds results/dye_tmix.csv points if they exist.
# dye "A" is drawn against both A models, whichever it lands nearer is a result.
import sys, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


pred = read(RES / "t_uniform.csv")
dye = read(RES / "dye_tmix.csv") if (RES / "dye_tmix.csv").exists() else []
MAP = {"A_planar": "A", "A_axisym": "A", "B_planar": "B"}
LAB = {"A_planar": "A (planar model)", "A_axisym": "A (axisym model)", "B_planar": "B (planar model)"}
C = {"A_planar": "#2a78d6", "A_axisym": "#1baf7a", "B_planar": "#eb6834"}

rows = []
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e5e1", "font.size": 9, "figure.dpi": 150, "savefig.bbox": "tight"})
fig, ax = plt.subplots(figsize=(6, 3.6))
for i, d in enumerate(LAB):
    sw = [float(r["t_uniform"]) for r in pred if r["design"] == d and r["N"] == "56"]
    cen = [float(r["t_uniform"]) for r in pred if r["design"] == d and r["N"] == "112"
           and r["layer"] == "none" and float(r["U0"]) == 0.05 and float(r["Dt"]) == 3e-5][0]
    ax.plot([i, i], [min(sw), max(sw)], color=C[d], lw=6, alpha=0.25, solid_capstyle="round")
    ax.plot(i, cen, "D", color=C[d], ms=7, label="predicted (central)" if i == 0 else None)
    meas = [float(r["t_uniform_s"]) for r in dye if r["design"] == MAP[d] and r["t_uniform_s"] not in ("", "nan")]
    if meas:
        ax.scatter(i + 0.18 + 0.03 * np.arange(len(meas)), meas, color="#2b2b2a", s=18, zorder=3,
                   label="measured (dye replicate)" if i == 0 else None)
    rows.append(dict(model=d, dye_design=MAP[d], pred_central=cen, pred_min=min(sw), pred_max=max(sw),
                     n_meas=len(meas), meas_mean=np.mean(meas) if meas else np.nan,
                     meas_sd=np.std(meas, ddof=1) if len(meas) > 1 else np.nan,
                     within_sweep=(min(sw) <= np.mean(meas) <= max(sw)) if meas else ""))
ax.set_xticks(range(3)); ax.set_xticklabels(LAB.values())
ax.set_yscale("log"); ax.set_ylabel("time to ±5% uniformity (s)")
ttl = "Measured vs predicted mixing time" if dye else "Predicted mixing time (no dye data yet)"
ax.set_title(ttl, loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8)
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "measured_vs_predicted.png"); plt.close(fig)
with open(RES / "measured_vs_predicted.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
