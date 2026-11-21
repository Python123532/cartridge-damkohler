# Commit order

One file per commit, in this order. Use lowercase filenames.

| # | file | layer |
|---|---|---|
| 1 | `README.md`, `requirements.txt`, `.gitignore` | setup |
| 2 | `docs/DERIVATION.md` | math |
| 3 | `p2/__init__.py`, `p2/params.py` | math |
| 4 | `p2/timescales.py` | math |
| 5 | `p2/fv_solver.py` | computation |
| 6 | `p2/designs.py` | computation |
| 7 | `tests/conftest.py`, `tests/test_verification.py` | computation |
| 8 | `p2/gci.py` | computation |
| 9 | `scripts/run_designs.py` | computation |
| 10 | `results/t_uniform.csv`, `results/grid_convergence.csv`, `figures/t_uniform_designs.png` | computation |
| 11 | `scripts/ladder_figure.py`, `figures/fig1_timescale_ladder.png` | math |
| **12** | **`docs/PREREGISTRATION.md`**, dated **before** any dye run | experiment |
| 13 | `docs/dye_protocol.md` | experiment |
| 14 | `analysis/dye_analysis.py`, `analysis/plot_measured_vs_predicted.py` | experiment |
| 15 | `docs/OPEN_QUESTION.md` | biology |
| 16 | `docs/RESULT.md` | summary |
| 17 | `COMMIT_ORDER.md` | — |

After the dye runs: `results/dye_tmix.csv`, the updated `figures/measured_vs_predicted.png`, and
the measured column in `docs/RESULT.md`.
