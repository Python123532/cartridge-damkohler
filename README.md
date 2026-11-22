# cartridge-damkohler

Does the shape of a bacteria-releasing cartridge control what reaches the plant in a small
hydroponic jar? A mechanical-engineering answer: compare how fast the jar mixes with how fast
the biology happens.

**Result:** the jar mixes in 6.7–302 s, at least 12× faster than the fastest biological step
($`\mathrm{Da} \le 0.084`$), so cartridge geometry is very unlikely to control dispersal at jar
scale. The one-page summary is in [`docs/RESULT.md`](docs/RESULT.md).

**Status:** the mechanical question is taken as far as it goes. What's left is biological: how
much flow shear changes bead release ($`\beta`$), which needs live bacteria and a flow rig. See
[`docs/OPEN_QUESTION.md`](docs/OPEN_QUESTION.md).

| layer | what | status |
|---|---|---|
| math | $`V_c`$, $`D_L`$, $`t_{\text{mix}}`$, Damköhler numbers, $`L^*`$ | done |
| computation | 2-D solver, designs A and B, 3 grids, tests | done |
| experiment | dye $`t_{\text{mix}}`$ for A and B | predicted; **measurement pending** |
| biology | the open question on $`\beta`$ | written up, not attempted |

## Run

```bash
pip install -r requirements.txt
python -m pytest -q                            # 6 tests: decay, mass, flow direction
python -m p2.timescales              # analytical ladder + Damköhler bracket
python scripts/run_designs.py        # A and B at 28/56/112 + sweep (~8 min)
python scripts/ladder_figure.py      # figure 1
python analysis/dye_analysis.py --self-test
```

After each dye video (protocol in `docs/dye_protocol.md`):

```bash
python analysis/dye_analysis.py --video ... --design A --replicate 1 --calibration ... --t-inject ... --roi ...
python analysis/plot_measured_vs_predicted.py
```

## Files

| path | what |
|---|---|
| `p2/params.py` | all inputs, each tagged [DESIGN] / [ASSUMED] / [LITERATURE] / [PROPERTY] |
| `p2/timescales.py` | correlations, $`\kappa`$, Damköhler numbers, $`L^*`$ |
| `p2/fv_solver.py` | finite-volume transport solver, planar + axisymmetric |
| `p2/designs.py` | flow fields for A and B |
| `p2/gci.py` | grid convergence index |
| `scripts/run_designs.py` | → `results/t_uniform.csv`, `results/grid_convergence.csv` |
| `analysis/` | dye video → $`t_{\text{uniform}}`$; measured vs predicted |
| `docs/` | derivation, pre-registration, dye protocol, open question, result |

## Assumed vs measured

Nothing is measured yet. Jar size, pump setting and source positions are design values from
CAD. Flow speed $`U_0`$ and eddy diffusivity $`D_t`$ are assumed and swept. The biological time
scales come from the literature, for other strains.

## Limitations

The flow is prescribed, not solved; 2-D only; the correlations are extrapolated below their
validated diameter; no roots, chemotaxis or 3-D effects.
