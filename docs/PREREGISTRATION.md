# Pre-registration: mixing time, Designs A and B

Predictions come from the 2-D solver (`results/t_uniform.csv`, made by `scripts/run_designs.py`).
They must not be edited after the first dye run. Corrections go in a dated addendum below.

## Predicted time to ±5% uniformity

Central case: $`U_0 = 0.05`$ m/s, $`D_t = 3\times10^{-5}`$ m²/s, finest grid ($`N = 112`$, GCI ≤ 2.3%).

| design | central | sweep ($`U_0 \times D_t`$, N = 56) | with stagnant layer |
|---|---|---|---|
| A, planar model | 24.6 s | 8.0–71 s | 30–302 s |
| A, axisymmetric model | 20.5 s | 6.7–58 s | 21–248 s |
| B, planar model | 22.9 s | 7.5–63 s | — |

## Predictions

1. **Range.** The measured $`t_{\text{uniform}}`$ for both designs falls inside **6.7–302 s**.
2. **No big geometry effect.** The mean measured $`t_{\text{uniform}}`$ for A and B differ by
   **less than 2×**. The model puts them within 20% of each other.
3. **Not predicted:** which design is faster. The model can't resolve it, so neither outcome
   confirms anything.

## Falsifier

A measured $`t_{\text{uniform}}`$ **above 600 s** for either design would give
$`\mathrm{Da} > 0.17`$ against 1 h root capture. In that case the claim that the jar is well
mixed relative to the biology is in trouble and the time-scale argument needs rework.

## Expected bias

The camera sees a depth-averaged image and the analysis averages over windows, so measured
values should come out **at or below** the solver's cell-wise values.

## Addenda

*(dated, after the fact only)*
