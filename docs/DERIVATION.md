# Why cartridge geometry cannot control jar-scale dispersal

Every step below can be derived by hand. Empirical correlations and assumed priors are marked
**[input]**. Nothing in this file counts as a finding until the dye data exist.

---

## 1. Governing equation and assumptions

```math
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf{u}\,c) = \nabla\cdot(D_t\nabla c) + S(\mathbf{x},t) - k\,c - q(\mathbf{x},c)
```

Here $`c`$ is cells per volume, $`S`$ is release from the cartridge, $`k`$ is first-order die-off,
and $`q`$ is root uptake. The walls have no flux, $`\mathbf{n}\cdot\nabla c = 0`$ and
$`\mathbf{u}\cdot\mathbf{n} = 0`$, and the liquid is incompressible, $`\nabla\cdot\mathbf{u} = 0`$.

**Assumptions:**

1. **Dilute.** Cells don't change the flow; $`c`$ is a passive scalar.
2. **No motility.** No swimming and no chemotaxis. Cells move only with the liquid.
3. **No settling.** Stokes settling of a 1 µm cell with $`\Delta\rho \approx 100`$ kg/m³ is
   $`v_s = \Delta\rho\,g\,d^2/(18\mu) \approx 6\times10^{-8}`$ m/s, which is $`5\times10^{-7}`$ of the circulation velocity
   (`timescales.settling_velocity`). The cell size and density are **[input]**.
4. **Eddy diffusivity.** Turbulence enters only through a scalar $`D_t`$, which is swept over
   $`10^{-5}`$–$`10^{-4}`$ m²/s **[input]**. No turbulence model is used.

**Where geometry enters:** only through $`\mathbf{u}`$ and through *where* $`S`$ is nonzero.

---

## 2. The Damköhler argument as an expansion

Scale $`\mathbf{x} = L\mathbf{x}^*`$, $`t = t_{\text{bio}}\,t^*`$ and $`\mathbf{u} = U\mathbf{u}^*`$.
Scale the reaction terms on the biological time scale, so $`S = (c_{\text{ref}}/t_{\text{bio}})S^*`$
and so on. Multiplying through by $`L^2/D_t`$ gives

```math
\varepsilon\,\partial_{t^*}c + \mathrm{Pe}\,\nabla^*\!\cdot(\mathbf{u}^*c) = \nabla^{*2}c + \varepsilon\left(S^* - k^*c - q^*\right),
\qquad
\varepsilon \equiv \frac{L^2/D_t}{t_{\text{bio}}} = \frac{t_{\text{mix}}}{t_{\text{bio}}},\quad
\mathrm{Pe} = \frac{UL}{D_t}.
```

So $`\varepsilon`$ is the Damköhler number. Expand $`c = c_0 + \varepsilon c_1 + \dots`$

**$`O(1)`$:** $`\nabla^{*2}c_0 = \mathrm{Pe}\,\nabla^*\!\cdot(\mathbf{u}^*c_0)`$, with no-flux walls.
Multiply by $`c_0`$ and integrate. Because $`\nabla\cdot\mathbf{u} = 0`$ and $`\mathbf{u}\cdot\mathbf{n} = 0`$,
the advective term is $`\tfrac{\mathrm{Pe}}{2}\int\nabla\cdot(\mathbf{u}^*c_0^2)\,dV = 0`$. That leaves

```math
\int |\nabla^* c_0|^2\,dV = 0 \;\;\Rightarrow\;\; c_0 = c_0(t^*)\ \text{(uniform in space)}.
```

**$`O(\varepsilon)`$:** $`\nabla^{*2}c_1 - \mathrm{Pe}\,\nabla^*\!\cdot(\mathbf{u}^*c_1) = \partial_{t^*}c_0 - S^* + k^*c_0 + q^*`$.
Integrating over the jar, the left side vanishes by the no-flux conditions. This solvability
condition is

```math
V\frac{dc_0}{dt} = \int_V S\,dV - kVc_0 - \int_V q\,dV.
```

**What this shows:**

- At leading order the jar is a well-mixed reservoir: a single ODE for $`c_0(t)`$ *follows* from
  the expansion; it isn't assumed.
- Geometry enters $`c_0`$ only through $`\int S\,dV`$, the **total** release rate. Where the source
  sits and what $`\mathbf{u}`$ looks like affect only $`c_1`$, which means spatial deviations of
  relative size $`O(\varepsilon)`$.
- So geometry can change the outcome only by changing how much is released. That is the shear
  question (§5), not the transport question.

**Advection can't slow the mixing (a rigorous bound).** Let $`c' = c - \bar c`$ and drop the slow
reaction terms. Then

```math
\frac{d}{dt}\int c'^2\,dV = -2D_t\int|\nabla c'|^2\,dV \le -2D_t\lambda_1\int c'^2\,dV,
```

where the advective term drops out by the same argument as above, and $`\lambda_1`$ is the first
nonzero Neumann eigenvalue of the domain (the Poincaré inequality). The variance therefore decays
at least as fast as $`e^{-2D_t\lambda_1 t}`$ for **any** divergence-free flow.

This bound uses the RMS deviation, not the ±5% maximum-deviation criterion used in the solver.
It uses the eddy diffusivity $`D_t`$, not the correlation $`D_L`$.

### The numbers (from `python -m p2.timescales`, which reads `results/t_uniform.csv`)

The 2-D solver range covers every design, every $`(U_0, D_t)`$ case and the stagnant-layer case:
$`t_{\text{mix}} = 6.7`$–$`302`$ s. Against the biological ranges:

| process | $`\mathrm{Da}`$ | worst-case separation |
|---|---|---|
| root capture (1–48 h) | $`3.8\times10^{-5}`$ – $`8.4\times10^{-2}`$ | 12× |
| die-off (1.3–20 d) | $`3.8\times10^{-6}`$ – $`2.7\times10^{-3}`$ | 372× |
| bead release (2–7 d) | $`1.1\times10^{-5}`$ – $`1.7\times10^{-3}`$ | 573× |

The worst case (12×) is the pessimistic stagnant layer, the slowest flow and the lowest $`D_t`$,
set against 1 h root capture. Even there $`\varepsilon \le 0.084`$. The biological ranges are
**[input]** from other strains.

---

## 3. The mixing-time prefactor $`\kappa`$ (one Fourier mode, 1-D)

This is a **1-D, single-mode estimate**, so read it as an order of magnitude only. For a
closed 1-D domain of length $`L`$ with an instantaneous point release at $`x_0`$:

```math
\frac{c(x,t)}{\bar c} = 1 + 2\sum_{n\ge1}\cos\frac{n\pi x_0}{L}\cos\frac{n\pi x}{L}\,e^{-n^2\pi^2 Dt/L^2}.
```

At late times the $`n = 1`$ mode dominates, and its worst-case amplitude over $`x`$ is
$`2|\cos(\pi x_0/L)|`$. Setting that equal to the tolerance $`\delta = 0.05`$ gives

```math
t_{\text{mix}} = \kappa\frac{L^2}{D},\qquad \kappa = \frac{1}{\pi^2}\ln\frac{2|\cos(\pi x_0/L)|}{\delta}.
```

| release | $`x_0/L`$ | amplitude | $`\kappa`$ |
|---|---|---|---|
| at an end (worst case) | 0 or 1 | 2.00 | 0.374 |
| **airlift outlets, 57 mm in a 73.7 mm column** | 0.773 | 1.51 | **0.346** |

**Caveats:** near $`x_0 = L/2`$ the first mode vanishes and the $`n = 2`$ mode sets
$`t_{\text{mix}}`$. Early-time higher modes are ignored. The real jar is 3-D.

---

## 4. Crossover length $`L^*`$: a method, not a channel prediction

Set the transport time equal to the fastest biological time, so $`\mathrm{Da} = 1`$:

```math
\kappa\frac{L^{*2}}{D_L} = t_{\text{bio}} \quad\Longrightarrow\quad L^* = \sqrt{\frac{D_L\,t_{\text{bio}}}{\kappa}}.
```

**Why advection makes this a lower bound.** If a net flow $`U`$ also carries material along the
system, the transport time is at most $`\min(\kappa L^2/D_L,\; L/U)`$, which is no larger than the
dispersive estimate. A shorter transport time pushes the $`\mathrm{Da} = 1`$ crossover to *larger*
$`L`$, so the dispersion-only $`L^*`$ can only underestimate it.

**What it does NOT give.** Plugging in the *jar's* $`D_L`$ (§6) and $`t_{\text{bio}} = 1`$ h gives
$`L^* = 1.4`$–$`4.9`$ m across the two correlations and three values of $`\kappa`$ (`python -m p2.timescales`). That is an illustration of the method, not a prediction for any real
system: a raft channel has its own cross-section, aeration and $`D_L`$. To make a channel claim,
compute $`D_L`$ for that channel first.

---

## 5. Where the mechanics stops

§2 shows that geometry can only change the outcome through the **total release rate**
$`\int S\,dV`$. In other words, through how bead release responds to flow shear. A common way
to write that is

```math
R(\tau) = R_0\left[1 + \beta\left(\frac{\tau}{\tau_{\text{ref}}}\right)^m\right],
```

where $`\beta`$ is the dimensionless shear sensitivity. Nobody has measured $`\beta`$ for alginate
beads releasing *Bacillus*, and the literature disagrees on whether it is large or near zero.
That's a biological measurement, not a transport calculation; see `docs/OPEN_QUESTION.md`.

---

## 6. Inputs: cited, not derived

| input | source | note |
|---|---|---|
| Circulation velocity $`V_c = 1.31(gDU_G)^{1/3}`$ | Joshi & Sharma 1979 | empirical |
| Axial dispersion $`D_L = 0.678\,D^{1.4}U_G^{0.3}`$ (cgs) | Deckwer et al. 1974 | empirical |
| Axial dispersion $`D_L = 0.35\,D^{4/3}(gU_G)^{1/3}`$ (SI) | Baird & Rice 1975 | control; differs 4.15× from Deckwer |
| Superficial gas velocity $`U_G = Q/(\pi D^2/4) = 9.6\times10^{-4}`$ m/s | jar CAD, pump setting | measured input |

All three correlations were fitted on columns wider than the 74.4 mm jar; the handoff gives
$`\gtrsim 0.1`$ m. **Confirm each paper's exact validated diameter and $`U_G`$ ranges before
printing them**, and state the extrapolation on the board. The V1 dye measurement is the
mitigation.

**The 2-D solver uses prescribed kinematics.** The flow is an assumed stream-function shape scaled
to a peak speed $`U_0`$. It does not solve the Navier–Stokes equations. It tests how source
location and circulation shape affect $`t_{\text{mix}}`$, not what the real flow is.

---

## 7. Limitations (named, not modelled)

- 3-D structure, including azimuthal flow and the real bubble plume.
- Roots: uptake geometry, exudates, biofilm.
- Chemotaxis and motility (assumption 2).
- Turbulence physics: $`D_t`$ is swept, not modelled.
