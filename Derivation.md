## 1. Governing equation

The bacteria are a passive scalar $`c(\mathbf{x},t)`$ in a turbulent, aerated liquid:

```math
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf{u}\,c) = \nabla\cdot(D_t\nabla c) + S(\mathbf{x},t) - k\,c - \text{(root uptake)}
```

$`S`$ is release from the cartridge (localised: geometry lives here), $`k`$ is die-off, and
root uptake is a sink at the roots. Geometry enters **only** through $`\mathbf{u}`$ and the
location of $`S`$.

## 2. Non-dimensionalise

Scale $`t`$ by the mixing time $`t_{\text{mix}}`$ and the reaction terms by their own time
scales $`t_i`$. Each biological process $`i`$ gets a Damköhler number

```math
\mathrm{Da}_i = \frac{t_{\text{mix}}}{t_i}.
```

If $`\mathrm{Da}_i \ll 1`$ for **every** $`i`$, the transport term relaxes $`c`$ to uniform long
before any biology acts. In that case the outcome depends on $`\int S\,dV`$ (the total released)
and not on *where* $`S`$ is or what $`\mathbf{u}`$ looks like. **Geometry can change the
outcome only through the total release rate, which is the shear question ($`\beta`$) and not the
transport question.**

## 3. Mixing time from bubble-column correlations

Superficial gas velocity:

```math
U_G = \frac{Q}{\pi D^2/4} = \frac{4.17\times10^{-6}}{4.35\times10^{-3}} = 9.6\times10^{-4}\ \text{m/s}
```

Axial dispersion (Deckwer 1974, cgs):

```math
D_L = 0.678\,D^{1.4}U_G^{0.3} = 5.57\times10^{-4}\ \text{m}^2/\text{s}
```

Second correlation, used as a control (Baird & Rice 1975, SI):

```math
D_L = 0.35\,D^{4/3}(gU_G)^{1/3} = 2.31\times10^{-3}\ \text{m}^2/\text{s}
```

The two correlations disagree by 4×, and both are extrapolated below their validated diameter
range. Report that; don't hide it.

```math
t_{\text{mix}} = \kappa\,\frac{H^2}{D_L}
```

**The prefactor $`\kappa`$ (new).** The handoff uses $`\kappa=1`$. For a 1-D closed domain with an
end release, the slowest mode is $`2\cos(\pi x/L)\,e^{-\pi^2 D t/L^2}`$. Requiring
$`2e^{-\pi^2Dt/L^2}\le 0.05`$:

```math
\kappa_{5\%} = \frac{\ln 40}{\pi^2} = 0.374
```

So $`t_{\text{mix}}`$ falls between 0.9 and 9.8 s from the correlations. The 2-D solver gives
longer times (up to roughly 300 s with the pessimistic stagnant layer), so the solver range is
the one to use for the conservative Damköhler numbers.

## 4. The crossover length $`L^*`$

Set $`\mathrm{Da}=1`$ against the fastest biological process (root capture, $`t_{\text{fast}}=1`$ h):

```math
\kappa\frac{L^{*2}}{D_L} = t_{\text{fast}} \quad\Rightarrow\quad L^* = \sqrt{\frac{D_L\,t_{\text{fast}}}{\kappa}}
```

| correlation | $`\kappa=1`$ | $`\kappa=0.374`$ |
|---|---|---|
| Deckwer | 1.42 m | 2.32 m |
| Baird–Rice | 2.88 m | 4.72 m |

This uses the *jar's* $`D_L`$ as a stand-in for a raft channel. A real channel
has its own cross-section and aeration, and any net recirculation velocity $`U`$ adds an
advective route with $`L^*_{\text{adv}} = U\,t_{\text{fast}}`$ that is far longer. So $`L^*`$ is a
**lower bound** on the scale where geometry can start to matter. The defensible statement is
"$`L^* \sim 1`$–$`5`$ m, i.e. raft-channel scale", not "1.2–2.9 m".

## 5. Shear: what the prefactor can and cannot do

Release model: $`R = R_0\left[1+\beta(\tau/\tau_{\text{ref}})^m\right]`$, with $`\tau_{\text{ref}}=\tau_B`$. Then

```math
\frac{C_A}{C_B} = \frac{1+\beta\,(\tau_A/\tau_B)^m}{1+\beta}.
```

For any shear law $`\tau = C\,U^n d^{-k}`$ the constant $`C`$ cancels:

```math
\frac{\tau_A}{\tau_B} = \left(\frac{U_A}{U_B}\right)^n\left(\frac{d_B}{d_A}\right)^k
```

So only the **flow regime around the bead** (the exponent $`n`$) matters:

| regime | law | $`n`$ | valid |
|---|---|---|---|
| Stokes (creeping) | $`\tau_{\max} = 3\mu U/d`$ | 1 | $`Re_d \lesssim 1`$–$`10`$ |
| Laminar boundary layer | $`\tau = 0.664\,\rho U^2/\sqrt{Re}`$ | 3/2 | $`Re_d \gtrsim 20`$–$`50`$ |

With the slip priors used here, Design B beads sit at median $`Re_d \approx 6`$: **98% of
draws fall below 20**, outside the boundary-layer law the handoff applied to them. The
regime-aware composite $`\tau=\max(\tau_{\text{Stokes}},\tau_{BL})`$ raises $`\tau_B`$, which lowers
the predicted A/B ratio and roughly **doubles** the $`\beta`$ needed for each threshold. This is
the most important correction in the repo.

## 6. Dose (well-mixed reservoir, justified by §2)

```math
\frac{dC}{dt} = \frac{m R_\infty}{V T}e^{-t/T} - kC
\;\;\Rightarrow\;\;
C(t) = \frac{m R_\infty}{V T}\,\frac{e^{-t/T}-e^{-kt}}{k-1/T}
```

with $`T=4.58`$ d from Balla's Pa curve and $`R_\infty = 8.1\times10^7`$ CFU/g from B25.
