# Result (one page)

**Question.** Does the shape of a bacteria-releasing cartridge (an airlift through the bead pack,
**A**, or a passive basket next to an offset air stone, **B**) control how bacteria spread through a
320 mL hydroponic jar?

**Answer (model; dye measurement pending).** Probably not. The jar mixes much faster than
anything biological happens, so where the cartridge releases bacteria stops mattering within
about a minute.

## 1. Analytical ladder

| quantity | value | source |
|---|---|---|
| superficial gas velocity $`U_G`$ | 0.96 mm/s | pump + jar size |
| circulation velocity $`V_c`$ | 0.116 m/s | Joshi & Sharma 1979 (simplified form, unverified) |
| axial dispersion $`D_L`$ | $`2.2\times10^{-3}`$ – $`2.3\times10^{-3}`$ m²/s | Deckwer 1974 / Baird & Rice 1975 (agree within 4%) |
| mixing time, correlations (order of magnitude) | 0.8–2.4 s | $`t_{\text{mix}} = \kappa H^2/D_L`$ |
| mixing time, 2-D solver | **6.7–302 s** | below |
| Damköhler number vs root capture (1–48 h) | **≤ 0.084** | worst-case separation 12× |
| crossover length $`L^*`$ (method, jar $`D_L`$) | 2.8–4.9 m | $`\sqrt{D_L t_{\text{bio}}/\kappa}`$ |
## 2. Verified solver

Time to ±5% uniformity, central case, three grids:

| design | N = 28 | 56 | 112 | GCI |
|---|---|---|---|---|
| A, planar | 24.7 | 24.5 | 24.6 s | 0.4% |
| A, axisymmetric | 20.7 | 20.3 | 20.5 s | 2.3% |
| B, planar | 22.5 | 22.8 | 22.9 s | 0.2% |

The designs differ by less than 20%. The tests check no-flow diffusion against the exact
$`e^{-\pi^2Dt/L^2}`$ decay, mass conservation, and the **flow direction**. Two earlier bugs in
`designs.py` (A flowing down its own riser, and B's plume in the wrong place) passed every mass test,
and only the direction tests caught them.

## 3. Dye measurement

| design | predicted (central) | measured |
|---|---|---|
| A | 20.5–24.6 s | *pending* |
| B | 22.9 s | *pending* |

The falsifier is anything above 600 s.

## 4. Where it stops

Transport isn't the lever at jar scale. The one remaining route is shear-dependent bead release,
$`\beta`$, and the literature is split between a strong effect (Rittmann 1982, $`\tau^{0.58}`$) and
none (Peyton & Characklis 1993, $`r^2 = 0.0004`$). Measuring it needs live bacteria, a flow
rig and approvals, which makes it a biological question. See `docs/OPEN_QUESTION.md`.

**Limitations:** the flow is prescribed, not solved; 2-D, not 3-D; the correlations are
outside their data range (jar diameter, gas rate and $`L/d`$ all below); jar values are from CAD, not measured; there are no
roots.
