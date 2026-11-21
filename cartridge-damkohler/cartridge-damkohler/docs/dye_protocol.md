# V1 dye protocol: measured mixing time for Designs A and B

Purpose: measure the time to ±5% uniformity, $`t_{\text{uniform}}`$, for each cartridge
geometry and compare it with the 2-D solver predictions in `PREREGISTRATION.md`. No bacteria,
no SRC and no biosafety cabinet are needed.

**Pre-registration first.** Commit `docs/PREREGISTRATION.md` (predicted mixing times)
before the first injection.

## 1. Equipment

| item | spec |
|---|---|
| Jar | The actual 320 mL DWC jar, filled to the design liquid height (73.7 mm). Mark the fill line. |
| Cartridges | Design A (airlift) and Design B (passive basket + offset stone), **with beads removed** (or blank alginate beads, the same volume) so the flow path matches |
| Air pump | The project pump. Set the flow to 0.25 L/min and check it with a rotameter or a bag-fill timing |
| Light source | A white LED panel (edge-lit, diffuse) behind the jar. A phone flashlight is not acceptable |
| Camera | Phone or webcam on a tripod, about 30 cm from the jar, perpendicular to the panel |
| Dye | Red food dye diluted 1:10 in water as the **stock**. Use one bottle for the whole study |
| Injection | A 1 mL syringe with a blunt needle (or tubing) long enough to reach the release point |
| Calibration | Five spare containers of the same glass and diameter as the jar, or the jar itself refilled |
| Other | Black card around the setup, timer, notebook, thermometer |

## 2. Fixed settings (write the values down and never change them)

- **Lighting:** room lights off. The LED panel is the only source, on the same brightness setting every time, and warmed up for 5 min.
- **Camera:** manual exposure, ISO, white balance and focus, all locked. No HDR, no auto-brightness. Use a frame rate of **30 fps** and the highest resolution that holds 30 fps.
- **Geometry:** tape marks for the tripod feet and the jar base, so the field of view is identical in every run.
- **Exposure check:** in clear water, no pixel in the jar ROI should be saturated (green channel < 240 of 255).
- **Temperature:** record the water temperature. Keep it within ±2 °C across all runs, because temperature changes the viscosity.

## 3. Calibration series (do it on each measurement day)

The camera's intensity is not linear in dye concentration, so calibrate it.

1. Make 6 standards in the jar's own liquid volume: **0, 0.25, 0.5, 1.0, 1.5 and 2.0×** the
   final well-mixed concentration of one injection. With a 0.20 mL injection into 320 mL,
   1.0× = 0.625 µL stock per mL. Make each standard by adding stock with the syringe to 320 mL
   and stirring by hand for 60 s.
2. Photograph each standard in the jar position, with the air on (bubbles affect what the
   camera sees), for 10 s of video or 10 stills.
3. Record `concentration_x, image_or_video_path` in `data/calibration_<date>.csv`, with
   concentration in multiples of the final mixed concentration.
4. The analysis fits mean absorbance against concentration and rejects the calibration if
   R² < 0.98. A peak reading above 2.0× can't be converted; if that happens, make a 3.0× standard.

## 4. Run procedure (one replicate)

1. Fill the jar to the mark with fresh tap water at room temperature. Install the cartridge.
2. Start the air. Wait **3 min** for the flow to settle.
3. Start recording. Record **10 s of clear water** (the blank).
4. Inject **0.20 mL stock** in about 1 s at the release point:
   - Design A: through the cartridge's outlet holes (about 57 mm height).
   - Design B: at the basket slots.
   Say "now" out loud, or tap the jar, to mark the injection frame.
5. Keep recording until the jar looks uniform, **plus 120 s**.
6. Stop recording. Drain, rinse twice, and refill.

**Replicates:** 5 per design, alternating A, B, A, B… so drift in the light or camera doesn't
bias one design. Log every run, including failed ones, in `data/run_log.csv`:
date, design, replicate, temperature, file name, notes.

## 5. Analysis

```bash
python analysis/dye_analysis.py --video data/A_rep1.mp4 --design A --replicate 1 \
    --calibration data/calibration_2026-12-01.csv --t-inject 10.3 --roi 120 860 300 700
```

- `--t-inject` is the injection time in seconds; t = 0 in the output.
- `--roi y0 y1 x0 x1` is the liquid region in pixels. Exclude the meniscus, the cartridge body and the walls' bright edges.
- The script converts intensity to absorbance to concentration (through the calibration), averages a 6×6 grid of windows, and finds the first time after which every window *stays* within ±5% of its final value. It appends the result to `results/dye_tmix.csv`.

Then run `python analysis/plot_measured_vs_predicted.py` to make the comparison figure.

## 6. Known limitations (state them on the board)

- The image is a **projection** through the jar's depth, so it averages along the line of sight. It will read as more uniform than the 3-D field. Treat measured $`t_{\text{uniform}}`$ as a lower bound on the true point-wise value.
- Window averaging also smooths small-scale structure. The solver uses cell-wise maxima.
- Food dye is a passive tracer. Bacteria are close to passive, but beads or alginate in the cartridge can hold dye back temporarily.
