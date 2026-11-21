# Dye video -> concentration -> t_uniform. Steps in docs/dye_protocol.md.
# python analysis/dye_analysis.py --video data/A_rep1.mp4 --design A --replicate 1 \
#     --calibration data/calibration.csv --t-inject 10.3 --roi y0 y1 x0 x1
# python analysis/dye_analysis.py --self-test
import sys, csv, argparse
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
RES = ROOT / "results"
TOL = 0.05


def load_frames(path, every=1):
    import imageio.v3 as iio
    path = str(path)
    if path.lower().endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff")):
        return np.asarray(iio.imread(path))[None], 1.0
    frames = np.asarray(iio.imread(path))[::every]
    try:
        fps = iio.immeta(path).get("fps", 30.0)
    except Exception:
        fps = 30.0
    return frames, fps / every


def absorbance(frames, I_blank, channel=1):
    # green channel for red dye. +1 avoids log(0)
    I = frames[..., channel].astype(float) + 1
    return -np.log10(I / (I_blank + 1))


def roi_mean(A, roi):
    y0, y1, x0, x1 = roi
    return A[..., y0:y1, x0:x1].mean(axis=(-1, -2))


def fit_calibration(conc, absb):
    # A = a1 c + a2 c^2. the c^2 term soaks up camera non-linearity
    c = np.asarray(conc, float)
    a = np.asarray(absb, float)
    X = np.column_stack([c, c ** 2])
    coef = np.linalg.lstsq(X, a, rcond=None)[0]
    pred = X @ coef
    r2 = 1 - np.sum((a - pred) ** 2) / np.sum((a - a.mean()) ** 2)
    return coef[0], coef[1], r2


def make_inverse(a1, a2, c_max):
    cg = np.linspace(0, 1.2 * c_max, 2000)
    Ag = a1 * cg + a2 * cg ** 2
    if np.any(np.diff(Ag) <= 0):
        raise ValueError("calibration curve turns over, add standards or check exposure")
    return lambda A: np.interp(A, Ag, cg, left=0.0, right=np.nan)


def window_means(F, roi, n=6):
    y0, y1, x0, x1 = roi
    sub = F[:, y0:y1, x0:x1]
    ys = np.linspace(0, sub.shape[1], n + 1).astype(int)
    xs = np.linspace(0, sub.shape[2], n + 1).astype(int)
    cols = []
    for i in range(n):
        for j in range(n):
            cols.append(sub[:, ys[i]:ys[i + 1], xs[j]:xs[j + 1]].mean(axis=(1, 2)))
    return np.stack(cols, axis=1)


def t_uniform(t, W, tol=TOL, plateau_frac=0.1):
    # first time after which every window STAYS within tol of its final value
    k = max(1, int(len(t) * plateau_frac))
    final = np.nanmean(W[-k:], axis=0)
    dev = np.nanmax(np.abs(W / final - 1), axis=1)
    bad = np.where(~(dev <= tol))[0]
    if len(bad) == 0:
        return float(t[0]), dev
    i = bad[-1] + 1
    if i >= len(t):
        return np.nan, dev
    return float(t[i]), dev


def calibrate(cal_csv, I_blank, roi):
    conc, absb = [], []
    with open(cal_csv) as f:
        for row in csv.DictReader(f):
            p = Path(row["image_or_video_path"])
            if not p.is_absolute():
                p = ROOT / p
            fr, _ = load_frames(p)
            conc.append(float(row["concentration_x"]))
            absb.append(float(roi_mean(absorbance(fr, I_blank), roi).mean()))
    return conc, absb


def analyse(frames, fps, t_inject, roi, cal, n=6, channel=1):
    t_all = np.arange(len(frames)) / fps
    before = t_all < t_inject - 0.2
    if before.sum() < 5:
        raise ValueError("need at least 5 clear-water frames before injection")
    I_blank = frames[before][..., channel].astype(float).mean(axis=0)

    a1, a2, r2 = fit_calibration(cal[0], cal[1])
    if r2 < 0.98:
        raise ValueError(f"calibration R^2 = {r2:.3f}, redo the standards")
    to_conc = make_inverse(a1, a2, max(cal[0]))

    after = t_all >= t_inject
    C = to_conc(absorbance(frames[after], I_blank, channel))
    W = window_means(C, roi, n)
    t = t_all[after] - t_inject
    tu, dev = t_uniform(t, W)
    out_of_range = float(np.isnan(C[:, roi[0]:roi[1], roi[2]:roi[3]]).mean())
    c_final = float(np.nanmean(W[-max(1, len(t) // 10):]))
    return dict(t=t, dev=dev, t_uniform=tu, cal_a1=a1, cal_a2=a2, cal_r2=r2,
                c_final=c_final, frac_out_of_cal=out_of_range)


def write_outputs(res, design, rep, video):
    RES.mkdir(exist_ok=True)
    with open(RES / f"dye_dev_{design}_{rep}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "dev"])
        for t, d in zip(res["t"], res["dev"]):
            w.writerow([round(t, 4), round(d, 5)])
    out = RES / "dye_tmix.csv"
    new = not out.exists()
    with open(out, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["design", "replicate", "video", "t_uniform_s", "cal_r2", "c_final_x", "frac_out_of_cal"])
        w.writerow([design, rep, video, f"{res['t_uniform']:.2f}", f"{res['cal_r2']:.4f}",
                    f"{res['c_final']:.3f}", f"{res['frac_out_of_cal']:.4f}"])


def fake_camera(c, rng=None):
    # Beer-Lambert then a gamma curve, so the calibration has real work to do
    I = 255 * (225 / 255 * 10 ** (-0.55 * c)) ** 0.8
    if rng is not None:
        I = I + rng.normal(0, 1.2, np.shape(I))
    return np.clip(I, 0, 255)


def self_test(design="B_planar", up=4):
    from p2.designs import DESIGNS
    from p2.fv_solver import run
    rng = np.random.default_rng(0)
    g, psi, c0 = DESIGNS[design](28, 0.05)
    fps = 10.0
    t_end = 90.0
    r = run(g, psi, 3e-5, c0, t_end=t_end, tol=-1, record_every=0.05,
            snapshot_times=list(np.arange(0, t_end, 1 / fps)))
    ts = np.array(sorted(r.snapshots))
    cs = np.array([r.snapshots[t] for t in ts])
    cbar = (c0 * g.vol).sum() / g.vol.sum()

    # solver field -> image (row 0 = top of jar), with 20 blank frames first
    field = np.transpose(cs / cbar, (0, 2, 1))[:, ::-1, :]
    field = np.kron(field, np.ones((1, up, up)))
    n_blank = 20
    field = np.concatenate([np.zeros((n_blank,) + field.shape[1:]), field])

    def to_rgb(F):
        return np.stack([np.full_like(F, 225.0), fake_camera(F, rng), np.full_like(F, 225.0)], -1)

    frames = to_rgb(field).astype(np.uint8)
    H, W = field.shape[1:]
    roi = (0, H, 0, W)
    standards = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    I_blank = to_rgb(np.zeros((1, H, W)))[..., 1].astype(float).mean(axis=0)
    cal_abs = []
    for c in standards:
        img = to_rgb(np.full((5, H, W), c)).astype(np.uint8)
        cal_abs.append(float(absorbance(img, I_blank).mean()))

    res = analyse(frames, fps, t_inject=n_blank / fps, roi=roi, cal=(standards, cal_abs), n=4)
    t_true, _ = t_uniform(ts, window_means(field[n_blank:], roi, 4))
    return res, t_true


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--video")
    ap.add_argument("--design")
    ap.add_argument("--replicate", default="1")
    ap.add_argument("--calibration")
    ap.add_argument("--t-inject", type=float)
    ap.add_argument("--roi", type=int, nargs=4)
    ap.add_argument("--windows", type=int, default=6)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        res, t_true = self_test()
        print(f"self-test: recovered t_uniform = {res['t_uniform']:.2f} s, truth = {t_true:.2f} s, "
              f"calibration R^2 = {res['cal_r2']:.4f}")
        return

    frames, fps = load_frames(a.video, a.every)
    before = np.arange(len(frames)) / fps < a.t_inject - 0.2
    I_blank = frames[before][..., 1].astype(float).mean(axis=0)
    cal = calibrate(a.calibration, I_blank, a.roi)
    res = analyse(frames, fps, a.t_inject, a.roi, cal, a.windows)
    write_outputs(res, a.design, a.replicate, a.video)
    print(f"{a.design} rep {a.replicate}: t_uniform = {res['t_uniform']:.1f} s  "
          f"(cal R^2 {res['cal_r2']:.3f}, {100*res['frac_out_of_cal']:.2f}% of pixels above calibration)")


if __name__ == "__main__":
    main()
