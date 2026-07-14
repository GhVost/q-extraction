#!/usr/bin/env python3
"""
Q-factor extraction from resonance curves in CSV files
(e.g. HFSS / COMSOL / Origin ASCII exports).

Methods (cross-checked against each other):
  1. Half-power bandwidth with interpolated crossings:
       magnitude |Y|      : level = peak / sqrt(2)
       conductance Re(Y)  : level = peak / 2        (use --conductance)
     Q = f0 / BW
  2. Butterworth-Van Dyke (BVD) least-squares fit around the peak:
       motional branch  Ym = (1/Rm) / (1 + 2jQ(f-f0)/f0)
       conductance mode : y = Re(e^{j phi} * Ym) + C   (phi = Fano skew)
       magnitude mode   : y = |Ym + B|, B = complex feedthrough ~ j*2pi*f0*C0
     Robust against asymmetric (Fano) peaks that bias the 3dB method;
     also yields Rm, Lm, Cm (and C0 / phi). f0 from this fit is the one
     to track for softening/hardening assessment.

Per-trace extras: peak height, baseline (median), asymmetry ratio
(f2-f0)/(f0-f1) (1.00 = symmetric Lorentzian).

Usage:
    python extract_q.py data.csv                      # |Y| magnitude data
    python extract_q.py data.csv --conductance        # G = Re(Y) data
    python extract_q.py data.csv --funit GHz --db
    python extract_q.py data.csv --fcol 0 --ycol 1 --ycol 3

Column A is treated as the frequency axis; all other numeric columns
are traces. Columns whose values fall inside the frequency range are
skipped as duplicate axes.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

UNIT_SCALE = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}


def load_data(path, fcol=None, ycols=None):
    df = pd.read_csv(path, sep=None, engine="python", comment="#")
    df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").dropna()
    if df.shape[1] < 2:
        sys.exit("Need at least two numeric columns (freq, trace).")
    if fcol is None:
        fcol = 0
    f = df.iloc[:, fcol].to_numpy(dtype=float)
    if ycols is None:
        ycols = [i for i in range(df.shape[1]) if i != fcol]
    traces = {str(df.columns[i]): df.iloc[:, i].to_numpy(dtype=float) for i in ycols}
    order = np.argsort(f)
    return f[order], {k: v[order] for k, v in traces.items()}


def crossings(f, y, level):
    s = y - level
    idx = np.nonzero(np.diff(np.signbit(s)))[0]
    return np.array([f[i] + (f[i+1]-f[i]) * (-s[i]) / (s[i+1]-s[i]) for i in idx])


def q_half_power(f, y, level_div):
    ipk = int(np.argmax(y))
    fpk, ypk = f[ipk], y[ipk]
    level = ypk / level_div
    c = crossings(f, y, level)
    left, right = c[c < fpk], c[c > fpk]
    if len(left) == 0 or len(right) == 0:
        return None
    f1, f2 = left.max(), right.min()
    bw = f2 - f1
    if 0 < ipk < len(f) - 1:  # parabolic refinement of f0 on log scale
        x = f[ipk-1:ipk+2]
        yl = np.log(np.clip(y[ipk-1:ipk+2], 1e-300, None))
        d = yl[0] - 2*yl[1] + yl[2]
        if d != 0:
            fpk = x[1] - 0.25*(x[2]-x[0])*(yl[2]-yl[0])/d
    return {"f0": fpk, "ypk": ypk, "bw": bw, "f1": f1, "f2": f2, "Q": fpk/bw, "level": level}


def bvd_fit(f, y, conductance, hp):
    """Butterworth-Van Dyke fit near resonance.

    Motional branch Ym = (1/Rm) / (1 + 2jQ(f-f0)/f0); the narrow band
    (BW/f0 ~ 1e-5) lets the feedthrough be a complex constant.
      conductance: y = Re(e^{j phi} Ym) + C     -> phi captures Fano skew
      magnitude:   y = |Ym + br + j bi|         -> bi ~ 2*pi*f0*C0
    With phi=0 / b=0 this reduces to the plain Lorentzian fit.
    """
    ipk = int(np.argmax(y))
    f0g = f[ipk]
    Qg, bw = (hp["Q"], hp["bw"]) if hp else (1e4, f0g/1e4)
    m = np.abs(f - f0g) < 6.0*bw
    if m.sum() < 8:
        m = slice(None)
    base = float(np.median(y))
    A0 = y[ipk] - base

    if conductance:
        def model(f, A, f0, Q, phi, C):
            z = A * np.exp(1j*phi) / (1.0 + 2j*Q*(f - f0)/f0)
            return z.real + C
        p0 = [A0, f0g, Qg, 0.0, base]
    else:
        def model(f, A, f0, Q, br, bi):
            z = A / (1.0 + 2j*Q*(f - f0)/f0)
            return np.abs(z + br + 1j*bi)
        p0 = [A0, f0g, Qg, base, 0.0]

    try:
        popt, pcov = curve_fit(model, f[m], y[m], p0=p0, maxfev=20000)
        A, f0, Q = abs(popt[0]), popt[1], abs(popt[2])
        Rm = 1.0 / A
        out = {"f0": f0, "Q": Q, "Rm": Rm,
               "Lm": Q*Rm/(2*np.pi*f0), "Cm": 1.0/(2*np.pi*f0*Q*Rm),
               "Q_err": np.sqrt(np.diag(pcov))[2]}
        if conductance:
            out["phi"] = popt[3]
        else:
            out["C0"] = popt[4] / (2*np.pi*f0)
        return out
    except Exception as e:
        print(f"  BVD fit failed: {e}", file=sys.stderr)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--conductance", action="store_true",
                    help="traces are motional conductance Re(Y); level = peak/2")
    ap.add_argument("--fcol", type=int, default=None)
    ap.add_argument("--ycol", type=int, action="append", default=None)
    ap.add_argument("--funit", default="Hz")
    ap.add_argument("--db", action="store_true", help="y data in dB")
    ap.add_argument("--out", default="q_extraction.png")
    args = ap.parse_args()

    level_div = 2.0 if args.conductance else np.sqrt(2.0)
    scale = UNIT_SCALE[args.funit.lower()]
    f, traces = load_data(args.csv, args.fcol, args.ycol)
    f = f * scale
    fmin, fmax = f.min(), f.max()

    fig, ax = plt.subplots(figsize=(9, 6))
    print(f"{'trace':<45}{'f0 (Hz)':>16}{'peak (S)':>12}{'BW (Hz)':>12}"
          f"{'Q (3dB)':>10}{'Q (BVD)':>10}{'Rm (Ohm)':>12}{'baseline':>12}{'asym':>8}")
    print("-" * 137)

    for label, y in traces.items():
        if args.db:
            y = 10 ** (y / 20.0)
        if y.min() >= fmin/scale*0.99 and y.max() <= fmax/scale*1.01:
            print(f"{label:<45}  looks like a frequency column - skipped")
            continue
        hp = q_half_power(f, y, level_div)
        bvd = bvd_fit(f, y, args.conductance, hp)
        f0 = hp["f0"] if hp else (bvd["f0"] if bvd else np.nan)
        ypk = y[np.argmax(y)]
        baseline = float(np.median(y))
        q3 = f"{hp['Q']:.0f}" if hp else "n/a"
        bw = f"{hp['bw']:.4g}" if hp else "n/a"
        qb = f"{bvd['Q']:.0f}" if bvd else "n/a"
        rm = f"{bvd['Rm']:.4g}" if bvd else "n/a"
        asym = (f"{(hp['f2']-hp['f0'])/(hp['f0']-hp['f1']):.2f}"
                if hp and hp["f0"] > hp["f1"] else "n/a")
        print(f"{label:<45}{f0:>16.6g}{ypk:>12.4g}{bw:>12}{q3:>10}{qb:>10}"
              f"{rm:>12}{baseline:>12.4g}{asym:>8}")
        if bvd:
            extra = (f"phi = {bvd['phi']:+.3f} rad" if args.conductance
                     else f"C0 = {bvd['C0']:.3g} F")
            print(f"    BVD:  f0 = {bvd['f0']:.10g} Hz   Q = {bvd['Q']:.0f}"
                  f"   Rm = {bvd['Rm']:.4g} Ohm   Lm = {bvd['Lm']:.4g} H"
                  f"   Cm = {bvd['Cm']:.4g} F   {extra}")

        ax.semilogy(f/scale, np.clip(y, 1e-300, None), lw=1, label=f"{label} (Q≈{q3})")
        if hp:
            ax.plot([hp["f1"]/scale, hp["f2"]/scale], [hp["level"]]*2, "r.-", ms=8, lw=1)

    ax.set_xlabel(f"freq ({args.funit})")
    ax.set_ylabel("G (S)" if args.conductance else "|Y| (S)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nAnnotated plot: {args.out}")


if __name__ == "__main__":
    main()
