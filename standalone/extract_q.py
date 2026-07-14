#!/usr/bin/env python3
"""
Q-factor extraction from resonance curves in CSV files
(e.g. HFSS / COMSOL / Origin ASCII exports).

Methods (cross-checked against each other):
  1. Half-power bandwidth with interpolated crossings:
       magnitude |Y|      : level = peak / sqrt(2)
       conductance Re(Y)  : level = peak / 2        (use --conductance)
     Q = f0 / BW
  2. Lorentzian least-squares fit around the peak.

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
    return {"f0": fpk, "bw": bw, "f1": f1, "f2": f2, "Q": fpk/bw, "level": level}


def lorentz_power(f, A, f0, Q, C):
    return A / (1.0 + (2.0*Q*(f - f0)/f0)**2) + C


def q_fit(f, y, conductance, hp):
    """Fit the power-like quantity: G directly, or |Y|^2."""
    p = y if conductance else y**2
    ipk = int(np.argmax(p))
    f0g = f[ipk]
    Qg, bw = (hp["Q"], hp["bw"]) if hp else (1e4, f0g/1e4)
    m = np.abs(f - f0g) < 6.0*bw
    if m.sum() < 6:
        m = slice(None)
    try:
        popt, pcov = curve_fit(lorentz_power, f[m], p[m],
                               p0=[p[ipk], f0g, Qg, np.median(p)], maxfev=20000)
        return {"f0": popt[1], "Q": abs(popt[2]),
                "Q_err": np.sqrt(np.diag(pcov))[2], "popt": popt}
    except Exception as e:
        print(f"  fit failed: {e}", file=sys.stderr)
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
    print(f"{'trace':<45}{'f0 (Hz)':>16}{'BW (Hz)':>12}{'Q (3dB)':>10}{'Q (fit)':>10}")
    print("-" * 93)

    for label, y in traces.items():
        if args.db:
            y = 10 ** (y / 20.0)
        if y.min() >= fmin/scale*0.99 and y.max() <= fmax/scale*1.01:
            print(f"{label:<45}  looks like a frequency column - skipped")
            continue
        hp = q_half_power(f, y, level_div)
        fit = q_fit(f, y, args.conductance, hp)
        f0 = hp["f0"] if hp else (fit["f0"] if fit else np.nan)
        q3 = f"{hp['Q']:.0f}" if hp else "n/a"
        bw = f"{hp['bw']:.4g}" if hp else "n/a"
        qf = f"{fit['Q']:.0f}" if fit else "n/a"
        print(f"{label:<45}{f0:>16.6g}{bw:>12}{q3:>10}{qf:>10}")

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
