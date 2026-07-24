#!/usr/bin/env python3
"""
Q-factor extraction from resonance curves in CSV files
(e.g. HFSS / COMSOL / Origin ASCII exports).

Mode is auto-detected per trace: a column of complex literals (trailing
"i", e.g. "4.656e-7+0.0141i") is real Y data -> conductance = Re(Y) is
used directly. A column of plain real numbers is assumed to already be
|Y| -> magnitude mode. --conductance/--magnitude override this for all
traces (needed since a real-only column can't self-identify as G vs |Y|).

Methods (cross-checked against each other):
  1. Half-power bandwidth with interpolated crossings:
       magnitude |Y|      : level = peak / sqrt(2)
       conductance Re(Y)  : level = peak / 2
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
    python extract_q.py data.csv                      # mode auto-detected
    python extract_q.py data.csv --conductance        # force all traces
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


def parse_column(vals):
    """Parse cells that are either plain reals or complex literals with a
    trailing i (e.g. "4.656e-7+0.0141i"). Returns (values, is_complex);
    values is a complex array if any cell parsed as complex, else float."""
    out = []
    is_complex = False
    for v in vals:
        s = str(v).strip()
        if s[-1:] in "ijIJ":
            try:
                out.append(complex(s[:-1] + "j"))
                is_complex = True
                continue
            except ValueError:
                pass
        try:
            out.append(float(s))
        except ValueError:
            out.append(float("nan"))
    if is_complex:
        values = np.array([x if isinstance(x, complex) else complex(x, 0) for x in out])
    else:
        values = np.array(out, dtype=float)
    return values, is_complex


def load_data(path, fcol=None, ycols=None):
    df = pd.read_csv(path, sep=None, engine="python", comment="#", dtype=str)
    if df.shape[1] < 2:
        sys.exit("Need at least two numeric columns (freq, trace).")
    if fcol is None:
        fcol = 0
    f, _ = parse_column(df.iloc[:, fcol])
    f = f.real if np.iscomplexobj(f) else f
    valid = np.isfinite(f)
    if ycols is None:
        ycols = [i for i in range(df.shape[1]) if i != fcol]
    traces = {}
    for i in ycols:
        y, is_complex = parse_column(df.iloc[:, i])
        valid &= np.isfinite(y)
        traces[str(df.columns[i])] = (y, is_complex)
    order = np.argsort(f[valid])
    f = f[valid][order]
    traces = {k: (y[valid][order], is_complex) for k, (y, is_complex) in traces.items()}
    return f, traces


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
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--conductance", action="store_true",
                       help="force all traces as motional conductance Re(Y); level = peak/2")
    mode.add_argument("--magnitude", action="store_true",
                       help="force all traces as admittance magnitude |Y|; level = peak/sqrt(2)")
    ap.add_argument("--fcol", type=int, default=None)
    ap.add_argument("--ycol", type=int, action="append", default=None)
    ap.add_argument("--funit", default="Hz")
    ap.add_argument("--db", action="store_true", help="y data in dB (real-valued traces only)")
    ap.add_argument("--out", default="q_extraction.png")
    args = ap.parse_args()

    scale = UNIT_SCALE[args.funit.lower()]
    f, traces = load_data(args.csv, args.fcol, args.ycol)
    f = f * scale
    fmin, fmax = f.min(), f.max()

    fig, ax = plt.subplots(figsize=(9, 6))
    print(f"{'trace':<45}{'f0 (Hz)':>16}{'peak (S)':>12}{'BW (Hz)':>12}"
          f"{'Q (3dB)':>10}{'Q (BVD)':>10}{'Rm (Ohm)':>12}{'baseline':>12}{'asym':>8}")
    print("-" * 137)

    seen_modes = set()
    for label, (raw, is_complex) in traces.items():
        if is_complex:
            conductance = not args.magnitude
            y = raw.real if conductance else np.abs(raw)
            tag = " [auto: complex Y -> G]" if conductance else " [auto: complex Y -> |Y|]"
        else:
            conductance = args.conductance
            y = raw
            if args.db:
                y = 10 ** (y / 20.0)
            tag = ""
        seen_modes.add(conductance)
        level_div = 2.0 if conductance else np.sqrt(2.0)
        if y.min() >= fmin/scale*0.99 and y.max() <= fmax/scale*1.01:
            print(f"{label:<45}  looks like a frequency column - skipped")
            continue
        hp = q_half_power(f, y, level_div)
        bvd = bvd_fit(f, y, conductance, hp)
        f0 = hp["f0"] if hp else (bvd["f0"] if bvd else np.nan)
        ypk = y[np.argmax(y)]
        baseline = float(np.median(y))
        q3 = f"{hp['Q']:.0f}" if hp else "n/a"
        bw = f"{hp['bw']:.4g}" if hp else "n/a"
        qb = f"{bvd['Q']:.0f}" if bvd else "n/a"
        rm = f"{bvd['Rm']:.4g}" if bvd else "n/a"
        asym = (f"{(hp['f2']-hp['f0'])/(hp['f0']-hp['f1']):.2f}"
                if hp and hp["f0"] > hp["f1"] else "n/a")
        print(f"{(label+tag):<45}{f0:>16.6g}{ypk:>12.4g}{bw:>12}{q3:>10}{qb:>10}"
              f"{rm:>12}{baseline:>12.4g}{asym:>8}")
        if bvd:
            extra = (f"phi = {bvd['phi']:+.3f} rad" if "phi" in bvd
                     else f"C0 = {bvd['C0']:.3g} F")
            print(f"    BVD:  f0 = {bvd['f0']:.10g} Hz   Q = {bvd['Q']:.0f}"
                  f"   Rm = {bvd['Rm']:.4g} Ohm   Lm = {bvd['Lm']:.4g} H"
                  f"   Cm = {bvd['Cm']:.4g} F   {extra}")

        ax.semilogy(f/scale, np.clip(y, 1e-300, None), lw=1, label=f"{label} (Q≈{q3})")
        if hp:
            ax.plot([hp["f1"]/scale, hp["f2"]/scale], [hp["level"]]*2, "r.-", ms=8, lw=1)

    ax.set_xlabel(f"freq ({args.funit})")
    ax.set_ylabel("G (S)" if seen_modes == {True} else "|Y| (S)" if seen_modes == {False} else "G or |Y| (S)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nAnnotated plot: {args.out}")


if __name__ == "__main__":
    main()
