#!/usr/bin/env python3
"""Generate synthetic Lorentzian data with known Q and verify extraction.

Cases:
  1. magnitude |Y|, symmetric        -> 3dB and BVD both recover Q
  2. conductance G, symmetric        -> 3dB and BVD both recover Q
  3. conductance G, Fano skew        -> BVD recovers Q (3dB is biased)
  4. magnitude |Y| with feedthrough  -> BVD recovers Q (3dB is biased)
  5. complex Y column, no flags      -> auto-detected as conductance
  6. complex Y column, --magnitude   -> forced to |Y| instead
BVD amplitude 4e-3 S -> Rm = 250 Ohm, checked in the conductance cases.

Run:  python tests/test_synthetic.py
"""
import subprocess
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
SCRIPT = HERE.parent / "standalone" / "extract_q.py"
F0, Q, RM = 3.312975e9, 52000, 250.0


def make_data(path, conductance, phi=0.0, feed=0.0, n=2001):
    f = np.linspace(F0 - 2e5, F0 + 2e5, n)
    z = (1.0 / RM) / (1.0 + 2j * Q * (f - F0) / F0)
    y = (z * np.exp(1j * phi)).real if conductance else np.abs(z + 1j * feed)
    pd.DataFrame({"freq_GHz": f / 1e9, "trace": y}).to_csv(path, index=False)


def make_complex_data(path, n=2001):
    """Full complex Y, written as 'a+bi' string literals like a real FEM export."""
    f = np.linspace(F0 - 2e5, F0 + 2e5, n)
    z = (1.0 / RM) / (1.0 + 2j * Q * (f - F0) / F0)
    cells = [f"{v.real:.17g}{v.imag:+.17g}i" for v in z]
    pd.DataFrame({"freq_GHz": f / 1e9, "trace": cells}).to_csv(path, index=False)


def run(csv, flags):
    out = subprocess.run([sys.executable, str(SCRIPT), str(csv), "--funit", "GHz"] + flags,
                         capture_output=True, text=True, cwd=HERE)
    print(out.stdout)
    r = {}
    for line in out.stdout.splitlines():
        tok = line.split()
        if len(tok) >= 9 and tok[0] == "trace" and tok[1] != "f0":
            # trace f0 peak BW Q3 QBVD Rm baseline asym
            r["Q3"] = tok[4]
        elif tok[:1] == ["BVD:"]:
            # BVD: f0 = v Hz Q = v Rm = v Ohm Lm = ...
            r.update(f0=float(tok[3]), QBVD=float(tok[7]), Rm=float(tok[10]))
    if "QBVD" not in r:
        sys.exit(f"could not parse output of {csv}:\n{out.stdout}\n{out.stderr}")
    return r


def check(name, got, want, tol=0.02):
    ok = abs(got - want) / want < tol
    print(f"{name}: got {got:.6g}, want {want:.6g} -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    ok = True

    make_data(HERE / "mag.csv", conductance=False)
    r = run(HERE / "mag.csv", [])
    ok &= check("magnitude   Q(3dB)", float(r["Q3"]), Q)
    ok &= check("magnitude   Q(BVD)", r["QBVD"], Q)

    make_data(HERE / "cond.csv", conductance=True)
    r = run(HERE / "cond.csv", ["--conductance"])
    ok &= check("conductance Q(3dB)", float(r["Q3"]), Q)
    ok &= check("conductance Q(BVD)", r["QBVD"], Q)
    ok &= check("conductance Rm    ", r["Rm"], RM)

    make_data(HERE / "fano.csv", conductance=True, phi=0.6)
    r = run(HERE / "fano.csv", ["--conductance"])
    ok &= check("Fano phi=.6 Q(BVD)", r["QBVD"], Q)
    ok &= check("Fano phi=.6 f0    ", r["f0"], F0, tol=1e-6)
    ok &= check("Fano phi=.6 Rm    ", r["Rm"], RM)

    make_data(HERE / "feed.csv", conductance=False, feed=1e-3)
    r = run(HERE / "feed.csv", [])
    ok &= check("feedthrough Q(BVD)", r["QBVD"], Q)

    make_complex_data(HERE / "complex.csv")
    r = run(HERE / "complex.csv", [])
    ok &= check("complex auto->G Q(BVD)", r["QBVD"], Q)
    ok &= check("complex auto->G Rm    ", r["Rm"], RM)

    r = run(HERE / "complex.csv", ["--magnitude"])
    ok &= check("complex ->|Y| Q(BVD)  ", r["QBVD"], Q)
    ok &= check("complex ->|Y| Rm      ", r["Rm"], RM)

    sys.exit(0 if ok else 1)
