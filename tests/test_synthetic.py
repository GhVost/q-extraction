#!/usr/bin/env python3
"""Generate synthetic Lorentzian data with known Q and verify extraction.
Run:  python tests/test_synthetic.py
"""
import subprocess
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
SCRIPT = HERE.parent / "standalone" / "extract_q.py"


def make_data(path, conductance, f0=3.312975e9, Q=52000, n=2001):
    f = np.linspace(f0 - 2e5, f0 + 2e5, n)
    power = 1.0 / (1.0 + (2*Q*(f - f0)/f0)**2)
    y = 4e-3 * power if conductance else 4e-3 * np.sqrt(power)
    pd.DataFrame({"freq_GHz": f/1e9, "trace": y}).to_csv(path, index=False)
    return f0, Q


def run(csv, flags):
    out = subprocess.run([sys.executable, str(SCRIPT), str(csv), "--funit", "GHz"] + flags,
                         capture_output=True, text=True, cwd=HERE)
    print(out.stdout)
    line = [l for l in out.stdout.splitlines() if l.startswith("trace") or "3.31" in l][-1]
    q3 = float(line.split()[-2])
    return q3


def check(name, got, want, tol=0.02):
    ok = abs(got - want) / want < tol
    print(f"{name}: extracted Q = {got:.0f}, true Q = {want} -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    ok = True
    f0, Q = make_data(HERE/"mag.csv", conductance=False)
    ok &= check("magnitude   ", run(HERE/"mag.csv", []), Q)
    f0, Q = make_data(HERE/"cond.csv", conductance=True)
    ok &= check("conductance ", run(HERE/"cond.csv", ["--conductance"]), Q)
    sys.exit(0 if ok else 1)
