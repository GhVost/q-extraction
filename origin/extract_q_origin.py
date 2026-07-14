# extract_q_origin.py - Q-factor extraction inside OriginLab (2021+, embedded Python)
# No external packages required (pure Python + originpro).
#
# DATA LAYOUT: column A = frequency; every other column = one resonance trace.
#              Column X/Y designations are ignored (imports often mislabel them).
#
# HOW TO RUN:
#   1. Click the data workbook so it is the active window.
#   2. Connectivity -> Open Untitled.py, paste this script, press Run (F5).
#   3. Results are printed and written to a new "QResults" sheet.

# ======================= CONFIGURATION =======================
# 'conductance' : traces are motional conductance G = Re(Y)
#                 -3 dB level = peak / 2        (G is power-like)
# 'magnitude'   : traces are admittance magnitude |Y|
#                 -3 dB level = peak / sqrt(2)
MODE = 'conductance'
# =============================================================

import math
import originpro as op

LEVEL_DIV = 2.0 if MODE == 'conductance' else math.sqrt(2.0)

wks = op.find_sheet('w')
if wks is None:
    raise SystemExit("No active worksheet - click your data workbook first!")

ncols = wks.cols
print(f"Mode: {MODE} (level = peak / {LEVEL_DIV:g})")
print(f"Processing sheet with {ncols} columns...")

longnames = wks.get_labels('L') or [''] * ncols
comments  = wks.get_labels('C') or [''] * ncols


def to_floats(vals):
    out = []
    for v in vals:
        try:
            x = float(v)
        except (TypeError, ValueError):
            x = float('nan')
        out.append(x)
    return out


def cross_left(f, y, ipk, lev):
    for k in range(ipk, 0, -1):
        if y[k-1] < lev <= y[k]:
            return f[k-1] + (f[k]-f[k-1]) * (lev - y[k-1]) / (y[k] - y[k-1])
    return None


def cross_right(f, y, ipk, lev):
    for k in range(ipk, len(y)-1):
        if y[k+1] < lev <= y[k]:
            return f[k] + (f[k+1]-f[k]) * (lev - y[k]) / (y[k+1] - y[k])
    return None


freq_all = to_floats(wks.to_list(0))
fvalid = [v for v in freq_all if math.isfinite(v)]
if len(fvalid) < 5:
    raise SystemExit("Column A does not look like a frequency column.")
fmin, fmax = min(fvalid), max(fvalid)
print(f"Frequency axis: col A, {len(fvalid)} points, {fmin:.6g} .. {fmax:.6g}")

names, f0s, bws, qs = [], [], [], []

for j in range(1, ncols):
    y_all = to_floats(wks.to_list(j))
    f, y = [], []
    for a, b in zip(freq_all, y_all):
        if math.isfinite(a) and math.isfinite(b):
            f.append(a)
            y.append(b)
    if len(f) < 5:
        continue

    if min(y) >= fmin * 0.99 and max(y) <= fmax * 1.01:
        print(f"  col {j+1}: looks like a frequency column - skipped")
        continue

    label = comments[j] or longnames[j] or f"col {j+1}"

    ypk = max(y)
    ipk = y.index(ypk)
    fpk = f[ipk]
    lev = ypk / LEVEL_DIV

    f1 = cross_left(f, y, ipk, lev)
    f2 = cross_right(f, y, ipk, lev)

    names.append(label)
    f0s.append(fpk)
    if f1 is None or f2 is None:
        print(f"  {label}: peak cut off by sweep edge - no Q")
        bws.append(None)
        qs.append(None)
    else:
        bw = f2 - f1
        q = fpk / bw
        bws.append(bw)
        qs.append(q)
        print(f"* {label}:  f0 = {fpk:.10g}   BW = {bw:.6g}   Q = {q:.0f}")

if names:
    res = op.new_sheet('w', 'QResults')
    res.cols = 4
    res.from_list(0, names, lname='Trace')
    res.from_list(1, f0s,   lname='f0',     units='(x units)')
    res.from_list(2, bws,   lname='BW 3dB', units='(x units)')
    res.from_list(3, qs,    lname='Q')
    print(f"\nDone: {len(names)} trace(s) -> QResults sheet.")
else:
    print("No traces processed - check the data layout.")
