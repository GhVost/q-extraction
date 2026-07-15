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
import re
import originpro as op

PARAM_RE = re.compile(r'^\s*([^=,]+?)\s*=\s*(.+?)\s*$')


def parse_params(label):
    """Pull name=value tokens out of a comma-separated trace label,
    e.g. 'R_l=0.1 Ohm, internal_ring=2.5 um, motional conductance'
    -> {'R_l': '0.1 Ohm', 'internal_ring': '2.5 um'}"""
    params = {}
    for part in label.split(','):
        m = PARAM_RE.match(part)
        if m:
            params[m.group(1)] = m.group(2)
    return params


LEVEL_DIV = 2.0 if MODE == 'conductance' else math.sqrt(2.0)

wks = op.find_sheet('w')
if wks is None:
    raise SystemExit("No active worksheet - click your data workbook first!")

ncols = wks.cols
print(f"Mode: {MODE} (level = peak / {LEVEL_DIV:g})")
print(f"Processing sheet with {ncols} columns...")

longnames = wks.get_labels('L') or [''] * ncols
comments  = wks.get_labels('C') or [''] * ncols
units     = wks.get_labels('U') or [''] * ncols
funit     = units[0] if units and units[0] else '(x units)'


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

names, params_list, f0s, ypks, bws, qs = [], [], [], [], [], []
rms, bases, asyms = [], [], []
yunit = '(y units)'

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
    if yunit == '(y units)' and j < len(units) and units[j]:
        yunit = units[j]

    ypk = max(y)
    ipk = y.index(ypk)
    fpk = f[ipk]
    lev = ypk / LEVEL_DIV

    f1 = cross_left(f, y, ipk, lev)
    f2 = cross_right(f, y, ipk, lev)

    base = sorted(y)[len(y) // 2]          # median = off-resonance floor
    # ponytail: Rm from peak height only (conductance mode); BVD fit lives
    # in standalone/extract_q.py (Origin's Python has no scipy)
    rm = 1.0 / (ypk - base) if MODE == 'conductance' and ypk > base else None

    names.append(label)
    params_list.append(parse_params(label))
    f0s.append(fpk)
    ypks.append(ypk)
    rms.append(rm)
    bases.append(base)
    if f1 is None or f2 is None:
        print(f"  {label}: peak cut off by sweep edge - no Q")
        bws.append(None)
        qs.append(None)
        asyms.append(None)
    else:
        bw = f2 - f1
        q = fpk / bw
        asym = (f2 - fpk) / (fpk - f1) if fpk > f1 else None
        bws.append(bw)
        qs.append(q)
        asyms.append(asym)
        astr = f"{asym:.2f}" if asym is not None else "n/a"
        print(f"* {label}:  f0 = {fpk:.10g}   peak = {ypk:.6g}   BW = {bw:.6g}"
              f"   Q = {q:.0f}   baseline = {base:.6g}   asym = {astr}"
              + (f"   Rm = {rm:.6g}" if rm is not None else ""))

if names:
    param_keys = []
    for p in params_list:
        for k in p:
            if k not in param_keys:
                param_keys.append(k)

    res = op.new_sheet('w', 'QResults')
    res.cols = 8 + len(param_keys)
    res.from_list(0, names, lname='Trace')
    col = 1
    for k in param_keys:
        res.from_list(col, [p.get(k, '') for p in params_list], lname=k)
        col += 1
    res.from_list(col,     f0s,   lname='f0',       units=funit); col += 1
    res.from_list(col,     ypks,  lname='peak',     units=yunit); col += 1
    res.from_list(col,     bws,   lname='BW 3dB',   units=funit); col += 1
    res.from_list(col,     qs,    lname='Q');                     col += 1
    res.from_list(col,     rms,   lname='Rm',       units='Ohm',
                  comments='1/(peak-baseline), conductance mode only'); col += 1
    res.from_list(col,     bases, lname='baseline', units=yunit); col += 1
    res.from_list(col,     asyms, lname='asym',
                  comments='(f2-f0)/(f0-f1), 1 = symmetric')
    print(f"\nDone: {len(names)} trace(s) -> QResults sheet.")
else:
    print("No traces processed - check the data layout.")
