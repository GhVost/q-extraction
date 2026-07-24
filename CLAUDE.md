# CLAUDE.md — project context for Claude Code

## What this project is
Q-factor extraction from resonator frequency sweeps (MEMS ring resonator,
f0 ≈ 3.313 GHz, Q ~ 5e4). Data comes from FEM parametric sweeps
(columns = geometry variants: internal_ring / external_ring widths),
analyzed either inside OriginLab or standalone.

## Physics rules (do not change without checking)
- Q = f0 / BW(-3dB), crossings interpolated between samples.
- Crossing level depends on the physical quantity:
  - motional conductance G = Re(Y): level = peak / 2   (G is power-like)
  - admittance magnitude |Y|:       level = peak / sqrt(2)
  Wrong level biases Q by ~1.55x. Default MODE = 'conductance'.
- Mode is auto-detected per trace in standalone/extract_q.py and
  origin/extract_q_origin.py: a column of complex literals with a
  trailing "i" (e.g. "4.656e-7+0.0141i") is real Y data, so G = Re(Y)
  is used directly (unambiguous — a real-only column can't self-report
  whether it's G or |Y|, so MODE/--conductance/--magnitude only matter
  as the fallback/override for those). LabTalk doesn't implement this
  detection (see Environment constraints).
- Q is unitless; frequency units (GHz vs Hz) don't matter.
- Asymmetric (Fano) peaks from feedthrough bias the 3dB method; the fix
  is the BVD fit in standalone/extract_q.py (rotated Lorentzian for G,
  complex feedthrough for |Y|). The asym column (1.00 = symmetric) flags
  when to trust BVD over 3dB. BVD f0 is the one to track for
  softening/hardening.

## Data format quirks (learned the hard way)
- Origin ASCII imports of multi-trace sweeps mislabel data columns as
  X2, X3, ... — column DESIGNATIONS CANNOT BE TRUSTED. Column A is the
  only frequency axis; every other numeric column is a trace.
- Columns whose values fall inside column A's range are duplicate
  frequency axes and must be skipped.
- Trace names live in Origin's Comments row (full parameter combo),
  Long Name only has partial info.
- Typical size: ~19k rows x ~27 columns.
- Trace labels are `name=value` lists (e.g. `R_l=0.1 Ohm, motional
  conductance`); origin/extract_q_origin.py parses these and adds one
  QResults column per distinct parameter name found (handles any number
  of swept parameters, not just one).
- Real units (GHz, S, ...) live in the worksheet's Units row
  (`wks.get_labels('U')`), not in column names — don't hardcode unit
  strings in QResults, read them from there.
- extract_q_origin.py writes QResults into a NEW book named
  "<source book> - QResults" (op.new_book, not new_sheet-into-source),
  so raw sweep data and derived results stay separate but traceable.

## Environment constraints
- origin/extract_q_origin.py runs INSIDE Origin (2025b) via
  Connectivity -> Open Untitled.py. It must stay dependency-free
  (no numpy — the user's Origin Python has no packages installed).
  originpro API only; column designations read via op.lt_exec/lt_int.
- labtalk/extract_q_diag.ogs: LabTalk is extremely fragile in this
  Origin install. Known constraints: no `break`/`continue` in loops,
  no comma-separated declarations (one per line), no C-style int
  declaration inside for(), errors are SILENT. The file keeps staged
  sections t1..t5 + main for bisecting failures. Does NOT auto-detect
  complex-vs-real traces (unlike the two Python front ends) — parsing
  "a+bi" string literals reliably isn't worth the risk in this
  interpreter; it always uses the plain 3dB level from its MODE setting.
- standalone/extract_q.py may use numpy/scipy/pandas/matplotlib freely.
- origin/app/QExtraction/ packages extract_q_origin.py as an Origin App
  (.opx) via Package Manager. Python isn't an officially listed App
  content type (X-Function/.ogs/.c/.h/DLL are); launch.ogs works around
  this with a documented run.python(path, 2) call. Package Manager
  bundles launch.ogs + the icon, NOT an install-time copy of the .py, so
  extract_q_origin.py must be manually placed in the User Files Folder
  (%Y) and re-copied there after edits. Never verified against a real
  Origin install — if the App silently does nothing, check that copy
  step first, then the %Y path expansion.

## Testing
- tests/test_synthetic.py generates Lorentzians with known Q=52000 and
  requires both modes to recover it (tolerance 2%). Run before commits:
      python tests/test_synthetic.py

## Style
- Keep the three front ends behaviorally identical (same level rules,
  same duplicate-axis skipping, same labeling fallbacks).
- Any new analysis method must be cross-checked against the 3dB method
  on synthetic data first.
