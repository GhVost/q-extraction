# q-extraction

Automated Q-factor extraction from resonator frequency sweeps
(FEM/EM simulation output: HFSS, COMSOL, ADS, ...), with three interchangeable
front ends:

| Path | File | Use when |
|---|---|---|
| **OriginLab, embedded Python** | `origin/extract_q_origin.py` | You work in Origin 2021+ (no extra packages needed) |
| **Standalone Python** | `standalone/extract_q.py` | You have CSV exports; adds a Lorentzian-fit cross-check |
| **Pure LabTalk** | `labtalk/extract_q_diag.ogs` | Origin without Python; staged diagnostic sections t1..t5 + main |

## Method

Half-power (-3 dB) bandwidth with linear interpolation between samples:

    Q = f0 / (f2 - f1)

where the crossing level depends on what the trace physically is:

- **motional conductance** G = Re(Y): level = **peak / 2**
  (G is already a power-like Lorentzian)
- **admittance magnitude** |Y|: level = **peak / sqrt(2)**

Using the wrong level biases Q by ~1.55x, so set the mode correctly:
`MODE = 'conductance' | 'magnitude'` at the top of the Origin script,
or the `--conductance` flag for the standalone script.

The standalone script additionally fits a Lorentzian to the power-like
quantity as an independent check; the two Q values should agree closely.

## Data layout

Column A = frequency; every other numeric column = one trace.
Column X/Y designations in Origin are ignored (multi-trace ASCII imports
often mislabel data columns as X2, X3, ...). Columns whose values fall
inside the frequency range are skipped as duplicate axes. Trace names are
taken from Origin's Comments row (falls back to Long Name).

## Usage

**Origin:** activate the data workbook, Connectivity -> Open Untitled.py,
paste `origin/extract_q_origin.py`, Run. Results print and land in a new
`QResults` sheet.

**Standalone:**

    python standalone/extract_q.py data.csv --funit GHz --conductance

**Tests** (synthetic Lorentzian with known Q; both modes must recover it):

    python tests/test_synthetic.py

## Accuracy notes

- Make sure the sweep has at least ~10 points inside the 3 dB bandwidth.
- A visibly asymmetric (Fano-like) peak from capacitive feedthrough biases
  the 3 dB method; fit a Butterworth-Van Dyke model instead in that case.
- Q is unitless: frequency units (GHz vs Hz) do not affect it.
