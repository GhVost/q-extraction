# q-extraction

Automated Q-factor extraction from resonator frequency sweeps
(FEM/EM simulation output: HFSS, COMSOL, ADS, ...), with three interchangeable
front ends:

| Path | File | Use when |
|---|---|---|
| **OriginLab, embedded Python** | `origin/extract_q_origin.py` | You work in Origin 2021+ (no extra packages needed) |
| **Standalone Python** | `standalone/extract_q.py` | You have CSV exports; adds a Butterworth-Van Dyke fit cross-check |
| **Pure LabTalk** | `labtalk/extract_q_diag.ogs` | Origin older than 2021, or any Origin without working Python |

## Method

Half-power (-3 dB) bandwidth with linear interpolation between samples:

    Q = f0 / (f2 - f1)

where the crossing level depends on what the trace physically is:

- **motional conductance** G = Re(Y): level = **peak / 2**
  (G is already a power-like Lorentzian)
- **admittance magnitude** |Y|: level = **peak / sqrt(2)**

Using the wrong level biases Q by ~1.55x, so the two Python front ends
auto-detect it per trace: a column of complex literals with a trailing
`i` (e.g. `4.656e-7+0.0141i`) is real Y data, so conductance = Re(Y) is
used directly — unambiguous, since that's exactly what conductance mode
already needs. A column of plain real numbers can't self-identify as G
vs |Y|, so it falls back to `MODE = 'conductance' | 'magnitude'` (top of
the Origin script) / `--conductance`/`--magnitude` (standalone script,
mutually exclusive; also force the mode for complex traces if given —
e.g. `--magnitude` on a complex column takes `abs(Y)` instead of
`Re(Y)`). The pure LabTalk front end doesn't implement this detection
(parsing `a+bi` reliably isn't worth the risk in that interpreter) —
it always uses the plain 3dB level from its `MODE` setting.

### Extracted parameters

| Parameter | 3dB method | BVD fit (standalone only) |
|---|---|---|
| f0 (resonance frequency) | peak sample + parabolic refinement | fit — use this one to track softening/hardening |
| peak value | raw maximum | — |
| BW, Q | interpolated -3 dB crossings | fit |
| baseline | median of trace (off-resonance floor) | fit offset |
| asymmetry | (f2-f0)/(f0-f1), 1.00 = symmetric | phi (conductance) / C0 (magnitude) |
| Rm, Lm, Cm | Rm = 1/(peak-baseline), conductance only | fit |

### Butterworth-Van Dyke fit (standalone script)

The standalone script also fits a BVD model near resonance:
motional branch `Ym = (1/Rm) / (1 + 2jQ(f-f0)/f0)`, plus a Fano skew
angle `phi` (conductance mode) or a complex feedthrough `~ j*2*pi*f0*C0`
(magnitude mode). On a symmetric peak it agrees with the 3dB method;
on an asymmetric (Fano) peak the 3dB Q is biased (the `asym` column
flags this — 1.00 means symmetric) and the BVD values are the ones to
trust. For softening/hardening assessment across drive levels, track
the BVD `f0`: it stays exact even when feedthrough shifts the raw peak
position by several kHz.

## Data layout

Column A = frequency; every other numeric column = one trace.
Column X/Y designations in Origin are ignored (multi-trace ASCII imports
often mislabel data columns as X2, X3, ...). Columns whose values fall
inside the frequency range are skipped as duplicate axes. Trace names are
taken from Origin's Comments row (falls back to Long Name).

### Sweep parameter columns (Origin script)

If a trace's Comments/Long Name is a `name=value` list — e.g.
`R_l=0.1 Ohm, motional conductance` — `origin/extract_q_origin.py` splits
it and adds one QResults column per distinct parameter found (`R_l` in
this example). Sweeping more than one parameter (e.g.
`internal_ring=2.5 um, external_ring=4 um, ...`) adds one column each,
so a multi-parameter sweep unstacks into separate, sortable/filterable
columns instead of one combined text label. Tokens without an `=` (like
the trailing `motional conductance`) are ignored. `f0`/`peak`/`BW`
columns get their units from the source worksheet's Units row (falls
back to `(x units)`/`(y units)` if that row is empty).

## Usage

### Origin 2021 and newer (embedded Python)

Origin 2021 introduced the embedded Python environment with the
`originpro` module — that is all `origin/extract_q_origin.py` needs
(deliberately no numpy/scipy, so a stock install works):

1. Activate the data workbook (click it so it is the active window).
2. **Connectivity -> Open Untitled.py**, paste `origin/extract_q_origin.py`.
3. Set `MODE = 'conductance' | 'magnitude'` at the top, press **Run** (F5).
4. Results print to the message log and land in a new book named
   `<source book> - QResults`, with units on f0/BW/peak/baseline copied
   from the source columns.

Tested on Origin 2025b; anything 2021+ should behave the same.

Pasting into Open Untitled.py every time is fine for occasional use but
gets old fast if you run this often. Two less copy-paste-y ways to launch
the same script, both documented by OriginLab:

- **One-click button**: bind a button (worksheet button, or
  **Format -> Toolbars** custom toolbar button) to a LabTalk command that
  runs the saved file directly (the trailing `2` tells `run.python()` its
  argument is a file path, not inline code — see
  [Run(object)](https://docs.originlab.com/labtalk/ref/run-obj/)):

      run.python("C:\path\to\origin\extract_q_origin.py", 2);

  No paste step; you edit the `.py` file on disk and the button always
  runs the current version.
- **Origin App (`origin/app/QExtraction/`)**: a real Apps Gallery icon,
  see below.

#### Building the Origin App

`origin/app/QExtraction/` in this repo has the two files Origin's App
packaging system needs — `launch.ogs` (entry point) and `AppIcon.png`
(32x32 icon) — following the
[Create and Update Apps](https://docs.originlab.com/tutorials/create-update-apps/)
convention (Python itself isn't one of the officially listed App content
types — X-Function/.ogs/.c/.h/DLL are — but `.ogs` calling
`run.python()` is a documented, working way to run a Python script from
one). To build and install it:

1. Copy `origin/extract_q_origin.py` into your Origin **User Files
   Folder** (`%Y` in LabTalk — find it via **Tools -> User Files
   Folder**). `launch.ogs` runs `%Yextract_q_origin.py`, so re-copy it
   there whenever the script changes; Package Manager bundles the
   launcher and icon, not an arbitrary file-copy-on-install step, so the
   `.py` can't ride inside the `.opx` itself.
2. In Origin, open **Code Builder**, find `origin/app/QExtraction/` (or
   copy it into your Apps workspace folder), right-click it, and choose
   **Generate...** to open **Package Manager**. Point it at `launch.ogs`
   as the entry script and `AppIcon.png` as the icon, set a name/version,
   and generate the `.opx`.
3. Double-click the resulting `.opx` to install; the icon appears in the
   Apps Gallery. With a data workbook active, click it to run — same
   output as the manual flow (`<source book> - QResults` book).

This last step (running Package Manager and testing the installed App)
has to happen inside Origin itself — I can't drive that GUI or verify it
end-to-end from here, so treat step 2/3 as the part to sanity-check
first if the icon doesn't behave as expected.

### Origin older than 2021 (no `originpro`)

Pre-2021 Origin has no embedded Python (`Connectivity` menu and
`originpro` do not exist; 2017-2020 only offer the external PyOrigin
bridge, which needs its own Python install). Two options, in order of
preference:

1. **Export + standalone script** (full feature set incl. BVD fit):
   **File -> Export -> ASCII**, comma or tab separated, then run
   `standalone/extract_q.py` on the exported file (see below).
2. **Pure LabTalk** (3dB method only, results in the Script Window):
   open the Script Window (**Window -> Script Window**) and run

       run.section(<path>\labtalk\extract_q_diag.ogs, main);

   with the data worksheet active. If `main` fails silently (old
   LabTalk interpreters reject different constructs), run the staged
   sections `t1`..`t5` one by one — the first section that prints its
   header but not its `== tN OK ==` line names the offending construct.
   Note: the LabTalk path uses the `peak/sqrt(2)` (magnitude) level and
   X/Y column designations; check both match your data.

### Standalone script

    python standalone/extract_q.py data.csv --funit GHz --conductance

**Tests** (synthetic Lorentzian with known Q; both modes must recover it):

    python tests/test_synthetic.py

## Accuracy notes

- Make sure the sweep has at least ~10 points inside the 3 dB bandwidth.
- A visibly asymmetric (Fano-like) peak from capacitive feedthrough biases
  the 3 dB method; watch the `asym` column and trust the BVD-fit values
  (standalone script) when it deviates from 1.00.
- Q is unitless: frequency units (GHz vs Hz) do not affect it.
