# Using BHRuler (`bhruler.py`)

A practical guide to running the Black Hole Ruler on your own source list—what to install, what columns to provide, what comes out, and how to interpret it.

---

## 1) What this is

`bhruler.py` is a single-file, command-line tool that computes a consistent set of black-hole scale quantities from a few observables such as mass, spin, Eddington ratio, and host-galaxy context.

It implements the core of the Black Hole Ruler framework:

- gravitational time and length scales
- spin-aware ISCO radius and frequency
- horizon magnetic field and Blandford–Znajek jet power under two accretion prescriptions
- environmental context via sphere of influence
- a simple TDE capture boundary

It also supports direct ingestion of the **Survey V2** CSV format, including columns such as:

- `Name`
- `Class`
- `Mass_Msun`
- `Spin_a`
- `Spin_Known`
- `Regime`

This lets you use the same single-file tool for both standard BHRuler inputs and the survey table without rewriting your CSV by hand.

---

## 2) Quick start

### Prereqs

- Python 3.8+
- `pandas`

Install:

```bash
pip install pandas
````

### Run with the built-in demo (no input file)

```bash
python bhruler.py --output trio_out.csv
```

This writes a 3-row CSV (`Cygnus X-1`, `Sgr A*`, `M87*`) using both accretion branches by default.

### Run on your own standard BHRuler CSV

```bash
python bhruler.py --input objects.csv --branch both --output derived.csv
```

### Run directly on Survey V2

```bash
python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv
```

### Branch options

* `core` — mass/spin/timing/ISCO/TDE quantities only
* `eta_bridge` — efficiency-bridge accretion prescription
* `adaf` — ADAF/RIAF accretion prescription
* `both` — computes both accretion branches

For **Survey V2**, use `--branch core` unless you have added accretion fields such as `lambda_Edd`, `eta_eff`, and `kappa`.

---

## 3) Input CSV schema

Provide a header row. Extra columns are ignored.

Column matching is tolerant to:

* case differences
* spaces
* underscores
* punctuation

That means `Mass_Msun`, `mass msun`, and `Mass [M_sun]` can all be resolved correctly if they match one of the accepted aliases.

### Standard BHRuler schema

| Column key   | Meaning                     | Units / notes                        |
| ------------ | --------------------------- | ------------------------------------ |
| `name`       | Identifier                  | string                               |
| `M_Msun`     | Black hole mass             | solar masses                         |
| `a_star`     | Dimensionless spin          | approximately `-1` to `+1`           |
| `lambda_Edd` | Eddington ratio             | `L/L_Edd`                            |
| `eta_eff`    | Radiative efficiency        | used by `eta_bridge`                 |
| `kappa`      | ADAF/RIAF constant          | used by `adaf`                       |
| `sigma_kms`  | Stellar velocity dispersion | km/s                                 |
| `Re_kpc`     | Effective radius            | kpc                                  |
| `branch`     | Optional per-row override   | `core`, `eta_bridge`, `adaf`, `both` |

### Accepted aliases

#### Name / identifier

* `name`
* `Name`
* `Object`
* `object`

#### Mass

* `M_Msun`
* `Mass_Msun`
* `M`
* `mass_Msun`
* `Mass [M_sun]`
* `M [M_sun]`
* `M [Msun]`

#### Spin

* `a_star`
* `Spin_a`
* `spin`
* `Spin a_* (prograde)`
* `a_* (assumed/est.)`

#### Eddington ratio

* `lambda_Edd`
* `lambda`
* `L/L_Edd`
* `L/L_Edd (assumed)`

#### Optional survey metadata

* `Class`
* `Spin_Known`
* `Regime`

### Per-row branch override

A `branch` column can override the CLI `--branch` for individual rows.

Accepted values:

* `core`
* `eta_bridge`
* `adaf`
* `both`

---

## 4) Example input files

### Example standard CSV

```csv
name,M_Msun,a_star,lambda_Edd,eta_eff,kappa,sigma_kms,Re_kpc,branch
MyAGN,1.5e8,0.7,0.01,0.1,0.1,220,5.4,both
QuiescentCore,5.0e9,0.5,3e-5,0.1,0.1,300,8.0,adaf
XRB-1,12.4,0.95,0.02,0.1,0.1,,,
```

### Example Survey V2-style CSV

```csv
Name,Class,Mass_Msun,Spin_a,Spin_Known,Regime
M87*,SMBH,6.5e9,0.9,Yes,Jet-dominated
Sgr A*,SMBH,4.3e6,0.5,Estimated,Quiescent
Cygnus X-1,XRB,21.2,0.95,Yes,Accreting
```

---

## 5) What the script computes

### Core ruler (mass & spin only)

* `t_g_s` — gravitational time, `t_g = GM/c^3` (seconds)
* `r_s_km` — Schwarzschild radius, `r_s = 2GM/c^2` (km)
* `f_ISCO_Schw_Hz` — Schwarzschild ISCO orbital frequency
* `fISCO_tg_invariant` — mass-invariant product `f_ISCO_Schw * t_g`
* `r_ISCO_rg` — Kerr ISCO radius in `r_g = GM/c^2` units (prograde)
* `f_ISCO_Kerr_Hz` — Kerr ISCO frequency (prograde)
* `f_ratio` — `f_ISCO_Kerr / f_ISCO_Schw`

The invariant `f_ISCO_Schw * t_g` is computed for every row as a built-in sanity check.

### Environment (optional inputs)

* `rinfl_pc` — sphere of influence, `GM/sigma^2` (pc)
* `rinfl_over_Re` — `r_infl / R_e`
* `t_fb_days` — canonical TDE fallback time (days; solar-type star)

### TDE capture flag

* `tde_possible` — boolean from a logistic capture boundary with tunable midpoint and width

### Accretion + jet outputs

When using `eta_bridge`, `adaf`, or `both`:

* `B_H_G_eta`
* `P_BZ_erg_s_eta`
* `B_H_G_adaf`
* `P_BZ_erg_s_adaf`

If a branch is selected but required inputs are missing or invalid, the script writes:

* `warn_eta_bridge`
* `warn_adaf`

and returns `NaN` for that branch’s derived quantities.

### Survey-compatible output aliases

When a survey-like CSV is detected, the output also includes survey-style aliases alongside the standard BHRuler fields:

* `Name`
* `Mass_Msun`
* `Spin_a`
* `Grav_Time_tg_s`
* `Schwarzschild_Rad_km`
* `Orbital_Period_s`
* `Invariant_fISCO_tg`

This makes it easier to compare the output directly with Survey V2 tables while preserving the canonical BHRuler schema for downstream use.

---

## 6) Accretion prescriptions

### `core`

Computes only the geometry, timing, spin, environment, and TDE quantities.

No accretion inputs are required.

Use this for:

* Survey V2 runs
* mass-and-spin-only catalogs
* quick baseline ruler checks

### `eta_bridge`

Assumes:

`dot m = lambda_Edd / eta_eff`

Horizon field scales as:

`B_H ∝ dot m^(1/2) M^(-1/2)`

Jet power scales as:

`P_BZ ∝ a_*^2 B_H^2 M^2`

Columns end in `_eta`.

### `adaf`

Assumes:

`lambda_Edd = kappa dot m^2`

so that:

`dot m = (lambda_Edd / kappa)^(1/2)`

Columns end in `_adaf`.

### `both`

Runs both accretion prescriptions and writes both sets of columns.

In practice, low-`lambda` SMBHs often motivate the ADAF branch, while brighter thin-disk systems may fit the efficiency bridge better. Running both gives you a simple model-dependence bracket.

---

## 7) TDE boundary controls

The tool tags whether disruptions are possible given `M` and `a_*` using a logistic capture factor `S(M; a_*)` with:

* midpoint mass: `--tde-mcrit`
  default: `3e7` solar masses
* transition width in dex: `--tde-width`
  default: `0.15`

Spin shifts the midpoint modestly, so prograde spin can slightly favor disruption at fixed mass.

The code reports:

`tde_possible = (S < 0.5)`

### Example

Make capture easier and the transition sharper:

```bash
python bhruler.py --input my.csv --tde-mcrit 2.0e7 --tde-width 0.10 --output out.csv
```

---

## 8) Worked examples

### A) Built-in demo

```bash
python bhruler.py --output trio_out.csv
```

### B) Standard CSV with both accretion branches

```bash
python bhruler.py --input objects.csv --branch both --output derived.csv
```

### C) Survey V2 direct ingestion

```bash
python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv
```

### D) Only ADAF branch

```bash
python bhruler.py --input objects.csv --branch adaf --output out_adaf.csv
```

### E) Only efficiency bridge

```bash
python bhruler.py --input objects.csv --branch eta_bridge --output out_eta.csv
```

### F) Per-row branch override

Include a `branch` column in your CSV:

```csv
name,M_Msun,a_star,lambda_Edd,eta_eff,kappa,branch
Target-1,8.0e8,0.9,3e-3,0.1,0.1,adaf
Target-2,5.0e7,0.5,0.02,0.1,0.1,eta_bridge
Target-3,4.3e6,0.3,0.0,,,core
```

Then run:

```bash
python bhruler.py --input objects.csv --output out.csv
```

Each row uses its own branch value if present.

### G) TDE sweep for a single object

Duplicate a row with different masses or spins and compare:

* `tde_possible`
* `t_fb_days`

The fallback time scales approximately as `M^(1/2)` for solar-type stars.

---

## 9) Interpreting key columns

### `f_ratio`

Approximately `1` near zero spin.
Rises above `1` for prograde spin because the ISCO moves inward and the orbital frequency increases.

### `fISCO_tg_invariant`

This is the Schwarzschild mass-invariant benchmark. It should be constant across objects apart from numerical precision.

### `B_H_G_*`

Scales up with accretion rate and down with mass as `M^(-1/2)`. Stellar-mass systems can therefore have very large horizon fields.

### `P_BZ_erg_s_*`

Scales with spin, field strength, and mass. Useful as an order-of-magnitude jet-power bracket rather than a precision prediction.

### `rinfl_pc`, `rinfl_over_Re`

Useful for judging whether the BH sphere of influence is likely to be resolved in host-galaxy measurements.

### `tde_possible`

A heuristic disruption/capture flag. It is useful for quick classification sweeps, not a full relativistic disruption calculation.

---

## 10) Programmatic use (optional)

Although designed as a CLI, the file can also be imported as a Python module.

```python
import pandas as pd
import bhruler as br

row = dict(
    name="Test",
    M_Msun=1e8,
    a_star=0.7,
    lambda_Edd=0.01,
    eta_eff=0.1,
    kappa=0.1,
    sigma_kms=200,
    Re_kpc=4.0
)

out = br.compute_row(row, branch="both", tde_mcrit=3e7, tde_width=0.15)
pd.DataFrame([out]).to_csv("one_row.csv", index=False)
```

Useful functions include:

* `compute_row`
* `t_g_seconds`
* `f_isco_schwarz_hz`
* `f_isco_kerr_hz`
* `B_H_eta`
* `B_H_adaf`
* `P_BZ_erg_s`

---

## 11) Troubleshooting

### Survey CSV fails to load

Use the survey-enabled script and run:

```bash
python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv
```

The updated parser accepts Survey V2 field names such as:

* `Name`
* `Mass_Msun`
* `Spin_a`

### NaNs in `_eta` or `_adaf`

Supply the required parameters:

* `eta_eff` for `eta_bridge`
* `kappa` for `adaf`

Or use:

```bash
python bhruler.py --input your.csv --branch core --output out.csv
```

if you only want the geometry/timing outputs.

### No `rinfl` values

Provide `sigma_kms`. Without it, the environment fields remain empty.

### Spins outside `(-1, 1)`

Values are softly clamped to `±0.9999` for the Kerr ISCO calculations.

### Units look wrong

Check that:

* mass is in solar masses
* velocity dispersion is in km/s
* `Re_kpc` is in kpc

### Windows CSV paths

Quote any path with spaces:

```bash
python bhruler.py --input "C:\path with spaces\in.csv" --output out.csv
```

---

## 12) Physics references embedded in the tool

* scale-invariant baseline: `t_g`, `r_s`, `f_ISCO`, and the invariant `f_ISCO * t_g`
* Kerr ISCO: Bardeen–Press–Teukolsky expressions for `r_ISCO(a_*)` and `f_ISCO(a_*)`
* accretion branches and BZ scaling: MAD-motivated horizon fields and standard `P_BZ` proportionality
* environment and TDEs: sphere of influence `GM/sigma^2` and a tunable logistic capture boundary

---

## 13) Reproducibility knobs

* Physical constants (`G`, `c`, `M_sun`) are defined at the top of the file.
* The BZ normalization is explicit in the source.
* The TDE boundary is controlled by:

  * `--tde-mcrit`
  * `--tde-width`

Adjust these only if you are intentionally testing alternative calibrations or prescriptions.

---

## 14) FAQ

### Which branch should I use?

* Use `core` for survey-style mass/spin catalogs or when accretion inputs are unavailable.
* Use `eta_bridge` for radiatively efficient thin-disk systems.
* Use `adaf` for dim, low-Eddington, LLAGN-like systems.
* Use `both` to bracket model dependence.

### Can I run Survey V2 directly?

Yes:

```bash
python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv
```

### Do I need `sigma_kms` and `Re_kpc`?

Only if you want:

* `rinfl_pc`
* `rinfl_over_Re`

Otherwise, leave them blank.

### Retrograde orbits?

The CLI computes prograde ISCO by default. Extending the CLI to expose retrograde output would require a small code edit where `f_isco_kerr_hz(..., prograde=True)` is called.

### How is `tde_possible` decided?

Via a tunable logistic capture boundary around a critical mass, with a modest spin shift. Adjust:

* `--tde-mcrit`
* `--tde-width`

to explore alternate assumptions.

---

## 15) Citation

If you use this tool or its outputs in a report, repository, or analysis note, cite the accompanying concept note and the code version you actually used.

Example concept-note citation:

> Nagy, S. L. (2025). *A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties.*

Example code reference:

> `bhruler.py` — Black Hole Ruler single-file implementation, survey-enabled update.

---

That’s it. Drop in a CSV, choose a branch, and you’ll get a tidy table of derived scales, spin-aware frequencies, jet-power brackets, and simple environment/TDE diagnostics—consistent across XRBs, SMBHs, and survey-style catalogs.
