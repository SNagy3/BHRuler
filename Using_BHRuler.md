# Using BHRuler (`bhruler.py`)

A practical guide to running the **Black Hole Ruler** on your own source list—what to install, what columns to provide, what comes out, and how to interpret it.

---

## 1) What this is

`bhruler.py` is a single-file, command-line tool that computes a consistent set of black-hole scale quantities from a few observables (mass, spin, Eddington ratio, etc.). It implements the core of the **Black Hole Ruler (V1)**: gravitational time/length scales, spin-aware ISCO frequency, horizon magnetic field and Blandford–Znajek jet power under two accretion prescriptions, plus environmental context (sphere of influence) and a simple TDE capture boundary.

---

## 2) Quick start

### Prereqs
- Python 3.8+  
- `pandas` (`pip install pandas`)

### Run with the built-in demo (no input file)
```bash
python bhruler.py --output trio_out.csv
```
This writes a 3-row CSV (Cygnus X-1, Sgr A*, M87*) using both accretion branches by default.

### Run on your own CSV
```bash
python bhruler.py --input objects.csv --branch both --output derived.csv
```
- `--branch` can be `eta_bridge`, `adaf`, or `both` (default `both`).

---

## 3) Input CSV schema (case-insensitive)

Provide a header row; extra columns are ignored. Recognized keys:

| Column key | Meaning | Units / notes |
|---|---|---|
| `name` | Identifier | string |
| `M_Msun` (or `M`, `mass_Msun`) | Black hole mass | in solar masses |
| `a_star` (or `spin`) | Dimensionless spin | −0.999… to +0.999… (code clamps internally) |
| `lambda_Edd` (or `lambda`) | Eddington ratio | \(L/L_{\rm Edd}\) |
| `eta_eff` | Radiative efficiency | used by **Efficiency Bridge** branch |
| `kappa` | ADAF/RIAF constant | used by **ADAF** branch via \(\lambda=\kappa\dot m^2\) |
| `sigma_kms` | Stellar velocity dispersion | km/s, for sphere of influence |
| `Re_kpc` | Effective radius | kpc; used with `sigma_kms` |
| `branch` | Optional per-row override | `eta_bridge`/`adaf`/`both` |

> The script accepts common synonyms where shown (e.g., `M` for `M_Msun`, `spin` for `a_star`). Per-row `branch` overrides the CLI `--branch`.

**Example minimal CSV**
```csv
name,M_Msun,a_star,lambda_Edd,eta_eff,kappa,sigma_kms,Re_kpc,branch
MyAGN,1.5e8,0.7,0.01,0.1,0.1,220,5.4,both
QuiescentCore,5.0e9,0.5,3e-5,0.1,0.1,300,8.0,adaf
XRB-1,12.4,0.95,0.02,0.1,0.1,,,
```

---

## 4) What the script computes

**Core ruler (mass & spin only)**
- `t_g_s` — gravitational time \(t_g=GM/c^3\) (seconds).  
- `r_s_km` — Schwarzschild radius \(r_s=2GM/c^2\) (km).  
- `f_ISCO_Schw_Hz` — Schwarzschild ISCO orbital frequency.  
- `r_ISCO_rg` — Kerr ISCO radius in \(r_g=GM/c^2\) units (prograde).  
- `f_ISCO_Kerr_Hz` — Kerr ISCO frequency (prograde).  
- `f_ratio` — \(f_{\rm ISCO,Kerr}/f_{\rm ISCO,Schw}\) (spin shift).  
(Formulas follow the Bardeen–Press–Teukolsky ISCO and the scale-invariant ISCO–\(t_g\) product discussed in the Ruler doc.)

**Environment (optional inputs)**
- `rinfl_pc` — sphere of influence, \(r_{\rm infl}=GM/\sigma^2\) (pc).  
- `rinfl_over_Re` — \(r_{\rm infl}/R_e\) (dimensionless).  
- `t_fb_days` — canonical TDE fallback time (days; solar-type star).

**TDE capture flag**
- `tde_possible` — boolean from a logistic capture boundary with tunable midpoint and width (see §6).

**Accretion + jet (per branch; columns are suffixed)**
- `B_H_G_<suffix>` — horizon-scale magnetic field estimate (Gauss).  
- `P_BZ_erg_s_<suffix>` — Blandford–Znajek jet-power scaling (erg/s).  
Suffix `<suffix>` is `_eta` (Efficiency Bridge) and/or `_adaf` (ADAF/RIAF).

> If a branch is selected but its required parameter is missing or non-positive, you’ll see `warn_eta_bridge` and/or `warn_adaf` and NaNs for that branch’s outputs.

---

## 5) Accretion prescriptions (choose one or both)

1) **Efficiency Bridge**  
   Assumes \( \dot m = \lambda_{\rm Edd} / \eta_{\rm eff} \).  
   Horizon field scales as \(B_H \propto \dot m^{1/2} M^{-1/2}\); jet \(P_{\rm BZ}\propto a_*^2 B_H^2 M^2\). Columns end in `_eta`.

2) **ADAF/RIAF**  
   Assumes \( \lambda_{\rm Edd}=\kappa \dot m^2 \Rightarrow \dot m=(\lambda/\kappa)^{1/2}\).  
   Same field/jet scaffolding; columns end in `_adaf`.

In practice, low-\(\lambda\) SMBHs prefer the ADAF branch; bright disks often fit the bridge. The README’s dual output lets you bracket model dependence.

---

## 6) TDE boundary controls

The tool tags whether disruptions are **possible** given \(M\) and \(a_*\) using a **logistic capture factor** \(S(M;a_*)\) with:
- Midpoint mass: `--tde-mcrit` (default \(3\times10^7\,M_\odot\), Schwarzschild, solar-type star)  
- Transition width (dex): `--tde-width` (default 0.15)  
- Spin shifts the midpoint (prograde spin modestly favors disruption, reducing capture at fixed mass).  
`tde_possible = (S<0.5)`. Tune the midpoint/width to test different stellar types or GR capture prescriptions.

Examples:
```bash
# Make capture easier (lower midpoint), sharper transition
python bhruler.py --input my.csv --tde-mcrit 2.0e7 --tde-width 0.10 --output out.csv
```

---

## 7) Worked examples

### A) Only ADAF branch
```bash
python bhruler.py --input objects.csv --branch adaf --output out_adaf.csv
```

### B) Per-row branch override
Include a `branch` column in your CSV:
```csv
name,M_Msun,a_star,lambda_Edd,eta_eff,kappa,branch
Target-1,8.0e8,0.9,3e-3,0.1,0.1,adaf
Target-2,5.0e7,0.5,0.02,0.1,0.1,eta_bridge
```
Run normally; each row uses its own branch.

### C) TDE sweep for a single object
Duplicate a line with different spins or pass multiple rows; compare `tde_possible` and `t_fb_days`. (The latter scales as \( \propto M^{1/2}\) for solar-type stars.)

---

## 8) Interpreting key columns (rules of thumb)

- **`f_ratio`**: \(\approx 1\) at \(a_*=0\); rises above 1 for prograde spin (higher ISCO frequency) and drops below 1 for retrograde (not computed by default).  
- **`B_H_G_*`**: Scales up with \(\dot m^{1/2}\) and down with \(M^{1/2}\); XRBs show large \(B_H\) despite tiny \(M\).  
- **`P_BZ_erg_s_*`**: \(\propto a_*^2 B_H^2 M^2\); grows with both spin and mass, bracketing observed LLAGN jets in ADAF-like states (e.g., M87*).  
- **`rinfl_pc`, `rinfl_over_Re`**: Useful to gauge whether stellar-dynamical measurements resolve the BH’s sphere of influence.

---

## 9) Programmatic use (optional)

Although designed as a CLI, you can import and call functions:

```python
import pandas as pd
import bhruler as br

row = dict(name="Test", M_Msun=1e8, a_star=0.7, lambda_Edd=0.01, eta_eff=0.1, kappa=0.1, sigma_kms=200, Re_kpc=4.0)
out = br.compute_row(row, branch="both", tde_mcrit=3e7, tde_width=0.15)
pd.DataFrame([out]).to_csv("one_row.csv", index=False)
```
See `compute_row`, `B_H_eta`, `B_H_adaf`, `P_BZ_erg_s`, etc., in the source.

---

## 10) Troubleshooting

- **NaNs in `_eta` or `_adaf`** → Supply the needed parameter (`eta_eff` for bridge; `kappa` for ADAF). The script sets `warn_eta_bridge` / `warn_adaf` when inputs are invalid.  
- **No `rinfl` values** → Provide `sigma_kms`; without it, environment fields are `None`.  
- **Spins outside (−1,1)** → Values are softly clamped to ±0.9999 for Kerr ISCO math.  
- **Units look off** → Check mass in **solar masses**; dispersion in **km/s**; `Re` in **kpc**.  
- **Windows CSV paths** → Quote paths with spaces: `--input "C:\path with spaces\in.csv"`.

---

## 11) Physics references embedded in the tool

- **Scale-invariant baseline**: \(t_g, r_s, f_{\rm ISCO}\) and the mass-invariant \(f_{\rm ISCO}t_g\) check.  
- **Kerr ISCO**: Bardeen–Press–Teukolsky expressions for \(r_{\rm ISCO}(a_*)\) and \(f_{\rm ISCO}(a_*)\) (prograde).  
- **Accretion branches & BZ scaling**: MAD-motivated horizon fields and a standard \(P_{\rm BZ}\) proportionality to bracket jet energetics.  
- **Environment & TDEs**: Sphere of influence \(GM/\sigma^2\); logistic capture boundary tunable by midpoint/width.

---

## 12) Reproducibility knobs

- Physical constants \(G, c, M_\odot\) are defined at the top of the file; adjust only if you need a different constants set.  
- The BZ normalization (e.g., \(10^{45}\) erg/s at chosen fiducials) is explicit—change with caution if you test alternative calibrations.

---

## 13) FAQ

- **Which branch should I use?**  
  If you expect a radiatively efficient thin disk, try **Efficiency Bridge**; if the source is dim/LLAGN-like, **ADAF/RIAF** is more appropriate. Use `both` to bracket.

- **Retrograde orbits?**  
  The CLI computes prograde by default. Extending to retrograde requires a small code edit where `f_isco_kerr_hz(..., prograde=True)` is called.

- **Do I need `sigma_kms` and `Re_kpc`?**  
  Only if you want `rinfl_pc` and `rinfl_over_Re`. Otherwise, leave them blank.

- **How is `tde_possible` decided?**  
  Via a tunable logistic around a critical mass, with a modest spin shift; tweak `--tde-mcrit`/`--tde-width` to explore scenarios.

---

## 14) Citation

If you use this tool or its outputs in a report or repository, please cite the accompanying concept note:

- **Nagy, S. L. (Sep 2025). _A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties_.** (Black Hole Ruler V1).

And reference the code:

- **`bhruler.py` — Black Hole Ruler (V1) single-file implementation.**

---

**That’s it!** Drop in a CSV, choose a branch (or both), and you’ll get a tidy table of derived scales, spin-aware frequencies, jet-power brackets, and simple environment/TDE diagnostics—consistent across XRBs, SMBHs, and beyond.
