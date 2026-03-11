# Black Hole Ruler
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18894798.svg)](https://doi.org/10.5281/zenodo.18894798)

A scale-invariant framework and reproducible dataset for comparing black holes across the mass spectrum (Nagy 2025).

**Stephen L. Nagy · Independent Researcher · September 2025 - Current**

*A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties*  
**License:** MIT

BHRuler implements the **Black Hole Ruler** framework for **scale-invariant**, **reproducible** comparison of black holes from stellar to ultra-massive regimes. It provides code and data to compute gravitational units and ISCO landmarks, spin-aware Kerr corrections, dual accretion prescriptions (η-bridge and ADAF/RIAF), Blandford–Znajek jet-power estimates, environmental metrics (σ, r_infl, r_infl/R_e), and a TDE module with a logistic capture boundary—plus a versioned catalog schema for cross-scale studies.

The single-file CLI, `bhruler.py`, now also supports direct ingestion of the **Survey V2** CSV format, including fields such as `Name`, `Class`, `Mass_Msun`, `Spin_a`, `Spin_Known`, and `Regime`, while preserving the standard BHRuler schema for pipeline use.

## What's New in V2.1.1-survey

- **Native Survey V2 CSV support** in `bhruler.py`
- **New `core` branch** for mass/spin/timing/ISCO/TDE calculations without requiring accretion inputs
- **Flexible CSV parsing** across case, spaces, underscores, and punctuation
- **Survey-compatible output aliases** written alongside canonical BHRuler columns:
  - `Name`
  - `Mass_Msun`
  - `Spin_a`
  - `Grav_Time_tg_s`
  - `Schwarzschild_Rad_km`
  - `Orbital_Period_s`
  - `Invariant_fISCO_tg`
- **Survey metadata pass-through** for `Class`, `Spin_Known`, and `Regime`
- **Direct Survey V2 workflow**:
  ```bash
  python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv
  ```

## What's New in V2

* **Five publication-quality figures** generated from computed data (cross-scale invariant, spin diagnostics, M87* sensitivity, 10-object atlas, worked uncertainty propagation)
* **Expanded reference list** (45 citations anchoring every equation to its theoretical source)
* **Worked uncertainty example** for M87* with both analytic error propagation and a 50,000-draw Monte Carlo, showing κ dominates the P_BZ error budget at 41% of variance
* **Strengthened case studies** with direct comparisons to observed jet powers (Cyg X-1), flare timescales (Sgr A*), and EHT constraints (M87*)
* **Complete equation-to-column appendix** mapping every `bhruler.py` output to its paper equation
* **Physical constants appendix** with CODATA/IAU sources
* Removed unimplemented TDE rate equation from the paper body (retained as future-work item)
* TDE spin coefficient explicitly framed as tunable heuristic with guidance to calibrate against Kesden (2012) and Servin & Kesden (2017)
* B_H normalization anchored to Tchekhovskoy+ (2011) and McKinney+ (2012) GRMHD MAD simulations

---

## Table of Contents

* [Why BHRuler?](#why-bhruler)
* [Core Equations](#core-equations)
* [Features](#features)
* [Repository Layout](#repository-layout)
* [Quick Start](#quick-start)
* [Data Products](#data-products)
* [Schema](#schema)
* [Reproducing Results](#reproducing-results)
* [Quality Control & Provenance](#quality-control--provenance)
* [Contributing](#contributing)
* [Citation](#citation)
* [License](#license)

---

## Why BHRuler?

* **Scale-invariant:** Everything is expressed in gravitational units so times, lengths, and frequencies map across 10+ orders of mass.
* **Universal invariant:** `f_ISCO × t_g = 1/(6^{3/2} 2π) ≈ 0.01083` provides a built-in sanity check, verified to machine precision across the mass scale.
* **Model-transparent:** Dual accretion prescriptions (η-bridge vs ADAF/RIAF) make **model dependence explicit** for jet-power predictions.
* **Cross-scale validation:** Unifies X-ray binaries, AGN, EHT imaging, and GW remnants on the same ruler.
* **Survey-ready:** The same single-file code now accepts both standard BHRuler inputs and Survey V2 tables directly.
* **Uncertainty-aware:** Analytic error propagation formulas plus a worked Monte Carlo example for M87*.
* **Reproducible:** Versioned catalog schema with per-row provenance and QC flags.

---

## Core Equations

**Gravitational units**

```
t_g = GM/c³,    r_g = GM/c²,    r_s = 2r_g
```

**Schwarzschild ISCO**

```
r_ISCO = 6 r_g
f_ISCO = c³ / (6^{3/2} 2π GM)
f_ISCO × t_g = 1 / (6^{3/2} 2π)    [mass-invariant]
```

**Kerr ISCO (equatorial, prograde/retrograde)**
Using Bardeen–Press–Teukolsky (1972):

```
f_ISCO(a*) = (c³ / 2π GM) × 1 / (r_ISCO^{3/2}(a*) ± a*)
```

**Accretion prescriptions**

* **η-bridge:** ṁ = λ_Edd / η_eff, B_H ∝ ṁ^{1/2} M^{-1/2}, P_BZ ∝ a*² (λ_Edd/η_eff) M
* **ADAF/RIAF:** λ_Edd = κ ṁ² ⟹ ṁ = (λ_Edd/κ)^{1/2}, B_H ∝ (λ_Edd/κ)^{1/4} M^{-1/2}, P_BZ ∝ a*² (λ_Edd/κ)^{1/2} M

**Blandford–Znajek (order-of-magnitude; normalization from Tchekhovskoy+ 2011 GRMHD MAD simulations)**

```
P_BZ ≈ 10⁴⁵ erg/s × (a*/0.9)² × (B_H/10⁴ G)² × (M/10⁹ M☉)²
```

**Environment & TDE**

```
r_infl = GM/σ²
t_fb ≈ 41 d × (M/10⁶ M☉)^{1/2} × (R*/R☉)^{3/2} × (M*/M☉)⁻¹
```

TDE capture boundary: logistic `S(M; a*)` with tunable midpoint and width.

---

## Features

* **Ruler core:** `t_g`, `r_s`, `f_ISCO`, plus invariant check
* **Spin-aware:** `r_ISCO(a*)`, `f_ISCO(a*)` via BPT
* **Accretion toggle:** `core`, η-bridge, ADAF/RIAF, or both
* **Jet power:** `B_H` estimates and `P_BZ` ranges (MAD-motivated scaling)
* **Environment:** `σ → r_infl`, `r_infl/R_e`
* **TDE:** disrupt-vs-swallow flag plus `t_fb`
* **Survey support:** direct ingestion of Survey V2 CSVs without manual column renaming
* **Flexible parser:** handles alternate naming conventions across cases, spaces, underscores, and punctuation
* **Schema + QC:** versioned, with provenance fields and automated checks

---

## Repository Layout

* `bhruler.py` — single-file CLI implementation of the Black Hole Ruler
* `Using_BHRuler.md` — full user guide, schema notes, examples, and troubleshooting
* `CHANGELOG.md` — version history and survey integration notes
* `Paper/` — manuscript source and figure-generation scripts
* `data/` — reproducible catalog products and supporting tables

---

## Quick Start

**Requirements:** Python ≥ 3.8, `pandas`

```bash
pip install pandas
```

### Run with built-in demo (no input file)

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

For survey-style catalogs that contain mass and spin but not accretion inputs, `--branch core` is the recommended mode.

See [Using_BHRuler.md](Using_BHRuler.md) for the full user guide, accepted aliases, worked examples, and troubleshooting.

### Programmatic use

```python
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
```

---

## Data Products

| File                              | Description                                                                         |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| `bh_ruler_check_2025.csv`         | Cross-scale invariant check: Gaia BH3, GW231123, ω Cen IMBH, CEERS-1019, J0529-4351 |
| `bh_spin_aware_ISCO_2025.csv`     | Spin-corrected ISCOs: Cyg X-1, LMC X-1, GW150914 remnant                            |
| `bh_trio_spin_field_jet_2025.csv` | Trio (Cyg X-1, Sgr A*, M87*) with `B_H` and `P_BZ` under dual branches              |
| `M87_RIAF_sensitivity_2025.csv`   | M87* jet-power sensitivity vs. accretion prescription (5 models)                    |
| `bh_ruler_env_SgrA_M87_V1.csv`    | Environment + TDE diagnostics for Sgr A* and M87*                                   |
| `50_BH_Survey_V2_verified.csv`    | Survey V2 source list supported directly by the updated parser                      |
| `survey_out.csv`                  | Example Survey V2-derived output using `--branch core`                              |

---

## Schema

### Standard CLI inputs

**Core inputs:**
`name, M_Msun, a_star, lambda_Edd, eta_eff, kappa, sigma_kms, Re_kpc, branch`

### Accepted aliases

**Name / identifier**

* `name`
* `Name`
* `Object`
* `object`

**Mass**

* `M_Msun`
* `Mass_Msun`
* `M`
* `mass_Msun`
* `Mass [M_sun]`
* `M [M_sun]`
* `M [Msun]`

**Spin**

* `a_star`
* `Spin_a`
* `spin`
* `Spin a_* (prograde)`
* `a_* (assumed/est.)`

**Eddington ratio**

* `lambda_Edd`
* `lambda`
* `L/L_Edd`
* `L/L_Edd (assumed)`

**Optional metadata**

* `Class`
* `Spin_Known`
* `Regime`

### Standard outputs

**Core:**
`t_g_s, r_s_km, f_ISCO_Schw_Hz, fISCO_tg_invariant`

**Spin-aware:**
`r_ISCO_rg, f_ISCO_Kerr_Hz, f_ratio`

**Environment:**
`rinfl_pc, rinfl_over_Re, t_fb_days`

**TDE:**
`tde_possible`

**Accretion/Jets:**
`B_H_G_eta, P_BZ_erg_s_eta, B_H_G_adaf, P_BZ_erg_s_adaf`

### Survey-compatible output aliases

When survey-style input is detected, `bhruler.py` also writes:

* `Name`
* `Mass_Msun`
* `Spin_a`
* `Grav_Time_tg_s`
* `Schwarzschild_Rad_km`
* `Orbital_Period_s`
* `Invariant_fISCO_tg`

This preserves compatibility with survey-style tables while keeping the canonical BHRuler columns in the same output file.

---

## Reproducing Results

### Regenerate all figures from data

```bash
cd Paper/figures
python fig1_invariant.py
python fig2_spin.py
python fig3_m87.py
python fig4_atlas.py
python fig5_uncertainty.py
```

Requires: `matplotlib`, `numpy`.

### Recompile the paper

```bash
cd Paper
pdflatex Black_Hole_Ruler_V2_2025.tex
pdflatex Black_Hole_Ruler_V2_2025.tex
```

### Quick data check

```python
import pandas as pd

# Verify invariant
check = pd.read_csv("data/bh_ruler_check_2025.csv")
print(check[["Object", "Mass [M_sun]", "f_ISCO * t_g (dimensionless)"]])

# M87* sensitivity
sens = pd.read_csv("data/M87_RIAF_sensitivity_2025.csv")
print(sens[["Model", "mdot (dotM/dotM_Edd)", "P_BZ [erg/s]"]])
```

---

## Quality Control & Provenance

* **Invariant verification:** `f_ISCO × t_g ≈ 0.01083` is checked automatically
* **Flexible schema handling:** alternate column names are normalized before parsing
* **Survey metadata preservation:** `Class`, `Spin_Known`, and `Regime` are retained when present
* **Accretion warnings:** missing `eta_eff` or `kappa` produce explicit warning columns for the affected branch
* Every table row can include provenance-style metadata for traceability

---

## Contributing

Issues and pull requests are welcome. Please:

1. Keep numeric changes **reproducible**.
2. Update documentation when adding or changing accepted input schema.
3. Run the invariant check (`f_ISCO × t_g ≈ 0.01083`) before submitting.
4. Keep `bhruler.py` self-contained, since the single-file CLI is the main public interface.

---

## Citation

If you use BHRuler or the accompanying datasets in your work, please cite:

> Nagy, S. L. (2025). *A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties.* Black Hole Ruler.

**BibTeX**

```bibtex
@software{Nagy2025BHRuler,
  author       = {Stephen L. Nagy},
  title        = {{A Scale-Invariant Ruler for Black Holes:
                   From Stellar-Mass to Ultra-Massive with Unified Uncertainties}},
  year         = {2025},
  version      = {2.1.1-survey},
  url          = {https://github.com/SNagy3/BHRuler},
  note         = {BHRuler software, data, and paper}
}
```

**Topics**

* black-holes
* astrophysics
* astronomy
* general-relativity
* gravitational-waves
* agn
* x-ray-binaries
* accretion
* relativistic-jets
* kerr-metric
* isco
* python
* scientific-computing
* scientific-python
* data-science
* simulation
* reproducible-research
* research
* space
* cosmology

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.
