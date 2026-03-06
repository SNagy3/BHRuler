# Black Hole Ruler
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18894798.svg)](https://doi.org/10.5281/zenodo.18894798)


A scale-invariant framework and reproducible dataset for comparing black holes across the mass spectrum (Nagy 2025).

**Stephen L. Nagy · Independent Researcher · September 2025**
*A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties*
**License:** MIT

BHRuler implements the **Black Hole Ruler** framework for **scale-invariant**, **reproducible** comparison of black holes from stellar to ultra-massive regimes. It provides code and data to compute gravitational units and ISCO landmarks, spin-aware Kerr corrections, dual accretion prescriptions (η-bridge and ADAF/RIAF), Blandford–Znajek jet-power estimates, environmental metrics (σ, r_infl, r_infl/R_e), and a TDE module with a logistic capture boundary—plus a versioned catalog schema for cross-scale studies.

---

## What's New in V2

- **Five publication-quality figures** generated from computed data (cross-scale invariant, spin diagnostics, M87* sensitivity, 10-object atlas, worked uncertainty propagation)
- **Expanded reference list** (45 citations anchoring every equation to its theoretical source)
- **Worked uncertainty example** for M87* with both analytic error propagation and a 50,000-draw Monte Carlo, showing κ dominates the P_BZ error budget at 41% of variance
- **Strengthened case studies** with direct comparisons to observed jet powers (Cyg X-1), flare timescales (Sgr A*), and EHT constraints (M87*)
- **Complete equation-to-column appendix** mapping every `bhruler.py` output to its paper equation
- **Physical constants appendix** with CODATA/IAU sources
- Removed unimplemented TDE rate equation from the paper body (retained as future-work item)
- TDE spin coefficient explicitly framed as tunable heuristic with guidance to calibrate against Kesden (2012) and Servin & Kesden (2017)
- B_H normalization anchored to Tchekhovskoy+ (2011) and McKinney+ (2012) GRMHD MAD simulations

---

## Table of Contents

- [Why BHRuler?](#why-bhruler)
- [Core Equations](#core-equations)
- [Features](#features)
- [Repository Layout](#repository-layout)
- [Quick Start](#quick-start)
- [Data Products](#data-products)
- [Schema (V1)](#schema-v1)
- [Reproducing Results](#reproducing-results)
- [Quality Control & Provenance](#quality-control--provenance)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Why BHRuler?

- **Scale-invariant:** Everything is expressed in gravitational units so times/lengths/frequencies map across 10+ orders of mass.
- **Universal invariant:** f_ISCO × t_g = 1/(6^{3/2} 2π) ≈ 0.01083 provides a built-in sanity check, verified to machine precision across 8.7 dex in mass.
- **Model-transparent:** Dual accretion prescriptions (η-bridge vs ADAF/RIAF) make **model dependence explicit** for jet-power predictions.
- **Cross-scale validation:** Unifies X-ray binaries, AGN, EHT imaging, and GW remnants on the same ruler.
- **Uncertainty-aware:** Analytic error propagation formulas + worked Monte Carlo example for M87*.
- **Reproducible:** Versioned catalog schema with per-row provenance and QC flags.

---

## Core Equations

**Gravitational units**

    t_g = GM/c³,    r_g = GM/c²,    r_s = 2r_g

**Schwarzschild ISCO**

    r_ISCO = 6 r_g
    f_ISCO = c³ / (6^{3/2} 2π GM)
    f_ISCO × t_g = 1 / (6^{3/2} 2π)    [mass-invariant]

**Kerr ISCO (equatorial, prograde/retrograde)**
Using Bardeen–Press–Teukolsky (1972):

    f_ISCO(a*) = (c³ / 2π GM) × 1 / (r_ISCO^{3/2}(a*) ± a*)

**Accretion prescriptions**

- **η-bridge:** ṁ = λ_Edd / η_eff, B_H ∝ ṁ^{1/2} M^{-1/2}, P_BZ ∝ a*² (λ_Edd/η_eff) M
- **ADAF/RIAF:** λ_Edd = κ ṁ² ⟹ ṁ = (λ_Edd/κ)^{1/2}, B_H ∝ (λ_Edd/κ)^{1/4} M^{-1/2}, P_BZ ∝ a*² (λ_Edd/κ)^{1/2} M

**Blandford–Znajek (order-of-magnitude; normalization from Tchekhovskoy+ 2011 GRMHD MAD simulations)**

    P_BZ ≈ 10⁴⁵ erg/s × (a*/0.9)² × (B_H/10⁴ G)² × (M/10⁹ M☉)²

**Environment & TDE**

    r_infl = GM/σ²
    t_fb ≈ 41 d × (M/10⁶ M☉)^{1/2} × (R*/R☉)^{3/2} × (M*/M☉)⁻¹

TDE capture boundary: logistic S(M; a*) with tunable midpoint and width.

---

## Features

- **Ruler core:** t_g, r_s, f_ISCO + invariant check
- **Spin-aware:** r_ISCO(a*), f_ISCO(a*) via BPT (prograde/retrograde)
- **Accretion toggle:** η-bridge vs ADAF/RIAF (user-selectable η_eff, κ)
- **Jet power:** B_H estimates and P_BZ ranges (MAD-motivated scaling)
- **Environment:** σ → r_infl, r_infl/R_e, morphology flags
- **TDE:** Disrupt vs. swallow flags; t_fb; tunable logistic boundary
- **Uncertainty:** Analytic log-space error formulas; worked MC example
- **Schema + QC:** Versioned, with provenance fields and automated checks

---

## Repository Layout

```
BHRuler/
├─ Paper/                           # LaTeX source and figures
│  ├─ Black_Hole_Ruler_V2_2025.tex  # Revised paper (V2.1)
│  ├─ Black_Hole_Ruler_V2_2025.pdf  # Compiled PDF
│  ├─ Black_Hole_Ruler_V1_2025.tex  # Original V1 (archived)
│  ├─ Equations.tex                 # Standalone equations reference
│  └─ figures/                      # Figure PDFs and generation scripts
│     ├─ fig1_invariant.{pdf,py}
│     ├─ fig2_spin.{pdf,py}
│     ├─ fig3_m87_sensitivity.{pdf,py}
│     ├─ fig4_atlas.{pdf,py}
│     └─ fig5_uncertainty.{pdf,py}
├─ data/                            # CSVs used in the paper & examples
│  ├─ bh_ruler_check_2025.csv
│  ├─ bh_spin_aware_ISCO_2025.csv
│  ├─ bh_trio_spin_field_jet_2025.csv
│  ├─ M87_RIAF_sensitivity_2025.csv
│  ├─ M87__jet-power_sensitivity_to_RIAF_prescription.csv
│  ├─ bh_ruler_env_SgrA_M87_V1.csv
│  └─ bh_ruler_env_v0_4_alpha_SgrA_M87.csv
├─ src/                             # Modular Python implementations
│  ├─ ruler.py
│  ├─ accretion.py
│  ├─ environment.py
│  └─ tde.py
├─ Notebooks/                       # Optional demos
├─ bhruler.py                       # Single-file CLI tool
├─ Using_BHRuler.md                 # User guide
├─ CITATION.cff
├─ CHANGELOG.md
├─ LICENSE
└─ README.md
```

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

This writes a 3-row CSV (Cygnus X-1, Sgr A*, M87*) using both accretion branches.

### Run on your own CSV

```bash
python bhruler.py --input objects.csv --branch both --output derived.csv
```

See [Using_BHRuler.md](Using_BHRuler.md) for the full user guide, input schema, worked examples, and troubleshooting.

### Programmatic use

```python
import bhruler as br

row = dict(name="Test", M_Msun=1e8, a_star=0.7, lambda_Edd=0.01,
           eta_eff=0.1, kappa=0.1, sigma_kms=200, Re_kpc=4.0)
out = br.compute_row(row, branch="both", tde_mcrit=3e7, tde_width=0.15)
```

---

## Data Products

| File | Description |
|------|-------------|
| `bh_ruler_check_2025.csv` | Cross-scale invariant check: Gaia BH3, GW231123, ω Cen IMBH, CEERS-1019, J0529-4351 |
| `bh_spin_aware_ISCO_2025.csv` | Spin-corrected ISCOs: Cyg X-1, LMC X-1, GW150914 remnant |
| `bh_trio_spin_field_jet_2025.csv` | Trio (Cyg X-1, Sgr A*, M87*) with B_H and P_BZ under dual branches |
| `M87_RIAF_sensitivity_2025.csv` | M87* jet-power sensitivity vs. accretion prescription (5 models) |
| `bh_ruler_env_SgrA_M87_V1.csv` | Environment + TDE diagnostics for Sgr A* and M87* |

---

## Schema (V1)

**Core:** `name, type, z, distance, M, M_err, anchor_method, mass_ref, a_star, a_ref, quality_grade, use_in_fits`

**Ruler:** `t_g, r_s, f_ISCO_schw, fISCO_tg_invariant`

**Spin-aware:** `r_ISCO_a, f_ISCO_a, f_ratio`

**Accretion:** `L_Edd, L_band, L_bol, lambda_Edd, eta_eff, branch{eta_bridge|ADAF}, kappa`

**Magnetics/Jets:** `B_H, phi_BH, P_BZ_pred, P_BZ_low, P_BZ_high, branch_recommendation, jet_power_ref`

**Environment:** `sigma, sigma_err, sigma_source, sigma_aperture, sigma_ref, k_factor, k_source, Re, Re_err, morph_type, bar_flag, morphology_warning, rinfl, rinfl_err, rinfl_over_Re`

**TDE:** `tde_flag{disrupt|swallow}, Mcrit_est, tfb_days`

**Versions:** `framework_version, catalog_version, row_version, last_updated, notes`

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
pdflatex Black_Hole_Ruler_V2_2025.tex   # second pass for cross-refs
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

- **M–σ outliers:** Flag when measured σ deviates from M–σ prediction by >3σ.
- **Aperture sanity:** Warn when σ aperture ≫ R_e.
- **Morphology:** `morphology_warning=true` for pseudo-bulges/bars; `k_source` recorded for dynamical fallbacks.
- Every table row includes `*_ref` fields for traceability.

---

## Contributing

Issues and pull requests are welcome. Please:

1. Keep numeric changes **reproducible** (include code or notebook).
2. Update `data/README` and `catalog_version` when adding/replacing measurements.
3. Run the invariant check (f_ISCO × t_g ≈ 0.01083) before submitting.

---

## Citation

If you use BHRuler or the accompanying datasets in your work, please cite:

> Nagy, S. L. (2025). *A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties.* BHRuler V2.

**BibTeX**

```bibtex
@software{Nagy2025BHRuler,
  author       = {Stephen L. Nagy},
  title        = {{A Scale-Invariant Ruler for Black Holes: From Stellar-Mass
                   to Ultra-Massive with Unified Uncertainties}},
  year         = {2025},
  version      = {2.1},
  url          = {https://github.com/SNagy3/BHRuler},
  note         = {BHRuler software, data, and paper}
}

topics:
black-holes
astrophysics
astronomy
general-relativity
gravitational-waves
agn
x-ray-binaries
accretion
relativistic-jets
kerr-metric
isco
python
scientific-computing
scientific-python
data-science
simulation
reproducible-research
research
space
cosmology

```

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.
