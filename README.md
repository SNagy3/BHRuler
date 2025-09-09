# BHRuler
A scale-invariant framework and reproducible dataset for comparing black holes across the mass spectrum (Nagy 2025).

Stephen L. Nagy · Independent Researcher · September 2025  
*A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties (V1)*  
**License:** MIT

BHRuler implements the **Black Hole Ruler** framework for **scale-invariant**, **reproducible** comparison of black holes from stellar to ultra-massive regimes. It provides code and data to compute gravitational units and ISCO landmarks, spin-aware Kerr corrections, dual accretion prescriptions (η-bridge and ADAF/RIAF), Blandford–Znajek jet-power estimates, environmental metrics (σ, \(r_{\rm infl}\), \(r_{\rm infl}/R_e\)), and a TDE module with a logistic capture boundary—plus a versioned catalog schema for cross-scale studies.

---

## Table of Contents
- [Why BHRuler?](#why-bhruler)
- [Core Equations](#core-equations)
- [Features](#features)
- [Repository Layout](#repository-layout)
- [Quick Start](#quick-start)
- [Data Products](#data-products)
- [Schema (V1)](#schema-v1)
- [Reproducing the Trio & Sensitivity Results](#reproducing-the-trio--sensitivity-results)
- [Quality Control & Provenance](#quality-control--provenance)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Why BHRuler?

- **Scale-invariant:** Everything is expressed in gravitational units so times/lengths/frequencies map across \(10^+\) orders of mass.  
- **Universal invariant:** \(f_{\rm ISCO}\,t_g = 1/(6^{3/2}2\pi)\approx 0.01083\) provides a built-in sanity check for Schwarzschild baselines.  
- **Model-transparent:** Dual accretion prescriptions (η-bridge vs ADAF/RIAF) make **model dependence explicit** for jet-power predictions.  
- **Cross-scale validation:** Unifies X-ray binaries, AGN, EHT imaging, and GW remnants on the same ruler.  
- **Reproducible:** Versioned catalog schema with per-row provenance and QC flags.

---

## Core Equations

**Gravitational units**
\[
t_g=\frac{GM}{c^3},\quad r_g=\frac{GM}{c^2},\quad r_s=2r_g
\]

**Schwarzschild ISCO**
\[
r_{\rm ISCO}=6r_g,\qquad f_{\rm ISCO}=\frac{c^3}{6^{3/2}2\pi GM},\qquad f_{\rm ISCO}\,t_g=\frac{1}{6^{3/2}2\pi}
\]

**Kerr ISCO (equatorial, prograde/retrograde)**  
Using Bardeen–Press–Teukolsky:
\[
f_{\rm ISCO}(a_*)=\frac{c^3}{2\pi GM}\,\frac{1}{r_{\rm ISCO}^{3/2}(a_*)\pm a_*}
\]

**Accretion prescriptions**

- **η-bridge:** \(\dot m=\lambda_{\rm Edd}/\eta_{\rm eff}\), \(B_H\propto \dot m^{1/2}M^{-1/2}\), \(P_{\rm BZ}\propto a_*^2(\lambda_{\rm Edd}/\eta_{\rm eff})M\).
- **ADAF/RIAF:** \(\lambda_{\rm Edd}=\kappa\dot m^2\Rightarrow \dot m=(\lambda_{\rm Edd}/\kappa)^{1/2}\), \(B_H\propto(\lambda_{\rm Edd}/\kappa)^{1/4}M^{-1/2}\), \(P_{\rm BZ}\propto a_*^2(\lambda_{\rm Edd}/\kappa)^{1/2}M\).

**Blandford–Znajek (order-of-magnitude)**
\[
P_{\rm BZ}\approx 10^{45}\ {\rm erg\,s^{-1}}\left(\frac{a_*}{0.9}\right)^2\left(\frac{B_H}{10^4\ {\rm G}}\right)^2\left(\frac{M}{10^9M_\odot}\right)^2
\]

**Environment & TDE**
\[
r_{\rm infl}=\frac{GM}{\sigma^2},\qquad
t_{\rm fb}\approx 41\,{\rm d}\left(\frac{M}{10^6M_\odot}\right)^{1/2}\!\left(\frac{R_*}{R_\odot}\right)^{3/2}\!\left(\frac{M_\odot}{m_*}\right)
\]
TDE capture boundary handled with a logistic truncation \(S(M;a_*)\).

---

## Features

- **Ruler core:** \(t_g, r_s, f_{\rm ISCO}\) + invariant check.  
- **Spin-aware:** \(r_{\rm ISCO}(a_*), f_{\rm ISCO}(a_*)\) (prograde/retrograde).  
- **Accretion toggle:** η-bridge vs ADAF/RIAF (user-selectable \(\eta_{\rm eff},\kappa\)).  
- **Jet power:** \(B_H\) estimates and \(P_{\rm BZ}\) ranges (MAD/SANE recommendation).  
- **Environment:** \(\sigma\rightarrow r_{\rm infl}\), \(r_{\rm infl}/R_e\), morphology flags.  
- **TDE:** Disrupt vs swallow flags; \(t_{\rm fb}\); per-galaxy rate model placeholders.  
- **Schema + QC:** Versioned, with provenance fields and automated checks.

---

## Repository Layout

```

BHRuler/
├─ data/                  # CSVs used in the paper & examples
│  ├─ bh\_ruler\_check\_2025.csv
│  ├─ bh\_spin\_aware\_ISCO\_2025.csv
│  ├─ bh\_trio\_spin\_field\_jet\_2025.csv
│  ├─ M87\_RIAF\_sensitivity\_2025.csv
│  └─ bh\_ruler\_env\_SgrA\_M87\_V1.csv
├─ src/                   # Python implementations
│  ├─ ruler.py            # core units, ISCO (Schwarz/Kerr), invariants
│  ├─ accretion.py        # η-bridge & ADAF/RIAF, B\_H, P\_BZ
│  ├─ environment.py      # σ→r\_infl, ratios, QC checks
│  └─ tde.py              # capture boundary S(M;a\*), t\_fb, rate stubs
├─ notebooks/             # optional demos reproducing figures/tables
├─ paper/                 # LaTeX source for V1
│  └─ Black\_Hole\_Ruler\_V1\_Nagy\_2025.tex
├─ README.md
├─ LICENSE
└─ CITATION.cff

---

## Quick Start

**Requirements:** Python ≥ 3.10

```bash
# (optional) create a virtual environment
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate

pip install numpy pandas matplotlib
````

**Minimal example (compute the ruler for a 10 M☉ Schwarzschild BH)**

```python
import numpy as np

G = 6.67430e-11
c = 2.99792458e8
M_sun = 1.98847e30

def tg_seconds(M_Msun): return G*(M_Msun*M_sun)/c**3
def rs_m(M_Msun):       return 2*G*(M_Msun*M_sun)/c**2
def fisco_hz(M_Msun):   return c**3/((6**1.5)*2*np.pi*G*(M_Msun*M_sun))

M = 10.0
print("t_g [s] =", tg_seconds(M))
print("r_s [km] =", rs_m(M)/1e3)
print("f_ISCO [Hz] =", fisco_hz(M))
print("Invariant f_ISCO * t_g =", fisco_hz(M)*tg_seconds(M))  # ~0.01083
```

**Kerr ISCO (prograde) helper**

```python
def r_isco_over_M(a):
    Z1 = 1 + (1 - a*a)**(1/3)*((1+a)**(1/3) + (1-a)**(1/3))
    Z2 = np.sqrt(3*a*a + Z1*Z1)
    return 3 + Z2 - np.sqrt((3 - Z1)*(3 + Z1 + 2*Z2))

def fisco_kerr_hz(M_Msun, a):
    rM = r_isco_over_M(a)
    denom = (rM**1.5 + a)
    return c**3/(2*np.pi*G*(M_Msun*M_sun)) * (1.0/denom)
```

---

## Data Products

* `bh_ruler_check_2025.csv` — cross-scale sanity check table (Gaia BH3, GW231123 remnant, ΩCen IMBH, CEERS-1019, J0529-4351).
* `bh_spin_aware_ISCO_2025.csv` — spin-corrected ISCOs for credible spins (e.g., Cyg X-1, LMC X-1, GW150914 remnant).
* `bh_trio_spin_field_jet_2025.csv` — **trio** (Cyg X-1, Sgr A\*, M87\*) with η-bridge/RIAF branches: $r_{\rm ISCO}(a_*), f_{\rm ISCO}(a_*), B_H, P_{\rm BZ}$.
* `M87_RIAF_sensitivity_2025.csv` — M87\* jet-power sensitivity vs. ADAF parameters (κ).
* `bh_ruler_env_SgrA_M87_V1.csv` — environment metrics and TDE flags for Sgr A\* and M87\*.

Each CSV includes a brief header and column units. Keep provenance in a companion `README` within `/data`.

---

## Schema (V1)

**Core:** `name, type, z, distance, M, M_err, anchor_method, mass_ref, a_star, a_ref, quality_grade, use_in_fits`
**Ruler:** `t_g, r_s, f_ISCO_schw, fISCO_tg_invariant`
**Spin-aware:** `r_ISCO_a, f_ISCO_a, f_ratio`
**Accretion:** `L_Edd, L_band, L_bol, lambda_Edd, eta_eff, branch{eta_bridge|ADAF}, kappa`
**Magnetics/Jets:** `B_H, phi_BH, P_BZ_pred, P_BZ_low, P_BZ_high, branch_recommendation, jet_power_ref`
**Environment:** `sigma, sigma_err, sigma_source, sigma_aperture, sigma_ref, k_factor, k_source, Re, Re_err, morph_type, bar_flag, morphology_warning, rinfl, rinfl_err, rinfl_over_Re`
**TDE:** `tde_flag{disrupt|swallow}, Mcrit_est, tfb_days, Gamma_gal, Gamma_err, rate_method, rate_params_ref`
**Imaging (if any):** `theta_shadow, image_band, baseline_req, polarization_frac, EHT_ref`
**Versions:** `framework_version, catalog_version, row_version, last_updated, notes`

---

## Reproducing the Trio & Sensitivity Results

```python
import pandas as pd

trio = pd.read_csv("data/bh_trio_spin_field_jet_2025.csv")
print(trio.filter(["Object","a_* (assumed/est.)","r_ISCO [r_g] (Kerr)","f_Kerr / f_Schwarz","B_H [G]","P_BZ [erg/s]"]))
```

**M87* ADAF sensitivity:*\*

```python
sens = pd.read_csv("data/M87_RIAF_sensitivity_2025.csv")
print(sens[["Model","mdot (dotM/dotM_Edd)","B_H [G]","P_BZ [erg/s]"]])
```

---

## Quality Control & Provenance

* **M–σ outliers:** flag when measured `sigma` deviates from M–σ prediction by >3σ.
* **Aperture sanity:** warn when `sigma_aperture` ≫ `Re`.
* **Morphology:** `morphology_warning=true` for pseudo-bulges/bars; `k_source={orbit_model|virial}` recorded for dynamical fallbacks.
* Every table row includes `*_ref` fields for traceability.

---

## Contributing

Issues and pull requests are welcome. Please:

1. Keep numeric changes **reproducible** (include code or notebook).
2. Update `data/README` and `catalog_version` when adding/replacing measurements.
3. Run simple invariant checks (e.g., $f_{\rm ISCO}t_g$) before submitting.

---

## Citation

If you use BHRuler or the accompanying datasets in your work, please cite:

> **Nagy, S. L. (2025)**, *A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties (V1)*, (arXiv\:XXXX.YYYYY).

**BibTeX**

```bibtex
@misc{Nagy2025BHRuler,
  author       = {Stephen L. Nagy},
  title        = {A Scale-Invariant Ruler for Black Holes: From Stellar-Mass to Ultra-Massive with Unified Uncertainties (V1)},
  year         = {2025},
  eprint       = {XXXX.YYYYY},
  archivePrefix= {arXiv},
  primaryClass = {astro-ph.HE},
  note         = {BHRuler software and data}
}
```

---

## License

This project is licensed under the **MIT License**. See `LICENSE` for details.
