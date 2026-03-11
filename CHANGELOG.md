# Changelog

All notable changes to the BHRuler project.

## [2.1.1-survey] - 2026-03 (Survey integration update)

### Repository

#### Added
- Native support for **Survey V2 CSV schema** in `bhruler.py`, including aliases for:
  - `Name`
  - `Class`
  - `Mass_Msun`
  - `Spin_a`
  - `Spin_Known`
  - `Regime`
- New **`core` branch** for mass/spin/timing/ISCO/TDE computations without requiring accretion inputs (`lambda_Edd`, `eta_eff`, `kappa`)
- Survey-compatible output aliases written alongside the standard BHRuler schema:
  - `Name`
  - `Mass_Msun`
  - `Spin_a`
  - `Grav_Time_tg_s`
  - `Schwarzschild_Rad_km`
  - `Orbital_Period_s`
  - `Invariant_fISCO_tg`
- Pass-through preservation of survey metadata fields (`Class`, `Spin_Known`, `Regime`) in output tables

#### Changed
- Input parsing made more tolerant across mixed CSV schemas by normalizing column matching over:
  - case differences
  - spaces
  - underscores
  - punctuation
- CLI usage and help text updated to document direct Survey V2 ingestion
- Recommended survey workflow updated to:
  - `python bhruler.py --input 50_BH_Survey_V2_verified.csv --branch core --output survey_out.csv`
- Output now includes both:
  - canonical BHRuler field names for pipeline compatibility
  - survey-style aliases for easier round-trip comparison with Survey V2 tables

#### Fixed
- Fixed failure to ingest Survey V2 files caused by schema mismatches between:
  - `Mass_Msun` vs `M_Msun`
  - `Spin_a` vs `a_star`
- Fixed brittle CSV parsing for alternate user/export naming conventions
- Fixed survey execution path so datasets without accretion parameters no longer produce unnecessary branch warnings when run with `--branch core`

#### Validated
- Replacement `bhruler.py` successfully executed on the full **49-row Survey V2 dataset**
- Generated survey output preserves the invariant check and produces the expected combined standard + survey-compatible columns

## [2.1] - 2025-09 (Paper revision)

### Paper (V2.1)

#### Added
- **5 publication-quality figures** generated from computed data:
  - Fig 1: Cross-scale invariant verification (dual-panel: f vs M + invariant bar chart)
  - Fig 2: Spin-aware ISCO radius and frequency ratio with BPT theory curves
  - Fig 3: M87* sensitivity analysis (B_H and P_BZ vs accretion prescription)
  - Fig 4: 10-object cross-scale atlas with detector bands (LIGO, LISA, PTA)
  - Fig 5: Worked uncertainty propagation (error budget + 50k-draw Monte Carlo)
- **Figure generation scripts** (Python/matplotlib) in `Paper/figures/`
- **Appendix A:** Complete equation-to-column mapping for `bhruler.py` output schema
- **Appendix B:** Physical constants table with CODATA/IAU sources
- **Worked uncertainty example** for M87* P_BZ (ADAF branch): analytic formula + MC validation showing ~0.5 dex 68% CI, dominated by κ (41% of variance)
- **45 references** (up from 6), including:
  - Bardeen, Press & Teukolsky (1972) — Kerr ISCO
  - Blandford & Znajek (1977) — jet power extraction
  - Shakura & Sunyaev (1973) — thin disk
  - Narayan & Yi (1994, 1995) — ADAF
  - Tchekhovskoy+ (2011), McKinney+ (2012) — GRMHD MAD simulations
  - Kesden (2012), Servin & Kesden (2017) — spin-dependent TDE disruption
  - Rees (1988), Stone & Metzger (2016), French+ (2020) — TDE theory/rates
  - Discovery papers for Gaia BH3, ω Cen IMBH, CEERS-1019, J0529-4351
  - Observational jet-power references for Cyg X-1 and M87*

#### Changed
- **B_H normalization** explicitly anchored to Tchekhovskoy+ (2011) GRMHD MAD saturated flux (φ_BH ~ 50)
- **TDE spin coefficient** (0.6) reframed as tunable heuristic with explicit guidance to calibrate against Kesden (2012)
- **Case studies** strengthened with direct observational comparisons:
  - Cyg X-1: P_BZ compared to Gallo+ (2005) and Russell+ (2007) jet-ISM estimates
  - Sgr A*: ISCO period linked to ~20-min QPO flare timescale (Genzel+ 2003, GRAVITY 2018)
  - M87*: Full 5-model sensitivity table with observed jet power bracketing
- **Uncertainty section** now honest about implementation status: formulas derived and demonstrated, code integration planned for V2
- **Introduction** rewritten with explicit design principles and stronger argument for pipeline integration value

#### Removed
- Unimplemented TDE rate equation (Γ_gal) removed from paper body; retained as future-work item

### Repository

#### Added
- `CHANGELOG.md` (this file)
- `Paper/figures/` directory with PDF figures and Python generation scripts

#### Changed
- `README.md` updated with V2 content, "What's New" section, figure reproduction instructions
- `CITATION.cff` updated to version 2.1

## [1.0] - 2025-09 (Initial release)

### Added
- Core framework: t_g, r_s, f_ISCO, Kerr ISCO via BPT
- Dual accretion prescriptions (η-bridge and ADAF/RIAF)
- B_H and P_BZ scaling (MAD-motivated)
- Environment module (r_infl, r_infl/R_e)
- TDE module with logistic capture boundary
- `bhruler.py` single-file CLI tool
- `Using_BHRuler.md` user guide
- 5 data CSVs (ruler check, spin-aware ISCO, trio, M87 sensitivity, environment)
- LaTeX paper (V1)
- Modular source in `src/`
