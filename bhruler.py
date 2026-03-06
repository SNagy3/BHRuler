#!/usr/bin/env python3
"""
BHRuler — single-file implementation of the Black Hole Ruler framework.

A scale-invariant framework for comparing black holes across the mass spectrum,
from stellar-mass X-ray binaries to ultra-massive black holes.

Usage:
  python bhruler.py --output out.csv
  python bhruler.py --input objects.csv --branch both --output derived.csv
  python bhruler.py --version

Input CSV columns (case-insensitive, extra columns ignored):
  name, M_Msun, a_star, lambda_Edd, eta_eff, kappa, sigma_kms, Re_kpc, branch

Branch options:
  - eta_bridge : uses dotm = lambda_Edd / eta_eff
  - adaf       : uses lambda_Edd = kappa * dotm^2  (=> dotm = sqrt(lambda_Edd / kappa))
  - both       : computes both and outputs columns with suffixes _eta and _adaf

Output columns (per branch):
  t_g_s, r_s_km, f_ISCO_Schw_Hz, fISCO_tg_invariant, r_ISCO_rg,
  f_ISCO_Kerr_Hz, f_ratio, B_H_G_{eta,adaf}, P_BZ_erg_s_{eta,adaf},
  rinfl_pc, rinfl_over_Re, t_fb_days, tde_possible

Notes:
  * The mass-invariant product f_ISCO_Schw * t_g = 1/(6^{3/2} 2 pi) ~ 0.01083
    is computed for every row as a built-in sanity check.
  * TDE capture boundary uses a logistic S(M; a*) with tunable midpoint and width.
    The spin coefficient (1 - 0.6*a*) is a heuristic; see Kesden (2012) and
    Servin & Kesden (2017) for detailed relativistic disruption calculations.
  * B_H normalization is MAD-motivated, calibrated against Tchekhovskoy+ (2011)
    GRMHD simulations. The scaling B_H ~ dot{m}^{1/2} M^{-1/2} is the key physics;
    the 4e4 G prefactor is order-of-magnitude.

References:
  Nagy, S. L. (2025). A Scale-Invariant Ruler for Black Holes: From Stellar-Mass
    to Ultra-Massive with Unified Uncertainties. BHRuler V2.
  https://github.com/SNagy3/BHRuler
"""

__version__ = "2.1.0"
__author__ = "Stephen L. Nagy"

import argparse
import math
import sys
from typing import Optional, Dict, Any, List, Tuple

import pandas as pd

# ---------------------------------------------------------------------------
# Physical constants (SI)
# Sources: G — CODATA 2018; c — exact (SI definition);
#          M_sun — IAU 2015 nominal; pc — IAU 2012
# ---------------------------------------------------------------------------
G     = 6.67430e-11              # m^3 kg^-1 s^-2
c     = 299_792_458.0            # m s^-1
M_sun = 1.98847e30               # kg
pc_m  = 3.085_677_581_491_367e16 # m

# Theoretical invariant: f_ISCO_Schw * t_g (Schwarzschild, any mass)
FISCO_TG_INVARIANT = 1.0 / (6**1.5 * 2.0 * math.pi)  # ≈ 0.010829122239...


# ===========================================================================
# Core ruler — gravitational scales (mass only)
# ===========================================================================

def t_g_seconds(M_Msun: float) -> float:
    """Gravitational time: t_g = GM/c^3  [seconds]."""
    return G * (M_Msun * M_sun) / c**3


def r_g_m(M_Msun: float) -> float:
    """Gravitational radius: r_g = GM/c^2  [metres]."""
    return G * (M_Msun * M_sun) / c**2


def r_s_m(M_Msun: float) -> float:
    """Schwarzschild radius: r_s = 2GM/c^2  [metres]."""
    return 2.0 * G * (M_Msun * M_sun) / c**2


def f_isco_schwarz_hz(M_Msun: float) -> float:
    """Schwarzschild ISCO orbital frequency: f = c^3 / (6^{3/2} 2π GM)  [Hz]."""
    return c**3 / (6**1.5 * 2.0 * math.pi * G * (M_Msun * M_sun))


# ===========================================================================
# Kerr spin corrections — Bardeen, Press & Teukolsky (1972)
# ===========================================================================

def kerr_isco_r_over_M(a_star: float, prograde: bool = True) -> float:
    """
    ISCO radius in units of r_g = GM/c^2, from the BPT expressions.

    Parameters
    ----------
    a_star : float
        Dimensionless spin, clamped internally to (-0.9999, +0.9999).
    prograde : bool
        If True, compute prograde (co-rotating) ISCO; else retrograde.

    Returns
    -------
    float
        r_ISCO / r_g.  Ranges from 1 (a* -> 1, prograde) to 9 (a* -> 1, retrograde).
    """
    a = max(-0.9999, min(0.9999, a_star))
    Z1 = 1.0 + (1.0 - a * a)**(1.0/3.0) * ((1.0 + a)**(1.0/3.0) + (1.0 - a)**(1.0/3.0))
    Z2 = math.sqrt(3.0 * a * a + Z1 * Z1)
    if prograde:
        return 3.0 + Z2 - math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2))
    else:
        return 3.0 + Z2 + math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2))


def f_isco_kerr_hz(M_Msun: float, a_star: float, prograde: bool = True) -> Tuple[float, float]:
    """
    Kerr ISCO orbital frequency and radius.

    Parameters
    ----------
    M_Msun : float
        Black hole mass in solar masses.
    a_star : float
        Dimensionless spin.
    prograde : bool
        Co-rotating (True) or counter-rotating (False) orbit.

    Returns
    -------
    (f_Hz, r_ISCO_rg) : tuple of float
        ISCO frequency in Hz and ISCO radius in r_g units.
    """
    r_over_M = kerr_isco_r_over_M(a_star, prograde=prograde)
    sign = a_star if prograde else -a_star
    denom = r_over_M**1.5 + sign
    f_Hz = c**3 / (2.0 * math.pi * G * (M_Msun * M_sun)) * (1.0 / denom)
    return f_Hz, r_over_M


# ===========================================================================
# Accretion prescriptions, horizon field, and jet power
# ===========================================================================

def B_H_eta(M_Msun: float, lambda_Edd: float, eta_eff: float) -> float:
    """
    Horizon magnetic field via the efficiency-bridge prescription.

    B_H ≈ 4×10^4 G × sqrt(mdot) × (10^9 M_sun / M)^{1/2}
    where mdot = lambda_Edd / eta_eff.

    Normalization: MAD-saturated flux from Tchekhovskoy+ (2011) GRMHD.
    """
    if eta_eff <= 0:
        return float("nan")
    mdot = max(0.0, lambda_Edd / eta_eff)
    return 4e4 * math.sqrt(mdot) * math.sqrt(1e9 / M_Msun)


def B_H_adaf(M_Msun: float, lambda_Edd: float, kappa: float) -> float:
    """
    Horizon magnetic field via the ADAF/RIAF prescription.

    lambda_Edd = kappa × mdot^2  =>  mdot = sqrt(lambda_Edd / kappa)
    Then same B_H scaling as the eta-bridge.

    The quadratic scaling arises because in a two-temperature ADAF, the
    radiative efficiency is itself proportional to mdot (Narayan & Yi 1995).
    """
    if kappa <= 0:
        return float("nan")
    mdot = math.sqrt(max(0.0, lambda_Edd / kappa))
    return 4e4 * math.sqrt(mdot) * math.sqrt(1e9 / M_Msun)


def P_BZ_erg_s(M_Msun: float, a_star: float, B_H: float) -> float:
    """
    Blandford–Znajek jet power (order-of-magnitude).

    P_BZ ≈ 10^45 erg/s × (a*/0.9)^2 × (B_H/10^4 G)^2 × (M/10^9 M_sun)^2

    Follows the standard P_BZ ∝ a*^2 Φ_BH^2 scaling (Blandford & Znajek 1977),
    confirmed in GRMHD simulations (Tchekhovskoy+ 2011, McKinney+ 2012).
    """
    return 1e45 * (a_star / 0.9)**2 * (B_H / 1e4)**2 * (M_Msun / 1e9)**2


# ===========================================================================
# Environment — sphere of influence and host coupling
# ===========================================================================

def r_infl_pc(M_Msun: float, sigma_kms: Optional[float]) -> Optional[float]:
    """Sphere of influence: r_infl = GM/σ^2  [parsecs].  Returns None if σ unavailable."""
    if sigma_kms is None or sigma_kms <= 0:
        return None
    M_kg = M_Msun * M_sun
    sigma_ms = sigma_kms * 1e3
    return (G * M_kg / sigma_ms**2) / pc_m


def r_infl_over_Re(rinfl_pc: Optional[float], Re_kpc: Optional[float]) -> Optional[float]:
    """Dimensionless ratio r_infl / R_e.  Returns None if inputs unavailable."""
    if rinfl_pc is None or Re_kpc in (None, 0):
        return None
    return (rinfl_pc / 1e3) / Re_kpc


def t_fb_days(M_Msun: float, Rstar_Rsun: float = 1.0, Mstar_Msun: float = 1.0) -> float:
    """
    Canonical TDE fallback time (Rees 1988, Phinney 1989).

    t_fb ≈ 41 d × (M/10^6 M_sun)^{1/2} × (R*/R_sun)^{3/2} × (M*/M_sun)^{-1}
    """
    return 41.0 * math.sqrt(M_Msun / 1e6) * (Rstar_Rsun**1.5) * (Mstar_Msun**-1.0)


# ===========================================================================
# TDE logistic capture boundary
# ===========================================================================

def logistic(x: float) -> float:
    """Standard logistic function, numerically stable for large |x|."""
    if x > 500:
        return 1.0
    if x < -500:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def tde_capture_factor(M_Msun: float, a_star: float,
                       Mcrit0: float = 3e7, width_dex: float = 0.15) -> float:
    """
    Logistic capture factor S(M; a*).

    S → 0 below M_crit (disruption dominates, TDE possible).
    S → 1 above M_crit (direct capture dominates, TDE suppressed).

    The spin-dependent midpoint M_crit = M_crit0 × (1 - 0.6 a*) is a tunable
    heuristic, NOT a fitted result.  For precision work, calibrate against
    Kesden (2012) and Servin & Kesden (2017).

    Parameters
    ----------
    M_Msun : float
        Black hole mass in solar masses (must be > 0).
    a_star : float
        Dimensionless spin (clamped to [0, 1] for the shift).
    Mcrit0 : float
        Schwarzschild midpoint mass [M_sun] (default 3e7, solar-type star).
    width_dex : float
        Logistic transition width in dex (default 0.15).
    """
    a_clamped = max(0.0, min(1.0, a_star))
    Mcrit = Mcrit0 * (1.0 - 0.6 * a_clamped)
    if Mcrit <= 0:
        Mcrit = 1e6
    x = (math.log10(M_Msun) - math.log10(Mcrit)) / width_dex
    return logistic(x)


def tde_possible_bool(M_Msun: float, a_star: float, Mcrit0: float = 3e7) -> bool:
    """Convenience: returns True if disruption dominates (S < 0.5)."""
    return tde_capture_factor(M_Msun, a_star, Mcrit0=Mcrit0) < 0.5


# ===========================================================================
# Input parsing helpers — flexible column-name resolution
# ===========================================================================

def _get_flex(obj: Dict[str, Any], *keys: str) -> Any:
    """Return the first non-None value found under any of the given keys."""
    for k in keys:
        v = obj.get(k)
        if v is not None and v != "":
            return v
    return None


def _float_flex(obj: Dict[str, Any], *keys: str, default: float = None) -> float:
    """Return float from the first matching key, or default."""
    v = _get_flex(obj, *keys)
    if v is None or v == "":
        if default is not None:
            return default
        raise ValueError(f"Required field missing; tried keys: {keys}")
    return float(v)


def _float_or_none(obj: Dict[str, Any], *keys: str) -> Optional[float]:
    """Return float from the first matching key, or None if all are absent/empty."""
    v = _get_flex(obj, *keys)
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


# ===========================================================================
# Row computation — the ruler pipeline
# ===========================================================================

def compute_row(obj: Dict[str, Any], branch: str,
                tde_mcrit: float, tde_width: float) -> Dict[str, Any]:
    """
    Apply the full Black Hole Ruler to one object.

    Parameters
    ----------
    obj : dict
        Input parameters.  Recognized keys (case-sensitive):
        name, M_Msun (or M, mass_Msun), a_star (or spin), lambda_Edd (or lambda),
        eta_eff, kappa, sigma_kms, Re_kpc.
    branch : str
        'eta_bridge', 'adaf', or 'both'.
    tde_mcrit : float
        Logistic capture midpoint [M_sun].
    tde_width : float
        Logistic transition width [dex].

    Returns
    -------
    dict
        Derived quantities.  See module docstring for column definitions.
    """
    # ---- Parse inputs with flexible synonym support ----
    # Handles column names from both the schema spec and the computed data CSVs.
    name  = _get_flex(obj, "name", "Name", "Object", "object")
    M     = _float_flex(obj, "M_Msun", "M", "mass_Msun", "Mass [M_sun]",
                        "M [M_sun]", "M [Msun]")
    a     = _float_flex(obj, "a_star", "spin", "Spin a_* (prograde)",
                        "a_* (assumed/est.)", default=0.0)
    lam   = _float_flex(obj, "lambda_Edd", "lambda", "L/L_Edd (assumed)",
                        default=0.0)
    eta   = _float_or_none(obj, "eta_eff")
    kap   = _float_or_none(obj, "kappa")
    sigma = _float_or_none(obj, "sigma_kms", "sigma [km/s]")
    Re    = _float_or_none(obj, "Re_kpc", "R_e [kpc]")

    # ---- Core ruler (mass only) ----
    tg  = t_g_seconds(M)
    f_s = f_isco_schwarz_hz(M)

    out = {
        "name":             name,
        "M_Msun":           M,
        "a_star":           a,
        "lambda_Edd":       lam,
        "sigma_kms":        sigma,
        "Re_kpc":           Re,
        "t_g_s":            tg,
        "r_s_km":           r_s_m(M) / 1e3,
        "f_ISCO_Schw_Hz":   f_s,
        "fISCO_tg_invariant": f_s * tg,   # should equal FISCO_TG_INVARIANT ≈ 0.01083
    }

    # ---- Spin-aware ISCO (Kerr) ----
    f_kerr, r_isco_rg = f_isco_kerr_hz(M, a, prograde=True)
    out.update({
        "r_ISCO_rg":        r_isco_rg,
        "f_ISCO_Kerr_Hz":   f_kerr,
        "f_ratio":          f_kerr / f_s,
    })

    # ---- Environment ----
    rinfl = r_infl_pc(M, sigma)
    out.update({
        "rinfl_pc":         rinfl,
        "rinfl_over_Re":    r_infl_over_Re(rinfl, Re),
        "t_fb_days":        t_fb_days(M),
    })

    # ---- TDE capture ----
    S = tde_capture_factor(M, a, Mcrit0=tde_mcrit, width_dex=tde_width)
    out["tde_possible"] = (S < 0.5)

    # ---- Accretion branches ----
    def add_branch(suffix: str, B: float):
        out[f"B_H_G{suffix}"]      = B
        out[f"P_BZ_erg_s{suffix}"]  = P_BZ_erg_s(M, a, B)

    if branch in ("eta_bridge", "both"):
        if eta is None or eta <= 0:
            out["warn_eta_bridge"] = "eta_eff missing or <=0"
            add_branch("_eta", float("nan"))
        else:
            add_branch("_eta", B_H_eta(M, lam, eta))

    if branch in ("adaf", "both"):
        if kap is None or kap <= 0:
            out["warn_adaf"] = "kappa missing or <=0"
            add_branch("_adaf", float("nan"))
        else:
            add_branch("_adaf", B_H_adaf(M, lam, kap))

    return out


# ===========================================================================
# CLI entry point
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="BHRuler — compute derived black-hole quantities across the mass spectrum.",
        epilog="See https://github.com/SNagy3/BHRuler for documentation and paper."
    )
    ap.add_argument("--version", action="version", version=f"BHRuler {__version__}")
    ap.add_argument(
        "--input",
        help="Input CSV (columns: name, M_Msun, a_star, lambda_Edd, eta_eff, kappa, "
             "sigma_kms, Re_kpc, branch).  If omitted, runs built-in demo trio.",
        default=None,
    )
    ap.add_argument(
        "--branch",
        choices=["eta_bridge", "adaf", "both"],
        default="both",
        help="Accretion prescription to compute (default: both).",
    )
    ap.add_argument("--output", default="bhruler_output.csv", help="Output CSV path.")
    ap.add_argument(
        "--tde-mcrit", type=float, default=3e7,
        help="TDE logistic capture midpoint mass in M_sun (default: 3e7).",
    )
    ap.add_argument(
        "--tde-width", type=float, default=0.15,
        help="TDE logistic transition width in dex (default: 0.15).",
    )

    args = ap.parse_args()

    # ---- Load or generate input ----
    if args.input is None:
        rows = [
            {"name": "Cygnus X-1", "M_Msun": 21.2,   "a_star": 0.95,
             "lambda_Edd": 2e-2,   "eta_eff": 0.1,    "kappa": 0.1,
             "sigma_kms": None,    "Re_kpc": None},
            {"name": "Sgr A*",     "M_Msun": 4.30e6,  "a_star": 0.5,
             "lambda_Edd": 1e-8,   "eta_eff": 0.1,    "kappa": 0.1,
             "sigma_kms": 120.0,   "Re_kpc": None},
            {"name": "M87*",       "M_Msun": 6.5e9,   "a_star": 0.9,
             "lambda_Edd": 3e-6,   "eta_eff": 0.01,   "kappa": 0.1,
             "sigma_kms": 350.0,   "Re_kpc": 6.3},
        ]
        print(f"BHRuler {__version__} — no input file; running built-in demo trio.")
    else:
        df_in = pd.read_csv(args.input)
        rows = df_in.to_dict(orient="records")
        print(f"BHRuler {__version__} — loaded {len(rows)} rows from {args.input}")

    # ---- Compute ----
    out_rows: List[Dict[str, Any]] = []
    for obj in rows:
        br = obj.get("branch", args.branch) or args.branch
        br = br.lower() if isinstance(br, str) else args.branch
        out_rows.append(compute_row(obj, br, args.tde_mcrit, args.tde_width))

    df_out = pd.DataFrame(out_rows)
    df_out.to_csv(args.output, index=False)

    # ---- Summary with invariant verification ----
    print(f"Wrote {args.output} with {len(df_out)} rows.")

    invariants = df_out["fISCO_tg_invariant"]
    max_dev = max(abs(v - FISCO_TG_INVARIANT) / FISCO_TG_INVARIANT for v in invariants)
    if max_dev < 1e-12:
        print(f"Invariant check PASSED: f_ISCO*t_g = {FISCO_TG_INVARIANT:.12f} "
              f"(max fractional deviation: {max_dev:.1e})")
    else:
        print(f"WARNING: Invariant check — max fractional deviation {max_dev:.3e} "
              f"exceeds 1e-12.  Possible unit or formula error.")


if __name__ == "__main__":
    main()
