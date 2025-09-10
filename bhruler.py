
#!/usr/bin/env python3
"""
BHRuler — single-file implementation of the Black Hole Ruler framework (V1)

Usage:
  python bhruler.py --output out.csv
  python bhruler.py --input objects.csv --branch both --output derived.csv

Input CSV columns (case-insensitive, extra columns ignored):
  name, M_Msun, a_star, lambda_Edd, eta_eff, kappa, sigma_kms, Re_kpc, branch

Branch options:
  - eta_bridge : uses dotm = lambda_Edd / eta_eff
  - adaf       : uses lambda_Edd = kappa * dotm^2  (=> dotm = sqrt(lambda_Edd / kappa))
  - both       : computes both and outputs columns with suffixes _eta and _adaf

Outputs (subset of V1 schema, per branch):
  t_g, r_s, r_ISCO_rg, f_ISCO_Schw, f_ISCO_Kerr, f_ratio, B_H, P_BZ,
  rinfl_pc, rinfl_over_Re, t_fb_days, tde_possible

Notes:
  * TDE capture boundary is modeled with a logistic S(M; a*) with default midpoint Mcrit=3e7 Msun (solar-type, Schwarzschild) and a spin scaling.
    This is configurable via --tde-mcrit and --tde-width parameters.
"""

import argparse
import math
import sys
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any, List
import pandas as pd

# Physical constants (SI)
G = 6.67430e-11
c = 299792458.0
M_sun = 1.98847e30
pc_m = 3.085677581491367e16

# ---------- Core ruler ----------
def t_g_seconds(M_Msun: float) -> float:
    return G*(M_Msun*M_sun)/c**3

def r_s_m(M_Msun: float) -> float:
    return 2*G*(M_Msun*M_sun)/c**2

def f_isco_schwarz_hz(M_Msun: float) -> float:
    return c**3/((6**1.5)*2*math.pi*G*(M_Msun*M_sun))

def kerr_isco_r_over_M(a_star: float, prograde: bool=True) -> float:
    """Bardeen-Press-Teukolsky ISCO radius (in r_g = GM/c^2 units)."""
    a = max(-0.9999, min(0.9999, a_star))
    Z1 = 1 + (1 - a*a)**(1/3) * ((1 + a)**(1/3) + (1 - a)**(1/3))
    Z2 = math.sqrt(3*a*a + Z1*Z1)
    if prograde:
        r_over_M = 3 + Z2 - math.sqrt((3 - Z1) * (3 + Z1 + 2*Z2))
    else:
        r_over_M = 3 + Z2 + math.sqrt((3 - Z1) * (3 + Z1 + 2*Z2))
    return r_over_M

def f_isco_kerr_hz(M_Msun: float, a_star: float, prograde: bool=True) -> float:
    r_over_M = kerr_isco_r_over_M(a_star, prograde=prograde)
    denom = (r_over_M**1.5 + (a_star if prograde else -a_star))
    return c**3/(2*math.pi*G*(M_Msun*M_sun)) * (1.0/denom), r_over_M

# ---------- Accretion, fields, jets ----------
def B_H_eta(M_Msun: float, lambda_Edd: float, eta_eff: float) -> float:
    """MAD-inspired horizon field: B ~ 4e4 G * sqrt(mdot) * (1e9/M)^{1/2}, with mdot = lambda_Edd / eta_eff."""
    if eta_eff <= 0: return float("nan")
    mdot = max(0.0, lambda_Edd/eta_eff)
    return 4e4 * math.sqrt(mdot) * math.sqrt(1e9 / M_Msun)

def B_H_adaf(M_Msun: float, lambda_Edd: float, kappa: float) -> float:
    """ADAF: lambda = kappa * mdot^2 => mdot = sqrt(lambda/kappa)."""
    if kappa <= 0: return float("nan")
    mdot = math.sqrt(max(0.0, lambda_Edd / kappa))
    return 4e4 * math.sqrt(mdot) * math.sqrt(1e9 / M_Msun)

def P_BZ_erg_s(M_Msun: float, a_star: float, B_H: float) -> float:
    return 1e45 * (a_star/0.9)**2 * (B_H/1e4)**2 * (M_Msun/1e9)**2

# ---------- Environment ----------
def r_infl_pc(M_Msun: float, sigma_kms: Optional[float]) -> Optional[float]:
    if sigma_kms is None or sigma_kms <= 0: return None
    M = M_Msun*M_sun
    sigma = sigma_kms*1e3
    r_m = G*M/sigma**2
    return r_m/pc_m

def r_infl_over_Re(rinfl_pc: Optional[float], Re_kpc: Optional[float]) -> Optional[float]:
    if rinfl_pc is None or Re_kpc in (None, 0): return None
    return (rinfl_pc/1e3)/Re_kpc

def t_fb_days(M_Msun: float, Rstar_Rsun: float=1.0, Mstar_Msun: float=1.0) -> float:
    return 41.0 * math.sqrt(M_Msun/1e6) * (Rstar_Rsun**1.5) * (Mstar_Msun**-1.0)

# ---------- TDE logistic boundary ----------
def logistic(x: float) -> float:
    return 1.0/(1.0+math.exp(-x))

def tde_capture_factor(M_Msun: float, a_star: float, Mcrit0: float=3e7, width_dex: float=0.15) -> float:
    """
    Logistic capture factor S(M; a*). Midpoint Mcrit = Mcrit0 * (1 - 0.6*a*) approx (spin raises disruption threshold).
    width_dex sets the transition breadth in log10(M).
    """
    # Spin scaling: prograde spin raises disruption threshold => lowers capture at fixed mass.
    # Here we shift the midpoint down modestly with a* to reflect enhanced disruption (fewer captures).
    Mcrit = Mcrit0 * (1.0 - 0.6*max(0.0, min(1.0, a_star)))  # simplistic heuristic
    if Mcrit <= 0: Mcrit = 1e6
    x = (math.log10(M_Msun) - math.log10(Mcrit)) / width_dex
    # S -> near 0 below Mcrit (disruption dominates), -> 1 above Mcrit (capture dominates)
    return logistic(x)

def tde_possible_bool(M_Msun: float, a_star: float, Mcrit0: float=3e7) -> bool:
    return tde_capture_factor(M_Msun, a_star, Mcrit0=Mcrit0) < 0.5

# ---------- Computation ----------
def compute_row(obj: Dict[str, Any], branch: str, tde_mcrit: float, tde_width: float) -> Dict[str, Any]:
    # Required parameters
    name = obj.get("name") or obj.get("Name")
    M = float(obj.get("M_Msun") or obj.get("M") or obj.get("mass_Msun"))
    a = float(obj.get("a_star") or obj.get("spin", 0.0))
    lam = float(obj.get("lambda_Edd") or obj.get("lambda", 0.0))
    eta = obj.get("eta_eff")
    eta = float(eta) if eta not in (None, "") else None
    kap = obj.get("kappa")
    kap = float(kap) if kap not in (None, "") else None
    sigma = obj.get("sigma_kms")
    sigma = float(sigma) if sigma not in (None, "") else None
    Re = obj.get("Re_kpc")
    Re = float(Re) if Re not in (None, "") else None

    out = {
        "name": name, "M_Msun": M, "a_star": a, "lambda_Edd": lam,
        "sigma_kms": sigma, "Re_kpc": Re,
        "t_g_s": t_g_seconds(M), "r_s_km": r_s_m(M)/1e3,
        "f_ISCO_Schw_Hz": f_isco_schwarz_hz(M)
    }
    f_kerr, r_isco_rg = f_isco_kerr_hz(M, a, prograde=True)
    out.update({
        "r_ISCO_rg": r_isco_rg,
        "f_ISCO_Kerr_Hz": f_kerr,
        "f_ratio": f_kerr / out["f_ISCO_Schw_Hz"]
    })
    # Environment
    rinfl = r_infl_pc(M, sigma)
    out.update({
        "rinfl_pc": rinfl,
        "rinfl_over_Re": r_infl_over_Re(rinfl, Re),
        "t_fb_days": t_fb_days(M)
    })
    # TDE
    S = tde_capture_factor(M, a, Mcrit0=tde_mcrit, width_dex=tde_width)
    out["tde_possible"] = (S < 0.5)

    def add_branch(suffix: str, B: float):
        out[f"B_H_G{suffix}"] = B
        out[f"P_BZ_erg_s{suffix}"] = P_BZ_erg_s(M, a, B)

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

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="BHRuler — compute derived quantities for black holes")
    ap.add_argument("--input", help="CSV with columns: name,M_Msun,a_star,lambda_Edd,eta_eff,kappa,sigma_kms,Re_kpc,branch", default=None)
    ap.add_argument("--branch", choices=["eta_bridge","adaf","both"], default="both", help="which accretion prescription to compute (default: both)")
    ap.add_argument("--output", help="output CSV path", default="bhruler_output.csv")
    ap.add_argument("--tde-mcrit", type=float, default=3e7, help="logistic capture midpoint mass in Msun (default 3e7)")
    ap.add_argument("--tde-width", type=float, default=0.15, help="logistic transition width in dex (default 0.15)")

    args = ap.parse_args()

    if args.input is None:
        # Use built-in trio as a demo
        rows = [
            {"name":"Cygnus X-1","M_Msun":21.2,"a_star":0.95,"lambda_Edd":2e-2,"eta_eff":0.1,"kappa":0.1,"sigma_kms":None,"Re_kpc":None},
            {"name":"Sgr A*","M_Msun":4.30e6,"a_star":0.5,"lambda_Edd":1e-8,"eta_eff":0.1,"kappa":0.1,"sigma_kms":120.0,"Re_kpc":None},
            {"name":"M87*","M_Msun":6.5e9,"a_star":0.9,"lambda_Edd":3e-6,"eta_eff":0.01,"kappa":0.1,"sigma_kms":350.0,"Re_kpc":6.3},
        ]
    else:
        df_in = pd.read_csv(args.input)
        rows = df_in.to_dict(orient="records")

    out_rows: List[Dict[str,Any]] = []
    for obj in rows:
        br = obj.get("branch", args.branch) or args.branch
        br = br.lower() if isinstance(br,str) else args.branch
        out_rows.append(compute_row(obj, br, args.tde_mcrit, args.tde_width))

    df_out = pd.DataFrame(out_rows)
    df_out.to_csv(args.output, index=False)
    print(f"Wrote {args.output} with {len(df_out)} rows.")

if __name__ == "__main__":
    main()
