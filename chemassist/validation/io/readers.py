from __future__ import annotations

import re
from typing import Dict


def parse_gaussian_log(text: str) -> Dict[str, float]:
    """Extract minimal metrics from Gaussian log text.

    Metrics:
    - scf_energy_hartree: final SCF energy (Hartree)
    - lowest_freq_cm1: lowest vibrational frequency (cm^-1)
    """
    metrics: Dict[str, float] = {}

    # SCF energy lines examples:
    #  SCF Done:  E(RHF) =  -76.0267654321     A.U. after  10 cycles
    m = re.search(r"SCF Done:\s+E\([A-Za-z0-9+\-]+\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)
    if m:
        metrics["scf_energy_hartree"] = float(m.group(1))

    # Frequencies -- 120.5  130.2 ...
    freqs = [float(x) for x in re.findall(r"Frequencies? --\s+([\d\.+\-Ee]+)", text)]
    if freqs:
        # the regex above captures the first number in each Frequencies line; collect all
        # refine: capture all numbers on those lines
        lines = re.findall(r"Frequencies? --\s+([^\n]+)", text)
        all_vals = []
        for line in lines:
            all_vals.extend([float(v) for v in re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)])
        if all_vals:
            metrics["lowest_freq_cm1"] = min(all_vals)

    return metrics


def parse_gromacs_mdlog(text: str) -> Dict[str, float | str]:
    """Extract minimal metrics from GROMACS md.log text.

    Metrics:
    - avg_temp_K: average temperature
    - pressure_bar: average pressure (if available)
    - ensemble: NVT/NPT/EM/MD (heuristic from thermostat/barostat/integrator)
    """
    metrics: Dict[str, float | str] = {}

    # Typical md.log contains something like:
    #   Average Temperature: 300.12 K
    m = re.search(r"Average\s+Temperature:\s*([-+]?[0-9]*\.?[0-9]+)\s*K", text, re.IGNORECASE)
    if m:
        metrics["avg_temp_K"] = float(m.group(1))

    m = re.search(r"Average\s+pressure:\s*([-+]?[0-9]*\.?[0-9]+)\s*bar", text, re.IGNORECASE)
    if m:
        metrics["pressure_bar"] = float(m.group(1))

    # Ensemble heuristic
    integrator = None
    thermostat = None
    barostat = None

    mi = re.search(r"integrator\s*=\s*([A-Za-z0-9_\-]+)", text)
    if mi:
        integrator = mi.group(1).lower()

    mt = re.search(r"(T-?coupl|thermostat)\s*=\s*([A-Za-z0-9_\-]+)", text, re.IGNORECASE)
    if mt:
        thermostat = mt.group(2).lower()

    mp = re.search(r"(P-?coupl|barostat)\s*=\s*([A-Za-z0-9_\-]+)", text, re.IGNORECASE)
    if mp:
        barostat = mp.group(2).lower()

    if integrator in {"md", "sd"} and thermostat and not barostat:
        metrics["ensemble"] = "NVT"
    elif integrator in {"md", "sd"} and thermostat and barostat:
        metrics["ensemble"] = "NPT"
    elif integrator in {"steep", "cg"} or (integrator and integrator.startswith("em")):
        metrics["ensemble"] = "EM"
    else:
        metrics["ensemble"] = "MD"

    return metrics


