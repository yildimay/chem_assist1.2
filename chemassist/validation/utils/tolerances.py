from __future__ import annotations

from typing import Dict


class Tolerances:
    @staticmethod
    def default_gaussian() -> Dict[str, float]:
        return {
            "scf_energy_hartree": 1e-6,
            "lowest_freq_cm1": 5.0,
        }

    @staticmethod
    def default_gromacs() -> Dict[str, float]:
        return {
            "avg_temp_K": 3.0,
            "pressure_bar": 50.0,
        }


