from __future__ import annotations

"""SMILES → DFT recommendation engine.

Pure‑Python heuristics, no self‑import. Safe to import from any module.
"""

from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

METALS = {\
    3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,\
    31,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,\
    *range(58,72),*range(72,81),*range(81,84),*range(84,88),*range(89,104),*range(104,119)\
}


@dataclass
class Suggestion:
    smiles: str
    method: str
    basis: str
    reason: str
    charge: int
    multiplicity: int
    xyz: str  # Cartesian coordinates block


# ───────────────────────────────────────────────────────────────
# Helper functions
# ───────────────────────────────────────────────────────────────

def _embed_xyz(mol: Chem.Mol) -> str:
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise ValueError("RDKit 3‑D embedding failed.")
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    lines = [
        f"{atom.GetSymbol():2} {conf.GetAtomPosition(i).x:>10.5f} {conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}"
        for i, atom in enumerate(mol.GetAtoms())
    ]
    return "\n".join(lines)


def _contains_metal(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in METALS for a in mol.GetAtoms())
