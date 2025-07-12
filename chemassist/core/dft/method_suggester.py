from __future__ import annotations

# chemassist/core/dft/method_suggester.py
# --------------------------------------
# Minimal, self-contained recommender: SMILES  →  (method, basis, reason, coords)

from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


@dataclass
class Suggestion:
    smiles: str
    method: str
    basis: str
    reason: str
    charge: int
    multiplicity: int
    xyz: str  # Cartesian coordinates


# ----- helpers ---------------------------------------------------------------


def _embed_xyz(mol: Chem.Mol) -> str:
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise ValueError("RDKit 3-D embedding failed.")
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    return "\n".join(
        f"{a.GetSymbol():2} {conf.GetAtomPosition(i).x:>10.5f}"
        f" {conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}"
        for i, a in enumerate(mol.GetAtoms())
    )


METALS = {
    3, 4, 11, 12, 13,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    55, 56, 57, *range(58, 72), *range(72, 81), *range(81, 84),
    *range(84, 88), *range(89, 104), *range(104, 119)
}
HALOGENS = {9, 17, 35, 53, 85}


def _contains_metal(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in METALS for a in mol.GetAtoms())


def _contains_halogen(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in HALOGENS for a in mol.GetAtoms())


# ----- public API ------------------------------------------------------------


def recommend(smiles: str) -> Suggestion:          #  ←  THIS is the symbol UI imports
    """
    Quick heuristic:
        • metal → PBE0 / def2-TZVP
        • >50 heavy atoms → ωB97X-D / def2-SVP
        • halogen → B3LYP-D3(BJ) / 6-311+G(d,p)
        • else → B3LYP / 6-31G(d)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    charge = Chem.GetFormalCharge(mol)
    multiplicity = 1 if (sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge) % 2 == 0 else 2

    if _contains_metal(mol):
        method, basis, reason = (
            "PBE0",
            "def2-TZVP",
            "Metal detected → hybrid GGA with def2 triple-ζ basis."
        )
    elif heavy > 50:
        method, basis, reason = (
            "ωB97X-D",
            "def2-SVP",
            "Large molecule (>50 heavy atoms) → range-separated functional, modest basis."
        )
    elif _contains_halogen(mol):
        method, basis, reason = (
            "B3LYP-D3(BJ)",
            "6-311+G(d,p)",
            "Halogen present → hybrid + diffuse & polarisation functions."
        )
    else:
        method, basis, reason = (
            "B3LYP",
            "6-31G(d)",
            "Medium organic molecule → classic hybrid and split-valence basis."
        )

    xyz = _embed_xyz(mol)

    return Suggestion(
        smiles=smiles,
        method=method,
        basis=basis,
        reason=reason,
        charge=charge,
        multiplicity=multiplicity,
        xyz=xyz,
    )
def _embed_with_translation(mol, spacing=10.0):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    coords_lines = []
    shift = 0.0
    for frag in frags:
        xyz = _embed_xyz(frag)              # current helper
        for line in xyz.splitlines():
            el, x, y, z = line.split()
            coords_lines.append(
                f\"{el} {float(x)+shift:>10.5f} {float(y):>10.5f} {float(z):>10.5f}\"
            )
        shift += spacing
    return \"\\n\".join(coords_lines)
