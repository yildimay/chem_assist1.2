from __future__ import annotations

"""SMILES → DFT recommendation engine.

Pure-Python heuristics, no self-import. Safe to import from any module.
"""

from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import AllChem, PeriodicTable, rdMolDescriptors

PT = PeriodicTable.GetPeriodicTable()


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
        raise ValueError("RDKit 3-D embedding failed.")
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    lines = [
        f"{atom.GetSymbol():2} {conf.GetAtomPosition(i).x:>10.5f} {conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}"
        for i, atom in enumerate(mol.GetAtoms())
    ]
    return "\n".join(lines)


def _contains_metal(mol: Chem.Mol) -> bool:
    return any(PT.GetElement(a.GetAtomicNum()).IsMetal for a in mol.GetAtoms())


def _contains_halogen(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in {9, 17, 35, 53, 85} for a in mol.GetAtoms())


def _electron_count(mol: Chem.Mol, charge: int) -> int:
    return sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge


# ───────────────────────────────────────────────────────────────
# Public API
# ───────────────────────────────────────────────────────────────

def recommend(smiles: str) -> Suggestion:  # noqa: D401
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    metal = _contains_metal(mol)
    hal   = _contains_halogen(mol)

    charge = Chem.GetFormalCharge(mol)
    mult   = 1 if _electron_count(mol, charge) % 2 == 0 else 2

    if metal:
        method, basis, reason = (
            "PBE0",
            "def2-TZVP",
            "Molecule contains a metal – hybrid GGA plus def2 basis is a robust starting combo.",
        )
    elif heavy > 50:
        method, basis, reason = (
            "ωB97X-D",
            "def2-SVP",
            "Large system (>50 heavy atoms); range-separated functional with moderate basis to control cost.",
        )
    elif hal:
        method, basis, reason = (
            "B3LYP-D3(BJ)",
            "6-311+G(d,p)",
            "Halogen atoms present; diffuse & polarization functions recommended.",
        )
    else:
        method, basis, reason = (
            "B3LYP",
            "6-31G(d)",
            "Medium organic molecule; classic hybrid functional + split-valence basis.",
        )

    xyz = _embed_xyz(mol)

    return Suggestion(
        smiles=smiles,
        method=method,
        basis=basis,
        reason=reason,
        charge=charge,
        multiplicity=mult,
        xyz=xyz,
    )
