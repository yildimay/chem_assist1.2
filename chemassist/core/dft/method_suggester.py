from __future__ import annotations

"""SMILES ➜ heuristic DFT recommendation engine.

*   Automatically embeds **each disconnected fragment** and translates
    subsequent fragments along +x so atoms never overlap in the output
    coordinates (useful for salts such as [Co(NH3)6]Cl3).
*   Returns a :class:`Suggestion` dataclass with method, basis, charge,
    multiplicity and a ready-to-paste Cartesian block.

If you later replace this with an ML/LLM model, keep the *public API*
(`recommend(smiles: str) -> Suggestion`) identical so UI code needn’t
change.
"""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# ───────────────────────────────────────────────────────────────
# Dataclass returned to the UI layer
# ───────────────────────────────────────────────────────────────


@dataclass
class Suggestion:
    smiles: str
    method: str
    basis: str
    reason: str
    charge: int
    multiplicity: int
    xyz: str  # Cartesian coordinates (Å)


# ───────────────────────────────────────────────────────────────
# Helper functions
# ───────────────────────────────────────────────────────────────


def _embed_xyz(mol: Chem.Mol) -> str:
    """RDKit ETKDG + UFF → XYZ block."""
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise ValueError("RDKit 3-D embedding failed.")
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    return "\n".join(
        f"{atom.GetSymbol():2} {conf.GetAtomPosition(i).x:>10.5f} "
        f"{conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}"
        for i, atom in enumerate(mol.GetAtoms())
    )


def _embed_with_translation(mol: Chem.Mol, spacing: float = 10.0) -> str:
    """Embed each fragment separately, translate along +x to avoid overlap."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    lines: list[str] = []
    offset = 0.0
    for frag in frags:
        for line in _embed_xyz(frag).splitlines():
            el, x, y, z = line.split()
            lines.append(
                f"{el:2} {float(x) + offset:>10.5f} {float(y):>10.5f} {float(z):>10.5f}"
            )
        offset += spacing
    return "\n".join(lines)


# Atomic numbers treated as *metals* (for heuristic branching)
_METALS: set[int] = {
    3, 4, 11, 12, 13,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    55, 56, 57,
}
# add lanthanides & beyond
_METALS.update(range(58, 72))  # Ce–Lu
_METALS.update(range(72, 81))  # Hf–Hg
_METALS.update({81, 82, 83})   # Tl, Pb, Bi
_METALS.update(range(84, 89))  # Po–Ac
_METALS.update(range(89, 104)) # Actinides up to Lr

_HALOGENS = {9, 17, 35, 53, 85}


def _contains(seq: set[int], mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in seq for a in mol.GetAtoms())


# ───────────────────────────────────────────────────────────────
# Main public API
# ───────────────────────────────────────────────────────────────


def recommend(smiles: str) -> Suggestion:  # noqa: D401 – simple facade
    """Return a heuristic *Suggestion* for the given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    has_metal = _contains(_METALS, mol)
    has_halogen = _contains(_HALOGENS, mol)

    charge = Chem.GetFormalCharge(mol)
    electron_count = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    multiplicity = 1 if electron_count % 2 == 0 else 2

    if has_metal:
        method, basis, reason = (
            "PBE0",
            "def2-TZVP",
            "Metal detected → hybrid GGA with triple-ζ def2 basis is a robust starting point.",
        )
    elif heavy_atoms > 50:
        method, basis, reason = (
            "ωB97X-D",
            "def2-SVP",
            "Large system (>50 heavy atoms) → range-separated functional with moderate basis to control cost.",
        )
    elif has_halogen:
        method, basis, reason = (
            "B3LYP-D3(BJ)",
            "6-311+G(d,p)",
            "Halogen present → hybrid functional with diffuse & polarisation functions recommended.",
        )
    else:
        method, basis, reason = (
            "B3LYP",
            "6-31G(d)",
            "Medium organic molecule → classic B3LYP and split-valence basis.",
        )

    xyz_block = _embed_with_translation(mol) if "." in smiles else _embed_xyz(mol)

    return Suggestion(
        smiles=smiles,
        method=method,
        basis=basis,
        reason=reason,
        charge=charge,
        multiplicity=multiplicity,
        xyz=xyz_block,
    )
