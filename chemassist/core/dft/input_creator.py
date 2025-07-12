from __future__ import annotations

"""SMILES ➜ heuristic DFT recommendation engine.

*  Embeds **each fragment** separately. If RDKit's ETKDG fails (common for
   metal complexes), falls back to CoordGen plus a simple ±1 Å z-lift to
   avoid planar collapse – the "old-fashioned" but robust route.
*  Subsequent fragments are translated along +x so dot-separated salts
   never overlap.
"""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdCoordGen

# ───────────────────────────────────────────────────────────────
# Dataclass handed to the UI layer
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
# Helper: single-fragment embedding
# ───────────────────────────────────────────────────────────────


def _embed_one(mol: Chem.Mol) -> str:
    """Try ETKDG+UFF → fallback to CoordGen with ±1 Å z jitter."""
    m = Chem.AddHs(mol)
    try:
        if AllChem.EmbedMolecule(m, AllChem.ETKDG()) == 0:
            AllChem.UFFOptimizeMolecule(m)
        else:
            raise ValueError("ETKDG failed")
    except Exception:  # noqa: BLE001
        # 2-D CoordGen then manual z-lift
        rdCoordGen.AddCoords(m)
        conf = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            pos.z = 1.0 if i % 2 else -1.0
            conf.SetAtomPosition(i, pos)
    conf = m.GetConformer()
    return "\n".join(
        f"{a.GetSymbol():2} {conf.GetAtomPosition(i).x:>10.5f} "
        f"{conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}"
        for i, a in enumerate(m.GetAtoms())
    )


def _embed_with_translation(mol: Chem.Mol, spacing: float = 10.0) -> str:
    """Embed each disconnected fragment and shift along x-axis."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    lines, offset = [], 0.0
    for frag in frags:
        for line in _embed_one(frag).splitlines():
            el, x, y, z = line.split()
            lines.append(
                f"{el:2} {float(x) + offset:>10.5f} {float(y):>10.5f} {float(z):>10.5f}"
            )
        offset += spacing
    return "\n".join(lines)


# ───────────────────────────────────────────────────────────────
# Heuristic method/basis picker helpers
# ───────────────────────────────────────────────────────────────

_METALS = {
    3, 4, 11, 12, 13,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    55, 56, 57, *range(58, 72), *range(72, 81), 81, 82, 83,
    *range(84, 89), *range(89, 104)
}
_HALOGENS = {9, 17, 35, 53, 85}


def _contains(seq: set[int], mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() in seq for a in mol.GetAtoms())


# ───────────────────────────────────────────────────────────────
# Public API
# ───────────────────────────────────────────────────────────────


def recommend(smiles: str) -> Suggestion:  # noqa: D401
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    has_metal = _contains(_METALS, mol)
    has_hal = _contains(_HALOGENS, mol)

    charge = Chem.GetFormalCharge(mol)
    e_cnt = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    mult = 1 if e_cnt % 2 == 0 else 2

    if has_metal:
        method, basis, reason = (
            "PBE0",
            "def2-TZVP",
            "Metal detected → PBE0 with def2 triple-ζ basis is a safe baseline.",
        )
    elif heavy > 50:
        method, basis, reason = (
            "ωB97X-D",
            "def2-SVP",
            "Large molecule (>50 heavy atoms) → range-separated functional with moderate basis.",
        )
    elif has_hal:
        method, basis, reason = (
            "B3LYP-D3(BJ)",
            "6-311+G(d,p)",
            "Halogen present → hybrid functional with diffuse & polarisation functions.",
        )
    else:
        method, basis, reason = (
            "B3LYP",
            "6-31G(d)",
            "Medium organic molecule → classic B3LYP and split-valence basis.",
        )

    xyz = _embed_with_translation(mol) if "." in smiles else _embed_one(mol)

    return Suggestion(
        smiles=smiles,
        method=method,
        basis=basis,
        reason=reason,
        charge=charge,
        multiplicity=mult,
        xyz=xyz,
    )
