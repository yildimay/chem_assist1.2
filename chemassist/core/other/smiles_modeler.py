from __future__ import annotations

"""SMILES to 3D Modeler - Convert SMILES strings to molecular structures.

Based on the original implementation from chemassist_platform.
Exports:
    • SmilesModeler   – main class for SMILES processing
    • generate_2d_model(smiles) -> str
    • generate_3d_model(smiles) -> str  
    • generate_xyz_coordinates(smiles) -> str
"""

import io
from typing import Optional, Tuple
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, MolToMolBlock
    from rdkit.Chem.rdmolfiles import MolToXYZBlock
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@dataclass
class MolecularStructure:
    """Container for molecular structure data."""
    smiles: str
    mol_2d: Optional[object] = None  # RDKit mol object
    mol_3d: Optional[object] = None  # RDKit mol object with 3D coords
    xyz_coords: str = ""
    mol_block: str = ""  # MOL format for 3D visualization
    error_message: str = ""


class SmilesModeler:
    """Main class for SMILES to 3D conversion."""
    
    def __init__(self):
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is required for SMILES modeling")
    
    def process_smiles(self, smiles: str) -> MolecularStructure:
        """Process SMILES string and generate all structure formats."""
        structure = MolecularStructure(smiles=smiles)
        
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                structure.error_message = "Invalid SMILES string"
                return structure
            
            # Generate 2D structure
            structure.mol_2d = mol
            
            # Generate 3D structure
            mol_3d = Chem.AddHs(mol)  # Add hydrogens
            AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol_3d)
            structure.mol_3d = mol_3d
            
            # Generate MOL block for 3D visualization
            structure.mol_block = MolToMolBlock(mol_3d)
            
            # Generate XYZ coordinates
            structure.xyz_coords = MolToXYZBlock(mol_3d)
            
        except Exception as e:
            structure.error_message = f"Error processing SMILES: {str(e)}"
        
        return structure
    
    def validate_smiles(self, smiles: str) -> Tuple[bool, str]:
        """Validate SMILES string."""
        if not smiles.strip():
            return False, "SMILES string cannot be empty"
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string format"
        
        return True, "Valid SMILES string"


# Convenience functions for direct use
def generate_2d_model(smiles: str) -> str:
    """Generate 2D model from SMILES."""
    modeler = SmilesModeler()
    structure = modeler.process_smiles(smiles)
    return structure.mol_2d if not structure.error_message else None

def generate_3d_model(smiles: str) -> str:
    """Generate 3D model from SMILES."""
    modeler = SmilesModeler()
    structure = modeler.process_smiles(smiles)
    return structure.mol_3d if not structure.error_message else None

def generate_xyz_coordinates(smiles: str) -> str:
    """Generate XYZ coordinates from SMILES."""
    modeler = SmilesModeler()
    structure = modeler.process_smiles(smiles)
    return structure.xyz_coords if not structure.error_message else ""


if __name__ == "__main__":  # pragma: no cover
    # Test with benzene
    test_smiles = "C1=CC=CC=C1"
    modeler = SmilesModeler()
    result = modeler.process_smiles(test_smiles)
    print(f"2D mol: {result.mol_2d is not None}")
    print(f"3D mol: {result.mol_3d is not None}")
    print(f"MOL block length: {len(result.mol_block)}")
    print(f"XYZ coords:\n{result.xyz_coords}")
