from __future__ import annotations

"""SMILES to 3D Modeler - Convert SMILES strings to molecular structures.

Exports:
    • SmilesModeler   – main class for SMILES processing
    • generate_2d_model(smiles) -> str
    • generate_3d_model(smiles) -> str  
    • generate_xyz_coordinates(smiles) -> str
"""

import io
import tempfile
import os
from typing import Optional, Tuple
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdMolDescriptors
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
    svg_2d: str = ""
    svg_3d: str = ""
    html_3d: str = ""  # Interactive 3D HTML
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
            structure.svg_2d = self._generate_2d_svg(mol)
            
            # Generate 3D structure
            mol_3d = Chem.AddHs(mol)  # Add hydrogens
            AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol_3d)
            structure.mol_3d = mol_3d
            
            # Generate 3D SVG
            structure.svg_3d = self._generate_3d_svg(mol_3d)
            
            # Generate interactive 3D HTML
            structure.html_3d = self._generate_3d_html(mol_3d)
            
            # Generate XYZ coordinates
            structure.xyz_coords = self._generate_xyz_coordinates(mol_3d)
            
        except Exception as e:
            structure.error_message = f"Error processing SMILES: {str(e)}"
        
        return structure
    
    def _generate_2d_svg(self, mol) -> str:
        """Generate 2D SVG representation."""
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    
    def _generate_3d_svg(self, mol) -> str:
        """Generate 3D SVG representation."""
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    
    def _generate_3d_html(self, mol) -> str:
        """Generate interactive 3D HTML using NGL Viewer."""
        try:
            # Convert to PDB format for 3D visualization
            pdb_content = Chem.MolToPDBBlock(mol)
            
            # Create interactive HTML with NGL Viewer
            html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>3D Molecular Structure</title>
    <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>
    <style>
        body { margin: 0; padding: 0; }
        #viewport { width: 100%; height: 400px; }
    </style>
</head>
<body>
    <div id="viewport"></div>
    <script>
        var stage = new NGL.Stage("viewport");
        stage.setParameters({
            backgroundColor: "white"
        });
        
        var pdbData = `{pdb_content}`;
        
        stage.loadFile(new Blob([pdbData], {{type: "text/plain"}}), {{ext: "pdb"}}).then(function (component) {{
            component.addRepresentation("ball+stick", {{
                color: "element"
            }});
            component.autoView();
        }});
    </script>
</body>
</html>
"""
            return html_template.format(pdb_content=pdb_content)
            
        except Exception as e:
            # Fallback to simple 3D SVG if HTML generation fails
            return self._generate_3d_svg(mol)
    
    def _generate_xyz_coordinates(self, mol) -> str:
        """Generate XYZ coordinate file content."""
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        # Get molecular formula for title
        formula = rdMolDescriptors.CalcMolFormula(mol)
        
        xyz_content = f"{num_atoms}\n{formula}\n"
        
        for i in range(num_atoms):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            xyz_content += f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n"
        
        return xyz_content
    
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
    """Generate 2D SVG model from SMILES."""
    modeler = SmilesModeler()
    structure = modeler.process_smiles(smiles)
    return structure.svg_2d if not structure.error_message else ""

def generate_3d_model(smiles: str) -> str:
    """Generate 3D SVG model from SMILES."""
    modeler = SmilesModeler()
    structure = modeler.process_smiles(smiles)
    return structure.svg_3d if not structure.error_message else ""

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
    print(f"2D SVG length: {len(result.svg_2d)}")
    print(f"3D SVG length: {len(result.svg_3d)}")
    print(f"3D HTML length: {len(result.html_3d)}")
    print(f"XYZ coords:\n{result.xyz_coords}")
