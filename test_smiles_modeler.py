#!/usr/bin/env python3
"""Test script for SMILES to 3D Modeler functionality."""

from chemassist.core.other.smiles_modeler import SmilesModeler

def test_smiles_modeler():
    """Test the SMILES modeler with various molecules."""
    
    # Test molecules
    test_molecules = [
        ("C1=CC=CC=C1", "Benzene"),
        ("CC(=O)O", "Acetic acid"),
        ("CCO", "Ethanol"),
        ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "Ibuprofen"),
        ("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "Testosterone")
    ]
    
    modeler = SmilesModeler()
    
    print("🧬 SMILES to 3D Modeler Test")
    print("=" * 50)
    
    for smiles, name in test_molecules:
        print(f"\n📋 Testing: {name}")
        print(f"SMILES: {smiles}")
        
        try:
            # Validate SMILES
            is_valid, validation_msg = modeler.validate_smiles(smiles)
            print(f"Validation: {'✅' if is_valid else '❌'} {validation_msg}")
            
            if is_valid:
                # Process SMILES
                structure = modeler.process_smiles(smiles)
                
                if structure.error_message:
                    print(f"❌ Error: {structure.error_message}")
                else:
                    print(f"✅ Success!")
                    print(f"   - 2D SVG: {len(structure.svg_2d)} characters")
                    print(f"   - 3D SVG: {len(structure.svg_3d)} characters")
                    print(f"   - XYZ coords: {len(structure.xyz_coords.split(chr(10)))} lines")
                    
                    # Show first few lines of XYZ coordinates
                    xyz_lines = structure.xyz_coords.split('\n')
                    if len(xyz_lines) >= 3:
                        print(f"   - XYZ preview: {xyz_lines[0]} atoms, {xyz_lines[1]}")
                        if len(xyz_lines) > 3:
                            print(f"     First atom: {xyz_lines[2]}")
            
        except Exception as e:
            print(f"❌ Exception: {e}")
    
    print("\n" + "=" * 50)
    print("✅ Test completed!")

if __name__ == "__main__":
    test_smiles_modeler() 