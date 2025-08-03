#!/usr/bin/env python3
"""
Test ZIP archive functionality
"""

from chemassist.core.md.input_creator import MDSpec, build_mdp, build_slurm
from chemassist.core.dft.input_creator import JobSpec, build_input
from chemassist.utils.file_io import create_md_zip_archive, create_dft_zip_archive

def test_md_zip():
    """Test MD ZIP archive creation."""
    print("🧪 Testing MD ZIP Archive...")
    
    spec = MDSpec(
        title="Test Protein MD",
        structure_file="protein.gro",
        topology_file="protein.top",
        output_prefix="test_md"
    )
    
    try:
        zip_data = create_md_zip_archive(spec, build_mdp, build_slurm)
        print(f"✅ MD ZIP created: {len(zip_data)} bytes")
        
        # Check if it's a valid ZIP
        import zipfile
        import io
        with zipfile.ZipFile(io.BytesIO(zip_data)) as zf:
            files = zf.namelist()
            print(f"   Contains files: {files}")
            
        return True
    except Exception as e:
        print(f"❌ MD ZIP failed: {e}")
        return False

def test_dft_zip():
    """Test DFT ZIP archive creation."""
    print("\n🧪 Testing DFT ZIP Archive...")
    
    spec = JobSpec(
        title="Test DFT Calculation",
        method="B3LYP",
        basis="6-31G(d)",
        charge=0,
        multiplicity=1,
        coords="C 0.0 0.0 0.0\nH 1.0 0.0 0.0"
    )
    
    try:
        zip_data = create_dft_zip_archive(spec, build_input)
        print(f"✅ DFT ZIP created: {len(zip_data)} bytes")
        
        # Check if it's a valid ZIP
        import zipfile
        import io
        with zipfile.ZipFile(io.BytesIO(zip_data)) as zf:
            files = zf.namelist()
            print(f"   Contains files: {files}")
            
        return True
    except Exception as e:
        print(f"❌ DFT ZIP failed: {e}")
        return False

def main():
    print("📦 Testing ZIP Archive Functionality")
    print("=" * 50)
    
    md_success = test_md_zip()
    dft_success = test_dft_zip()
    
    print("\n" + "=" * 50)
    if md_success and dft_success:
        print("🎉 All ZIP tests passed!")
    else:
        print("❌ Some ZIP tests failed")
    
    return md_success and dft_success

if __name__ == "__main__":
    main() 