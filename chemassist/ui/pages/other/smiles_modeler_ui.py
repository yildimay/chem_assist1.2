from __future__ import annotations

"""SMILES to 3D Modeler UI - Convert SMILES strings to molecular structures.

Features:
    • 2D molecular visualization (inline display)
    • 3D molecular visualization (inline display + interactive viewer)
    • XYZ coordinate file generation
    • Download functionality for essential formats only
"""

import io
import streamlit as st
from typing import Optional

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from chemassist.core.other.smiles_modeler import SmilesModeler, MolecularStructure

# ──────────────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────────────

def _create_download_button(content: str, filename: str, mime_type: str, button_text: str) -> None:
    """Create a download button for file content."""
    if content:
        st.download_button(
            label=button_text,
            data=content,
            file_name=filename,
            mime=mime_type,
            key=f"download_{filename}"
        )

def _display_svg(svg_content: str, title: str) -> None:
    """Display SVG content in Streamlit."""
    if svg_content:
        st.markdown(f"**{title}**")
        st.markdown(svg_content, unsafe_allow_html=True)
    else:
        st.warning(f"No {title.lower()} available")

def _display_3d_interactive(html_content: str, title: str) -> None:
    """Display interactive 3D HTML content in Streamlit."""
    if html_content:
        st.markdown(f"**{title}**")
        st.components.html(html_content, height=450, scrolling=False)
    else:
        st.warning(f"No interactive {title.lower()} available")

# ──────────────────────────────────────────────────────────────────────────────
# Main page
# ──────────────────────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🧬 SMILES → 3D Modeler")
    st.markdown(
        "Convert SMILES strings into molecular structures. Generate 2D and 3D "
        "visualizations, plus XYZ coordinate files for computational chemistry."
    )
    
    # Initialize session state
    if "smiles_structure" not in st.session_state:
        st.session_state.smiles_structure = None
    
    # Input section
    st.subheader("📝 Input")
    col1, col2 = st.columns([3, 1])
    
    with col1:
        smiles_input = st.text_input(
            "SMILES String",
            placeholder="C1=CC=CC=C1 (benzene), CC(=O)O (acetic acid), etc.",
            help="Enter a valid SMILES string to generate molecular structures"
        )
    
    with col2:
        process_button = st.button("Generate Structures", type="primary")
    
    # Example SMILES
    with st.expander("💡 Example SMILES"):
        st.markdown("""
        **Common molecules:**
        - `C1=CC=CC=C1` - Benzene
        - `CC(=O)O` - Acetic acid
        - `CCO` - Ethanol
        - `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` - Ibuprofen
        - `CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C` - Testosterone
        
        **Tips:**
        - Use standard SMILES notation
        - Include hydrogens implicitly (RDKit will add them for 3D)
        - Check the validation message below
        """)
    
    # Process SMILES
    if process_button and smiles_input:
        try:
            modeler = SmilesModeler()
            
            # Validate SMILES
            is_valid, validation_msg = modeler.validate_smiles(smiles_input)
            
            if is_valid:
                with st.spinner("Generating molecular structures..."):
                    structure = modeler.process_smiles(smiles_input)
                    st.session_state.smiles_structure = structure
                
                if structure.error_message:
                    st.error(f"Processing error: {structure.error_message}")
                else:
                    st.success("✅ Structures generated successfully!")
            else:
                st.error(f"❌ {validation_msg}")
                
        except ImportError:
            st.error("❌ RDKit is not available. Please install RDKit to use this feature.")
        except Exception as e:
            st.error(f"❌ Unexpected error: {str(e)}")
    
    # Display results
    if st.session_state.smiles_structure and not st.session_state.smiles_structure.error_message:
        structure = st.session_state.smiles_structure
        
        # Molecular info
        st.subheader("📊 Molecular Information")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("SMILES", structure.smiles)
        
        with col2:
            if structure.mol_2d:
                num_atoms = structure.mol_2d.GetNumAtoms()
                st.metric("Atoms", num_atoms)
        
        with col3:
            if structure.mol_2d:
                num_bonds = structure.mol_2d.GetNumBonds()
                st.metric("Bonds", num_bonds)
        
        # Three sections for different outputs
        tab_2d, tab_3d, tab_xyz = st.tabs(["🖼️ 2D Structure", "🎯 3D Structure", "📄 XYZ Coordinates"])
        
        # Tab 1: 2D Structure
        with tab_2d:
            st.markdown("**2D molecular visualization**")
            
            if structure.svg_2d:
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    _display_svg(structure.svg_2d, "2D Molecular Structure")
                
                with col2:
                    st.markdown("**Download Options**")
                    _create_download_button(
                        structure.svg_2d,
                        "molecule_2d.svg",
                        "image/svg+xml",
                        "📥 Download 2D SVG"
                    )
            else:
                st.warning("No 2D structure available")
        
        # Tab 2: 3D Structure
        with tab_3d:
            st.markdown("**3D molecular visualization**")
            
            if structure.html_3d:
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    _display_3d_interactive(structure.html_3d, "Interactive 3D Molecular Structure")
                    st.markdown("*💡 **Tip:** Drag to rotate, scroll to zoom, right-click for more options*")
                
                with col2:
                    st.markdown("**Download Options**")
                    _create_download_button(
                        structure.svg_3d,
                        "molecule_3d.svg",
                        "image/svg+xml",
                        "📥 Download 3D SVG"
                    )
                    
                    # Also provide PDB file for 3D structure
                    if structure.mol_3d:
                        pdb_content = Chem.MolToPDBBlock(structure.mol_3d)
                        _create_download_button(
                            pdb_content,
                            "molecule_3d.pdb",
                            "text/plain",
                            "📥 Download PDB File"
                        )
            else:
                st.warning("No 3D structure available")
        
        # Tab 3: XYZ Coordinates
        with tab_xyz:
            st.markdown("**XYZ coordinate file**")
            
            if structure.xyz_coords:
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    st.markdown("**XYZ File Content**")
                    st.code(structure.xyz_coords, language="text")
                
                with col2:
                    st.markdown("**Download Options**")
                    _create_download_button(
                        structure.xyz_coords,
                        "molecule.xyz",
                        "text/plain",
                        "📥 Download XYZ File"
                    )
            else:
                st.warning("No XYZ coordinates available")
        
        # Additional information
        with st.expander("ℹ️ About the generated structures"):
            st.markdown("""
            **2D Structure:**
            - Planar representation of the molecule
            - Shows connectivity and stereochemistry
            - Suitable for publication and documentation
            
            **3D Structure:**
            - Interactive three-dimensional visualization
            - Generated using RDKit's 3D embedding
            - Optimized using MMFF force field
            - Includes all hydrogen atoms
            - **Interactive features:** Rotate, zoom, and explore the molecule
            
            **XYZ Coordinates:**
            - Cartesian coordinates for all atoms
            - Compatible with most computational chemistry software
            - Can be used for further calculations
            """)

def _create_gaussian_input(structure: MolecularStructure) -> str:
    """Create a basic Gaussian input file from the structure."""
    if not structure.xyz_coords:
        return ""
    
    # Extract number of atoms and coordinates
    lines = structure.xyz_coords.strip().split('\n')
    if len(lines) < 3:
        return ""
    
    num_atoms = int(lines[0])
    title = lines[1]
    coords = '\n'.join(lines[2:2+num_atoms])
    
    # Basic Gaussian input template
    gjf_content = f"""%NProcShared=4
%Mem=4GB
# B3LYP/6-31G(d) opt

{title}

0 1
{coords}

"""
    return gjf_content
