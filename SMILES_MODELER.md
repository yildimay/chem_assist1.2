# SMILES to 3D Modeler

## Overview

The SMILES to 3D Modeler is a powerful tool that converts SMILES (Simplified Molecular Input Line Entry System) strings into various molecular structure formats. This tool is designed to make molecular modeling easier by providing multiple output formats from a single SMILES input.

## Features

### 🖼️ 2D Structure Visualization
- **Planar molecular representation**
- **SVG format** for high-quality graphics
- **Shows connectivity and stereochemistry**
- **Suitable for publications and documentation**

### 🎯 3D Structure Visualization  
- **Three-dimensional molecular conformation**
- **Generated using RDKit's 3D embedding**
- **Optimized using MMFF force field**
- **Includes all hydrogen atoms**
- **SVG format for visualization**

### 📄 XYZ Coordinate Files
- **Cartesian coordinates for all atoms**
- **Compatible with most computational chemistry software**
- **Can be used for further calculations**
- **Includes molecular formula in header**

## Usage

### Input
1. Enter a valid SMILES string in the input field
2. Click "Generate Structures" button
3. Wait for processing to complete

### Output Sections

#### 1. Molecular Information
- **SMILES string** - The input SMILES
- **Number of atoms** - Total atom count
- **Number of bonds** - Total bond count

#### 2. 2D Structure Tab
- **2D molecular visualization** in SVG format
- **Download options:**
  - 2D SVG file
  - SMILES text file

#### 3. 3D Structure Tab  
- **3D molecular visualization** in SVG format
- **Download options:**
  - 3D SVG file
  - 3D SMILES text file

#### 4. XYZ Coordinates Tab
- **XYZ coordinate file content** displayed as text
- **Download options:**
  - XYZ coordinate file (.xyz)
  - Gaussian input file (.gjf)

## Example SMILES

### Simple Molecules
- `C1=CC=CC=C1` - Benzene
- `CC(=O)O` - Acetic acid
- `CCO` - Ethanol

### Complex Molecules
- `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` - Ibuprofen
- `CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C` - Testosterone

## Technical Details

### Dependencies
- **RDKit** - Core molecular modeling library
- **Streamlit** - Web interface framework

### Processing Pipeline
1. **SMILES Validation** - Check if input is valid
2. **2D Structure Generation** - Create planar representation
3. **3D Structure Generation** - Add hydrogens and embed in 3D
4. **Force Field Optimization** - MMFF optimization for realistic geometry
5. **Coordinate Extraction** - Generate XYZ coordinates
6. **SVG Rendering** - Create visual representations

### File Formats

#### XYZ Format
```
<number_of_atoms>
<molecular_formula>
<atom_symbol> <x> <y> <z>
...
```

#### Gaussian Input Format
```
%NProcShared=4
%Mem=4GB
# B3LYP/6-31G(d) opt

<molecular_formula>

0 1
<xyz_coordinates>
```

## Error Handling

The tool includes comprehensive error handling for:
- **Invalid SMILES strings**
- **RDKit import errors**
- **Processing failures**
- **Memory issues with large molecules**

## Tips for Best Results

1. **Use standard SMILES notation**
2. **Include hydrogens implicitly** (RDKit will add them for 3D)
3. **Check validation messages** for input errors
4. **For large molecules**, processing may take longer
5. **Download files** for use in other software

## Integration

This tool integrates seamlessly with the ChemAssist framework and can be used in combination with:
- **DFT Input Creator** - Use generated XYZ coordinates
- **Error Fixer** - Process problematic molecular structures
- **Other computational chemistry tools**

## Testing

Run the test script to verify functionality:
```bash
python3 test_smiles_modeler.py
```

This will test various molecules and display processing results. 