#!/usr/bin/env python3
"""Test current setup to verify all fixes work properly."""

import os
import sys

def test_imports():
    """Test all critical imports."""
    print("🧪 Testing imports...")
    
    try:
        import streamlit as st
        print(f"✅ Streamlit: {st.__version__}")
    except Exception as e:
        print(f"❌ Streamlit import failed: {e}")
        return False
    
    try:
        import groq
        print(f"✅ Groq: {groq.__version__}")
    except Exception as e:
        print(f"❌ Groq import failed: {e}")
        return False
    
    try:
        from chemassist.models.llm_router import call_llm
        print("✅ LLM router import successful")
    except Exception as e:
        print(f"❌ LLM router import failed: {e}")
        return False
    
    try:
        from chemassist.ui.pages.other.smiles_modeler_ui import show_page
        print("✅ SMILES modeler UI import successful")
    except Exception as e:
        print(f"❌ SMILES modeler UI import failed: {e}")
        return False
    
    return True

def test_groq_api():
    """Test Groq API functionality."""
    print("\n🔑 Testing Groq API...")
    
    groq_key = os.environ.get("GROQ_API_KEY")
    if not groq_key:
        print("⚠️ No GROQ_API_KEY found in environment")
        return False
    
    try:
        import groq
        client = groq.Groq(api_key=groq_key)
        print("✅ Groq client created successfully")
        
        # Test a simple API call
        response = client.chat.completions.create(
            model="llama3-70b-8192",
            messages=[{"role": "user", "content": "Say 'Hello from ChemAssist!'"}],
            max_tokens=10
        )
        print(f"✅ Groq API call successful: {response.choices[0].message.content}")
        return True
        
    except Exception as e:
        print(f"❌ Groq API test failed: {e}")
        return False

def test_smiles_modeler():
    """Test SMILES modeler functionality."""
    print("\n🧬 Testing SMILES modeler...")
    
    try:
        from chemassist.core.other.smiles_modeler import SmilesModeler
        
        modeler = SmilesModeler()
        result = modeler.process_smiles("C1=CC=CC=C1")  # Benzene
        
        if result.error_message:
            print(f"❌ SMILES modeler failed: {result.error_message}")
            return False
        
        print(f"✅ SMILES modeler successful:")
        print(f"   - 2D SVG: {len(result.svg_2d)} characters")
        print(f"   - 3D SVG: {len(result.svg_3d)} characters")
        print(f"   - 3D HTML: {len(result.html_3d)} characters")
        print(f"   - XYZ coords: {len(result.xyz_coords.split(chr(10)))} lines")
        
        return True
        
    except Exception as e:
        print(f"❌ SMILES modeler test failed: {e}")
        return False

def test_llm_router():
    """Test LLM router functionality."""
    print("\n🤖 Testing LLM router...")
    
    groq_key = os.environ.get("GROQ_API_KEY")
    if not groq_key:
        print("⚠️ No GROQ_API_KEY found - skipping LLM router test")
        return True
    
    try:
        from chemassist.models.llm_router import call_llm
        
        response = call_llm(
            messages=[{"role": "user", "content": "Say 'Hello from ChemAssist!'"}],
            model="llama3-70b-8192"
        )
        
        print(f"✅ LLM router successful: {response}")
        return True
        
    except Exception as e:
        print(f"❌ LLM router test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🔧 ChemAssist Setup Test")
    print("=" * 50)
    
    # Test imports
    if not test_imports():
        print("\n❌ Import tests failed. Please check your installation.")
        return
    
    # Test SMILES modeler
    if not test_smiles_modeler():
        print("\n❌ SMILES modeler test failed.")
        return
    
    # Test Groq API
    if not test_groq_api():
        print("\n⚠️ Groq API test failed. Check your API key.")
    
    # Test LLM router
    if not test_llm_router():
        print("\n⚠️ LLM router test failed. Check your API key.")
    
    print("\n" + "=" * 50)
    print("✅ All tests completed!")
    print("\nIf you're still having issues:")
    print("1. Make sure you've pulled the latest changes from GitHub")
    print("2. Restart your deployment/application")
    print("3. Check that your API key is properly set")

if __name__ == "__main__":
    main() 