#!/usr/bin/env python3
"""Test Streamlit components functionality."""

import streamlit as st

def test_streamlit_components():
    """Test if Streamlit components work correctly."""
    
    st.title("🧪 Streamlit Components Test")
    
    # Test basic HTML component
    st.subheader("HTML Component Test")
    
    test_html = """
    <div style="background-color: #f0f0f0; padding: 20px; border-radius: 10px;">
        <h3>Test HTML Content</h3>
        <p>This is a test of the HTML component functionality.</p>
        <p>If you can see this styled content, the components are working!</p>
    </div>
    """
    
    try:
        st.components.html(test_html, height=200, scrolling=False)
        st.success("✅ HTML component working!")
    except Exception as e:
        st.error(f"❌ HTML component failed: {e}")
    
    # Test if we can import the SMILES modeler UI
    st.subheader("SMILES Modeler UI Test")
    
    try:
        from chemassist.ui.pages.other.smiles_modeler_ui import show_page
        st.success("✅ SMILES modeler UI imports successfully!")
    except Exception as e:
        st.error(f"❌ SMILES modeler UI import failed: {e}")
    
    # Test environment variables
    st.subheader("Environment Variables Test")
    
    groq_key = st.secrets.get("GROQ_API_KEY", "NOT_FOUND")
    openai_key = st.secrets.get("OPENAI_API_KEY", "NOT_FOUND")
    
    st.info(f"GROQ_API_KEY: {'Found' if groq_key != 'NOT_FOUND' else 'Not found'}")
    st.info(f"OPENAI_API_KEY: {'Found' if openai_key != 'NOT_FOUND' else 'Not found'}")
    
    if groq_key == "NOT_FOUND" and openai_key == "NOT_FOUND":
        st.warning("⚠️ No API keys found in secrets. Add them to .streamlit/secrets.toml for full functionality.")

if __name__ == "__main__":
    test_streamlit_components() 