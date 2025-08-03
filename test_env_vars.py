#!/usr/bin/env python3
"""Test environment variable handling for Streamlit deployment."""

import os
import streamlit as st

def test_environment_variables():
    """Test how environment variables are accessed."""
    
    st.title("🔧 Environment Variable Test")
    
    # Test different ways to access environment variables
    st.subheader("Environment Variable Access Methods")
    
    # Method 1: Direct os.environ
    groq_key_direct = os.environ.get("GROQ_API_KEY", "NOT_FOUND")
    openai_key_direct = os.environ.get("OPENAI_API_KEY", "NOT_FOUND")
    
    # Method 2: Streamlit secrets (for Streamlit Cloud)
    try:
        groq_key_secrets = st.secrets.get("GROQ_API_KEY", "NOT_FOUND")
        openai_key_secrets = st.secrets.get("OPENAI_API_KEY", "NOT_FOUND")
    except Exception as e:
        groq_key_secrets = f"ERROR: {e}"
        openai_key_secrets = f"ERROR: {e}"
    
    # Display results
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Direct os.environ:**")
        st.code(f"GROQ_API_KEY: {groq_key_direct[:10]}..." if len(groq_key_direct) > 10 else f"GROQ_API_KEY: {groq_key_direct}")
        st.code(f"OPENAI_API_KEY: {openai_key_direct[:10]}..." if len(openai_key_direct) > 10 else f"OPENAI_API_KEY: {openai_key_direct}")
    
    with col2:
        st.markdown("**Streamlit secrets:**")
        st.code(f"GROQ_API_KEY: {groq_key_secrets[:10]}..." if len(str(groq_key_secrets)) > 10 else f"GROQ_API_KEY: {groq_key_secrets}")
        st.code(f"OPENAI_API_KEY: {openai_key_secrets[:10]}..." if len(str(openai_key_secrets)) > 10 else f"OPENAI_API_KEY: {openai_key_secrets}")
    
    # Test LLM functionality
    st.subheader("LLM Functionality Test")
    
    if groq_key_direct != "NOT_FOUND" or groq_key_secrets != "NOT_FOUND":
        st.success("✅ API key found!")
        
        # Test the actual LLM call
        try:
            from chemassist.models.llm_router import call_llm
            
            response = call_llm(
                messages=[{"role": "user", "content": "Say 'Hello from ChemAssist!'"}],
                model="llama3-70b-8192"
            )
            st.success("✅ LLM call successful!")
            st.info(f"Response: {response}")
            
        except Exception as e:
            st.error(f"❌ LLM call failed: {e}")
    else:
        st.warning("⚠️ No API key found. Set GROQ_API_KEY or OPENAI_API_KEY in Streamlit secrets.")
    
    # Show all environment variables (for debugging)
    with st.expander("🔍 All Environment Variables"):
        env_vars = {k: v for k, v in os.environ.items() if "API" in k or "KEY" in k}
        if env_vars:
            for key, value in env_vars.items():
                display_value = value[:10] + "..." if len(value) > 10 else value
                st.code(f"{key}: {display_value}")
        else:
            st.info("No API-related environment variables found.")

if __name__ == "__main__":
    test_environment_variables() 