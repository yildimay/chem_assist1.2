from __future__ import annotations

"""DFT • Chat-based Error Fixer UI

A chat interface for Gaussian error fixing that uses the LLM router.
Based on the original gaussian_fixer_ui.py but with proper API handling.
"""

import streamlit as st
from PIL import Image
import pytesseract
from typing import Optional

from chemassist.models.llm_router import call_llm

def show_page() -> None:  # noqa: D401
    st.header("🧠 Gaussian Error Fixer Chat")

    st.markdown("""
    Welcome to the **Gaussian Error Fixer Chat**.
    - Describe your Gaussian issue in the chatbox below.
    - Optionally, upload a screenshot of the error.
    - Optionally, upload Gaussian input (`.gjf`) and/or output (`.log` / `.out`) files.

    **Note**: Text input is required. Other inputs are optional.
    """)

    # Initialize chat history
    if "chat_messages" not in st.session_state:
        st.session_state.chat_messages = []

    # Display chat history
    for msg in st.session_state.chat_messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    # User input and uploads
    with st.form("chat_form"):
        user_input = st.text_area(
            "Describe your Gaussian error:", 
            height=100, 
            placeholder="e.g. Link 9999 or SCF not converging..."
        )
        uploaded_image = st.file_uploader(
            "Optional: Upload error screenshot", 
            type=["png", "jpg", "jpeg"]
        )
        uploaded_input_file = st.file_uploader(
            "Optional: Upload Gaussian Input File (.gjf)", 
            type=["gjf"]
        )
        uploaded_output_file = st.file_uploader(
            "Optional: Upload Gaussian Output File (.log, .out)", 
            type=["log", "out"]
        )
        submitted = st.form_submit_button("Send")

    if submitted:
        if user_input.strip() == "":
            st.warning("Text input is required.")
        else:
            # Save user message
            st.session_state.chat_messages.append({"role": "user", "content": user_input})
            with st.chat_message("user"):
                st.markdown(user_input)

            # Process image if available
            image_text = ""
            if uploaded_image is not None:
                try:
                    image = Image.open(uploaded_image)
                    image_text = pytesseract.image_to_string(image)
                    st.markdown("*Extracted from image:* \n" + image_text)
                except Exception as e:
                    st.error(f"Failed to process image: {e}")

            # Read file contents
            input_file_text = ""
            output_file_text = ""
            
            if uploaded_input_file is not None:
                try:
                    input_file_text = uploaded_input_file.read().decode("utf-8", errors="ignore")
                except Exception as e:
                    st.error(f"Failed to read input file: {e}")
                    
            if uploaded_output_file is not None:
                try:
                    output_file_text = uploaded_output_file.read().decode("utf-8", errors="ignore")
                except Exception as e:
                    st.error(f"Failed to read output file: {e}")

            # Combine and process all inputs
            final_input = f"""USER MESSAGE:
{user_input}

IMAGE TEXT:
{image_text}

INPUT FILE:
{input_file_text}

OUTPUT FILE:
{output_file_text}"""

            # Call LLM using the router
            try:
                ai_response = call_llm(
                    messages=[{"role": "user", "content": final_input}],
                    system_prompt="""You are a helpful assistant for Gaussian error fixing. 
                    
                    Your role is to:
                    1. Analyze Gaussian error messages and input files
                    2. Provide clear explanations of what went wrong
                    3. Suggest specific fixes and solutions
                    4. Help users understand computational chemistry concepts
                    
                    If the question is not related to Gaussian software or Gaussian errors, 
                    reply with: 'This tool is specifically for Gaussian error fixing. 
                    Please ask a Gaussian-related question.'
                    
                    Be helpful, clear, and provide actionable advice."""
                )
            except Exception as e:
                ai_response = f"❌ Error while contacting AI service: {e}"
                st.error("Make sure you have set up your API key correctly.")

            # Save and display AI response
            st.session_state.chat_messages.append({"role": "assistant", "content": ai_response})
            with st.chat_message("assistant"):
                st.markdown(ai_response)

    # Clear chat button
    if st.session_state.chat_messages:
        if st.button("🗑️ Clear Chat History"):
            st.session_state.chat_messages = []
            st.rerun() 