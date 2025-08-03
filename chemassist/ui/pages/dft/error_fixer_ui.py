from __future__ import annotations

"""DFT • Error Fixer UI

Upload:
  • original input file (e.g., .gjf, .com)
  • matching .log / .out tail (or full file)

Or use the chat interface for interactive error fixing.
"""

import streamlit as st
import requests
import os
from PIL import Image
import pytesseract

from chemassist.core.dft.error_fixer import fix_input

# ───────────────────────────────────────────────────────────────
# Streamlit page
# ───────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🚑 DFT • Error Fixer")
    
    # Mode selection
    mode = st.radio("Choose Mode", ["File Upload", "Chat Interface"], horizontal=True)
    
    if mode == "File Upload":
        show_file_upload_mode()
    else:
        show_chat_mode()

def show_file_upload_mode():
    """Show the original file upload mode."""
    st.write("Upload your broken Gaussian (or other) input and the last part of the log file; ChemAssist will diagnose and propose a fix.")

    with st.expander("Instructions", expanded=False):
        st.markdown(
            "* **Input file** – original .gjf/.com/.inp as you submitted.  \n"
            "* **Log tail** – drag in the full .log/.out or just the last ~300 lines.  \n"
            "* The tool sends only those texts (no filenames) to the LLM backend.  \n"
            "* Supports Gaussian for now; ORCA/NWChem coming soon."
        )

    col_inp, col_log = st.columns(2)
    with col_inp:
        input_file = st.file_uploader("Input file (.gjf, .com, .inp)", type=["gjf", "com", "inp", "txt"])
    with col_log:
        log_file = st.file_uploader("Log / out file", type=["log", "out", "txt"])

    program = st.selectbox("Program", ["Gaussian"], index=0, disabled=True)

    col1, col2 = st.columns([1, 1])
    with col1:
        diagnose_button = st.button("🩹 Diagnose & Fix", disabled=not (input_file and log_file))
    with col2:
        retry_button = st.button("🔄 Retry with different format", disabled=not (input_file and log_file))
    
    if diagnose_button or retry_button:
        with st.spinner("🤖 AI is analyzing your error..."):
            try:
                # Reset file pointers
                input_file.seek(0)
                log_file.seek(0)
                
                result = fix_input(
                    input_text=input_file.read().decode("utf-8", "replace"),
                    log_text=log_file.read().decode("utf-8", "replace"),
                    program="gaussian",
                )
                
                st.subheader("✅ Diagnosis")
                st.info(result["diagnosis"] or "(LLM returned no explanation)")

                st.subheader("✅ Corrected input")
                st.code(result["fixed_input"], language="text")
                st.download_button(
                    label="💾 Download fixed .gjf",
                    data=result["fixed_input"].encode(),
                    file_name="fixed_input.gjf",
                    mime="text/plain",
                )
                
                st.success("🎉 Error analysis complete!")
                
            except Exception as exc:  # noqa: BLE001
                st.error(f"❌ Failed to fix input: {exc}")
                
                # Show helpful troubleshooting tips
                with st.expander("🔧 Troubleshooting Tips", expanded=True):
                    st.markdown("""
                    **Common issues and solutions:**
                    
                    1. **"LLM did not return a FIXED_INPUT block"**
                       - Try the "Retry with different format" button
                       - Check that your input file is a valid Gaussian input
                       - Make sure your log file contains error messages
                    
                    2. **"Invalid API Key"**
                       - Make sure you have set your GROQ_API_KEY environment variable
                       - Get a free API key from: https://console.groq.com/
                    
                    3. **"No response from LLM"**
                       - Check your internet connection
                       - Try again in a few moments
                    
                    4. **Still having issues?**
                       - Try uploading a smaller log file (just the error part)
                       - Make sure your input file is complete
                    """)

def show_chat_mode():
    """Show the chat interface mode."""
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

            # Call Groq API directly like the original code
            try:
                # Get API key from environment or Streamlit secrets
                api_key = os.getenv('GROQ_API_KEY')
                if not api_key:
                    try:
                        api_key = st.secrets.get("GROQ_API_KEY")
                    except:
                        pass
                
                if not api_key:
                    ai_response = "❌ No API key found. Please set GROQ_API_KEY environment variable or add it to Streamlit secrets."
                else:
                    headers = {
                        "Authorization": f"Bearer {api_key}",
                        "Content-Type": "application/json"
                    }

                    data = {
                        "model": "llama3-70b-8192",
                        "messages": [
                            {"role": "system", "content": "You are a helpful assistant for Gaussian error fixing. If the question is not related to Gaussian software or Gaussian errors, reply with: 'This tool is specifically for Gaussian error fixing. Please ask a Gaussian-related question.'"},
                            {"role": "user", "content": final_input}
                        ]
                    }

                    response = requests.post("https://api.groq.com/openai/v1/chat/completions", headers=headers, json=data)
                    response.raise_for_status()
                    ai_response = response.json()['choices'][0]['message']['content']
                    
            except Exception as e:
                ai_response = f"❌ Error while contacting Groq API: {e}"

            # Save and display AI response
            st.session_state.chat_messages.append({"role": "assistant", "content": ai_response})
            with st.chat_message("assistant"):
                st.markdown(ai_response)

    # Clear chat button
    if st.session_state.chat_messages:
        if st.button("🗑️ Clear Chat History"):
            st.session_state.chat_messages = []
            st.rerun()
