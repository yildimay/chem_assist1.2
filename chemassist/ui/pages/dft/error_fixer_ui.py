from __future__ import annotations

"""DFT • Error Fixer UI

Upload:
  • original input file (e.g., .gjf, .com)
  • matching .log / .out tail (or full file)

The page calls :pyfunc:`chemassist.core.dft.error_fixer.fix_input` which
invokes the LLM router under the hood.  On success it displays the
problem diagnosis and offers the corrected input for download.
"""

import streamlit as st

from chemassist.core.dft.error_fixer import fix_input

# ───────────────────────────────────────────────────────────────
# Streamlit page
# ───────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🚑 DFT • Error Fixer")
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
                       - Make sure you're running with: `./launch.sh YOUR_API_KEY`
                       - Get a free API key from: https://console.groq.com/
                    
                    3. **"No response from LLM"**
                       - Check your internet connection
                       - Try again in a few moments
                    
                    4. **Still having issues?**
                       - Try uploading a smaller log file (just the error part)
                       - Make sure your input file is complete
                    """)
