from __future__ import annotations

"""MD • Error Fixer UI

Upload:
  • original .mdp file
  • matching .log / .out file

The page calls :pyfunc:`chemassist.core.md.error_fixer.fix_md_input` which
invokes the LLM router under the hood.  On success it displays the
problem diagnosis and offers the corrected .mdp for download.
"""

import streamlit as st
from chemassist.core.md.error_fixer import fix_md_input, analyze_md_log
from chemassist.utils.file_io import create_md_error_fix_zip_archive

# ───────────────────────────────────────────────────────────────
# Streamlit page
# ───────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🚑 MD • Error Fixer")
    st.write("Upload your broken GROMACS .mdp file and the corresponding .log file; ChemAssist will diagnose and propose a fix.")

    with st.expander("📖 Instructions", expanded=False):
        st.markdown(
            """
            **This tool helps with:**
            - GROMACS parameter errors
            - Simulation crashes
            - Performance issues
            - File format problems
            
            **Upload files:**
            - **MDP file** – your .mdp parameter file
            - **Log file** – the .log/.out file with error messages
            
            **Supported stages:**
            - Energy minimization (em)
            - NVT equilibration (nvt)
            - NPT equilibration (npt)
            - Production MD (md)
            
            **The tool will:**
            1. Analyze the error in your log file
            2. Identify the problem in your .mdp file
            3. Generate a corrected .mdp file
            4. Provide detailed diagnosis and suggestions
            """
        )

    # ───────────────────────────────────────────────────────────────
    # File upload section
    # ───────────────────────────────────────────────────────────────

    st.subheader("📁 Upload Files")
    
    col_mdp, col_log = st.columns(2)
    with col_mdp:
        mdp_file = st.file_uploader("MDP file (.mdp)", type=["mdp", "txt"])
    with col_log:
        log_file = st.file_uploader("Log file (.log, .out)", type=["log", "out", "txt"])

    # ───────────────────────────────────────────────────────────────
    # Analysis options
    # ───────────────────────────────────────────────────────────────

    st.subheader("🔧 Analysis Options")
    
    col1, col2 = st.columns(2)
    with col1:
        stage = st.selectbox(
            "Simulation Stage",
            ["md", "em", "nvt", "npt"],
            index=0,
            help="Select the stage that failed"
        )
    
    with col2:
        analysis_type = st.radio(
            "Analysis Type",
            ["Fix MDP File", "Analyze Log Only"],
            index=0,
            help="Choose whether to fix the MDP file or just analyze the log"
        )

    # ───────────────────────────────────────────────────────────────
    # Analysis buttons
    # ───────────────────────────────────────────────────────────────

    col1, col2 = st.columns([1, 1])
    with col1:
        fix_button = st.button(
            "🩹 Fix MDP File", 
            disabled=not (mdp_file and log_file),
            help="Analyze and fix the MDP file"
        )
    with col2:
        analyze_button = st.button(
            "📊 Analyze Log Only",
            disabled=not log_file,
            help="Just analyze the log file for issues"
        )

    # ───────────────────────────────────────────────────────────────
    # Process analysis
    # ───────────────────────────────────────────────────────────────

    if fix_button and mdp_file and log_file:
        with st.spinner("🤖 AI is analyzing your GROMACS error..."):
            try:
                # Reset file pointers
                mdp_file.seek(0)
                log_file.seek(0)
                
                result = fix_md_input(
                    mdp_text=mdp_file.read().decode("utf-8", "replace"),
                    log_text=log_file.read().decode("utf-8", "replace"),
                    stage=stage,
                )
                
                st.success("✅ Analysis complete!")
                
                # Display results
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("🔍 Diagnosis")
                    st.info(result["diagnosis"] or "(No diagnosis provided)")
                
                with col2:
                    st.subheader("📊 Simulation Stage")
                    st.info(f"Analyzed as: **{stage.upper()}** stage")
                
                # Show corrected MDP
                st.subheader("✅ Corrected MDP File")
                st.code(result["fixed_mdp"], language="text")
                
                # Download buttons
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.download_button(
                        label="💾 Download fixed .mdp",
                        data=result["fixed_mdp"].encode(),
                        file_name=f"fixed_{stage}.mdp",
                        mime="text/plain",
                    )
                
                with col2:
                    # Create a simple comparison
                    mdp_file.seek(0)
                    original_mdp = mdp_file.read().decode("utf-8", "replace")
                    st.download_button(
                        label="📋 Download comparison",
                        data=f"ORIGINAL:\n{original_mdp}\n\nFIXED:\n{result['fixed_mdp']}".encode(),
                        file_name=f"{stage}_comparison.txt",
                        mime="text/plain",
                    )
                
                with col3:
                    try:
                        zip_data = create_md_error_fix_zip_archive(
                            original_mdp, result["fixed_mdp"], result["diagnosis"], stage
                        )
                        zip_filename = f"md_error_fix_{stage}.zip"
                        
                        st.download_button(
                            label="📦 Download ZIP",
                            data=zip_data,
                            file_name=zip_filename,
                            mime="application/zip",
                            help="Download all error fix files as ZIP"
                        )
                    except Exception as e:
                        st.error(f"Failed to create ZIP: {e}")
                
                st.success("🎉 MDP file fixed successfully!")
                
            except Exception as exc:  # noqa: BLE001
                st.error(f"❌ Failed to fix MDP file: {exc}")
                
                # Show troubleshooting tips
                with st.expander("🔧 Troubleshooting Tips", expanded=True):
                    st.markdown("""
                    **Common issues and solutions:**
                    
                    1. **"LLM did not return a FIXED_MDP block"**
                       - Try the "Analyze Log Only" option first
                       - Check that your .mdp file is valid GROMACS format
                       - Make sure your .log file contains error messages
                    
                    2. **"Invalid API Key"**
                       - Make sure you're running with: `./launch.sh YOUR_API_KEY`
                       - Get a free API key from: https://console.groq.com/
                    
                    3. **"No response from LLM"**
                       - Check your internet connection
                       - Try again in a few moments
                    
                    4. **Still having issues?**
                       - Try uploading a smaller log file (just the error part)
                       - Make sure your .mdp file is complete
                       - Check that the simulation stage is correct
                    """)

    elif analyze_button and log_file:
        with st.spinner("🤖 AI is analyzing your log file..."):
            try:
                # Reset file pointer
                log_file.seek(0)
                
                result = analyze_md_log(
                    log_text=log_file.read().decode("utf-8", "replace")
                )
                
                st.success("✅ Log analysis complete!")
                
                # Display results
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("📊 Analysis")
                    st.info(result["analysis"] or "(No analysis provided)")
                
                with col2:
                    st.subheader("💡 Suggestions")
                    st.info(result["suggestions"] or "(No suggestions provided)")
                
                # Show log preview
                log_file.seek(0)
                log_content = log_file.read().decode("utf-8", "replace")
                st.subheader("📄 Log File Preview (Last 50 lines)")
                st.code(log_content.splitlines()[-50:], language="text")
                
            except Exception as exc:  # noqa: BLE001
                st.error(f"❌ Failed to analyze log: {exc}")

    # ───────────────────────────────────────────────────────────────
    # Help section
    # ───────────────────────────────────────────────────────────────

    if not (mdp_file or log_file):
        st.info("👆 Upload your .mdp and .log files above to get started")
        
        with st.expander("📚 Common GROMACS Errors", expanded=False):
            st.markdown("""
            **Frequent GROMACS issues:**
            
            **Parameter Errors:**
            - Negative temperatures or pressures
            - Invalid cutoff values
            - Missing required parameters
            - Incompatible parameter combinations
            
            **File Format Issues:**
            - Missing semicolons in comments
            - Incorrect parameter names
            - Extra spaces or characters
            
            **Performance Problems:**
            - Memory allocation failures
            - Domain decomposition issues
            - Neighbor list problems
            
            **System Issues:**
            - Box size problems
            - Atom overlap
            - Missing topology information
            """)
