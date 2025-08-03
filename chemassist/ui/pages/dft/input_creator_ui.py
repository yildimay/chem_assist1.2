from __future__ import annotations

"""DFT › Input Creator – two-path workflow.

* Tab 1  📑  PDF extraction → pick the most frequent method/basis.
* Tab 2  🧬  SMILES suggestion → fully automatic recommendation.

Both paths converge to `build_input()` which returns a .gjf ready for
Gaussian (other engines later).
"""

import io
import streamlit as st

from chemassist.core.dft.input_creator import JobSpec, build_input
from chemassist.core.dft.pdf_methods import extract as extract_methods
from chemassist.core.dft.method_suggester import recommend
from chemassist.utils.file_io import create_dft_zip_archive

# ──────────────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────────────


def _to_dataframe(hits):  # lazy import pandas only when needed
    import pandas as pd

    return pd.DataFrame([h.__dict__ for h in hits])


# ──────────────────────────────────────────────────────────────────────────────
# Main page
# ──────────────────────────────────────────────────────────────────────────────


def show_page() -> None:  # noqa: D401
    st.header("🧮 DFT • Input Creator")
    tab_pdf, tab_smiles = st.tabs(["📑 From research paper", "🧬 From SMILES"])

    # hold generated spec between tabs via session state
    if "draft_spec" not in st.session_state:
        st.session_state.draft_spec = None

    # ------------------------------------------------------------------
    # TAB 1 – PDF extraction
    # ------------------------------------------------------------------
    with tab_pdf:
        pdf_file = st.file_uploader("Upload article PDF", type="pdf")

        if st.button("Extract method/basis", disabled=not pdf_file):
            hits = extract_methods(pdf_file.read())
            if not hits:
                st.warning("Could not find any recognised functional/basis.")
            else:
                st.success("Candidates found – choose one below.")
                df = _to_dataframe(hits)
                idx = st.radio(
                    "Select combo",
                    df.index,
                    format_func=lambda i: f"{df.loc[i,'method']} / {df.loc[i,'basis']} ({df.loc[i,'count']})",
                )
                chosen_method = df.loc[idx, "method"]
                chosen_basis = df.loc[idx, "basis"]

                # Minimal placeholder coordinates; user edits later.
                coords_placeholder = "<paste coordinates here>"
                st.session_state.draft_spec = JobSpec(
                    title="From PDF extraction",
                    method=chosen_method,
                    basis=chosen_basis,
                    coords=coords_placeholder,
                )
                st.info("Scroll to *Finalize* section below to tweak and download.")

    # ------------------------------------------------------------------
    # TAB 2 – SMILES suggestion
    # ------------------------------------------------------------------
    with tab_smiles:
        smiles = st.text_input("SMILES", placeholder="C1=CC=CC=C1 (benzene)")
        if st.button("Suggest and build", disabled=not smiles):
            try:
                sug = recommend(smiles)
                st.success(
                    f"{sug.method} / {sug.basis} selected.  \n"
                    f"*Reason:* {sug.reason}"
                )
                st.session_state.draft_spec = JobSpec(
                    title=f"Auto-{smiles}",
                    charge=sug.charge,
                    multiplicity=sug.multiplicity,
                    method=sug.method,
                    basis=sug.basis,
                    coords=sug.xyz,
                )
                st.info("Scroll to *Finalize* section below to tweak and download.")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Failed to generate suggestion: {exc}")

    # ------------------------------------------------------------------
    # Finalisation (common)
    # ------------------------------------------------------------------
    st.subheader("Finalize & download")
    spec: JobSpec | None = st.session_state.draft_spec

    if spec is None:
        st.info("Use one of the tabs above to populate a draft input first.")
        return

    with st.form("final_form"):
        title = st.text_input("Job title", value=spec.title)
        col1, col2, col3 = st.columns(3)
        charge = col1.number_input("Charge", value=spec.charge, step=1)
        multiplicity = col2.number_input("Multiplicity", value=spec.multiplicity, step=1)
        method = col1.text_input("Method", value=spec.method)
        basis = col2.text_input("Basis", value=spec.basis)
        coords = st.text_area("Cartesian coordinates", value=spec.coords, height=200)
        extra = st.text_area("Extra options", value=spec.extra_options or "")
        ready = st.form_submit_button("Create .gjf")

    if ready:
        if not coords.strip() or "<paste" in coords:
            st.error("Please supply a coordinates block before generating the file.")
            return
        final = JobSpec(
            title=title,
            charge=int(charge),
            multiplicity=int(multiplicity),
            method=method,
            basis=basis,
            coords=coords.strip(),
            extra_options=extra.strip() or None,
        )
        gjf = build_input(final)
        st.code(gjf)
        
        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                "💾 Download .gjf",
                data=gjf.encode(),
                file_name=f"{title.replace(' ', '_')}.gjf",
                mime="text/plain",
            )
        
        with col2:
            try:
                zip_data = create_dft_zip_archive(final, build_input)
                zip_filename = f"{title.replace(' ', '_')}_dft_files.zip"
                
                st.download_button(
                    "📦 Download ZIP",
                    data=zip_data,
                    file_name=zip_filename,
                    mime="application/zip",
                    help="Download input file with README"
                )
            except Exception as e:
                st.error(f"Failed to create ZIP: {e}")
