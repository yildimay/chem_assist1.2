from __future__ import annotations

import streamlit as st

from chemassist.core.dft.input_creator import JobSpec, build_input
from chemassist.core.dft.method_suggester import recommend

# ───────────────────────────────────────────────────────────────────────────────
# Streamlit page – fully automatic SMILES → .gjf
# ───────────────────────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🧬 SMILES-only DFT Input Creator")
    st.write("Enter a SMILES string; ChemAssist picks a method, basis, and builds a ready-to-run Gaussian input file.")

    smiles = st.text_input("SMILES", placeholder="C1=CC=CC=C1 (benzene)")

    if st.button("Generate .gjf", disabled=not smiles):
        try:
            suggestion = recommend(smiles)
            st.success(
                f"**Chosen:** {suggestion.method} / {suggestion.basis}  \n"
                f"*Reason:* {suggestion.reason}  \n"
                f"*Charge:* {suggestion.charge}, *Multiplicity:* {suggestion.multiplicity}"
            )

            spec = JobSpec(
                title=f"Auto-generated for {smiles}",
                charge=suggestion.charge,
                multiplicity=suggestion.multiplicity,
                method=suggestion.method,
                basis=suggestion.basis,
                coords=suggestion.xyz,
            )
            gjf_text = build_input(spec)
            st.code(gjf_text, language="text")
            st.download_button(
                label="💾 Download .gjf",
                data=gjf_text.encode(),
                file_name=f"{smiles.replace('/', '_')}.gjf",
                mime="text/plain",
            )
        except Exception as exc:  # noqa: BLE001
            st.error(f"Failed to generate input: {exc}")
