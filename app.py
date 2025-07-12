import importlib
import streamlit as st

# ──────────────────────────────────────────────────────────────────────────────
# App configuration
# ──────────────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ChemAssist v2",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ──────────────────────────────────────────────────────────────────────────────
# Sidebar – high-level navigation
# ──────────────────────────────────────────────────────────────────────────────

st.sidebar.title("⚗️ ChemAssist")

CATEGORY = st.sidebar.radio(
    "Tool category",
    ("DFT", "Molecular Dynamics", "Other"),
    key="category",
)

# Mapping: (category, tool) → dotted import path for the corresponding Streamlit
# UI function.  Each UI file must expose a `show_page()` entry point.
PAGE_REGISTRY: dict[tuple[str, str], str] = {
    ("DFT", "Input Creator"): "chemassist.ui.dft.input_creator_ui",
    ("DFT", "Error Fixer"): "chemassist.ui.dft.error_fixer_ui",
    ("Molecular Dynamics", "Input Creator"): "chemassist.ui.md.input_creator_ui",
    ("Molecular Dynamics", "Error Fixer"): "chemassist.ui.md.error_fixer_ui",
    ("Other", "Environmental Risk Calculator"): "chemassist.ui.other.env_risk_ui",
    ("Other", "SMILES → 3D Modeler"): "chemassist.ui.other.smiles_modeler_ui",
}

if CATEGORY == "DFT":
    TOOL = st.sidebar.radio(
        "DFT tools",
        ("Input Creator", "Error Fixer"),
        key="dft_tool",
    )

elif CATEGORY == "Molecular Dynamics":
    TOOL = st.sidebar.radio(
        "MD tools",
        ("Input Creator", "Error Fixer"),
        key="md_tool",
    )

else:  # Other
    TOOL = st.sidebar.radio(
        "Misc tools",
        ("Environmental Risk Calculator", "SMILES → 3D Modeler"),
        key="misc_tool",
    )

# ──────────────────────────────────────────────────────────────────────────────
# Dynamic import & execution
# ──────────────────────────────────────────────────────────────────────────────

module_path = PAGE_REGISTRY.get((CATEGORY, TOOL))

if module_path is None:
    st.error("Selected page is not registered (yet).")
else:
    try:
        module = importlib.import_module(module_path)
        # All UI modules must implement `show_page()`.
        if hasattr(module, "show_page"):
            module.show_page()
        else:
            st.warning(
                f"`{module_path}` does not expose a `show_page()` function. "
                "Please implement it."
            )
    except ModuleNotFoundError as err:
        st.error(
            f"Module {module_path!r} cannot be imported. Have you created the file?\n{err}"
        )
    except Exception as err:
        st.exception(err)

# ──────────────────────────────────────────────────────────────────────────────
# Footer
# ──────────────────────────────────────────────────────────────────────────────

st.sidebar.markdown("---")
st.sidebar.caption("ChemAssist v2 – AI-powered computational chemistry toolkit (free)")
