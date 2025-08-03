import importlib
import streamlit as st
from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()

# ─────────────────────────────────────────────────────────────
# ChemAssist v2 – main entry
# This version assumes UI files live under:
#   chemassist/ui/pages/<category>/<page>_ui.py
# ─────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ChemAssist v2",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────────────────────────────────────────
# Sidebar – navigation
# ─────────────────────────────────────────────────────────────

st.sidebar.title("⚗️ ChemAssist")

# Check for API keys and show helpful message
if not os.environ.get("GROQ_API_KEY") and not os.environ.get("OPENAI_API_KEY"):
    st.sidebar.warning(
        "⚠️ **LLM features disabled**\n\n"
        "To use error fixing and AI features, add your API key:\n\n"
        "**Option 1:** Edit `.env` file\n"
        "**Option 2:** Set environment variable:\n"
        "```bash\nexport GROQ_API_KEY='your_key_here'\n```\n\n"
        "Get free API key: https://console.groq.com/"
    )

CATEGORY = st.sidebar.radio(
    "Tool category",
    ("DFT", "Molecular Dynamics", "Other"),
    key="category",
)

# Mapping: (category, tool) → dotted import path
PAGE_REGISTRY: dict[tuple[str, str], str] = {
    ("DFT", "Input Creator"): "chemassist.ui.pages.dft.input_creator_ui",
    ("DFT", "Error Fixer"):   "chemassist.ui.pages.dft.error_fixer_ui",

    ("Molecular Dynamics", "Input Creator"): "chemassist.ui.pages.md.input_creator_ui",
    ("Molecular Dynamics", "Error Fixer"):   "chemassist.ui.pages.md.error_fixer_ui",

    ("Other", "Environmental Risk Calculator"): "chemassist.ui.pages.other.env_risk_ui",
    ("Other", "SMILES → 3D Modeler"):           "chemassist.ui.pages.other.smiles_modeler_ui",
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

# ─────────────────────────────────────────────────────────────
# Dynamic import & execution
# ─────────────────────────────────────────────────────────────

module_path = PAGE_REGISTRY.get((CATEGORY, TOOL))

if module_path is None:
    st.error("Selected page is not registered yet.")
else:
    try:
        module = importlib.import_module(module_path)
        if hasattr(module, "show_page"):
            module.show_page()
        else:
            st.warning(
                f"`{module_path}` does not expose a `show_page()` function. "
                "Please implement it."
            )
    except ModuleNotFoundError as err:
        st.error(
            f"Module {module_path!r} cannot be imported. Did you create the file?\n{err}"
        )
    except Exception as err:
        st.exception(err)

# ─────────────────────────────────────────────────────────────
# Footer
# ─────────────────────────────────────────────────────────────

st.sidebar.markdown("---")
st.sidebar.caption("ChemAssist v2 – AI-powered computational chemistry toolkit (free)")
