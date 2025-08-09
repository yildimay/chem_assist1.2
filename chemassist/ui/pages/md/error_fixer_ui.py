from __future__ import annotations

"""MD • Error Fixer UI

Interactive chat interface for GROMACS error fixing, matching the DFT chat UX.
"""

import streamlit as st
import requests
import os
from PIL import Image
import pytesseract

# ───────────────────────────────────────────────────────────────
# Streamlit page
# ───────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🚑 MD • Error Fixer")
    show_chat_mode()


def show_chat_mode() -> None:
    """Show the chat interface mode for GROMACS MD error fixing."""
    st.markdown(
        "Welcome to the **GROMACS Error Fixer Chat**.\n"
        "- Describe your MD issue in the chatbox below.\n"
        "- Optionally, upload a screenshot of the error.\n"
        "- Optionally, upload a GROMACS `.mdp` and/or `.log`/`.out` file.\n\n"
        "**Note**: Text input is required. Other inputs are optional."
    )

    # Initialize chat history (MD-specific)
    if "md_chat_messages" not in st.session_state:
        st.session_state.md_chat_messages = []

    # Display chat history
    for msg in st.session_state.md_chat_messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    # User input and uploads
    with st.form("md_chat_form"):
        user_input = st.text_area(
            "Describe your GROMACS MD error:",
            height=100,
            placeholder="e.g. LINCS warnings, PME errors, segmentation fault...",
        )
        stage = st.selectbox(
            "Simulation Stage (optional)",
            ["md", "em", "nvt", "npt"],
            index=0,
            help="Select the stage related to the issue",
        )
        uploaded_image = st.file_uploader(
            "Optional: Upload error screenshot",
            type=["png", "jpg", "jpeg"],
        )
        uploaded_mdp_file = st.file_uploader(
            "Optional: Upload GROMACS MDP File (.mdp)",
            type=["mdp"],
        )
        uploaded_log_file = st.file_uploader(
            "Optional: Upload GROMACS Log File (.log, .out)",
            type=["log", "out"],
        )
        submitted = st.form_submit_button("Send")

    if submitted:
        if user_input.strip() == "":
            st.warning("Text input is required.")
        else:
            # Save user message
            st.session_state.md_chat_messages.append({"role": "user", "content": user_input})
            with st.chat_message("user"):
                st.markdown(user_input)

            # Process image if available
            image_text = ""
            if uploaded_image is not None:
                try:
                    image = Image.open(uploaded_image)
                    image_text = pytesseract.image_to_string(image)
                    st.markdown("*Extracted from image:* \n" + image_text)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Failed to process image: {e}")

            # Read file contents
            mdp_file_text = ""
            log_file_text = ""

            if uploaded_mdp_file is not None:
                try:
                    mdp_file_text = uploaded_mdp_file.read().decode("utf-8", errors="ignore")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Failed to read MDP file: {e}")

            if uploaded_log_file is not None:
                try:
                    log_file_text = uploaded_log_file.read().decode("utf-8", errors="ignore")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Failed to read log file: {e}")

            # Combine and process all inputs
            final_input = f"""USER MESSAGE:
{user_input}

IMAGE TEXT:
{image_text}

MDP FILE:
{mdp_file_text}

LOG FILE:
{log_file_text}

STAGE:
{stage}"""

            # Call Groq API (same approach as DFT chat)
            try:
                api_key = os.getenv("GROQ_API_KEY")
                if not api_key:
                    try:
                        api_key = st.secrets.get("GROQ_API_KEY")
                    except Exception:  # noqa: BLE001
                        pass

                if not api_key:
                    ai_response = (
                        "❌ No API key found. Please set GROQ_API_KEY environment variable or add it to Streamlit secrets."
                    )
                else:
                    headers = {
                        "Authorization": f"Bearer {api_key}",
                        "Content-Type": "application/json",
                    }

                    data = {
                        "model": "llama3-70b-8192",
                        "messages": [
                            {
                                "role": "system",
                                "content": (
                                    "You are a helpful assistant for GROMACS error fixing. "
                                    "If the question is not related to GROMACS or MD errors, reply with: "
                                    "'This tool is specifically for GROMACS error fixing. Please ask a GROMACS-related question.'"
                                ),
                            },
                            {"role": "user", "content": final_input},
                        ],
                    }

                    response = requests.post(
                        "https://api.groq.com/openai/v1/chat/completions", headers=headers, json=data
                    )
                    response.raise_for_status()
                    ai_response = response.json()["choices"][0]["message"]["content"]

            except Exception as e:  # noqa: BLE001
                ai_response = f"❌ Error while contacting Groq API: {e}"

            # Save and display AI response
            st.session_state.md_chat_messages.append({"role": "assistant", "content": ai_response})
            with st.chat_message("assistant"):
                st.markdown(ai_response)

    # Clear chat button
    if st.session_state.md_chat_messages:
        if st.button("🗑️ Clear Chat History"):
            st.session_state.md_chat_messages = []
            st.rerun()
