from __future__ import annotations

"""Minimal LLM wrapper used across ChemAssist.

* If the env var ``GROQ_API_KEY`` is present → call Groq.
* Else if ``OPENAI_API_KEY`` is present → call OpenAI.
* Otherwise raise RuntimeError – prompts the deployer to add a key.

All modules should **only** import and call ``call_llm`` from here, never
from the vendor SDK directly.  This keeps vendor-specific changes local.
"""

import os
from typing import List, Dict, Any

MODEL_FALLBACK = "gpt-4o-mini"  # sensible default for OpenAI key

# Optional dependencies – import lazily
try:
    import groq
except ImportError:  # pragma: no cover
    groq = None  # type: ignore

try:
    import openai
except ImportError:  # pragma: no cover
    openai = None  # type: ignore


def _call_groq(model: str, messages: List[Dict[str, str]], system_prompt: str | None) -> str:
    if groq is None:
        raise RuntimeError("groq python package is not installed – add it to requirements.txt")

    try:
        # Use the current Groq client API
        client = groq.Groq(api_key=os.environ["GROQ_API_KEY"])
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "system", "content": system_prompt or ""}] + messages,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        raise RuntimeError(f"Groq API call failed: {e}")


def _call_openai(model: str, messages: List[Dict[str, str]], system_prompt: str | None) -> str:  # noqa: D401
    if openai is None:
        raise RuntimeError("openai python package is not installed.")

    openai.api_key = os.environ["OPENAI_API_KEY"]
    response = openai.chat.completions.create(
        model=model,
        messages=[{"role": "system", "content": system_prompt or ""}] + messages,
    )
    return response.choices[0].message.content.strip()


def call_llm(
    *,
    model: str | None = None,
    messages: List[Dict[str, str]],
    system_prompt: str | None = None,
) -> str:  # noqa: D401
    """Dispatch to Groq or OpenAI depending on available API keys."""

    # Try to get API keys from environment variables or Streamlit secrets
    groq_key = os.environ.get("GROQ_API_KEY")
    openai_key = os.environ.get("OPENAI_API_KEY")
    
    # Fallback to Streamlit secrets if not in environment
    if not groq_key or not openai_key:
        try:
            import streamlit as st
            if not groq_key:
                groq_key = st.secrets.get("GROQ_API_KEY")
            if not openai_key:
                openai_key = st.secrets.get("OPENAI_API_KEY")
        except Exception as e:
            # Streamlit secrets not available or not configured
            # Try alternative paths for deployment platforms
            try:
                import toml
                # Try common secrets file paths
                secrets_paths = [
                    ".streamlit/secrets.toml",
                    "/opt/render/.streamlit/secrets.toml",
                    "/opt/render/project/src/.streamlit/secrets.toml"
                ]
                
                for path in secrets_paths:
                    try:
                        if os.path.exists(path):
                            with open(path, 'r') as f:
                                secrets = toml.load(f)
                                if not groq_key:
                                    groq_key = secrets.get("GROQ_API_KEY")
                                if not openai_key:
                                    openai_key = secrets.get("OPENAI_API_KEY")
                                break
                    except Exception:
                        continue
            except ImportError:
                pass  # toml not available
    
    if groq_key:
        # Temporarily set the environment variable for the API call
        original_groq_key = os.environ.get("GROQ_API_KEY")
        os.environ["GROQ_API_KEY"] = groq_key
        try:
            return _call_groq(model or "llama3-70b-8192", messages, system_prompt)
        finally:
            # Restore original value
            if original_groq_key:
                os.environ["GROQ_API_KEY"] = original_groq_key
            else:
                os.environ.pop("GROQ_API_KEY", None)
    
    if openai_key:
        # Temporarily set the environment variable for the API call
        original_openai_key = os.environ.get("OPENAI_API_KEY")
        os.environ["OPENAI_API_KEY"] = openai_key
        try:
            return _call_openai(model or MODEL_FALLBACK, messages, system_prompt)
        finally:
            # Restore original value
            if original_openai_key:
                os.environ["OPENAI_API_KEY"] = original_openai_key
            else:
                os.environ.pop("OPENAI_API_KEY", None)

    raise RuntimeError("No LLM API key found – set GROQ_API_KEY or OPENAI_API_KEY in environment or Streamlit secrets")

