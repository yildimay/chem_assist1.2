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

    client = groq.Groq(api_key=os.environ["GROQ_API_KEY"])
    response = client.chat.completions.create(
        model=model,
        messages=[{"role": "system", "content": system_prompt or ""}] + messages,
    )
    return response.choices[0].message.content.strip()


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

    if "GROQ_API_KEY" in os.environ:
        return _call_groq(model or "llama3-70b-8192", messages, system_prompt)
    if "OPENAI_API_KEY" in os.environ:
        return _call_openai(model or MODEL_FALLBACK, messages, system_prompt)

    raise RuntimeError("No LLM API key found – set GROQ_API_KEY or OPENAI_API_KEY")

