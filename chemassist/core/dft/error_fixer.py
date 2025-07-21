# chemassist/core/dft/error_fixer.py
from __future__ import annotations

"""
Gaussian-style error fixer.

Call:
    from chemassist.core.dft.error_fixer import fix_input
    result = fix_input(input_text, log_text)
    # result = {"diagnosis": "...", "fixed_input": "..."}
"""

from typing import Dict

from chemassist import call_llm, get_logger

LOGGER = get_logger(__name__)

_SYSTEM_PROMPT = (
    "You are an expert in computational chemistry software. "
    "A user provides a broken input file plus the tail of the log/out.\n"
    "Tasks:\n"
    "1. Diagnose the root cause (<100 words).\n"
    "2. Provide a corrected input that should run.\n"
    "Return exactly two fenced blocks:\n"
    "```DIAGNOSIS\n<one paragraph>\n```\n"
    "```FIXED_INPUT\n<full corrected input>\n```"
)


def _tail(text: str, n: int = 200) -> str:
    """Return last *n* lines of *text*."""
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def fix_input(
    input_text: str,
    log_text: str,
    *,
    program: str = "gaussian",
) -> Dict[str, str]:
    """
    Diagnose failure and return patched input.

    Returns: {"diagnosis": str, "fixed_input": str}
    Raises  : RuntimeError if LLM response malformed.
    """
    snippet = _tail(log_text, 200)
    user_msg = f"<INPUT FILE>\n{input_text.strip()}\n<LOG TAIL>\n{snippet.strip()}"

    response = call_llm(
        system_prompt=_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_msg}],
    )

    diagnosis, fixed = "", ""
    mode = None
    for line in response.splitlines():
        if line.startswith("```DIAGNOSIS"):
            mode = "d"
            continue
        if line.startswith("```FIXED_INPUT"):
            mode = "f"
            continue
        if line.startswith("```"):
            mode = None
            continue
        if mode == "d":
            diagnosis += line + "\n"
        elif mode == "f":
            fixed += line + "\n"

    if not fixed.strip():
        raise RuntimeError("LLM did not return a FIXED_INPUT block.")

    return {"diagnosis": diagnosis.strip(), "fixed_input": fixed.rstrip() + "\n"}
