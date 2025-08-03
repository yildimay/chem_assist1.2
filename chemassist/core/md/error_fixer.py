from __future__ import annotations

"""
GROMACS MD error fixer.

Call:
    from chemassist.core.md.error_fixer import fix_md_input
    result = fix_md_input(mdp_text, log_text)
    # result = {"diagnosis": "...", "fixed_mdp": "..."}
"""

from typing import Dict

from chemassist import call_llm, get_logger

LOGGER = get_logger(__name__)

_SYSTEM_PROMPT = (
    "You are an expert in GROMACS molecular dynamics simulations. "
    "A user provides a broken .mdp file plus the tail of the .log/.out file.\n"
    "Tasks:\n"
    "1. Diagnose the root cause (<100 words).\n"
    "2. Provide a corrected .mdp file that should run.\n"
    "Return exactly two fenced blocks:\n"
    "```DIAGNOSIS\n<one paragraph>\n```\n"
    "```FIXED_MDP\n<complete corrected .mdp file>\n```\n"
    "Common GROMACS issues to watch for:\n"
    "- Incorrect parameter values (negative temperatures, invalid cutoffs)\n"
    "- Missing required parameters\n"
    "- Incompatible parameter combinations\n"
    "- File format issues\n"
    "- Memory/performance problems"
)


def _tail(text: str, n: int = 200) -> str:
    """Return last *n* lines of *text*."""
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def fix_md_input(
    mdp_text: str,
    log_text: str,
    *,
    stage: str = "md",
) -> Dict[str, str]:
    """
    Diagnose GROMACS failure and return patched .mdp file.

    Args:
        mdp_text: Content of the .mdp file
        log_text: Content of the .log/.out file
        stage: Simulation stage ("em", "nvt", "npt", "md")

    Returns: {"diagnosis": str, "fixed_mdp": str}
    Raises  : RuntimeError if LLM response malformed.
    """
    snippet = _tail(log_text, 200)
    user_msg = f"<MDP FILE>\n{mdp_text.strip()}\n<LOG TAIL>\n{snippet.strip()}\n<STAGE>\n{stage}"

    response = call_llm(
        system_prompt=_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_msg}],
    )

    diagnosis, fixed = "", ""
    mode = None
    
    # Debug logging
    LOGGER.info(f"Parsing LLM response for {stage} stage (length: {len(response)})")
    
    for line in response.splitlines():
        if line.startswith("```DIAGNOSIS"):
            mode = "d"
            continue
        if line.startswith("```FIXED_MDP"):
            mode = "f"
            continue
        if line.startswith("```"):
            mode = None
            continue
        if mode == "d":
            diagnosis += line + "\n"
        elif mode == "f":
            fixed += line + "\n"

    # More detailed error handling
    if not fixed.strip():
        LOGGER.error(f"Failed to parse LLM response. Response was:\n{response}")
        
        # Try fallback parsing for different formats
        if "```" in response:
            # Try to extract any code blocks
            import re
            code_blocks = re.findall(r'```(?:.*?)\n(.*?)```', response, re.DOTALL)
            if code_blocks:
                LOGGER.info(f"Found {len(code_blocks)} code blocks, using last one as FIXED_MDP")
                fixed = code_blocks[-1]
                diagnosis = "Auto-extracted from LLM response (format parsing failed)"
            else:
                raise RuntimeError(
                    f"LLM did not return a FIXED_MDP block. "
                    f"Response length: {len(response)} characters. "
                    f"Please try again or check the input format."
                )
        else:
            raise RuntimeError(
                f"LLM did not return a FIXED_MDP block. "
                f"Response length: {len(response)} characters. "
                f"Please try again or check the input format."
            )

    return {"diagnosis": diagnosis.strip(), "fixed_mdp": fixed.rstrip() + "\n"}


def analyze_md_log(log_text: str) -> Dict[str, str]:
    """
    Analyze GROMACS log file for common issues.
    
    Args:
        log_text: Content of the .log file
        
    Returns: {"analysis": str, "suggestions": str}
    """
    analysis_prompt = (
        "You are a GROMACS expert. Analyze this log file and provide:\n"
        "1. A brief analysis of what happened\n"
        "2. Specific suggestions to fix the issue\n"
        "Return in this format:\n"
        "```ANALYSIS\n<analysis>\n```\n"
        "```SUGGESTIONS\n<suggestions>\n```"
    )
    
    try:
        response = call_llm(
            system_prompt=analysis_prompt,
            messages=[{"role": "user", "content": log_text}],
        )
        
        analysis, suggestions = "", ""
        mode = None
        
        for line in response.splitlines():
            if line.startswith("```ANALYSIS"):
                mode = "a"
                continue
            if line.startswith("```SUGGESTIONS"):
                mode = "s"
                continue
            if line.startswith("```"):
                mode = None
                continue
            if mode == "a":
                analysis += line + "\n"
            elif mode == "s":
                suggestions += line + "\n"
        
        return {
            "analysis": analysis.strip(),
            "suggestions": suggestions.strip()
        }
        
    except Exception as e:
        LOGGER.error(f"Failed to analyze log: {e}")
        return {
            "analysis": "Could not analyze log file automatically.",
            "suggestions": "Please check the log file manually for error messages."
        }

