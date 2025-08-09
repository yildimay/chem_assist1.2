from __future__ import annotations

import re


def normalize_text_for_diff(text: str) -> str:
    """Normalize text for robust diff comparison.

    - Strip trailing spaces
    - Remove consecutive blank lines
    - Remove comments (#... for Gaussian, ; ... for mdp)
    - Normalize key order for mdp-like key=value lines by sorting lines within blocks
    """
    # Remove comments
    lines = []
    for raw in text.splitlines():
        line = raw
        # Gaussian comments often after ! or full-line; keep it simple
        line = re.sub(r"\s*!.*$", "", line)
        # mdp comments start with ';'
        line = re.sub(r"\s*;.*$", "", line)
        lines.append(line.rstrip())

    # Collapse multiple blank lines
    cleaned = []
    blank = False
    for l in lines:
        if l.strip() == "":
            if not blank:
                cleaned.append("")
            blank = True
        else:
            cleaned.append(l)
            blank = False

    # For mdp-like lines (contains '='), sort lines to reduce ordering noise
    if any("=" in l for l in cleaned):
        kv = [l for l in cleaned if "=" in l]
        other = [l for l in cleaned if "=" not in l]
        kv_sorted = sorted(kv, key=lambda s: s.split("=")[0].strip())
        cleaned = other + kv_sorted

    return "\n".join(cleaned).strip() + "\n"


