from __future__ import annotations

"""Extracts DFT methods & basis sets from a research-paper PDF.

Exports a single public helper:
    >>> hits = extract(pdf_bytes)
    >>> for row in hits: print(row.method, row.basis, row.count)

If *no* recognised pair is found, falls back to an LLM call—and returns a
single candidate suggested by Groq (via chemassist.call_llm).
"""

from dataclasses import dataclass
from collections import Counter
from io import BytesIO
from pathlib import Path
from typing import Iterable, List

import fitz  # PyMuPDF

from chemassist import call_llm, get_logger

LOGGER = get_logger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Lookup tables – EXTEND as needed
# Sources: Gaussian manual + popular literature reviews (B3LYP, ωB97X-D…)
# ──────────────────────────────────────────────────────────────────────────────

METHODS = [
    "B3LYP",
    "M06-2X",
    "ωB97X-D",
    "PBE0",
    "CAM-B3LYP",
    "TPSSh",
    "B97-D3",
    "BLYP",
    "HF",
]

BASIS_SETS = [
    "6-31G",
    "6-31G(d)",
    "6-31G(d,p)",
    "6-311+G(d,p)",
    "def2-SVP",
    "def2-TZVP",
    "cc-pVDZ",
    "cc-pVTZ",
    "LANL2DZ",
]

# Pre-compile lowercase variants for fast membership tests.
_METHODS_LC = {m.lower(): m for m in METHODS}
_BASIS_LC = {b.lower(): b for b in BASIS_SETS}


@dataclass
class Hit:
    method: str
    basis: str
    count: int


# ──────────────────────────────────────────────────────────────────────────────
# Core logic
# ──────────────────────────────────────────────────────────────────────────────

def _extract_text(pdf_bytes: bytes) -> str:
    """Return concatenated text from every page using PyMuPDF."""
    with fitz.open(stream=pdf_bytes, filetype="pdf") as doc:
        return "\n".join(page.get_text("text") for page in doc)


def _find_candidates(text: str) -> List[Hit]:
    tok = text.lower()

    method_hits: Counter[str] = Counter()
    basis_hits: Counter[str] = Counter()

    for m_lc, m in _METHODS_LC.items():
        method_hits[m] = tok.count(m_lc)

    for b_lc, b in _BASIS_LC.items():
        basis_hits[b] = tok.count(b_lc)

    # Any zero-count entry is irrelevant.
    method_hits = Counter({k: v for k, v in method_hits.items() if v})
    basis_hits = Counter({k: v for k, v in basis_hits.items() if v})

    hits: list[Hit] = []
    for m, m_cnt in method_hits.most_common():
        for b, b_cnt in basis_hits.most_common():
            hits.append(Hit(method=m, basis=b, count=min(m_cnt, b_cnt)))

    # Sort overall by combined count (descending)
    hits.sort(key=lambda h: h.count, reverse=True)
    return hits


def _llm_fallback(text: str) -> Hit | None:  # pragma: no cover
    """Call Groq/OpenAI to guess a method/basis if regex failed."""
    prompt = (
        "You are an expert in computational chemistry. A user gave you the "
        "following excerpt of a research paper and wants to know which DFT "
        "functional and basis set were most likely used. "
        "If nothing is obvious, suggest a sensible B3LYP-tier combo.\n\n"
        f"---\n{text[:5000]}\n---\n\n"
        "Respond with exactly one line: <METHOD> ; <BASIS> ; <RATIONALE>."
    )

    try:
        response = call_llm(
            model="gpt-4o",  # or llama3-70b
            system_prompt="You are a careful assistant.",
            messages=[{"role": "user", "content": prompt}],
        )
        parts = [p.strip() for p in response.split(";")]
        if len(parts) >= 2:
            return Hit(method=parts[0], basis=parts[1], count=1)
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("LLM fallback failed: %s", exc)
    return None


def extract(pdf_bytes: bytes, *, use_llm_fallback: bool = True) -> List[Hit]:
    """Return a ranked list of (method, basis, count) discovered in *pdf_bytes*.

    If *use_llm_fallback* and regex finds nothing, try an LLM suggestion.
    """

    text = _extract_text(pdf_bytes)
    hits = _find_candidates(text)

    if hits or not use_llm_fallback:
        return hits

    llm_hit = _llm_fallback(text)
    return [llm_hit] if llm_hit else []


# ──────────────────────────────────────────────────────────────────────────────
# Debug helper
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":  # pragma: no cover
    sample_pdf = Path("example.pdf").read_bytes()
    for hit in extract(sample_pdf):
        print(hit)
