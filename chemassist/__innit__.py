from __future__ import annotations

import logging
from importlib.metadata import version, PackageNotFoundError

# ── Version string ───────────────────────────────────────────
try:
    __version__ = version("chemassist")
except PackageNotFoundError:  # when running from source
    __version__ = "0.0.dev"

# ── Logger helper ────────────────────────────────────────────
def get_logger(name: str | None = None, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name or "chemassist")
    if not logger.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s"))
        logger.addHandler(h)
    logger.setLevel(level)
    return logger

# ── LLM router re-export ─────────────────────────────────────
try:
    from .models.llm_router import call_llm  # noqa: F401  (re-export)
except ModuleNotFoundError:  # file not yet created
    def call_llm(*_a, **_kw):  # type: ignore
        raise RuntimeError(
            "LLM router not available. Add chemassist/models/llm_router.py "
            "or install vendor SDKs."
        )

__all__ = ["__version__", "get_logger", "call_llm"]
