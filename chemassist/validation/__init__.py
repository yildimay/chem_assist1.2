from __future__ import annotations

"""ChemAssist validation suite.

Modules:
- golden: Golden-set validators for Gaussian and GROMACS
- bench:  Error-fix benchmark runners
- io:     Lightweight log readers/parsers
- utils:  Tolerances and normalization helpers
- metrics:Report writers (JSON/JUnit/Markdown)
"""

__all__ = [
    "golden",
    "bench",
    "io",
    "metrics",
    "utils",
]


