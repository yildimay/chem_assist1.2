from __future__ import annotations

"""Minimal, program-agnostic DFT input builder.

Right now only **Gaussian** is implemented so we can ship a working demo.
Other engines (ORCA, NWChem, CP2K…) can be added by dropping a Jinja2
Template into `_TEMPLATES` and extending `build_input()`.
"""

from dataclasses import dataclass, field
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict

# ──────────────────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class JobSpec:
    """Container for misc. input-file parameters."""

    title: str = "Untitled job"
    charge: int = 0
    multiplicity: int = 1
    method: str = "B3LYP"
    basis: str = "6-31G(d)"
    coords: str = field(
        default_factory=lambda: dedent(
            """
            O    0.0000    0.0000    0.0000
            H    0.0000    0.0000    0.9697
            H    0.0000    0.7572   -0.4848
            """
        ).strip()
    )
    extra_options: str | None = None  # raw text appended after the card

    def to_dict(self) -> Dict[str, Any]:  # convenience for Jinja2
        return self.__dict__


# ──────────────────────────────────────────────────────────────────────────────
# Templates
# ──────────────────────────────────────────────────────────────────────────────

_GAUSSIAN_TEMPLATE = dedent(
    """
    %NProcShared=8
    %Mem=8GB
    # {{ method }} / {{ basis }} {{ ' '.join(opts) }}

    {{ title }}

    {{ charge }} {{ multiplicity }}
    {{ coords }}

    {{ extra_options or '' }}
    """
).lstrip()

_TEMPLATES: dict[str, str] = {
    "gaussian": _GAUSSIAN_TEMPLATE,
    # "orca": "…", "cp2k": "…" can be added later.
}


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def build_input(spec: JobSpec, program: str = "gaussian") -> str:  # noqa: D401
    """Return a ready-to-save input string for *program* using *spec* settings."""

    program = program.lower()

    if program not in _TEMPLATES:
        raise ValueError(f"Program '{program}' not supported yet.")

    template = _TEMPLATES[program]

    # Quick inline rendering without Jinja2 to avoid template engine overhead.
    # For more complex files we can swap to real Jinja2 later.
    rendered = template.replace("{{ title }}", spec.title)
    rendered = rendered.replace("{{ charge }}", str(spec.charge))
    rendered = rendered.replace("{{ multiplicity }}", str(spec.multiplicity))
    rendered = rendered.replace("{{ method }}", spec.method)
    rendered = rendered.replace("{{ basis }}", spec.basis)

    opts = []
    if spec.extra_options:
        opts.append(spec.extra_options.strip())
    rendered = rendered.replace("{{ ' '.join(opts) }}", " ".join(opts))
    rendered = rendered.replace("{{ coords }}", spec.coords.strip())
    rendered = rendered.replace("{{ extra_options or '' }}", spec.extra_options or "")

    # Clean double newlines introduced by missing extras
    return "\n".join(line.rstrip() for line in rendered.splitlines()).strip() + "\n"


# ──────────────────────────────────────────────────────────────────────────────
# Quick CLI / debug helper
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":  # pragma: no cover
    example = JobSpec(title="Water optimisation", method="B3LYP", basis="6-31G(d,p)")
    print(build_input(example))

