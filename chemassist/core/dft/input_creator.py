from __future__ import annotations

"""Generic DFT input-file builder (Gaussian first).

Exports:
    • JobSpec   – dataclass container for form values
    • build_input(spec, program="gaussian") -> str
"""

from dataclasses import dataclass, field
from textwrap import dedent
from typing import Any, Dict

# ───────────────────────────────────────────────────────────────
# Public dataclass
# ───────────────────────────────────────────────────────────────


@dataclass
class JobSpec:
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
    extra_options: str | None = None  # appended verbatim after geometry

    def as_dict(self) -> Dict[str, Any]:
        return self.__dict__


# ───────────────────────────────────────────────────────────────
# Templates – can add ORCA, CP2K later
# ───────────────────────────────────────────────────────────────

_GAUSSIAN_TEMPLATE = dedent(
    """
    %NProcShared=8
    %Mem=8GB
    # {{method}} / {{basis}}{{opts}}

    {{title}}

    {{charge}} {{multiplicity}}
    {{coords}}

    {{extra}}
    """
).lstrip()

_TEMPLATES = {
    "gaussian": _GAUSSIAN_TEMPLATE,
}


# ───────────────────────────────────────────────────────────────
# Core render function (quick string replace to avoid heavy deps)
# ───────────────────────────────────────────────────────────────

def build_input(spec: JobSpec, program: str = "gaussian") -> str:  # noqa: D401
    prog = program.lower()
    if prog not in _TEMPLATES:
        raise ValueError(f"Unsupported program: {program}")
    tpl = _TEMPLATES[prog]

    rendered = (
        tpl.replace("{{title}}", spec.title)
        .replace("{{charge}}", str(spec.charge))
        .replace("{{multiplicity}}", str(spec.multiplicity))
        .replace("{{method}}", spec.method)
        .replace("{{basis}}", spec.basis)
        .replace("{{coords}}", spec.coords.strip())
    )
    opts = f" {spec.extra_options.strip()}" if spec.extra_options else ""
    rendered = rendered.replace("{{opts}}", opts)
    rendered = rendered.replace("{{extra}}", spec.extra_options or "")

    # strip trailing whitespace lines
    return "\n".join(l.rstrip() for l in rendered.splitlines()).strip() + "\n"


# ───────────────────────────────────────────────────────────────
# Debug CLI helper
# ───────────────────────────────────────────────────────────────

if __name__ == "__main__":  # pragma: no cover
    print(build_input(JobSpec()))
