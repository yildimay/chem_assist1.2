from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


def write_summary_json(result: Dict[str, Any], path: Path) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    except Exception:
        pass


def write_summary_md(result: Dict[str, Any], path: Path) -> None:
    try:
        suite = result.get("suite", "")
        num_cases = result.get("num_cases", 0)
        pr = result.get("pass_rate")
        fsr = result.get("fix_success_rate")
        degr = result.get("degradation_rate")

        lines = [f"# Validation Summary ({suite})", ""]
        if pr is not None:
            lines.append(f"- Golden pass rate: {pr:.2%} ({num_cases} cases)")
        if fsr is not None:
            lines.append(f"- Fix success rate: {fsr:.2%}")
        if degr is not None:
            lines.append(f"- Degradation rate: {degr:.2%}")

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception:
        pass


