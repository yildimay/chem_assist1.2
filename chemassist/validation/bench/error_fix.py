from __future__ import annotations

import difflib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from ...core.dft.error_fixer import fix_input as gaussian_fix
from ...core.md.error_fixer import fix_md_input as gromacs_fix
from ..utils.normalize import normalize_text_for_diff


@dataclass
class BenchCaseResult:
    case_id: str
    similarity: float
    passed: bool
    degradation_detected: bool
    details: Dict[str, str | float | bool]


def _read_case(case_dir: Path, suite: str) -> tuple[str, str]:
    if suite == "gaussian":
        broken = (case_dir / "broken.gjf").read_text(encoding="utf-8")
        context = (case_dir / "context.log").read_text(encoding="utf-8")
        return broken, context
    else:
        broken = (case_dir / "broken.mdp").read_text(encoding="utf-8")
        context = (case_dir / "context.log").read_text(encoding="utf-8")
        return broken, context


def _expected_patch(case_dir: Path) -> str:
    return (case_dir / "expected_patch.diff").read_text(encoding="utf-8")


def _compute_similarity(a: str, b: str) -> float:
    seq = difflib.SequenceMatcher(a=a.splitlines(), b=b.splitlines())
    return seq.ratio()


def _apply_fix(suite: str, broken: str, context: str) -> str:
    if suite == "gaussian":
        res = gaussian_fix(input_text=broken, log_text=context)
        fixed = res.get("fixed_input", "")
    else:
        res = gromacs_fix(mdp_text=broken, log_text=context, stage="md")
        fixed = res.get("fixed_mdp", "")
    return fixed


def run_bench_errorfix(data_root: Path, *, suite: str) -> dict:
    cases_dir = Path(data_root)
    results: List[BenchCaseResult] = []
    for case_dir in sorted(cases_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        case_id = case_dir.name
        broken, context = _read_case(case_dir, suite)
        expected = _expected_patch(case_dir)

        fixed = _apply_fix(suite, broken, context)
        norm_fixed = normalize_text_for_diff(fixed)
        norm_expected = normalize_text_for_diff(expected)

        similarity = _compute_similarity(norm_fixed, norm_expected)

        # basic degradation heuristic: fixed must not remove route line (Gaussian) or change integrator (GROMACS)
        degradation = False
        if suite == "gaussian":
            if "#" in broken and "#" not in fixed:
                degradation = True
        else:
            import re
            m_old = re.search(r"^\s*integrator\s*=\s*(\S+)", broken, re.MULTILINE)
            m_new = re.search(r"^\s*integrator\s*=\s*(\S+)", fixed, re.MULTILINE)
            if m_old and m_new and m_old.group(1) != m_new.group(1):
                degradation = True

        passed = similarity >= 0.90 and not degradation

        results.append(
            BenchCaseResult(
                case_id=case_id,
                similarity=similarity,
                passed=passed,
                degradation_detected=degradation,
                details={
                    "similarity": round(similarity, 4),
                    "degradation": degradation,
                },
            )
        )

    success_rate = sum(1 for r in results if r.passed) / max(1, len(results))
    degradation_rate = sum(1 for r in results if r.degradation_detected) / max(1, len(results))
    return {
        "suite": suite,
        "num_cases": len(results),
        "fix_success_rate": success_rate,
        "degradation_rate": degradation_rate,
        "cases": [r.__dict__ for r in results],
    }


