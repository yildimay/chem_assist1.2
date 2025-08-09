from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from ..io.readers import parse_gromacs_mdlog
from ..utils.tolerances import Tolerances


@dataclass
class GoldenCaseResult:
    case_id: str
    passed: bool
    metrics: Dict[str, float | str]
    expected: Dict[str, float | str]
    tolerances: Dict[str, float]
    diffs: Dict[str, float]


def _load_expected(case_dir: Path) -> dict:
    data = json.loads((case_dir / "expected.json").read_text(encoding="utf-8"))
    return data


def _compare(metrics: dict, expected: dict, tolerances: dict) -> tuple[bool, dict]:
    diffs: Dict[str, float] = {}
    ok = True
    for key, exp_val in expected.items():
        if key == "ensemble":
            ok = ok and (str(metrics.get(key, "")).strip() == str(exp_val).strip())
            continue
        if key not in metrics:
            ok = False
            diffs[key] = float("inf")
            continue
        mval = float(metrics[key])
        tol = tolerances.get(key, Tolerances.default_gromacs().get(key, 0.0))
        diff = abs(mval - float(exp_val))
        diffs[key] = diff
        if diff > tol:
            ok = False
    return ok, diffs


def run_golden_gromacs(data_root: Path) -> dict:
    cases_dir = Path(data_root) / "cases"
    results: List[GoldenCaseResult] = []
    for case_dir in sorted(cases_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        exp = _load_expected(case_dir)
        expected_metrics = exp.get("metrics", {})
        tolerances = exp.get("tolerances", {})

        log_path = case_dir / "md.log"
        metrics = parse_gromacs_mdlog(log_path.read_text(encoding="utf-8"))

        passed, diffs = _compare(metrics, expected_metrics, tolerances)
        results.append(
            GoldenCaseResult(
                case_id=exp.get("id", case_dir.name),
                passed=passed,
                metrics=metrics,
                expected=expected_metrics,
                tolerances={**Tolerances.default_gromacs(), **tolerances},
                diffs=diffs,
            )
        )

    pass_rate = sum(1 for r in results if r.passed) / max(1, len(results))
    return {
        "suite": "gromacs",
        "num_cases": len(results),
        "pass_rate": pass_rate,
        "cases": [r.__dict__ for r in results],
    }


