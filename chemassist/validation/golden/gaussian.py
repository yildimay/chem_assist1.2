from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from ..io.readers import parse_gaussian_log
from ..utils.tolerances import Tolerances


@dataclass
class GoldenCaseResult:
    case_id: str
    passed: bool
    metrics: Dict[str, float]
    expected: Dict[str, float]
    tolerances: Dict[str, float]
    diffs: Dict[str, float]


def _load_expected(case_dir: Path) -> dict:
    data = json.loads((case_dir / "expected.json").read_text(encoding="utf-8"))
    return data


def _compare(metrics: dict, expected: dict, tolerances: dict) -> tuple[bool, dict]:
    diffs: Dict[str, float] = {}
    ok = True
    for key, exp_val in expected.items():
        if key not in metrics:
            ok = False
            diffs[key] = float("inf")
            continue
        mval = metrics[key]
        tol = tolerances.get(key, Tolerances.default_gaussian().get(key, 0.0))
        diff = abs(mval - exp_val)
        diffs[key] = diff
        if diff > tol:
            ok = False
    return ok, diffs


def run_golden_gaussian(data_root: Path) -> dict:
    cases_dir = Path(data_root) / "cases"
    results: List[GoldenCaseResult] = []
    for case_dir in sorted(cases_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        exp = _load_expected(case_dir)
        expected_metrics = exp.get("metrics", {})
        tolerances = exp.get("tolerances", {})

        log_path = case_dir / "output.log"
        metrics = parse_gaussian_log(log_path.read_text(encoding="utf-8"))

        passed, diffs = _compare(metrics, expected_metrics, tolerances)
        results.append(
            GoldenCaseResult(
                case_id=exp.get("id", case_dir.name),
                passed=passed,
                metrics=metrics,
                expected=expected_metrics,
                tolerances={**Tolerances.default_gaussian(), **tolerances},
                diffs=diffs,
            )
        )

    pass_rate = sum(1 for r in results if r.passed) / max(1, len(results))
    return {
        "suite": "gaussian",
        "num_cases": len(results),
        "pass_rate": pass_rate,
        "cases": [r.__dict__ for r in results],
    }


