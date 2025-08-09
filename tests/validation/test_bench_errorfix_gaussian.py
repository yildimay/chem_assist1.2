from __future__ import annotations

from pathlib import Path

from chemassist.validation.bench.error_fix import run_bench_errorfix


def test_bench_errorfix_gaussian(tmp_path: Path):
    data_root = Path("validation_data/bench/gaussian_fixes")
    result = run_bench_errorfix(data_root, suite="gaussian")
    out = tmp_path / "reports/bench/gaussian.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("{}")
    assert result["suite"] == "gaussian"
    assert "fix_success_rate" in result


