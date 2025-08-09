from __future__ import annotations

from pathlib import Path

from chemassist.validation.bench.error_fix import run_bench_errorfix


def test_bench_errorfix_gromacs(tmp_path: Path):
    data_root = Path("validation_data/bench/gromacs_fixes")
    result = run_bench_errorfix(data_root, suite="gromacs")
    out = tmp_path / "reports/bench/gromacs.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("{}")
    assert result["suite"] == "gromacs"
    assert "fix_success_rate" in result


