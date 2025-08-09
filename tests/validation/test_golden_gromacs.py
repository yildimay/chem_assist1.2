from __future__ import annotations

from pathlib import Path

from chemassist.validation.golden.gromacs import run_golden_gromacs


def test_golden_gromacs(tmp_path: Path):
    data_root = Path("validation_data/golden/gromacs")
    result = run_golden_gromacs(data_root)
    out = tmp_path / "reports/golden/gromacs.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("{}")
    assert result["suite"] == "gromacs"
    assert "pass_rate" in result


