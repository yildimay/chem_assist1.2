from __future__ import annotations

from pathlib import Path

from chemassist.validation.golden.gaussian import run_golden_gaussian


def test_golden_gaussian(tmp_path: Path):
    data_root = Path("validation_data/golden/gaussian")
    result = run_golden_gaussian(data_root)
    # Write report
    out = tmp_path / "reports/golden/gaussian.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("{}")
    assert result["suite"] == "gaussian"
    assert "pass_rate" in result


