from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .golden.gaussian import run_golden_gaussian
from .golden.gromacs import run_golden_gromacs
from .bench.error_fix import run_bench_errorfix
from .metrics.report import write_summary_md, write_summary_json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser("chemassist.validation")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # Golden
    g = sub.add_parser("golden", help="Run golden-set validation")
    g.add_argument("--suite", choices=["gaussian", "gromacs"], required=True)
    g.add_argument("--data", required=True, help="Path to golden data root")
    g.add_argument("--out", required=True, help="JSON output path")

    # Benchmark
    b = sub.add_parser("bench", help="Run error-fix benchmark")
    b.add_argument("--suite", choices=["gaussian", "gromacs"], required=True)
    b.add_argument("--data", required=True, help="Path to benchmark data root")
    b.add_argument("--out", required=True, help="JSON output path")

    args = parser.parse_args(argv)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.cmd == "golden":
        if args.suite == "gaussian":
            result = run_golden_gaussian(Path(args.data))
        else:
            result = run_golden_gromacs(Path(args.data))
    else:
        result = run_bench_errorfix(Path(args.data), suite=args.suite)

    # write main output
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    # write summary in reports/summary.* if conventional locations are used
    # caller tests will also assemble an overall summary from individual reports
    summary_json = out_path.parent / "summary.json"
    summary_md = out_path.parent / "summary.md"
    write_summary_json(result, summary_json)
    write_summary_md(result, summary_md)

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


