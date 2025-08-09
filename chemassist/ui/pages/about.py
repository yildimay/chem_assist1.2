from __future__ import annotations

import json
from pathlib import Path
import streamlit as st


def _try_load_metrics() -> dict:
    metrics: dict[str, float] = {}
    try:
        # Golden pass rate: average of available golden reports
        golden_rates = []
        for p in [
            Path("reports/golden/gaussian.json"),
            Path("reports/golden/gromacs.json"),
        ]:
            if p.exists():
                data = json.loads(p.read_text(encoding="utf-8"))
                if "pass_rate" in data:
                    golden_rates.append(float(data["pass_rate"]))
        if golden_rates:
            metrics["golden_pass_rate"] = sum(golden_rates) / len(golden_rates)

        # Bench metrics: average across available bench reports
        bench_success = []
        bench_degr = []
        for p in [
            Path("reports/bench/gaussian.json"),
            Path("reports/bench/gromacs.json"),
        ]:
            if p.exists():
                data = json.loads(p.read_text(encoding="utf-8"))
                if "fix_success_rate" in data:
                    bench_success.append(float(data["fix_success_rate"]))
                if "degradation_rate" in data:
                    bench_degr.append(float(data["degradation_rate"]))
        if bench_success:
            metrics["fix_success_rate"] = sum(bench_success) / len(bench_success)
        if bench_degr:
            metrics["degradation_rate"] = sum(bench_degr) / len(bench_degr)
    except Exception:
        pass
    return metrics


def show_page() -> None:
    st.title("About")

    # Optional metrics line
    m = _try_load_metrics()
    if {
        "golden_pass_rate",
        "fix_success_rate",
        "degradation_rate",
    }.issubset(m.keys()):
        st.markdown(
            f"Golden pass rate: {m['golden_pass_rate']*100:.0f}% • "
            f"Fix success: {m['fix_success_rate']*100:.0f}% • "
            f"Degradation: {m['degradation_rate']*100:.0f}%"
        )

    # Sections
    st.markdown(
        """
        <h2>About ChemAssist</h2>
        <p>
        ChemAssist is a local, scriptable toolkit for computational chemistry workflows.
        It helps prepare inputs, analyze logs, and diagnose issues for DFT (Gaussian) and
        Molecular Dynamics (GROMACS) without requiring internet access for validation.
        </p>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <h2>Why ChemAssist is Trustworthy</h2>
        <p>
        ChemAssist includes a built-in validation suite that parses reference outputs and
        compares key metrics against case-specific tolerances. Error-fix benchmarks evaluate
        normalized patches and check for unintended changes. Results can be run in CI and are
        reported as JSON and Markdown for full transparency.
        </p>
        """,
        unsafe_allow_html=True,
    )


