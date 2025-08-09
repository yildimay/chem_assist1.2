## ChemAssist Validation Suite

Run golden-set verification and error-fix benchmarks for Gaussian and GROMACS.

### CLI

```
python -m chemassist.validation golden --suite gaussian --data validation_data/golden/gaussian --out reports/golden/gaussian.json
python -m chemassist.validation golden --suite gromacs  --data validation_data/golden/gromacs  --out reports/golden/gromacs.json

python -m chemassist.validation bench  --suite gaussian --data validation_data/bench/gaussian_fixes --out reports/bench/gaussian.json
python -m chemassist.validation bench  --suite gromacs  --data validation_data/bench/gromacs_fixes  --out reports/bench/gromacs.json
```

### Data layout

See `validation_data/` for examples. Case-level `expected.json` can override default tolerances.

### Reports

- JSON results at the `--out` path
- Markdown summary alongside (`summary.md`)


