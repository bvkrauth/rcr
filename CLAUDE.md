# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Implements the estimator from Krauth (2016), "Bounding a linear causal effect using relative
correlation restrictions" (*Journal of Econometric Methods*), in three languages that all share
one computational core:

- **`python/`** — the `rcrbounds` module (`rcrbounds.py`). This is the canonical implementation;
  everything else calls into it.
- **`stata/`** — the Stata package `rcr` (`rcr.ado` + helpers). At runtime it shells out to either
  the Python module (`rcrbounds.py`, preferred on Stata ≥16 with Python configured) or a legacy
  Fortran executable (`rcr.exe` / `rcr`), communicating via plain text input/output files. See
  "Cross-language architecture" below.
- **`r/rcrbounds/`** — the R package `rcrbounds`. It calls the Python module through `reticulate`
  (see `r/rcrbounds/R/rcr.R`); it has no independent computation of its own.
- **`src/`** — deprecated Fortran 95 source for the legacy `rcr`/`rcr.exe` executable, kept only
  for Stata versions without Python support. Not under active development.
- **`bin/`** — prebuilt deprecated Fortran executables (`rcr.exe` for Windows, `rcr` for Linux; no
  macOS build).

Each folder has its own `README.md` (or `readme.md`) with installation/usage details — read the
relevant one before working in that folder.

## Cross-language architecture

- The Python `rcrbounds.py` module serves double duty: as an importable module (`RCR` →
  `.fit()` → `RCRResults`, with `.summary()` / `.rcrplot()`), and as a standalone script invoked
  by Stata (entry point `stata_exe()` near the bottom of the file). When run as a script it reads
  a moment-vector input file and writes result/log/detail files — see `read_data`,
  `write_results`, `write_details` in `rcrbounds.py`.
- `stata/rcr.ado` decides at runtime whether to call Python (`python script rcrbounds.py ...`) or
  the Fortran executable (`winexec`/`shell`), based on Stata version and Python availability. Use
  `rcr_config` (Stata command, `stata/rcr_config.ado`) to diagnose which backend will be used.
- The R package never reimplements the math — `install_rcrpy()` installs the Python `rcrbounds`
  package via `reticulate`, and `rcr()` calls into it directly.
- Because of this, a change to the core estimation logic almost always belongs in
  `python/rcrbounds.py`; the Stata and R layers are thin wrappers plus argument/result marshaling.

## Branching/release workflow (see `local-notes.md`)

- Day-to-day work happens on `devel`; `master` is the stable/release branch.
- Before pushing to `devel`: run the full local test and lint suites (Python).
- Before merging `devel` → `master`: rely on the remote GitHub Actions test/lint runs.
- Releases are cut per-language (see "Creating distributions" in `local-notes.md` for the full
  Python PyPI / Stata install procedure) — these are manual, human-driven steps, not something to
  automate unprompted.

## Common commands

### Python (`python/`)

```
cd python
pip install --editable .          # editable install of rcrbounds
python -m pytest                  # run the full test suite (also runs from repo root)
flake8                            # lint (repo's .flake8 rules used in CI: --max-complexity=10 --max-line-length=127)
pylint rcrbounds.py rcr_example.py # must be run from python/ with explicit filenames
```

Run a single test file with `pytest testing/test_rc.py`, or a single test with
`pytest testing/test_rc.py::test_name`. Tests live in `python/testing/`; test file names are
abbreviations of the function/feature under test (e.g. `test_em.py` → `estimate_model`,
`test_rc.py` → `read_data`/`rc_*`). Some tests require network access to fetch
`rcr_example.dta` from `sfu.ca`.

CI (`.github/workflows/test-python-package.yml`, `pylint.yml`) runs pytest + flake8 on
Windows/macOS/Linux across Python 3.9 and 3.14, and pylint on Ubuntu/Python 3.13.

### R (`r/rcrbounds/`)

Standard R package conventions (devtools/testthat). Tests are in
`r/rcrbounds/tests/testthat/test-rcr.R`. Requires the Python `rcrbounds` module to be available via
`reticulate` (run `rcrbounds::install_rcrpy()` once per environment) — most examples are guarded
with `@examplesIf reticulate::py_module_available("rcrbounds")`.

### Stata (`stata/`)

There is no automated Stata test runner; testing is interactive, from within Stata:

```stata
do rcr_certification      // in stata/testing — full package certification suite
do rcr_example            // in stata/ — basic smoke test
```

Individual `stata/testing/test_*.do` files reproduce specific bugs/features (see
`stata/testing/README.md` for the current list of known-bug repros, e.g. `test_bug012.do`,
`test_bug013.do`, `test_bug014.do`).

### Fortran (`src/`, deprecated)

Windows: `compile.bat` (Intel Fortran Compiler + MKL). Linux: `compile.linux` (gfortran +
OpenBLAS/Netlib LAPACK). Do not extend this code for new features — it's maintained only for
Stata versions that lack Python support. See `src/README.md` for compiler/library setup details.
