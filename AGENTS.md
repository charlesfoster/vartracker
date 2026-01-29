# Repository Guidelines

## Project Structure & Module Organization
Core code lives in `vartracker/`, with domain modules such as `analysis.py`, `vcf_processing.py`, and the CLI entrypoint `main.py`. Reference assets and demo inputs sit under `vartracker/data/` and `vartracker/test_data/`. Pytest suites mirror the package layout in `tests/`, while helper linters live in `scripts/`. Generated artefacts such as wheels (`dist/`) and run outputs (`results/`, `vartracker_results/`) should stay out of feature branches unless you are intentionally updating them. Use `examples/` for runnable scenario notebooks or sample CSVs.

## Build, Test, and Development Commands
Install developer tooling with `pip install -e '.[dev]'` (or `uv pip install -e '.[dev]'` if you prefer `uv`). Run a local smoke test via `vartracker --test -o results/dev_run` to validate the CLI. Execute unit tests with `pytest`, and add coverage when touching analysis code using `pytest --cov=vartracker`. Keep the linters green with `pre-commit run --all-files`, which wraps Black, Flake8, and the repo-specific scripts in `scripts/`.

## Coding Style & Naming Conventions
Python files follow 4-space indentation, Black formatting (line length 88), and Flake8 with `max-line-length` 120; configure editors to respect both. Use `snake_case` for functions and variables, `PascalCase` for classes, and keep module names lowercase (matching the current package). Maintain explicit type hints in new public APIs and prefer f-strings for formatting. Run `black .` and `flake8` before opening a PR; the same hooks are enforced via pre-commit.

## Testing Guidelines
Place new tests under `tests/` with names like `test_module_feature.py`, and target behaviours rather than implementation details. Use Pytest markers sparingly—`--strict-markers` is enabled—so declare any custom marker in `pyproject.toml`. When extending CLI workflows, add end-to-end coverage beside `tests/test_cli_integration.py`. Fail fast by running `pytest tests/test_core.py -k your_case` during iteration, then finish with `pytest --cov=vartracker`.

## Commit & Pull Request Guidelines
Adopt concise, imperative commit subjects (e.g., “Add variant frequency filter”), mirroring the existing history. Each PR should link relevant GitHub issues, describe user-facing impacts, and note any new dependencies or data files. Include a “Testing” section in the description (commands run, sample outputs) and attach screenshots or tables when updating reporting artefacts. Ensure documentation changes accompany CLI or analysis adjustments.

## External Tools & Data
Several workflows rely on system binaries—install:
- `vartracker vcf`: `bcftools` and `tabix`
- `vartracker bam`: `bcftools`, `tabix`, `samtools`, `lofreq`
- `vartracker e2e`: `bcftools`, `tabix`, `samtools`, `lofreq`, `bwa`, `fastp`
Place large intermediate files outside the repo or add them to `.gitignore`, and cite any new reference datasets in `README.md` to keep provenance clear.
