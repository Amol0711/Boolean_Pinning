# Boolean Pinning Computations

Reproducible simulations, biological benchmarks and independent numerical checks
for permanent pinning in Boolean networks. The package contains five computational
workflows, their inputs and retained observations, solver certificates, statistical
summaries, and native LaTeX/TikZ/PGFPlots figure sources.

## Quick start

Use Python 3.12 and the dependencies in `requirements.txt`. From the repository root:

```sh
python -m pip install -r requirements.txt
python tools/reproduce.py list
python tools/reproduce.py verify
python tools/reproduce.py doctor scc-validation
python tools/reproduce.py run scc-validation --output runs/scc-check-001
```

Every run requires a new output directory. The runner copies the selected workflow
into `OUTPUT/workspace/`, records command output in `OUTPUT/logs/`, and writes
execution and integrity results to `OUTPUT/run_summary.json`. In-repository outputs
belong under `runs/`. The distributed observations remain available for comparison.
Keep Python assertions enabled; do not use `-O`, `-OO`, or a nonzero `PYTHONOPTIMIZE`.

The pinned numerical environment uses NumPy 2.3.5 and SciPy 1.17.0. A C++17-capable
`g++` enables accelerated simulations. Figure builds additionally need TeX and
Poppler tools. Task-specific `doctor` checks report required programs and versions.

## Workflows

| Workflow | Computational scope | Starting task |
|---|---|---|
| [SCC solver](workflows/scc/README.md) | Permanent-pinning solver and independent finite-state checks | `scc-validation` |
| [Biological benchmarks](workflows/biological/README.md) | Small coded models and CD4/macrophage target analysis | `biological-small-checks` |
| [Bernoulli study](workflows/bernoulli/README.md) | 250,000 observations and complete-sample cost statistics | `bernoulli-quick` |
| [Joint-rule study](workflows/joint_rules/README.md) | 18,000 observations in nine cells and exact finite checks | `joint_rules-exact` |
| [Independent J1 sample](workflows/joint_replication/README.md) | Separate 20,000-observation sample and 2,000-observation comparison cell | `replication-preflight` |

Use [REPRODUCTION.md](docs/REPRODUCTION.md) for execution instructions and
[RESULTS_INDEX.md](docs/RESULTS_INDEX.md) to locate observations and validation
reports. The two J1 samples have separate denominators. The Bernoulli 2,500-draw
subset is contained in its 250,000-draw sample and supplies no additional observations.
Unresolved costs retain their bounds, stop reasons and feasible pin sets.

## Figures

```sh
python tools/verify_figures.py
python tools/build_figures.py --output runs/computational-figures-001
```

These commands check plotting inputs and compile native figure sources in an
isolated directory. They do not rerun the simulations. See the
[figure guide](figures/README.md) for the available plots and dependencies.

## Integrity and source credits

[MANIFEST.json](MANIFEST.json) records distributed file sizes and SHA256 hashes.
[VALIDATION.md](docs/VALIDATION.md) distinguishes saved computation reports from
checks executed on a new run. [RELEASE.md](docs/RELEASE.md) describes the package
layout and provenance records.

Preserve the source credits and rights notices associated with biological inputs.
[LICENSE_STATUS.md](LICENSE_STATUS.md) records the applicable redistribution
questions; [THIRD_PARTY_SOURCES.md](docs/THIRD_PARTY_SOURCES.md) identifies pinned
sources, original publications and model-specific qualifications.
