# Computational reproduction

## Environment and task selection

Use Python 3.12 with the pinned packages in `requirements.txt`. Numerical production
uses NumPy 2.3.5 and SciPy 1.17.0. A C++17-capable `g++` enables the accelerated
backend. SCC checks use the Python standard library. Figure tasks require the TeX
and Poppler programs reported by their environment check.

```sh
python -m pip install -r requirements.txt
python tools/reproduce.py list
python tools/reproduce.py doctor
python tools/reproduce.py verify
```

`list` reads the task registry in `config/tasks.json`. `doctor TASK` checks the
programs required by that task. `verify` checks file sizes and hashes without
executing a scientific workflow. Keep assertions enabled in all Python commands.

## Run a task

```sh
python tools/reproduce.py doctor scc-validation
python tools/reproduce.py run scc-validation --output runs/scc-check-001
```

The destination must not already exist. Use a new directory below `runs/` or a
new external directory. The runner verifies the package, copies the selected
workflow into `OUTPUT/workspace/`, executes the registered commands there, and
verifies package integrity again. It stops when a command fails and preserves
available logs and outputs. Retry with a new destination after resolving the cause.

`OUTPUT/run_summary.json` records command completion and integrity checks.
`OUTPUT/logs/` contains stdout and stderr. Workflow-specific results are under
`OUTPUT/workspace/`. Inspect those results for unresolved bounds, guards and exact
coverage; a completed command does not mean that every search resolved exactly.

Add `--workers N` for tasks with parallel work. The declared draw identities and
replicate order remain fixed. Use the seed and identity comparisons in the workflow
reports to assess a replay; runtime and resource measurements are machine-dependent.

## Available computations

| Task family | Purpose | Key records |
|---|---|---|
| `scc-validation` | Examples, independent state-space checks, component-budget checks and reduction checks | `results/solver_validation.json`, `results/targeted_cap_check.json`, `results/reduction_validation.json` |
| `biological-small-checks` | Small-model certificates, source-export comparisons and engine checks | `biological/results/`, `biological/provenance/` |
| `biological-census` | Complete context-target census and independent certificate verification | `results/all_targets.json`, `results/independent_biological_validation.json`, `results/async_memory_check.json` |
| `biological-original-audit`, `biological-original-replay` | Check or regenerate the capped 11-context benchmark | `larger_networks/results/`, `larger_networks/verification/` |
| `bernoulli-quick` | Small-network oracle and plotting-data consistency checks | `verification/`, `qa/bernoulli_figures/` |
| `bernoulli-statistics`, `bernoulli-boundary` | Complete-sample statistics and special-function boundary checks | `results/bernoulli_statistics/`, `results/independent_statistics.json` |
| `bernoulli-pilot`, `bernoulli-reproduce` | Reproduce the 2,500-draw subset or all 250,000 observations | `simulation/` |
| `joint_rules-exact`, `joint_rules-statistics`, `joint_rules-reproduce` | Exact finite comparisons, nine-cell statistics or full 18,000-draw replay | `validation/`, `qa/` |
| `replication-preflight`, `replication-statistics`, `replication-reproduce` | Independent J1 finite checks, statistics or 20,000-draw replay | `workflows/joint_replication/` |
| `biological-figures`, `bernoulli-figures` | Compile workflow-specific plots from retained inputs | Workflow `figures/` directories |

Paths in the last column are relative to the selected workflow unless shown from
the repository root. The registry is the executable source for command arguments.
Each workflow README supplies its schemas, update conventions, search budgets,
statistical definitions and optional full verification commands.

## Data conventions

Permanent pins fix selected coordinates to their target values from the initial
update onward. Baseline clamps and additional pins are recorded separately.
Repeated input slots read the same physical coordinate during state updates.
Independent checks distinguish nominal-slot support from physical-variable support.

Cost intervals represent incomplete optimization under recorded budgets. Their
endpoints and feasible upper witnesses remain part of the data. Sampling intervals
are separate quantities, computed using all declared observations. Poisson
comparison parameters are calculated from the declared sampling laws and are not
estimated by fitting the observed costs.

The Bernoulli subset overlaps the full sample. The independent J1 sample and its
comparison cell are separate batches and are not pooled. Preserve seeds, draw
identities, declared probabilities and comparison denominators when rerunning or
summarizing the studies.

## Figures and interpretation of checks

```sh
python tools/verify_figures.py
python tools/build_figures.py --output runs/computational-figures-001
```

The root figure commands validate retained plotting data, displayed coordinate
maps and generated sources, then compile figures in a new workspace. See
[figures/README.md](../figures/README.md). These commands do not execute a network
sample or repeat a biological target census.

[RESULTS_INDEX.md](RESULTS_INDEX.md) locates the scientific records.
[VALIDATION.md](VALIDATION.md) explains the distinct roles of integrity checks,
production replay, independent dynamical checks and statistical reaggregation.
