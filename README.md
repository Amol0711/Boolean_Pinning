# Boolean Pinning and Poisson Limits

Reproducible code, frozen observations, exact biological benchmarks, and
independent computational checks for permanent pinning in Boolean networks.

Redistribution terms for the project material and biological model inputs
require clarification; see [LICENSE_STATUS.md](LICENSE_STATUS.md) and
[source credits](docs/THIRD_PARTY_SOURCES.md).

## Start here

Use Python 3.12 and the pinned numerical dependencies, NumPy 2.3.5 and SciPy
1.17.0. From the repository root, in an activated Python 3.12 environment:

```sh
python -m pip install -r requirements.txt
python tools/reproduce.py list
python tools/reproduce.py verify
python tools/reproduce.py doctor scc-validation
python tools/reproduce.py run scc-validation --output runs/scc-check-001
```

Choose a **new** output directory for every run. An in-repository destination
must be below `runs/`; an external new directory is also accepted. The runner
checks the distributed files, copies the selected workflow into
`OUTPUT/workspace/`, and executes there. Logs are in `OUTPUT/logs/` and execution
and integrity status are in `OUTPUT/run_summary.json`. Frozen repository files
remain the comparison baseline. A successful command exit is interpreted
together with its scientific report and stated verification coverage.

A C++17-capable `g++` supports accelerated numerical production. Optional figure
tasks need the TeX/plotting tools listed by their task-specific `doctor` check.
`doctor` reports availability and versions; use the pinned installation above
for reproducible numerical runs. Keep Python assertions enabled: do not use
`-O`, `-OO`, or a nonzero `PYTHONOPTIMIZE`.

## Workflows

All commands below use `python tools/reproduce.py run TASK --output NEW_PATH`.
Add `--workers N` for tasks that support parallel execution.

| Workflow | Frozen scope | Starting task |
|---|---|---|
| [SCC solver](workflows/scc/README.md) | Deterministic solver and finite independent checks | `scc-validation` |
| [Biological benchmarks](workflows/biological/README.md) | Two small coded models; CD4/macrophage context and target analysis | `biological-small-checks` |
| [Bernoulli study](workflows/bernoulli/README.md) | 250,000 original draws and SCC-based cost refinement | `bernoulli-quick` |
| [Joint-rule study](workflows/joint_rules/README.md) | 18,000 draws across nine cells, plus exact finite checks | `joint_rules-exact` |
| [J1 follow-up](workflows/j1_followup/README.md) | Separate 20,000-draw follow-up and original-cell comparison | `j1-preflight` |

See [the reproduction guide](docs/REPRODUCTION.md) for task scope, cost, outputs,
and full replay commands. [The results index](docs/RESULTS_INDEX.md) locates the
raw records, certificates, summaries, and independent checks.

## Evidence and integrity

The five workflows contain fixed scientific inputs and observations. Saved
reports describe the computations they record. Repository checks and their
coverage are documented in
[VALIDATION.md](docs/VALIDATION.md) and
[REPOSITORY_VALIDATION.json](provenance/REPOSITORY_VALIDATION.json); do not infer a fresh
full replay merely from an archived `PASS`.

The top-level [MANIFEST.json](MANIFEST.json) covers the distributed tree.
[Source lineage](provenance/SOURCE_LINEAGE.json) distinguishes original and
current file hashes; [commitments](provenance/COMMITMENTS.json) retain recorded
source/protocol/gate identities. Names and explanatory metadata have been
adapted while scientific numerical data are preserved. Hash verification
establishes integrity against a specified inventory, not authorship or a
digital signature.

Protocols, seed identities, generation order, budgets, and theoretical
parameters remain fixed. Unresolved costs retain certified bounds and feasible
witnesses; they are not omitted or replaced. No resampling or parameter refit is
authorized by diagnostic discrepancies. Sampling intervals and finite examples
do not establish an exact finite-size Poisson law or an asymptotic rate.

[Repository structure and integrity](docs/RELEASE.md) explains the layout,
reference documentation and hash representations.
