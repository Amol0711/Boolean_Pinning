# Repository structure and integrity

| Location | Contents |
|---|---|
| `workflows/scc/` | Permanent-pinning solver, examples and finite independent checks |
| `workflows/biological/` | Model rules, targets, certificates, context protocols and verification |
| `workflows/bernoulli/` | Fixed observations, cost bounds, statistical analysis and validation |
| `workflows/joint_rules/` | Nine-cell joint-rule study and finite exact comparisons |
| `workflows/joint_replication/` | Independent J1 sample and separate comparison cell |
| `config/tasks.json` | Declared task commands and execution requirements |
| `tools/` | Isolated execution, integrity checks and figure compilation |
| `figures/` | Computational plots, native sources, data and consistency checks |
| `provenance/` | Current computational file inventory and hashes |
| `MANIFEST.json` | Distributed-file sizes and SHA256 hashes |

The top-level manifest describes the distributed files.
`provenance/COMPUTATIONAL_INVENTORY.json` identifies current computational inputs,
code and results. Individual workflow records identify declared sample identities
and computation inputs. Compare hashes only within the scope of their named record.

Run tasks through the [reproduction guide](REPRODUCTION.md). The runner copies the
selected workflow to a new workspace and checks repository integrity before and
after execution. Numerical observations, model rules, target assignments, seeds,
search budgets and statistical definitions determine the scientific comparison.
Runtime measurements and execution timestamps can vary between machines.

Keep cost intervals, feasible upper witnesses and stop reasons together. Keep the
two J1 samples separate, and count the Bernoulli subset only within its parent
sample. [VALIDATION.md](VALIDATION.md) describes the scope of the available checks.

[Redistribution terms](../LICENSE_STATUS.md) and
[source credits](THIRD_PARTY_SOURCES.md) apply to source files and to model
representations in derived records.
