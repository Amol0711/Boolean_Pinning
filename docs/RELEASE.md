# Repository structure and integrity

The five workflows contain solver code, model inputs, protocols, fixed seeds,
observations, certificates, statistical summaries and independent checks.

| Location | Purpose |
|---|---|
| `workflows/scc/` | SCC solver, fixtures and finite independent checks |
| `workflows/biological/` | Model rules, target census, certificates and verification |
| `workflows/bernoulli/` | Bernoulli observations, successive bounds and statistics |
| `workflows/joint_rules/` | Joint-rule study and exact validation |
| `workflows/j1_followup/` | Separate J1 follow-up and baseline-cell comparison |
| `config/tasks.json`, `tools/` | Declared commands and isolated execution |
| `docs/computational_details/` | Detailed methods, tables and source qualifications |
| `provenance/` | Source lineage, recorded commitments and verification coverage |
| `MANIFEST.json` | Distributed-file sizes and SHA256 hashes |

Names, explanatory metadata and documentation have been adapted to this layout.
Scientific numerical observations, target values, model rules, seeds and
statistical definitions are preserved. This is not a claim that every file has
the bytes of its source package: renamed paths and rewritten metadata change
file hashes, and figure metadata can change without changing rendered data.

The [top manifest](../MANIFEST.json) verifies the distributed representation.
[Source lineage](../provenance/SOURCE_LINEAGE.json) records original and current
file hashes. [Commitments](../provenance/COMMITMENTS.json) retain the recorded
source/protocol/gate commitments and distinguish them from current file hashes.
A commitment to original bytes must not be silently reinterpreted as a
commitment to rewritten bytes. Hash checks establish integrity against a
specified record; they are not signatures or evidence of authorship.

[VALIDATION.md](VALIDATION.md) identifies checks actually performed and their
coverage. Saved validation reports and a new execution are distinct evidence.
Production replay, independent dynamical verification, statistical reaggregation
and integrity checking each establish different properties.

Use the [runner](REPRODUCTION.md) with a new output path. It copies the selected
workflow into an isolated workspace, retains stdout/stderr and checks the
repository before and after execution. Runtime measurements can vary while
scientific records agree. Keep all observations, caps, stop reasons, cost
intervals and feasible witnesses. Do not resample, discard difficult draws,
refit limiting parameters or pool the J1 batches.

[Redistribution terms](../LICENSE_STATUS.md) and
[source credits](THIRD_PARTY_SOURCES.md) apply to the complete collection,
including model representations in derived records.
