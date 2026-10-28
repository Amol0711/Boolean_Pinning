# Validation coverage

The package provides four complementary kinds of checks.

| Check | What it establishes |
|---|---|
| File integrity | Distributed sizes and hashes agree with `MANIFEST.json` |
| Production replay | Declared seeds and inputs regenerate the compared scientific outputs |
| Independent computation | Separate state-space, graph, model or statistical implementations agree on the tested cases |
| Figure verification | Plotting records, displayed maps and generated sources agree with retained numerical inputs |

Use `python tools/reproduce.py verify` for integrity and
`python tools/reproduce.py list` to select a scientific task. New executions record
their commands, exit status and before/after integrity checks in
`OUTPUT/run_summary.json`. The invoked scientific reports state their detailed
coverage, unresolved cases and numerical comparisons.

Saved reports document the computations that produced them. Their source hashes,
protocols and input identities define that scope. A saved PASS does not by itself
record execution on another machine or with changed sources. Consult
workflow-specific reports for recorded execution details and associated
comparisons. `provenance/COMPUTATIONAL_INVENTORY.json` identifies the current
computational files and their hashes.

The SCC tests enumerate compatible physical states and pin sets for finite cases,
check independent support and component calculations, and examine incomplete search
budgets. Biological checks compare source encodings, target lists and certificates.
Numerical checks preserve draw identities and recompute complete-sample summaries.
Unresolved costs remain intervals and participate in the endpoint statistics.

Root figure checks use retained records and finite displayed maps. They do not
repeat the network studies or biological census. Native builds use a separate
workspace and report their own comparisons.

Numerical environments, resource measurements and archived output formats are
documented with the workflows. Reproduction of scientific payloads and reproduction
of elapsed times are different comparisons. Hash verification establishes file
integrity against an inventory; it does not authenticate authorship.

See [REPRODUCTION.md](REPRODUCTION.md) for execution instructions and
[RESULTS_INDEX.md](RESULTS_INDEX.md) for retained records. Biological source credits
and rights notices are retained in [THIRD_PARTY_SOURCES.md](THIRD_PARTY_SOURCES.md)
and [LICENSE_STATUS.md](../LICENSE_STATUS.md).
