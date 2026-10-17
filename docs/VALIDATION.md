# Repository verification

All 18 declared tasks passed using the distributed code and command contracts. Repository integrity checks passed before and after each run. The [machine-readable record](../provenance/REPOSITORY_VALIDATION.json) contains task results, executed source hashes and links to deterministic replay comparisons.

## Scientific replay and preservation

- The Bernoulli producer replayed all 250,000 frozen draws. All 50 compressed raw/detail streams and both operational summaries matched their recorded scientific payloads. All 792 unresolved intervals remain included.
- The joint-rule producer replayed all 18,000 frozen draws and the bounded exact enumeration. All 23 compared scientific files matched.
- The separate J1 follow-up replayed all 20,000 frozen draws. All five compared scientific files matched; ordered seed/network identities and completion, integrity and provenance checks passed.
- Biological checks reproduced the small-model certificates, complete 50-target census and original 11-context capped benchmark. Independent census and original-candidate audits, selected fair-coordinate checks and finite SCC solver/oracle checks passed.
- Statistical reaggregation, special-function boundary checks, the no-refit regression and both figure-building tasks passed.

The replays use the same frozen samples; they are not new observations. The complete independent numerical network/dynamical search suites were **not rerun**. Their archived reports retain their original scope, and optional rerun commands are documented in [REPRODUCTION.md](REPRODUCTION.md).

All 172 compressed observation streams remain byte-identical to the input package. All 455,056 numeric cells across 89 CSV files and 599,352 numeric, Boolean or null fields across 65 scientific JSON files were preserved under explicit metadata-key renaming. Four Boolean rule files and four upstream Python source files remained byte-identical. Scientific model payloads, targets, seeds, budgets and statistical definitions were preserved.

Descriptive names and explanatory metadata replace internal labels. Three figure PDFs had metadata-only edits; their rendered pixels matched the originals exactly. The remaining five PDFs retained their original bytes. The complete file scan covered relative names, text, decompressed data, PDF text, metadata and hidden string objects. No unapproved internal labels or version-status wording remained. Genuine scientific variables, cited source identifiers and standard URIs retain their meanings.

## Integrity and execution scope

Sixteen runner/manifest integration tests passed, and all 73 Python sources parsed successfully. The final command definitions and executed code hashes match the distributed files. [Source lineage](../provenance/SOURCE_LINEAGE.json), [reference lineage](../provenance/REFERENCE_LINEAGE.json) and [recorded commitments](../provenance/COMMITMENTS.json) distinguish original byte commitments from current descriptive representations. No rewritten file is claimed to have existed at the historical freeze time.

Validation used Linux x86_64, Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and g++ 13.3.0. Optional figures used installed TeX/Poppler tools. This records the tested installed environment; no clean dependency installation, hosted GitHub Actions run or cross-platform coverage is claimed. Resource timings are operational observations, not comparative performance measurements.

The manifest and provenance maps establish integrity against a specified inventory, not authorship or a digital signature. [Redistribution terms](../LICENSE_STATUS.md) and [source credits](THIRD_PARTY_SOURCES.md) apply.
