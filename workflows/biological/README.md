# Exact biological benchmarks

This workflow retains two small coded Boolean implementations and the frozen
CD4/macrophage models, contexts, targets, certificates, and independent checks.
Upstream credits and unresolved redistribution terms are documented in
[THIRD_PARTY_SOURCES.md](../../docs/THIRD_PARTY_SOURCES.md).

Run tasks from the repository root:

| Task | Scope |
|---|---|
| `biological-small-checks` | Small-model certificates, source-export audit, and larger-engine tests |
| `biological-census` | Regenerate and independently verify the census; check selected fair-update certificates |
| `biological-original-audit` | Audit the earlier larger-model selected-target certificates |
| `biological-original-replay` | Reproduce the original larger-context computation and independently audit its fresh records |
| `biological-figures` | Rebuild figures from the archived validated plotting inputs |

```sh
python tools/reproduce.py doctor biological-small-checks
python tools/reproduce.py run biological-small-checks --output runs/biological-small-001
```

Use `list` and [REPRODUCTION.md](../../docs/REPRODUCTION.md) for exact task
commands, output files, and any optional figure task. The original
[README.txt](README.txt), [README_Biological.txt](README_Biological.txt), and
[small-model guide](biological/README.txt) retain the scientific definitions.

The model rules, full-state targets, external baselines, eligible coordinates,
and update semantics are part of each problem. Permanent pins take target
values. Main cost searches use synchronous physical dynamics; the additional
fair-coordinate-update analysis applies only to the specified selected OFF/APC
targets and interventions. It supplies no schedule-uniform update-time bound
or random-ensemble conclusion.

`results/all_targets.json` is the complete later census, including isolated
truth tables, all local minimizers, and restricted eligibility. Its 46 CD4 and
four macrophage targets must be distinguished from the earlier orbit-selected
targets in `larger_networks/results/`. The later record supplies complete APC
minimum families and the exhaustive macrophage BACTERIA no-fixed-point result.
The orbit-selected and all-target records retain their distinct scopes.

Keep the T-LGL Fas and source-target qualifications and the Drosophila CIA/CIR
variant distinction. Do not reinterpret a full-state target as a biological
phenotype. `biological/provenance/build_inventory.py` constructs source inventories.
Reproduction tasks use the included model exports and pinned source files.

See [RESULTS_INDEX.md](../../docs/RESULTS_INDEX.md) for evidence locations and
[VALIDATION.md](../../docs/VALIDATION.md) for repository verification coverage. The root manifest describes the
distributed tree; source and current hashes are distinguished in provenance/.
