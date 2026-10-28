# Results index

The following paths locate retained computational inputs and outputs. Paths in
each table are relative to its named workflow. Execution instructions are in
[REPRODUCTION.md](REPRODUCTION.md); each workflow README documents the record fields.

## SCC solver

Root `workflows/scc/`.

| Records | Location |
|---|---|
| Input schema and examples | `input_schema.json`, `example_network.json`, `example_two_slot.json` |
| Independent physical-state, support and graph checks | `results/solver_validation.json` |
| Independent finite reduction and eligibility checks | `results/reduction_validation.json` |
| Component-budget checks | `results/targeted_cap_check.json` |
| Oracle definitions and coverage | `verification/README.txt` |

Costs count additional pins beyond the declared baseline. Unresolved exact costs
are null; their bounds, feasible pin sets and stopping reasons remain available.

## Biological benchmarks

Root `workflows/biological/`.

| Records | Location |
|---|---|
| Small-model rules, node order and targets | `biological/models/` |
| Small-model certificates and validation | `biological/results/` |
| Source copies, comparisons and model qualifications | `biological/provenance/` |
| CD4/macrophage BNET models and source metadata | `larger_networks/data/` |
| Eleven-context assignments and computation protocol | `larger_networks/protocol/` |
| Capped benchmark certificates and resource records | `larger_networks/results/` |
| Complete fixed-target, cost and minimum-pin-set records | `results/all_targets.json` |
| Context and target summaries | `results/context_summary.csv`, `results/all_target_summary.csv` |
| Independent census verification | `results/independent_biological_validation.json` |
| Selected fair-coordinate-update checks | `results/async_memory_check.json` |

The complete census contains 46 CD4 targets and four macrophage targets. The
macrophage BACTERIA context has no fixed point. Keep full-state target enumeration
separate from orbit-based target discovery. Fair-coordinate-update checks cover
the specified OFF/APC cases and intervention sets.

## Bernoulli study

Root `workflows/bernoulli/`.

| Records | Location |
|---|---|
| Sampling and solver protocol | `simulation/Protocol_Bernoulli.json` |
| Per-draw cost bounds and component details | `simulation/results/main/raw/`, `simulation/results/main/details/` |
| Complete-sample statistical protocol | `verification/Bernoulli_Analysis_Protocol.json` |
| Cell summaries and pointwise intervals | `results/bernoulli_statistics/` |
| Independent network checks | `results/independent_bernoulli_validation.json`, `results/independent_bernoulli_targeted.json` |
| Independent statistical checks | `results/independent_statistics.json` |
| Plotting data and consistency reports | `figures/data/`, `qa/bernoulli_figures/` |

All 250,000 observations contribute to the summaries, including 792 unresolved
cost intervals. The 2,500-draw subset belongs to this sample. Earlier bound records
and their generators preserve comparisons for the same draws; they do not supply
additional independent observations.

## Joint-rule study

Root `workflows/joint_rules/`.

| Records | Location |
|---|---|
| Probabilities, seeds, budgets and update conventions | `validation/Protocol_JointRules.json` |
| Exact finite checks and exports | `qa/independent_exact_checks.json`, `validation/results/exact/` |
| 18,000 observations and component details | `validation/results/main/raw/`, `validation/results/main/details/` |
| Nine-cell summaries | `validation/results/main/summary.json`, `validation/results/main/summary.csv` |
| Independent export, network and statistical checks | `qa/independent_exact_exports.json`, `qa/independent_joint_replay.json`, `qa/independent_joint_statistics.json` |

## Independent J1 sample

Root `workflows/joint_replication/`. Its README identifies the fixed protocol,
complete seed list, generation code and verification commands.

| Records | Location |
|---|---|
| Separate 2,000-draw comparison cell | `original/raw/`, `original/details/` |
| 20,000 independent observations | `results/main/raw/`, `results/main/details/` |
| Cost summaries and rule counts | `results/main/summary.json`, `results/main/summary.csv`, `results/main/rule_counts.json` |
| Batch comparisons and rule-frequency diagnostics | `results/main/statistical_summary.csv`, `results/main/rule_frequency_summary.csv` |
| Independent network and statistical checks | `qa/independent_network_replay.json`, `qa/statistical_result.json` |
| Invariance under diagnostic rule-frequency changes | `qa/production_no_refit.json` |

The batches have separate sample sizes and identities. The comparison-cell copy
is an excerpt of the joint-rule study, whose accompanying summaries can cover
more cells than the excerpt. Rule-frequency diagnostics do not alter declared
sampling probabilities or comparison parameters.

## Computational figures

[figures/README.md](../figures/README.md) indexes the root plot sources, PDFs and
data checks. Workflow-local figure folders supply their corresponding computation
views. Figure verification checks those retained inputs without repeating the
full simulations.
