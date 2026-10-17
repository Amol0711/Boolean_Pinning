# Results index

These paths locate the fixed scientific evidence distributed with the five
workflows. Recorded evidence and repository verification are distinguished in [VALIDATION.md](VALIDATION.md) and
[the repository verification record](../provenance/REPOSITORY_VALIDATION.json). Task definitions and new
output locations are documented in [REPRODUCTION.md](REPRODUCTION.md).

## SCC solver

Paths are relative to [`workflows/scc/`](../workflows/scc/).

| Evidence | Location |
|---|---|
| Solver input schema and two worked inputs | `input_schema.json`, `example_network.json`, `example_two_slot.json` |
| Independent physical-state/support/graph validation | `results/solver_validation.json` |
| Independent finite theorem and eligibility checks | `results/independent_theorem_check.json` |
| Targeted component-budget checks | `results/targeted_cap_check.json` |
| Oracle coverage | `verification/README.txt` |

Solver `lower`/`upper` bounds concern **additional** pins beyond a declared
baseline. An unresolved `exact_cost` is null; its feasible upper witness and
stop reason remain available. These finite checks are separate from the large
synthetic studies.

## Biological benchmarks

Paths are relative to [`workflows/biological/`](../workflows/biological/).

| Evidence | Location |
|---|---|
| Small-model rules, node orders, and full targets | `biological/models/` |
| Small-model certificates and implementation validation | `biological/results/` |
| Frozen source copies, exact source comparisons, and model qualifications | `biological/provenance/` |
| CD4/macrophage source BNET and metadata | `larger_networks/data/` |
| Eleven-context protocol, assignments, and freeze record | `larger_networks/protocol/` |
| Earlier orbit-selected certificates and resource records | `larger_networks/results/CONTEXT/certificate.json`, `larger_networks/results/CONTEXT/resources.json` |
| Complete fixed-target/cost/component/minimizer records | `results/all_targets.json` |
| All context totals and one row per full target | `results/context_summary.csv`, `results/all_target_summary.csv` |
| Independent all-target audit | `results/independent_biological_validation.json` |
| Selected OFF/APC fair-coordinate-update and memory checks | `results/async_memory_check.json` |

The all-target record contains 46 CD4 targets and four macrophage targets; the
macrophage BACTERIA context has no fixed point. The JSON retains all local
minimum pin sets and isolated truth tables; the CSV summary is not a substitute
for that evidence. It also distinguishes the original zero-initialized target
selection from the later exhaustive census.

Earlier APC records established the optimum without a complete minimum-set
family; the later `all_targets.json` supplies the complete family. Earlier
cycling discovery is distinct from the later exhaustive no-fixed-point result.
Keep these scopes separate. The fair-coordinate-update result concerns the
specified OFF/APC targets and interventions, not every census target or every
update convention. Eligibility sensitivity is a separate full-state problem:
prohibiting the proliferation pin leaves seven CD4 targets feasible and 39
infeasible under the stated model and baseline conventions.

## Bernoulli study

Paths are relative to [`workflows/bernoulli/`](../workflows/bernoulli/).

| Evidence | Location |
|---|---|
| SCC replay protocol and fixed caps | `simulation/Protocol_Bernoulli.json` |
| Original sample, earlier bounds, generator, and independent checks | `simulation/legacy/` |
| Current frozen per-draw costs and residual/component evidence | `simulation/results/main/raw/*.csv.gz`, `simulation/results/main/details/*.jsonl.gz` |
| Retained pilot | `simulation/results/pilot_10m/` |
| Complete-sample statistical protocol | `verification/Bernoulli_Analysis_Protocol.json` |
| All 25 cell summaries and pointwise intervals | `results/bernoulli_statistics/summary.json`, `results/bernoulli_statistics/cell_summary.csv` |
| Full independent network audit and targeted guard resolution | `results/independent_bernoulli_validation.json`, `results/independent_bernoulli_targeted.json` |
| Independent statistics | `results/independent_statistics.json` |
| Pilot-gate limitation | `results/pilot_gate_preservation.json`, `results/pilot_gate_reconstruction_check.json` |
| Plot data and consistency checks | `figures/data/`, `figures/bernoulli_plot_metadata.json`, `qa/bernoulli_figures/` |

The 250,000 retained observations comprise 249,208 exact costs and 792 certified
intervals in the frozen final summaries. All draws contribute to the reported
endpoint statistics. The 2,500-draw pilot is part of the original sample, not a
separate sample to add to its denominator. The original pilot-gate report bytes
were overwritten by a later extended audit; the preserved limitation record
must accompany any statement about that gate's recorded hash.

## General joint-rule study

Paths are relative to [`workflows/joint_rules/`](../workflows/joint_rules/).

| Evidence | Location |
|---|---|
| Joint corner probabilities, seeds, caps, and conventions | `validation/Protocol_JointRules.json` |
| Preserved independent exact pre-main gate | `qa/independent_exact_checks.json` |
| Production finite exact realizations and support patterns | `validation/results/exact/` |
| 18,000 raw observations and component details | `validation/results/main/raw/`, `validation/results/main/details/` |
| Nine-cell complete-sample summaries | `validation/results/main/summary.json`, `validation/results/main/summary.csv` |
| Independent exact-export, network, and statistics checks | `qa/independent_exact_exports.json`, `qa/independent_joint_replay.json`, `qa/independent_joint_statistics.json` |
| Saved interpretation and coverage summary | `validation/qa/main_results_report.md`, `qa/joint_rules_independent_validation_summary.txt` |

The frozen study resolves every cost under its declared budgets. Its original
J1 largest-size second-moment departure remains in the data and summaries;
the later follow-up does not overwrite it.

## J1 follow-up

Paths are relative to [`workflows/j1_followup/`](../workflows/j1_followup/).

| Evidence | Location |
|---|---|
| Fixed design, complete seed identities, and local freeze receipt | `Protocol_J1_Followup.json`, `seed_manifest.csv`, `freeze_receipt.json` |
| Original 2,000-draw J1 comparison cell | `original/raw/`, `original/details/` |
| Preserved production and independent preflights | `qa/production_preflight.json`, `qa/independent_preflight.json` |
| Separate 20,000-draw follow-up | `results/main/raw/`, `results/main/details/` |
| Production summaries and observed rule counts | `results/main/summary.json`, `results/main/summary.csv`, `results/main/rule_counts.json` |
| Independent batch comparison and rule diagnostics | `results/main/statistical_summary.csv`, `results/main/rule_frequency_summary.csv` |
| Independent replay, statistics, and no-refit check | `qa/independent_network_replay.json`, `qa/statistical_result.json`, `qa/production_no_refit.json` |
| Statistical qualifications | `qa/statistical_interpretation.txt` |

The frozen follow-up has 20,000 exact costs and is kept separate from its
2,000-draw comparison cell. Its `original/` directory is an excerpt of the
larger joint-rule study: baseline summaries there can cover more
cells than the excerpt contains. The local freeze receipt is not a public
preregistration. Between-batch contrasts are exploratory, and observed rule
frequencies do not refit the theoretical parameters.

## Reading a new run

For `--output PATH`, begin with `PATH/run_summary.json` and `PATH/logs/`, then
inspect the workflow-specific report in `PATH/workspace/`. The runner records
command completion and baseline integrity; exact scientific checks, remaining
guards, comparison fields, and limitations are recorded by the invoked tools.
Runtime/resource timestamps need not reproduce bytewise. Preserve the raw
records, computational bounds, and sampling uncertainty as distinct quantities.
