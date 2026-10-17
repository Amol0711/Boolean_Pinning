# Reproduction guide

Run the declared tasks through `tools/reproduce.py` from the repository root. The runner verifies the repository's listed files, copies the selected workflow into a new writable workspace, executes its commands there, and checks the repository again. Preserved source, protocol, seed and result files remain the reference inputs.

The executable command arrays, working directories and requirements are in [`config/tasks.json`](../config/tasks.json). The task descriptions below explain their scientific scope; the registry specifies the exact invocation.

## Environment and first checks

The full numerical environment uses Python 3.11 or newer, NumPy **2.3.5** and SciPy **1.17.0**. The audited environment was Linux, Python 3.12.14 and g++ 13.3.0. Standalone standard-library scripts use Python 3.10-compatible syntax, but the pinned NumPy/SciPy versions require Python 3.11 or newer. J1 production checks the frozen numerical package versions explicitly.

For an isolated numerical environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python tools/reproduce.py list
python tools/reproduce.py doctor bernoulli-pilot
python tools/reproduce.py verify
```

The repository's `requirements.txt` is the integration dependency list. The retained Bernoulli and the joint-rule study requirements files pin the same two versions above. Biological's current census needs NumPy; the independent biological validators and SCC checks use the standard library. The runner itself uses the standard library.

Numerical production tasks require `g++` with C++17 support. Machine-specific shared libraries are excluded from this repository; the retained C++ source builds the accelerator inside the run workspace. The archived engine has a Python fallback, but the declared production tasks require the compiler so that the intended accelerated path is available. A bounded Bernoulli pilot was independently checked after removing the copied compiled binary and rebuilding from source.

The original biological benchmark and its engine tests import the POSIX `resource` module. Use Linux, WSL or another compatible POSIX environment for `biological-small-checks` and `biological-original-replay`; native Windows support is not claimed. Optional figure tasks require `pdflatex`, `pdftoppm`, TikZ/PGFPlots and Latin Modern. `doctor` checks executable/package availability; it does not test every TeX package or guarantee compiler compatibility.

**Do not use Python `-O` or `-OO`.** Leave `PYTHONOPTIMIZE` unset, empty or `0`. Many scientific checks use assertions; the runner rejects optimized execution. Registry commands use `-B` to suppress ordinary Python bytecode writes. No source download or network access is needed after obtaining this repository and installing the dependencies.

## Running a task and reading its result

Choose a destination that does not exist:

```sh
python tools/reproduce.py run scc-validation --output runs/scc-check-01
python tools/reproduce.py run biological-census --output runs/biology-census-01
python tools/reproduce.py run bernoulli-pilot --workers 4 --output runs/bernoulli-pilot-01
```

`--output` is required. An output inside the repository must be a new directory below `runs/`; a new external directory is also accepted. Existing destinations and symlink paths are rejected. Use a different destination for each invocation. `--workers` defaults to **1** and is passed only to commands with a worker argument; it does not change sample counts. J1 production accepts 1–8 workers.

Each destination contains:

- `run_summary.json`: task, resolved command arguments, environment, timing, exit status and before/after baseline checks.
- `logs/`: numbered stdout and stderr logs for every executed command.
- `workspace/`: the complete copied workflow and its generated reports, data and build products.

Output locations in the tables below are relative to this `workspace/`. A runner status of `passed` means the declared commands completed successfully and the protected baseline remained intact. Inspect the workflow's scientific reports as well: a finite validation is not a proof of the general theorem, and a valid capped result can retain an unresolved cost interval. Commands stop at the first failure. Partial results and logs remain available; retry in a new destination after correcting the cause.

## Task catalogue

The times below are orientation values from the recorded verification or preserved recorded runs. They exclude some copying, integrity checks, compiler startup and comparison overhead. They are not performance guarantees; some recorded checks ran concurrently. Full-study times marked **recorded** are not claims that the central runner has completed those studies.

### SCC solver and biological computations

| Task | Computation and scope | Main outputs in `workspace/` | Cost guidance |
|---|---|---|---|
| `scc-validation` | Both examples; exhaustive small solver checks, cap checks and independent finite reduction/SCC checks | Example JSON in command logs; `results/solver_validation.json`, `results/targeted_cap_check.json`, `results/independent_theorem_check.json` | About 3.3 s in the recorded verification; standard library |
| `biological-small-checks` | Regenerate both small-model certificates, validate their engine, independently audit the source exports and test the larger-network engine | Five files in `reproduced_small/`; `biological/provenance/independent_export_audit.json`; engine test JSON in logs | Under 1 s in the recorded verification; POSIX `resource` required |
| `biological-census` | Regenerate the 50-target census over 11 contexts, independently validate the fresh census, then check the selected fair-schedule examples | `reproduced_census/all_targets.json`, two CSV summaries, `independent_biological_validation.json`, `async_memory_check.json` | Producer about 1.3 s, independent validator 6.3 s, async check 0.3 s; NumPy for production |
| `biological-original-audit` | Independently check the preserved original capped selected-target benchmark | `reproduced_original_audit.json` | Recorded audit about 85 s; 4,210 candidate sets and 48,038,912 compatible-state checks |
| `biological-original-replay` | Replay the frozen original benchmark budgets and independently audit its fresh records | `reproduced_original_benchmark/`; `reproduced_original_audit.json` | Recorded production analysis about 65 s plus independent audit; POSIX `resource` required |

The current census contains 46 CD4 targets and four macrophage targets. Its fixed-point enumeration tests 805,380 assignments. The current producer retains the original selected-target certificates as agreement inputs, so `larger_networks/results/` is a required dependency. The async check uses selected original OFF/APC targets; it is not an asynchronous verification of all 50 census targets.

The original capped benchmark is a separate selected-target calculation. Its protocol allows residual size 18, additional pin count 6, 30,000,000 charged work units, 2,000 target-discovery updates and nominal analysis time 120 seconds per context. These are algorithmic budgets, not a hard whole-process runtime bound. Resource measurements vary by hardware. An earlier wall cutoff can change the reached search bound while leaving it valid; inspect such a change rather than treating it as exact reproduction. The independent audit report also contains elapsed time and is not expected to be byte-identical as a whole.

The older selected-target APC record retains an incomplete minimizer-list status, while the later census completes that list. Likewise, the original macrophage discovery orbit and the later exhaustive no-fixed-point conclusion have different scopes. Replaying the older protocol does not replace the current census.

### Bernoulli study

| Task | Computation and scope | Main outputs in `workspace/` | Cost guidance |
|---|---|---|---|
| `bernoulli-quick` | Independent tiny-network oracle plus checks of archived numerical plot data | Oracle JSON in logs; `qa/bernoulli_figures/plot_data_validation.json` | About 0.25 s; standard library |
| `bernoulli-statistics` | Regenerate six statistical exports in the copied canonical location and independently audit them against original/refined/SCC records | `results/bernoulli_statistics/`; `results/reproduced_independent_statistics.json` | About 3 s for exports and 14 s for independent audit; NumPy/SciPy |
| `bernoulli-boundary` | Regenerate all 10,001 binomial success counts and 200 chi-square boundary cases, then independently validate them with the archived statistics | `results/bernoulli_statistics/numerical_boundary_cases.json`; `results/reproduced_boundary_validation.json` | About 23 s for regeneration plus 14 s for independent audit; NumPy/SciPy |
| `bernoulli-pilot` | Replay the first 100 frozen identities in each of 25 cells and compare the scientific files | `simulation/results/reproduced_pilot/`, including `reproduction_comparison.json` | About 3.2 s with four workers in the recorded verification; NumPy and C++ compiler |
| `bernoulli-reproduce` | Replay all 250,000 frozen identities and compare scientific files | `simulation/results/reproduced_full/`, including `reproduction_comparison.json` | Recorded main computation about 228 s with four workers; NumPy and C++ compiler |

Bernoulli has five mixtures, five sizes `(100, 400, 1600, 6400, 25600)` and 10,000 draws per cell. The original, refined and SCC records are three bound calculations of the same **250,000** networks. The 2,500 pilot identities overlap those networks; they are not additional observations. Only the `pilot_10m` reference directory is shipped, although the simulation README also discusses an earlier pilot.

`bernoulli-statistics` has two ordered operations. The first replaces six exports only in the workspace's `results/bernoulli_statistics/`: `summary.json`, `cell_summary.csv`, `cdf_intervals.csv`, `pmf_intervals.csv`, `residual_survival.csv` and `local_slopes.csv`. The second independently checks that same canonical location, so it validates the fresh exports. Similarly, `bernoulli-boundary` replaces the copied boundary-case export before independent verification. The reference files under the repository's `workflows/bernoulli/` remain unchanged. The underlying independent script has no option to select an arbitrary alternate statistical-export directory.

The independent statistics audit covers 250,000 unique draws across 750,000 analysis rows. It does not regenerate every source network. The production replay wrappers compare raw/detail scientific streams and operational summaries, excluding timing and backend metadata. Independent dynamical revalidation is a separate optional operation below.

Retain all 792 unresolved intervals in the full study. The statistical exports distinguish computational ranges from sampling uncertainty. The Poisson values are fixed asymptotic references; no task refits them. `bernoulli-quick` checks numerical plot records only: the standalone Bernoulli package omits biological figure assets that appear in baseline figure metadata.

### Joint-rule study and the J1 follow-up

| Task | Computation and scope | Main outputs in `workspace/` | Cost guidance |
|---|---|---|---|
| `joint_rules-exact` | Independent bounded exact checks, validation of archived exact exports and exact-only reproduction/comparison | `qa/reproduced_exact_checks.json`, `qa/independent_exact_exports.json`, `validation/results/reproduced_exact/` | About 2.4 s combined; NumPy/SciPy are imported by the production module even in exact-only mode |
| `joint_rules-statistics` | Independently verify all nine archived joint-rule cells | `qa/reproduced_joint_statistics.json` | About 0.7 s; NumPy/SciPy |
| `joint_rules-reproduce` | Reproduce exact exports and all 18,000 frozen draws, then compare scientific files | `validation/results/reproduced_full/exact/`, `main/`, and `reproduction_comparison.json` | Recorded main computation about 7.7 s with four workers, excluding exact/comparison overhead; NumPy/SciPy and C++ compiler |
| `j1-preflight` | Independent finite preflight plus the production no-refit regression on preserved data | `qa/reproduced_preflight.json`, `qa/reproduced_no_refit.json` | About 2.4 s combined; NumPy/SciPy |
| `j1-statistics` | Independently verify the separate follow-up and regenerate its statistical CSVs | `qa/reproduced_statistics.json`; `results/main/statistical_summary.csv`, `rule_frequency_summary.csv` | About 2.3 s; NumPy/SciPy |
| `j1-reproduce` | Reproduce the 20,000 frozen follow-up identities, then compare scientific files, completion/provenance and ordered identities | `results/reproduced_full/`, including `reproduction_comparison.json` | Recorded production about 49.6 s with four workers, excluding comparison; exact frozen NumPy/SciPy versions and C++ compiler |

the joint-rule study has three rule-law families, sizes `(100, 1600, 25600)` and 2,000 draws per cell: **18,000** observations. J1 is a separate **20,000**-draw follow-up at size 25,600. Its copied `original/` data contain the 2,000-row comparison cell from the joint-rule study. Those 2,000 rows are not another dataset to add to the joint-rule study or pool with the follow-up.

The joint-rule study's reproduction wrapper compares the exact outputs and frozen simulated scientific records. J1's additional `verify_reproduction.py` compares five scientific files, both compressed and decoded stream payloads, all 20,000 ordered seed/network identities, and completion/integrity/provenance fields. Timing, shards and non-scientific run metadata are not byte-compared. Neither comparison performs the full independent dynamical optimization again.

The the joint-rule study and J1 statistics scripts read fixed canonical main-result locations. J1 statistics rewrites two CSVs in its copied `results/main/`, even when the JSON report is redirected. This is one reason to use the central isolated runner. These statistics tasks audit the supplied reference results; they do not automatically find a prior run's fresh production directory.

### Optional figures

| Task | Scope | Outputs |
|---|---|---|
| `biological-figures` | Regenerate the current census figure and original benchmark figure from archived validated records | PDFs and plotting data under `figures/`; build/check reports and previews under `qa/` |
| `bernoulli-figures` | Regenerate the three numerical figures from canonical statistical exports | `figures/fig1_cost_distribution.pdf`, `fig2_certificate_failure.pdf`, `fig3_residual_and_moments.pdf`, plotting data and `qa/bernoulli_figures/` reports/previews |

Both require the TeX/Poppler dependencies described above. They use their own copied reference inputs, not outputs from another run directory. The builders can replace plotting data and remove obsolete exports inside the workspace. Rebuilt PDFs can depend on the installed TeX environment; scientific data comparison and rendered-figure inspection are distinct checks.

## Optional independent network revalidation

The registry deliberately separates production replay from the more extensive independent network verifiers. After a successful full production task, run the relevant commands **inside that run's `workspace/`**, never under `workflows/`. The commands below are additional manual operations and are not included in the already-written central `run_summary.json`. Their own JSON reports record the scientific results.

For Bernoulli, using the output example `runs/bernoulli-full-01`:

```sh
cd runs/bernoulli-full-01/workspace
python -B verification/validate_bernoulli_networks.py --results simulation/results/reproduced_full --workers 4 --output results/independent_bernoulli_validation.json
python -B verification/targeted_component_check.py --results simulation/results/reproduced_full --output results/reproduced_targeted_validation.json
```

The first command intentionally replaces only the copied canonical independent report. The targeted check reads that fixed path even when its own output is redirected, so this ordering links the targeted completion to the fresh full verification. The full verifier uses caps of 18 coordinates, 100,000 pin candidates and 20 million compatible states; the targeted 17-node component uses a 100-million-state cap. For the frozen pilot, `--expected-rows 100` is available with the pilot result directory; it is not an interface for arbitrary new sampling.

For the joint-rule study, from its full run workspace:

```sh
python -B verification/validate_joint_simulations.py --results validation/results/reproduced_full/main --workers 4 --output qa/reproduced_network_validation.json
```

This validates all nine cells with the frozen protocol and gate; its default whole-residual cap is 12. `--expected-rows` does not slice the released 2,000-row streams: a smaller value requires a corresponding valid complete nine-cell dataset.

For J1, from its full run workspace:

```sh
python -B verification/validate_followup_simulations.py --results results/reproduced_full --workers 4 --output qa/reproduced_network_validation.json
python -B verification/test_production_no_refit.py --results results/reproduced_full --output qa/reproduced_full_no_refit.json
```

The verifier requires all 20,000 frozen rows, receipt/source/seed linkage and the frozen independent preflight. Its component search caps are 18 coordinates, 100,000 candidates and 100 million compatible states, with whole-residual cap 12. `--chunk-size` changes batching, not sample size. Keep the no-refit input directory inside the copied package because the script records package-relative paths.

To rerun a joint-rule/J1 statistics checker on newly produced records, first stage the complete fresh result set at the checker's canonical main-result location **in another disposable complete workflow copy**, preserving protocol, source, gate and baseline linkage. Do not move or overwrite the reference records in this repository. The direct statistics scripts do not provide a general fresh-follow-up `--results` option.

## Frozen design, budgets and provenance

A successful replay reruns the same seeded networks. It is not a new independent sample or a fresh validation cohort. The frozen production CLIs do not offer arbitrary new seeds, counts or model families. A new study needs a separate protocol identity and seed namespace; reducing worker count or selecting the frozen pilot does not create one.

Bernoulli, the joint-rule study and J1 production use per-complex-component caps of size 18, pin count 18 and charged work 10,000,000. These are not global runtime or memory bounds. Changing worker count, compiler, platform or numerical package versions can change elapsed time and resources; fixed scientific payloads and identities are the reproduction targets. Use recorded status/bounds and declared comparison exclusions instead of requiring timing records to match.

The top-level [`MANIFEST.json`](../MANIFEST.json) describes the distributed tree.
`verify` checks file sizes and SHA256 hashes. [Source lineage](../provenance/SOURCE_LINEAGE.json)
and [commitments](../provenance/COMMITMENTS.json) distinguish source records from
current file representations. Names, locators and explanatory metadata have been
adapted, so source and current hashes may differ while scientific numerical data
agree. A recorded protocol or gate commitment is not silently replaced by a
current representation hash. These checks are not signatures or proof of authorship.

The raw scientific records and their frozen provenance must not be edited to make a check pass. The runner copies only manifest-covered workflow files, verifies the copies and compares the repository baseline before/after execution. Keep package layouts and vendored engine/solver copies intact: changing a shared import or consolidating byte-identical helpers can still break source-hash and freeze obligations.

Bernoulli retains a disclosed pilot-gate hash mismatch: the original pre-main report bytes were overwritten and were not recovered. The preservation reports document the discrepancy. Later broader validation does not restore that missing byte-level provenance; neither this repository manifest nor a fresh replay should be described as doing so.

## Specialist utilities outside the task registry

These files remain provenance or specialist utilities, not normal reproduction tasks:

| File within its workflow | Reason |
|---|---|
| Biological `biological/provenance/build_inventory.py` | Source-inventory construction utility; use its declared inputs rather than a reproduction task |
| Biological `biological/src/extract_models.py` | Source extraction utility; the supplied JSON/BNET models and pinned source files provide reproduction inputs |
| Biological `biological/provenance/source_files/*.py` | Preserved learning/topology source inputs; read relevant rules via AST rather than executing or importing these scripts |
| Bernoulli `qa/bernoulli_compute_theorem_constants.py` | Directed-arithmetic evaluation of the explicit finite-size theorem constants |
| Bernoulli `simulation/legacy/simulation/pilot_timing.py` | Old calibration procedure, not the frozen scientific outcome replay |
| Bernoulli `simulation/legacy/simulation/reproduce.py` | Earlier study workflow; not the current SCC replay |
| J1 `src/run_followup.py freeze` | Creates a new freeze receipt/seed manifest; reproduction uses `run`. Never delete or recreate the preserved receipt |
| J1 `src/joint_validation.py` as a standalone CLI | Retained the joint-rule study module imported by J1; its old standalone CLI expects an absent the joint-rule study protocol at a different root |

Some specialist scripts have fixed output paths, and provenance records retain source acquisition locators. Do not discover tasks by executing every Python file or by assuming every script is a reproduction entry point. In particular, the original Biological audit needs `--protocol`; the registry supplies it. The SCC theorem checker has no `--output` option, and its validator always regenerates `verification/fixtures.json`. Isolation contains these writes without changing the frozen scripts.

Use [`RESULTS_INDEX.md`](RESULTS_INDEX.md) for result locations and [`computational_details/`](computational_details/) for the retained detailed protocols, uncertainty definitions, biological tables and source qualifications.
