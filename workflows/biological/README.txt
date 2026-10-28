BIOLOGICAL BOOLEAN-NETWORK COMPUTATION PACKAGE
==================================================

This standalone package contains exact Boolean-model computations, frozen
source rules, certificates, verification code and editable figures. 

Python 3.10+ and NumPy are needed for the all-target producer.
The independent all-target verifier, asynchronous verifier, SCC vendor solver
and original small-model certificates use only Python's standard library.
Optional figure builds require pdflatex, TikZ, PGFPlots, Latin Modern and
pdftoppm. No network connection is needed after extraction.

Target-census reproduction
-------------------
From this extracted package root, run sequentially:

  python src/run_all_targets.py
  python verification/verify_all_targets.py
  python verification/async_memory_check.py
  python figures/fig6_build.py

The first three commands regenerate results/all_targets.json,
context_summary.csv, all_target_summary.csv, independent_biological_validation.json,
and async_memory_check.json. The figure command verifies and exports data and
rebuilds figures/fig6_cd4_target_mechanisms.pdf. To compile Figure 6 using only
its included plotting data, run:

  python figures/fig6_build.py --compile-only

For a separate production output directory without replacing the included
reference results, use:

  python src/run_all_targets.py --output new_all_target_results
  python verification/verify_all_targets.py --production new_all_target_results/all_targets.json --output new_all_target_results/independent_biological_validation.json
  python verification/async_memory_check.py --output new_all_target_results/async_memory_check.json

Figure 6's full exporter deliberately reads the fixed results/ reference
paths; regenerate those first when rebuilding a changed analysis. Scientific
outputs are deterministic. Execution logs and diagnostic previews are not
used as scientific evidence and are not included in the manifest's reproducible
output comparison. Read README_Biological.txt and verification/README.txt for full
methods, input semantics, declared search guards and verification coverage.

Target-census results
-----------
All eleven original input contexts are preserved. The complete fixed-point
enumeration gives 46 CD4 targets plus four macrophage targets. CD4 counts in
OFF, ON, APC, ECM, IL2 and GALPHAQ are 6, 6, 15, 6, 7 and 6. The actual IL2
context is CGC=IL2RB=IL2_e=1, with all other external inputs zero.

The macrophage BACTERIA context has no fixed point: every one of the 524288
assignments left after forced-value propagation violates a fixed-point
equation. This strengthens the preserved earlier record that the selected
zero-state discovery orbit is periodic. It does not modify that orbit or its
source rules. Across all eleven conditions, 805380 assignments are examined.

Every full target has its exact physical residual, essential functional
dependencies, SCCs, isolated truth tables, local minimum pin sets and global
feasible witness recorded. All ten original selected-target optima remain
unchanged. Raw marker bits are reported without assigning phenotype labels.
When the modeled proliferation coordinate is ineligible for control, 39 CD4
full-state target problems are infeasible; seven remain feasible. The separate
asynchronous result concerns the selected OFF and APC targets and convergence
under coordinate-fair one-coordinate schedules, with exact costs 3 and 4.
It supplies no uniform update-count bound under unbounded fair delays and no
random-ensemble conclusion.

Selected-target computation
---------------------------
The capped orbit-selected searches and exhaustive target census are separate
calculations. The capped APC search establishes cost four without enumerating
all minimum sets; results/all_targets.json provides the complete family. The
macrophage discovery orbit records a cycle; exhaustive enumeration establishes
the absence of a fixed point. Original budgets and numerical results are retained.

Reproduce both original small-model certificates and their validation with:

  python biological/src/reproduce.py --output /absolute/path/to/new_small_results

Run the inexpensive original larger-model implementation checks with:

  python larger_networks/src/test_runner.py

The original larger-context full searches remain available but are not needed
to reproduce the target census. To execute them explicitly in a new directory:

  python larger_networks/src/benchmark.py --protocol larger_networks/protocol/contexts.json --output /absolute/path/to/new_large_results

Independently audit the preserved original larger-model certificates:

  python larger_networks/verification/verify_results.py --protocol larger_networks/protocol/contexts.json --results larger_networks/results --output /absolute/path/to/new_original_audit.json

The original Figure 5 sources and exact cycle data are also retained. Its
figure generator can be run with:

  python scripts/build_larger_figure.py

Source and result identity
--------------------------
The repository MANIFEST.json verifies distributed files. Source lineage and
COMMITMENTS.json under the repository provenance/ directory distinguish recorded
source/protocol identities from hashes of this representation. Documentation
and locators have been adapted; scientific numerical records are preserved.
Use the repository-level runner for isolated execution and actual check reports.

The archived model assertions concern the specified logical rules, inputs,
targets, eligibility and update semantics. No experimental biological
measurements, phenotype-control problem, 250000-network synthetic replay or
database-wide model scan is added by this package.
