Independent J1 replication simulation package
=================================

This package contains the fixed 20,000-draw replication for the J1 planted
two-input configuration model at n = 25,600. It also retains the original
2,000-draw J1 cell at that size for a separate comparison. No observations
from these two batches are pooled or replaced.

The protocol and seed list specify every scheduled draw. All observations are retained. Model probabilities remain fixed throughout the analysis. integrity_record.json records the current source, protocol and validation-file checksums.

MODEL AND FIXED PARAMETERS

Vertices 0..12,799 have output degree one; vertices 12,800..25,599 have output
degree three. A uniform permutation of all output stubs supplies two labeled
input slots per vertex. Self-loops and repeated physical parents are retained.
The target is zero, updates are synchronous, every vertex is eligible for a
permanent zero pin, and success requires convergence from every compatible
physical initial state.

Corner code c = B1 + 2*B2 + 4*B12, with B1 = g(1,0), B2 = g(0,1), and
B12 = g(1,1). Code order 0..7 is

  ZERO, x1 AND NOT x2, NOT x1 AND x2, XOR, AND, COPY1, COPY2, OR.

Degree-one probabilities are [30,0,0,0,15,45,10,0]/100.
Degree-three probabilities are [60,0,0,0,10,10,20,0]/100.
These laws are unchanged from the original J1 study. The fixed references are

  lambda = 29/80, rho = 19/40, mu = -log(51/80),
  E_reference[C^2] = mu + mu^2, P_reference(C=0) = 51/80.

They come from the declared probabilities, not the empirical rule counts.
The raw files retain all eight corner-code counts in each degree class. The
six displayed rule categories are ZERO, COPY1, COPY2, AND, OR, and XOR.

REPLICATION IDENTITIES AND SOLVER CAPS

The NumPy stream is Generator(PCG64(SeedSequence(
  [2026101601, 3101, 25600, replicate]
))). Replicate is exactly 0..19,999, independently of worker count or task
order. The original joint-rule study stream used master seed 2026101504 and family ID 1.
Protocol_J1_Replication.json gives the exact random-call sequence and all assumptions.

The component solver permits at most 18 component coordinates,
18 pins, and 10,000,000 deterministic work units per noncopy cyclic SCC.
Each unresolved draw would retain its certified cost interval, stop reason,
and feasible upper witness. None is discarded or imputed. The independent
oracle separately caps a component at 18 coordinates, 100,000 candidate pin
sets, and 100,000,000 compatible states; full-residual enumeration is limited
to dimension 12. Any oracle guard is explicitly recorded.

CONTENTS

  Protocol_J1_Replication.json       Fixed replication design and statistical conventions.
  seed_manifest.csv       All 20,000 ordered seed identities.
  integrity_record.json     Current protocol, source, validation and version records.

  original/               Selected baseline JointRules inputs.
    Protocol_JointRules.json       Original full JointRules protocol, for provenance.
    raw/J1_n25600.csv.gz   Original 2,000 observations used in the comparison.
    details/...           Original residual and certificate records.
    summary.json          Baseline joint-rule summaries, including its other cells.
    manifest.json         Baseline joint-rule run metadata.
  The original directory is an excerpt for this replication, not the complete
  18,000-draw JointRules archive. Its comparison manifest and summaries have that wider scope.

  src/run_replication.py      Fixed replication production entrypoint.
  src/joint_validation.py  Imported JointRules generator, analyzer, and
                           summary formulas. Its standalone CLI is
                           not the entrypoint for this package.
  vendor/                 Engine and physical SCC solver sources.

  results/main/raw/J1_n25600.csv.gz
                           One ordered record per new draw, including seed,
                           network hash, rule counts, support and SCC sizes,
                           G/G_prime, cost bounds, work, stops, and pins.
  results/main/details/J1_n25600.jsonl.gz
                           Complete nominal residual rows and SCC details;
                           physical repeated parents remain identifiable.
  results/main/summary.json and summary.csv
                           Production all-draw estimates and uncertainty.
  results/main/rule_counts.json
                           Degree-conditioned aggregate rule counts.
  results/main/statistical_summary.csv
                           Independently recomputed original/replication
                           comparisons, kept separate.
  results/main/rule_frequency_summary.csv
                           Independently recomputed rule diagnostics.
  results/main/io_integrity.json, manifest.json
                           Completed-record counts, stream/source hashes,
                           protocol linkage, and run status.
  results/main/resources.json, run_started.json
                           Execution timing and provenance. Timing and
                           execution timestamps need not reproduce exactly.

  verification/independent_network_oracle.py
                           Independent generator reconstruction, synchronous
                           support iteration, Tarjan SCCs, and physical
                           functional-graph indegree-peeling oracle.
  verification/preflight_replication_oracle.py
                           Exhaustive tiny-network and edge/guard fixtures.
  verification/validate_replication_simulations.py
                           Full independent network/certificate replay.
  verification/independent_replication_statistics.py
                           Independent statistics with exact integer sums,
                           rational parameter derivation, and rule diagnostics.
  verification/test_production_no_refit.py
                           Separate mutation check of the unchanged production
                           summarizer using deep copies of rule counts.

  qa/independent_preflight.json, production_preflight.json
                           Validation reports.
  qa/independent_network_replay.json
                           Completed full-network independent verification.
  qa/statistical_result.json, statistical_interpretation.txt
                           Independent numerical audit and interpretation.
  qa/production_no_refit.json
                           Production-summary invariance under altered
                           diagnostic counts, with original inputs preserved.
  qa/independent_validator_port_check.json
                           Port check on 100 baseline joint-rule J1 n=100 records.
                           This finite check is not new replication data.

DEPENDENCIES

The integrity record records Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
g++ 13.3.0, and x86_64 Linux with glibc 2.39. The production entrypoint
checks the exact NumPy and SciPy versions against the integrity record.

In a suitable Python 3.12 environment, install the numerical dependencies
with

  python -m pip install numpy==2.3.5 scipy==1.17.0

A C++17-capable g++ compiler reproduces the recorded accelerated backend.
vendor/engine.py builds engine_core.cpp locally on demand into vendor/_build.
The engine also has a Python fallback, which is slower. The independent
network replay itself does not use the production engine or C++ accelerator.
No external datasets, accounts, or network access are needed after installing
the dependencies. Run commands below from this package's root directory.

RECHECK THE ARCHIVED REPLICATION

  python verification/preflight_replication_oracle.py --output qa/rechecked_preflight.json
  python verification/validate_replication_simulations.py --workers 4 --preflight qa/rechecked_preflight.json --output qa/rechecked_network_replay.json
  python verification/independent_replication_statistics.py --output qa/rechecked_statistics.json
  python verification/test_production_no_refit.py --output qa/rechecked_production_no_refit.json

The first command checks all 1,024 ordered two-node zero-target networks,
all 4,096 physical pin sets, and 9,216 compatible initial states. It also
checks seven edge fixtures, four bounded-exit fixtures, and the empty residual.
Its deterministic report must match the distributed validation-report hash. The second command
regenerates and checks all 20,000 complete networks, not merely their summaries.

The independent statistics command recomputes its JSON report and regenerates
the two derived CSVs statistical_summary.csv and rule_frequency_summary.csv
under results/main. It does not modify raw observations, details, the fixed
protocol, or production summaries. The no-refit test alters only deep copies
in memory and verifies that the original inputs and files remain identical.

REPEAT THE PRODUCTION COMPUTATION

Use a fresh output directory. Do not delete or replace results/main, edit the
protocol or distributed sources.

  python src/run_replication.py run --workers 4 --output results/reproduced
  python verification/validate_replication_simulations.py --results results/reproduced --workers 4 --output qa/reproduced_network_replay.json
  python verification/test_production_no_refit.py --results results/reproduced --output qa/reproduced_production_no_refit.json

The production entrypoint refuses an existing output directory. It schedules
the same 20,000 identities, writes temporary shards, and combines them in
replicate order. Changing the worker count does not change any identity.
The fresh combined records and production summaries can be compared directly
with the archived ones, for example on a system with the cmp utility

  cmp results/main/raw/J1_n25600.csv.gz results/reproduced/raw/J1_n25600.csv.gz
  cmp results/main/details/J1_n25600.jsonl.gz results/reproduced/details/J1_n25600.jsonl.gz
  cmp results/main/summary.json results/reproduced/summary.json

Compression timestamps are fixed. Runtime, run-start timestamps, and resource
records are execution metadata and should not be required to match bytewise.
The statistics command targets the archived results/main by design. When the
fresh raw/detail/summary streams match, its statistical inputs are identical.

OBSERVED VALIDATION AND INTERPRETATION

The new sample contains 20,000 exact costs and zero unresolved intervals.
G holds in 19,999 draws; G_prime holds in all 20,000. Independent replay checked
512,000,000 source coordinates, 40,000 degree-conditioned rule-count vectors,
16,008 residual rows, and 13,611 SCCs. The SCCs comprise 9,040 copy cycles and
4,571 acyclic components. The sole G-failure has a three-coordinate residual;
independent full physical-state enumeration confirms its optimum of one.

There is no noncopy cyclic SCC in this new sample. Consequently the new data
do not measure performance on difficult noncopy SCCs. Their handling remains
covered by the preflight/edge checks and prior JointRules validation; do not describe
an absent case as an empirical stress test. Full networks include 14,948
repeated-parent rows; none remains in the new nominal residuals.

Moment intervals are pointwise approximate normal Monte Carlo intervals for
finite-size quantities. Probability intervals use Clopper-Pearson bounds.
The Poisson quantities are fixed limiting references, not an exact finite-n
null distribution or a quantitative finite-size error bound. Because the J1
replication was selected after inspecting the original batch, between-batch
contrasts are exploratory. Rule-frequency diagnostics do not estimate a new
lambda, rho, or Poisson mean. Sampling intervals and rule-frequency diagnostics retain their distinct statistical meanings.

CLEAN DISTRIBUTION

The source distribution omits generated __pycache__ directories, .pyc files,
native .so binaries, vendor/_build, transient .crdownload files, and redundant
results/main/shards. The combined streams, manifests, and all source files
needed for fresh generation and independent replay are retained. Resource
metadata may still list recorded shard hashes; those records document the
original run and do not make the omitted temporary shards an input dependency.
Fresh reruns recreate their own shards and local accelerator build.
