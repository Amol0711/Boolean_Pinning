# Focused validation of general joint two-input rule laws

This is a separate 18,000-network validation sample for the joint-rule model. It does not regenerate or modify the Bernoulli sample. The three laws, seed
namespace, nine cells, deterministic component caps, and statistical summaries
were fixed in `Protocol_JointRules.json` before the main run. The independent exact gate
is `../qa/independent_exact_checks.json`.

Each local rule is sampled jointly as an integer code
`c = B1 + 2*B2 + 4*B12`, where `B1=g(1,0)`, `B2=g(0,1)`, `B12=g(1,1)`, and
`g(0,0)=0`. Codes 0, 3, 4, 5, 6, and 7 are respectively ZERO, XOR, AND, COPY1,
COPY2, and OR. The eight code probabilities are exact integer weights divided
by 100. Individual truth-table corners are not sampled independently.

| Family | Rule-law feature | lambda | rho | max(lambda,rho) |
|---|---|---:|---:|---:|
| J1 | Asymmetric copy/AND laws associated with degrees 1 and 3 | 29/80 | 19/40 | 19/40 |
| J2 | Degree 2 with correlated OR/XOR/AND law | 1/2 | 9/20 | 1/2 |
| J3 | Degree-3 type has zero singleton probabilities but positive AND mass | 1/5 | 5/8 | 5/8 |

Here lambda is the degree-weighted sum of singleton marginals and rho is the
degree-weighted probability that a row is nonzero. Each family has mean output degree exactly two. All type
probabilities, fractions, and derived parameters are retained as rational
numbers in the protocol. The benchmark Poisson mean is fixed at
`mu=-log(1-lambda)`.

## Commands and dependencies

Use Python 3.11 or newer with the versions in `requirements.txt` (the recorded
run uses Python 3.12.14). A C++17 compiler enables the shared
certificate accelerator; the original engine also supports a Python fallback.
No online data are required. From this directory:

```sh
python src/joint_validation.py exact --output results/replayed_exact
python src/joint_validation.py main --workers 4 --gate ../qa/independent_exact_checks.json --output results/replayed_main
```

To repeat both calculations in a fresh location and compare all 23 scientific
files (including all 18 compressed main streams) with the saved results:

```sh
python reproduce.py --workers 4 --out results/verified_replay
```

Omitting `--out` uses a retained temporary directory. `--exact-only` performs
only the three-file exhaustive comparison. Compressed bytes and decoded
payloads are compared; timing and manifest metadata are not compared. The
default gate location assumes the archive retains the sibling `qa/` directory;
`--gate` can explicitly select the same preserved independent gate file.

Output directories must not already exist. The gate must report `PASS`. Main
generation uses PCG64 with `SeedSequence([2026101504,family_id,n,replicate])`.
It first permutes the output source-label multiset, then calls `rng.choice`
once per type block in declared order to sample the joint corner code. Sizes
are 100, 1600, and 25600, with replicates 0 through 1999 in every cell.
Self-loops and repeated parents are retained.

The three vendor files implement the shared generator/certificate
engine and the verified SCC solver. New random generation is implemented in
`src/joint_validation.py`; the old Bernoulli generator is not called.

## Exact checks

`results/exact/` exhaustively enumerates all distinct matching source sequences
and positive-probability rule assignments for six two-node laws. Source
sequences have equal labeled-stub multiplicity and therefore equal probability.
All 228 realized networks have exact rational weights, marked loop and two-cycle
counts, residual certificates, physical optimum costs, and pin witnesses.
The six laws comprise J1/J2/J3, a ZERO/AND mixture with lambda zero, constant
ZERO, and constant AND. All eight deterministic corner codes also have their
nominal active-slot support patterns recorded.

Exact identities cover expected marked loops, marked two-cycles, and the
second factorial moment of marked loops. The latter includes
`P(B1=B2=1)` for overlapping slots within a row, rather than a product of
singleton marginals. Independent verification checks the exported physical
networks and costs without importing production routines.

Constant AND has lambda zero and rho one, so it lies outside strict support
contraction. Its all-one state is a nonzero fixed point and a positive control
cost is required although all singleton-marked cycle counts vanish. The
ZERO/AND mixture has lambda zero and rho below one; its tiny finite networks
can still have positive control costs. 

## Main computation and records

The conservative nominal-slot residual and original certificate `G` are
computed by the original deterministic engine. Every residual row must have
exactly one residual nominal slot for `G`. The weaker certificate `G_prime`
allows arbitrary acyclic singleton rows and requires exactly one internal
marked slot per node of every cyclic SCC. Repeated physical-parent slots
remain distinct for these certificates.

For each non-`G` draw the SCC solver receives both nominal slots of
every residual row; exterior slots map to a constant-zero sentinel. Isolated
component dynamics identify repeated slots with one physical coordinate.
The solver sums component costs, with caps of 18 nodes, 18 pins, and 10,000,000
charged work units per complex SCC. A work unit is a local-row evaluation or
controlled successor-edge traversal. Every candidate pin set is tested
directly. All intervals and feasible upper witnesses are retained on a stop.

`results/main/raw/FAMILY_nSIZE.csv.gz` has one row per draw. Its columns record
family/size/replicate, full-network SHA256, residual size, two-slot row count,
G/G-prime flags, cycle count when certified, cyclic-component count, largest
complex SCC, cost lower/upper bounds, exact/stop reason, charged solver work,
and semicolon-separated original pin IDs. The network digest hashes parent
indices as little-endian signed int64 row-major bytes followed by corner
triples as uint8 row-major bytes.

The matching `details/*.jsonl.gz` files retain residual rows as
`[original_vertex,parent_slot0,parent_slot1,B1,B2,B12]`. On non-G draws they
also retain every SCC component record and `local_to_original`, with a final
null entry for the constant-zero sentinel. Component nodes and pins use this
local indexing; failure-state bits use the component node order. These fields
determine each exterior-target-restricted physical subsystem. On G draws,
components are omitted, while full residual rows and the exact cycle count
remain available.

Streams are written to excluded temporary `*.crdownload` paths, closed and
verified through their gzip trailers, atomically published, and checked again
by count/record IDs/hash before a completed manifest is emitted. Gzip headers
and scientific records are deterministic across worker counts; runtime is
stored separately in `resources.json`. `io_integrity.json` records all 18
scientific streams and their hashes.

## Complete-sample summaries

`summary.json` and `summary.csv` contain all nine cells. Mean and raw second
moment lower/upper endpoints average all 2,000 draw bounds. Endpoint SEs are
computed from the unbiased sample variance of the endpoint variables divided
by sample size. The pointwise normal envelope uses
`z=Phi^(-1)(0.975)=1.959963984540054`, with its lower end truncated at zero.
These approximate sampling envelopes are distinct from computational intervals.

P0 uses definite-zero and possible-zero counts. Its conservative pointwise 95%
Clopper–Pearson envelope takes the lower limit from the definite-zero count and
the upper limit from the possible-zero count. G and G-prime failure proportions
have pointwise two-sided 95% Clopper–Pearson intervals, including zero counts.
The reference mean, raw second moment, and P0 are respectively mu,
mu+mu squared, and 1-lambda. The reference values are fixed by the model probabilities. Sampling intervals describe finite-size quantities.
