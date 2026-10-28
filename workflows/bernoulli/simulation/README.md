# Bernoulli replay of the original configuration-model instances

This directory refines the control-cost calculation for the **same 250,000
networks** in the initial sample: five degree/rule mixtures, five sizes, and
10,000 replicates per cell. Every draw is regenerated from its original seed.
No new statistical sample is generated, and no original exact cost is replaced
by a different value. Bernoulli is a post hoc computational refinement; its SCC search
caps were selected after the initial sampling and frozen before the Bernoulli main
replay. See `Protocol_Bernoulli.json` for the complete computational protocol and the
separate `../verification/Bernoulli_Analysis_Protocol.json` for statistical analysis.

## Contents

- `baseline/` contains the initial protocol, initial and refined records, code,
  and figures.
- `vendor/engine.py` and `vendor/engine_core.cpp` are generator
  and certificate sources. The C++ file accelerates deterministic calculations;
  NumPy alone generates the random networks.
- `vendor/scc_control.py` is the independently verified SCC solver.
- `src/replay_scc.py` joins the original seeded generation with SCC localization.
- `results/main/` contains the complete Bernoulli main replay; `results/pilot_10m/` contains the first 100 original replicates per cell.
  References to another pilot identify a separate calculation, not extra data.
- `qa/environment.json` records the numerical execution environment.

The solver reads only files below this directory. Verification and statistical
programs live in the sibling `verification/` directory of the full archive.
Temporary C++ compilation occurs under `vendor/_build`, never under `baseline/`.

## Reproduce

Python 3.11 or newer and NumPy are required for the pinned requirements. The recorded execution used Python
3.12.14 and NumPy 2.3.5. A C++17 compiler enables the original deterministic
accelerator; the engine has a Python fallback if compilation fails.
No network access is needed.

The pinned runtime dependency is in `requirements.txt`; the recorded compiler
was g++ 13.3.0. Statistical analysis and figure generation have their own
dependencies documented in the full archive.

From this directory, run:

```sh
python src/replay_scc.py --workers 4 --output results/replayed_main
```

The program refuses to overwrite an existing output directory. The output path
can be anywhere, while all source inputs remain relative to this directory.
Worker count and task completion order do not change the scientific records.
The default work cap is 10,000,000 per complex SCC; the non-pilot command checks
that its caps agree with the frozen protocol. No wall-time stopping rule is used.

For a complete replay **and comparison** with all 50 saved scientific streams:

```sh
python reproduce.py --workers 4 --out results/verified_replay
```

Omitting `--out` creates a retained temporary output directory. The script
validates all reference and regenerated gzip trailers, record counts, and
consecutive IDs, then compares compressed bytes, decoded payloads, and both
operational summaries. It writes `reproduction_comparison.json` in the new
output directory and exits unsuccessfully on any difference. Runtime and
backend metadata are excluded from the bytewise comparison. Add `--pilot` for
the first 100 original replicates per cell (2,500 total) and comparison with
the saved 10M pilot. A full replay took about four minutes with four workers
on the recorded machine; exact time depends on hardware.

The declared timing pilot can be repeated with:

```sh
python src/replay_scc.py --pilot --workers 4 --work-cap 10000000 --output results/replayed_pilot
```

The pilot reuses the first 100 original replicates in each cell. All these
instances also remain in the main replay. The 10M pilot took 2.39 seconds with
four workers on the recorded machine and passed independent full seeded
reconstruction and component-optimum checks before the main replay began.

## Mathematical conventions

Vertices retain their original zero-based indices. The planted target is zero,
all vertices are eligible for permanent zero pins, and cost is the number of
pins needed for eventual convergence from every compatible initial state under
synchronous update.

The residual is the original conservative **nominal-slot** support residual.
Two slots with the same source remain two slots. The original certificate `G`
requires exactly one residual nominal slot in every residual row. The weaker
`G_prime` certificate requires every cyclic SCC to have exactly one internal
slot at each node, with singleton error mark one. Acyclic singleton SCCs need
no such condition. On `G_prime`, the number of cyclic SCCs is the original
marked-slot cycle count and the exact cost. The old certificate implies the
weaker certificate; failure of either certificate is not failure of control.

For every non-`G` draw, the full SCC solver receives both original slots of each
residual row. An exterior source maps to an added constant-zero sentinel. This
preserves the nominal graph and isolated physical maps. In a component's exact
dynamics, repeated slots use the same physical bit and every exterior error is
held at zero. Exact component costs add by SCC localization. These are isolated
component calculations, not an assumption that all components start at their
exterior target values in the original dynamics.

Every complex SCC is searched with caps of 18 vertices, 18 pins, and 10,000,000
charged work units. One unit is one local-row evaluation in successor-table
construction or one controlled successor-edge traversal during candidate
testing. Preprocessing, allocations, and bookkeeping are not charged. The cap
is per component, not per draw. An 18-node successor table requires 4,718,592
row evaluations. Pin sets are tested in increasing cardinality, and each set is
tested directly; success is not assumed monotone under adding pins.

On a stop, the solver retains a lower bound from completed rejected
cardinalities and a feasible upper witness (all component vertices if no
smaller witness was established). Search status is recorded per component.
The final draw interval is the intersection of the SCC interval and the
**final refined Baseline interval**. A returned global pin set witnesses the selected
upper bound. Every original exact value must be preserved by assertion. No
unresolved draw is discarded or assigned an invented point cost.

The graph algorithm underlying the certificate is linear in the input graph.
The portable SCC implementation sorts some output/index lists for reproducible
ordering, so this implementation does not claim a strict linear running-time
bound. Exact component enumeration is exponential and capped as stated above.

## Scientific output schema

Each cell has `raw/A_n100.csv.gz` and `details/A_n100.jsonl.gz`, with analogous
names for all 25 cells. There are exactly 10,000 records per file in increasing
replicate order. Gzip timestamps are zero and filenames are omitted from gzip
headers. Bytewise compressed reproducibility assumes the same compression
library; the uncompressed scientific content is deterministic across worker
counts. Timing is recorded separately in `resources.json`.

Streams are written to temporary `*.crdownload` paths, closed, fully decoded
through their gzip CRC/trailers, and checked for the expected consecutive record
IDs before atomic publication. All 50 final streams are checked again against
their recorded hashes before a completed manifest is published.
`io_integrity.json` records that final gate. This prevents a partially written
stream from being exposed as a completed scientific file in a synchronized
workspace.

CSV columns:

| Columns | Meaning |
|---|---|
| `mixture,n,replicate` | Original sample identity; mixture IDs A–E are 1–5. |
| `network_sha256` | Digest of the full regenerated parent/rule arrays; encoding is specified in the protocol. |
| `residual_size,double_rows` | Original nominal-slot residual size and rows with two residual slots. |
| `G,G_prime` | Boolean certificate indicators encoded as 0 or 1. |
| `certificate_cycle_count` | Exact marked-slot cycle count on `G_prime`; empty otherwise. |
| `cyclic_components,largest_complex` | Number of cyclic SCCs, and largest cyclic SCC that fails the marked-slot-cycle condition (zero when absent). |
| `scc_lower,scc_upper` | New SCC/certificate cost interval before intersection. |
| `original_lower,original_upper` | Final refined Baseline cost interval. |
| `cost_lower,cost_upper` | Final intersection interval, always present. Equality means exact. |
| `exact_reason` | One of the six classifications below. |
| `stop_summary` | JSON counts of unresolved component search reasons; these can remain even if the final intersection is exact. |
| `solver_work` | Sum of charged work across complex components. |
| `pins` | Semicolon-separated original vertex IDs witnessing `cost_upper`; empty when cost zero. |
| `bound_witness_source` | `baseline_feasible_witness` or `SCC_or_G_feasible_witness`. |

The six `exact_reason` values are `original_G`, `weaker_G_prime`, `SCC_exact`,
`preserved_original_exact`, `intersected_bounds_exact`, and
`bounded_unresolved`. Certificate status and computational resolution are
separate quantities.

Each JSON detail record retains the identity, digest, both certificate flags,
final bounds and witness, complete original refined CSV record, and:

- `residual_rows`: `[original_vertex, parent_slot_0, parent_slot_1, B1, B2, B12]`
  for every residual vertex. The zero corner is identically zero. All original
  parent IDs are retained, including exterior sources and repeated slots.
- `local_to_original`: residual vertices in increasing original order, followed
  by `null` for the constant-zero sentinel. It is `null` on the old-`G` fast path.
- `components`: full SCC records on every non-`G` draw. Component `nodes`, `pins`,
  slot endpoints, and failed-candidate indices are in the declared local
  residual index system and are translated with `local_to_original`.
  The component node order also defines the bits of recorded failure states.
  Exterior-target-restricted physical truth tables are uniquely recoverable
  from `residual_rows`, component membership, and the zero exterior convention.
  Components are omitted (`null`) on the original-`G` fast path, where the
  complete residual rows and certificate cycle count remain available.

Complex-component records include exact cost or bounds, a local feasible pin
set, caps, work split into row evaluations and successor traversals, rejected
candidate counts by cardinality, last rejected candidate/failure state, and a
precise stop reason. Symbolic acyclic/marked-cycle components have zero charged
search work. A whole-network upper witness can come from the original record
even when some new component searches remain bounded.

`operational_summary.json` and `.csv` aggregate certificate failures, exact
resolution, maximum sizes, work, and stop reasons without excluding any draw.
`run_manifest.json` identifies source/protocol/raw-input hashes and confirms
completion. Read the main directory only once that manifest says `completed`.
Statistical intervals and figures are generated separately from these complete
records. Finite-size Poisson diagnostics are comparisons with the asymptotic
benchmark, not exact finite-sample null statements or convergence-rate claims.
