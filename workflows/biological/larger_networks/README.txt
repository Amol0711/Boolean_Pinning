LARGER BOOLEAN-NETWORK BENCHMARKS

This package contains frozen source networks, eleven predeclared input contexts,
an exact support/residual runner, complete scientific certificates and search
ledgers, and isolated-process resource measurements. Python 3.10 or later and
the standard library are sufficient.

REPRODUCTION

From this directory run

    python src/test_runner.py
    python src/benchmark.py --protocol protocol/contexts.json --output /absolute/path/to/a/new/output_directory

The output directory must not exist. The runner checks the exact frozen protocol
SHA256 before execution and verifies each source BNET hash against that protocol.
It executes all eleven contexts in their declared order, with one fresh worker
process per context. qa/environment.json records the execution environment. No outcome is omitted. A rerun may take approximately one
minute on an environment similar to the reference run; this is not a performance
guarantee. Resource measurements are intentionally separate from scientific
certificate files and vary between reruns.

FROZEN DESIGN AND SOURCE SEMANTICS

protocol/contexts.json is the authority for context assignments, order, and caps.
Its current SHA256 is

    2d8f64156c5643cd0559a149b5ccde1cc7bed26d06bf0f517963a7db188aa47d

The original protocol commitment is retained separately in
protocol/freeze_record.json. Context assignments and budgets are unchanged.

The sources are BBM id050 CD4 T-cell signaling (154 regulated coordinates plus
34 free inputs, 188 total) and id001 signaling in macrophage activation (302
regulated coordinates plus 19 free inputs, 321 total), pinned to BBM commit
314b2fed94689c38d3aaaf8101ed622b3501c7f7. Source BNET files and metadata are in data/.
The complete node list comes from BNET update left-hand sides and referenced
variables. The metadata variable-names field is not a complete list in every
model and is not used to count nodes.

The BNET files omit updates for free inputs. The runner explicitly completes
each omitted input rule as an identity update. Dynamics are synchronous. Input
contexts fix all input coordinates from time zero. Additional controls are
permanent pins at the selected target values. Every regulated coordinate is
eligible, including output and phenotype coordinates such as proliferation.
This is not a gene-only intervention model or a clinical treatment claim.

Within every certificate, node_order is sorted regulated BNET names followed by
sorted input names. All numerical indices are zero based in this order. The
original names are retained. Each local truth table uses distinct physical
variables, so repeated references are evaluated consistently.

Four endpoint contexts were known from the preceding endpoint computation. Seven
new focal-stimulus contexts were selected from source coverage and explicit
input rules before examining their new outcomes. The frozen protocol records
that distinction. Focal-stimulus assignments are computational test conditions,
not assertions that a Boolean context or discovered full-state target is a
validated physiological state. Figure contexts were also selected before the
seven new outcomes were computed.

TARGET DISCOVERY AND VERIFICATION

For each context, the sole target-discovery procedure starts all regulated
coordinates at zero, fixes every input to its declared value, and follows the
synchronous map for at most 2,000 updates. A fixed candidate is verified both
against the physical truth tables and directly against the source-expression
AST. Discovery of a target is not itself a convergence certificate.

If that orbit enters a nontrivial cycle, the result records the entire cycle and
the initial state. It reports that this declared orbit did not find a fixed
target. It does not assert that the model has no other fixed point, and no
alternative target-discovery method or favorable replacement context is used.

For a verified target, exact local-cube support reduction provides an invariant
residual. A nonempty residual does not alone imply unsuccessful control. When
the residual contains at most 18 nodes, the runner constructs its exact physical
successor table and tests pin subsets in increasing cardinality. Each candidate
is classified over every compatible residual state. General physical success
is not assumed to be monotone under additional pins.

The frozen per-context caps are 30,000,000 work units, 120 seconds of analysis
wall time, and candidate cardinality at most six. One work unit is reserved per
residual successor evaluated and per compatible state in a candidate test.
Interrupted work is conservatively included in the reserved count. Completed
unsuccessful cardinality levels certify lower bounds. Pinning the complete
original residual always supplies a feasible upper bound if a cap is reached.
A pure-copy certificate
can give an exact cost directly if a residual exceeds the enumeration size cap.

The baseline inputs are identity coordinates and are individually mandatory if
all initial states are allowed before imposing the desired context. Therefore
total cost, when those clamps are charged, is baseline-input count plus the
additional optimum or bounds. If the environment is considered fixed for free,
only the additional cost is counted. Both conventions are reported explicitly.

SCIENTIFIC OUTPUTS

results/<context_id>/certificate.json contains the full fixed target or discovery
failure evidence, support trace, original-to-residual mapping, residual successor
table when computed, every completed candidate test, failing cycles, successful
transient witnesses, lower/upper cost bounds, chosen feasible intervention, and
the completeness status of any minimum-set list.

Ten of the eleven declared contexts produced a verified fixed target. All six
CD4 contexts have exact additional costs: three in OFF, ON, ECM, IL2, and GALPHAQ;
four in APC. All four target-producing macrophage contexts have exact additional
cost zero. MAC_BACTERIA follows a period-four orbit after an eleven-update
transient from the declared starting state; no fixed-target control result is
claimed for that context.

CD4_APC exhausts all candidate cardinalities zero through three and finds two
successful cardinality-four candidates before the work cap. Thus its optimum
four is exact, while its complete collection of optimum sets remains unknown.
The complete candidate ledger and the selected feasible set are retained. In
summary.csv, minimum_sets_count is blank for an incomplete collection, and
successful_candidates_found reports verified successes actually found. The
selected set for a cap-limited run is the last successful candidate encountered
before the cap; fully completed optimum levels select the first successful set.
This selection is deterministic and does not modify the fixed context or budget.

Uniform convergence bounds for the full network are proved upper bounds, not
claims of attained worst-case transient. When available, the runner reports the
minimum of two valid bounds: empty full support after the selected pins, and
baseline-support capture time plus exact residual hitting time. No 2^188 or
2^321 full state space is enumerated.

RESOURCE MEASUREMENTS

results/<context_id>/resources.json records one isolated worker execution.
Analysis wall/CPU times include parsing, truth-table compilation, target
discovery, support reduction, exact search, and certificate-object construction.
They exclude interpreter startup, imports, and file serialization. A separate
parent-measured process wall time includes startup and export. Peak worker RSS
includes the worker's entire process lifetime and interpreter/import overhead;
it is not an incremental algorithm-memory measurement. These single execution
measurements are a reproducibility diagnostic, not a scaling experiment or an
asymptotic complexity claim.

VERIFICATION AND HASHES

src/test_runner.py independently checks the finite-state classifier on all 64
two-coordinate zero-fixed successor maps and all four pin sets, checks the
nonmonotone-success example, tests all search cap statuses, verifies a pure-copy
fallback, checks periodic target-discovery status, and rejects unsupported or
unsafe formula syntax. The general local-cube engine is copied from the
previously validated biological certificate package.

An independent parser and certificate-ledger auditor are in verification/. They are
maintained separately from the runner. See the verification report for the exact
checks completed on this workflow. results/run_manifest.json records the frozen
protocol and source hashes and scientific result hashes. Resources and summary
timing fields are expected to vary when reproduced.

This package contains computational sources, model data, and verification
outputs.
