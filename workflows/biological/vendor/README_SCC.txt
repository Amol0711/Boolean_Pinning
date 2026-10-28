SCC-localized permanent-pinning solver
====================================

Python 3.10 or later; standard library only. No network access, compiled
extension or external experiment files are required.

Quick start
-----------
  python scc_control.py example_network.json -o example_result.json
  python scc_control.py example_two_slot.json

The first example has one complex cyclic component {0,1}, one acyclic
singleton {2}, and exact additional optimum 1 (pin 0). The second is a
two-slot certified copy example with exact additional optimum 1 (pin 0).

Independent verification
------------------------
  python verification/validate_scc_solver.py --output results/solver_validation.json

The verifier supplies its own full-state dynamics, support and graph oracles;
it does not call production helpers to establish the reference answers.
See verification/README.txt for cases, limits and the saved deterministic
verification record. This bounded regression is separate from a replay of
the original 250000-draw synthetic experiment.

JSON input
----------
input_schema.json documents the structural input schema. The CLI validates
the additional cross-field requirements without a schema-library dependency.
Required fields:
  parents: one ordered list of zero-based source indices per coordinate.
  tables:  one Boolean truth table per coordinate; a k-slot table has 2**k
           entries. Input slot j contributes value * 2**j to the table index.
Optional fields:
  target:  one physical target bit per coordinate; default all zero.
  baseline: distinct coordinate indices fixed permanently to target values;
            default empty. All other coordinates are eligible with unit cost.
  options: support_mode ("slot" default, or "physical"),
           max_component_size (18 default), max_pins (8 default),
           max_work (2000000 default).

Caps must be nonnegative integers, or null for an explicitly unlimited cap.
Use the finite defaults for routine runs. Work caps apply separately to EACH
complex component. They are not global wall-clock or memory limits. The
component size cap controls successor-table memory, proportional to 2**s.
All input tables must exist in memory; general-arity preprocessing can require
O(sum_i k_i*2**k_i) work. No linear-time assertion is made for arbitrary-arity
truth-table input. Do not disable caps on untrusted or very large inputs.

The structural graph passes use worklist deletion and iterative SCC traversal.
This reference implementation additionally sorts dependency lists, components
and pin witnesses for stable, reviewable output. Those presentation-ordering
steps can add logarithmic overhead; the Python program is not claimed to be
strictly linear on every unrestricted-degree two-slot input. The bounded-degree
linear-time certificate uses adjacency lists and traversal order without
requiring sorted output. Exact component optimization remains exponential.

The constructor verifies that the physical target is an uncontrolled fixed
point. Zero-input rows are supported. Repeated parent slots read the SAME
physical coordinate, during both local and full dynamics. Targets other than
zero are conjugated exactly into error coordinates for support and search.

Python API
----------
  from scc_control import BooleanNetwork, solve_network
  net = BooleanNetwork(parents, tables, target=None)
  result = solve_network(net, baseline=(), support_mode="slot",
                         max_component_size=18, max_pins=8,
                         max_work=2000000)

  evaluate_state(net, state, pins=())
      Returns one PHYSICAL synchronous successor integer. State bit i is
      coordinate i. Pinned inputs and outputs are set to physical target bits.

  support_residual(net, baseline=(), mode="slot")
      Returns residual, mode, evaluations and deletion_count. The evaluations
      count worklist row-support queries; preprocessing is separate. No
      synchronous time-to-capture value is inferred from this worklist count.

  solve_component(net, nodes, max_component_size=18, max_pins=8,
                  max_work=2000000)
      Searches the isolated subsystem with all exterior errors fixed to zero.
      It is also usable on any declared subset as an isolated map. Component
      state indices in its evidence use the sorted global nodes order.

Original two-input simulation adapter
------------------------------------
  from scc_control import from_two_input_arrays, solve_network
  net = from_two_input_arrays(parents, rules_nonzero_corners, target=None)
  answer = solve_network(net)

This accepts the original Baseline engine's n-by-2 parents and n-by-3 error-corner
arrays (NumPy arrays or ordinary lists). Rule columns are g(1,0), g(0,1),
g(1,1), with g(0,0)=0. With a nonzero optional target the adapter converts
these ERROR rules back into physical rules before construction. It copies
inputs and never rewrites existing engine code or scientific records. The
current package performs no full replay of the 250000 original networks.

Support and graph conventions
-----------------------------
"slot" support lets free nominal input slots vary independently. It matches
the original two-input conservative support even when parents repeat.
"physical" support synchronizes repeated references before testing support.
Both compute a greatest support that captures physical dynamics. SCCs use
ALL nominal dependencies within that support. Extra syntactic dependencies
can enlarge components but do not invalidate localization.

G is true if each residual row has exactly one nominal slot from the residual.
G_prime is true if each CYCLIC residual SCC has exactly one internal nominal
slot per row and the conjugated singleton value of that slot is 1. Acyclic
singleton SCCs do not contribute. Parallel slots are counted separately, so
two identical self-parent slots do not qualify as a single-slot cycle.

For the unpinned original two-slot ensemble with slot support, these fields
are the original G and the weaker G_prime. In baseline/general-arity settings
they are analogous deterministic certificates; do not apply the original
ensemble's C_n identity or probability law automatically. The flag
original_two_slot_unpinned_scope records the syntactic scope, not whether a
network was actually sampled from that random ensemble.

Outputs and computational evidence
----------------------------------
  exact_cost: exact ADDITIONAL minimum, or null if unresolved.
  lower, upper: rigorous bounds including every component.
  feasible_pins: a verified-by-construction additional pin witness of size upper.
  total_feasible_pins: that witness together with baseline clamps.
  components: all SCCs, including acyclic singleton components.
  certificate_cycle_count: cyclic SCC count only when G_prime holds; otherwise
      null. The solver does not enumerate all marked cycles of a complex graph.
  total_work: sum of charged exact-search work, excluding support, SCC and
      symbolic-certificate preprocessing.
  unresolved_components: count of components with unresolved exact cost.

Every component record includes nodes, exact_cost, lower, upper, pins, reason,
work, cyclic, marked_slot_cycle and the complete internal nominal-slot list.
Complex-search records additionally include row_evaluations,
successor_traversals, completed_cardinality, next_cardinality,
candidate_counts, rejected_candidates, and the last rejected pin set plus an
initial error-state witness. This is deterministic enumeration evidence, not
a stored witness for every rejected subset. Independent verification should
replay the listed small component searches.

Each candidate's global convergence is tested on ALL compatible local states
using its controlled successor graph. Pin sets are visited in increasing
cardinality. No monotonicity of Boolean pinning success is assumed. Success
and failure of a candidate are always evaluated directly. All smaller
cardinalities are exhausted before a minimum witness is returned.

Work is one local row evaluation during successor construction or one
controlled successor-edge traversal during candidate testing. Construction
requires s*2**s row evaluations; if this exceeds the finite work cap, preflight
returns work_cap_successor_construction with zero charged work and no large
allocation. Otherwise the same cap covers construction and all candidates.

Stopping reasons distinguish component_size_cap,
work_cap_successor_construction, work_cap_candidate_search and pin_count_cap.
A stopped search retains its proved lower bound and the feasible all-pinned
upper witness. The aggregate never discards unresolved components. Cases
where completed failures force lower=upper are exact_from_bounds. Symbolic
acyclic and marked-cycle certificates need no exponential state search and
remain available even with zero search budget.

Scope
-----
Synchronous finite Boolean dynamics, permanent target-valued pins from time
zero, all compatible initial states, every nonbaseline coordinate eligible,
unit cost. This code does not infer biological intervention feasibility or
implement asynchronous semantics, phenotype targets, weighted costs or random
ensemble guarantees. External inputs fixed in an application should be passed
as baseline clamps and their cost reported separately.
