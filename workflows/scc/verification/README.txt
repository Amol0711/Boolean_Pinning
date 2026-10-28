Independent solver validation
=============================
From workflows/scc/ run:

  python verification/validate_scc_solver.py --output results/solver_validation.json

Python standard library only; no external network, model files, external implementation,
or imports outside this package are used. The command regenerates fixtures.json
and a deterministic scientific validation record. Wall time and environment
metadata are deliberately excluded from that record.

independent_oracle.py has independent implementations of physical state updates,
full compatible-state reverse reachability, every pin-subset enumeration,
nominal-slot and physical-variable support elimination, SCCs using mutual
reachability, and vertex-simple directed marked-slot cycle enumeration.
Parallel nominal slots count separately. The production evaluator is itself
cross-checked against this separate physical map.

Completed coverage
------------------
* All 1024 zero-target two-node/two-slot cases: 16 ordered wiring choices and
  64 planted truth-table pairs. All four pin sets are checked in each case.
* 300 seeded arbitrary-rule cases: 50 each at n=3,4,5,6,7,8, with 0–4 input
  slots, arbitrary planted mixed targets, repeated physical parents, and
  randomly selected baseline clamps. Every additional pin set is examined
  by the full-state oracle; the solver is compared from each resulting
  baseline in BOTH support modes.
* Eight explicit boundary fixtures, including repeated self-XOR cancellation,
  repeated self-OR parallel cycles, AND synergy without marked cycles,
  downstream nonlinear rows with G false/G' true, ignored self-dependencies,
  nonzero target/baseline clamps, and a row with eight nominal slots.
* Four incomplete-budget records: component size cap, pin-count cap, work cap
  before successor construction, and work cap during candidate search.
  Every returned interval contains the independent optimum and every upper
  bound has a verified feasible full-network witness.

Totals (the boundary fixtures can overlap the exhaustive family)
-------------------------------------------------------------
1332 network instances, 17916 oracle pin-subset checks, 218580 compatible-state
successor checks, 29342 raw physical states, and 31740 solver calls, including
four incomplete-budget calls. All assertions passed. These are finite tests,
not exhaustive verification for arbitrary network sizes.

Certificate scope
-----------------
For nominal-slot support, every marked cycle of the original baseline-clamped
slot graph survives reduction. On G' the solver cost and certificate count are
checked against the independently counted original marked-slot cycles.
Physical-variable support can delete a repeated-slot cancellation cycle; its
G' certificate therefore counts marked cycles of that reduced residual graph,
not necessarily the original nominal-slot graph. The validator preserves this
distinction. Success under additional pins is never assumed monotone.

Independent targeted cap/representation check
--------------------------------------------
  python verification/targeted_cap_check.py --output results/targeted_cap_check.json

This preserves the separate 720-call check. It tests five
small fixtures (repeated XOR, repeated OR, a complex OR component, nonmonotone
pinning success, and a mixed-target baseline), both support modes, and the
Cartesian product of component-size, pin-count and work caps listed in the
JSON output. Its oracle directly follows original physical truth-table
trajectories from every compatible initial state and does not call production
evaluation/support/SCC/convergence helpers. Every reported interval, exact
answer and full feasible witness is checked; finite per-component work caps
are checked as well. The record includes stopping-reason counts and cases
stopped after a completed rejected cardinality. Output excludes timings and
is deterministic. These targeted finite checks do not replay the synthetic
network experiment or extend its scientific results.

Independent reduction and eligibility checks
--------------------------------------------
  python verification/validate_reduction.py

This separate implementation checks state-space reduction, component-local
pinning costs, eligibility restrictions and nominal-slot cycle certificates.
It uses exhaustive two-variable cases and seeded three-/four-variable cases,
without importing production solver or graph routines. Its deterministic output
is results/reduction_validation.json.
