Independent biological verification
======================================
From workflows/biological run:

  python verification/verify_all_targets.py

Optional explicit paths:
  python verification/verify_all_targets.py \
    --production results/all_targets.json \
    --output results/independent_biological_validation.json

This verifier requires only Python's standard library. It imports the packaged
independent_core.py parser/oracle; it imports neither the production solver nor
production BNET/parser routines. All referenced model, protocol and prior
certificate files are inside larger_networks/ in this package.

Independent methods
-------------------
1. Parse original BNET rules independently and build complete local truth tables.
2. Propagate forced values by matching every compatible row-table corner.
3. Evaluate every regulated fixed-point equation over every remaining assignment
   using arbitrary-precision Python integer bitsets. These bitsets encode the
   exhaustive Boolean assignments; they are not sampled states or numerical
   approximations. Every recovered full target is rechecked by scalar evaluation.
4. Recompute exact physical residuals. Enumerate each restricted row table and
   remove exactly the parents whose bit flips never alter that row. Compute SCCs
   by independent mutual reachability, preserving self-dependencies.
5. For each isolated component, enumerate all pin sets at every cardinality up
   through the minimum, checking every compatible local state by reverse
   reachability. Enumerate all minimum pin sets at that cardinality.
6. Verify each production global witness on every compatible residual state and
   verify invariance of all outside-residual rows. Check proliferation-ineligible
   component feasibility and every available restricted global witness.
7. Compare full targets, frozen contexts, rows, SCCs, local minima and costs to
   the production JSON. Verify original selected-target results against the
   unchanged prior certificates. No success monotonicity assumption is used.

Completed deterministic coverage
--------------------------------
11 contexts; 50 fixed targets (46 CD4 and 4 macrophage); 805380 fixed-point
candidate assignments; 478 component comparisons; 1155 independently tested
local pin candidates over 292247 compatible local states; 40 additional
restricted-eligibility local pin candidates over 80 compatible states.

All 50 production global witnesses were verified over 281092 compatible
residual states. The seven feasible proliferation-ineligible CD4 targets were
additionally verified over 49152 compatible states. The other 39 CD4 target
problems are infeasible under that eligibility restriction. This count concerns
specified full-state targets, not biological feasibility.

CD4 fixed-point counts by OFF, ON, APC, ECM, IL2, GALPHAQ are 6, 6, 15, 6, 7, 6.
The actual IL2 context has CGC=IL2RB=IL2_e=1 and all other inputs zero.
Macrophage BACTERIA retains 19 unforced coordinates; every one of its 524288
assignments violates at least one fixed-point equation. Its frozen input
condition therefore has no fixed point under the stated original rule map.

Every scientific assertion passed. The output contains deterministic expected
records, all source/production hashes and exact verification counts. Timings,
absolute host paths and environment metadata are not included.

Scope
-----
These checks establish properties of the archived rules, target values,
permanent pins and synchronous semantics. They do not implement phenotype
control, provide biological lineage labels, validate the random-ensemble
Poisson law, or rerun the 250000-network synthetic study.
