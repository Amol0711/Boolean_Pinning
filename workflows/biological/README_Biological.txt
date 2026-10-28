Exact fixed-target census for biological Boolean models
==================================

Run from this directory:
  python src/run_all_targets.py
or preserve a fresh reproduction separately:
  python src/run_all_targets.py --output reproduction_results

Python 3.10+ and NumPy are required. The computation is entirely local. The
SCC solver in vendor/ uses the standard library only. Production
outputs are deterministic; they contain no runtime timestamps or author data.
The code writes only to its designated output directory. It never edits the
frozen source models, original protocol, original certificates or SCC solver.

Outputs
-------
results/all_targets.json is the complete machine-readable record.
results/context_summary.csv summarizes every frozen condition, including the
condition with no fixed point. results/all_target_summary.csv contains one
row per fixed target with exact additional cost and raw marker values.

The JSON stores source hashes, protocol hash, repository commit, explicit
input assignments, model node orders, forced coordinates, all enumerated full
targets, original-zero-orbit selection flags, residuals, exact functional
dependencies, isolated SCC truth tables, all local minimum pin sets, global
feasible witnesses and restricted-eligibility results.

Source and protocol
-------------------
The study reuses the frozen eleven-context protocol and its
BNET files under larger_networks/. The immutable Biodivine repository commit
and BNET SHA256 values are checked against protocol/contexts.json. Inputs
missing update rows in the source BNET are explicitly completed by identity
rules and permanently clamped to the protocol assignments. Regulated rows
are compiled exactly as supplied using a restricted Boolean-expression AST.

This is an exhaustive fixed-target census under the specified conditions. The original context prespecification does not constitute
prespecification of the fixed-target census. No contexts are replaced,
no rule repairs are made, and the original orbit-selected targets remain
identified explicitly. In particular, CD4_IL2 means CGC=IL2RB=IL2_e=1 with
every other external input zero, exactly as frozen earlier.

Method
------
1. Propagate forced physical values from each complete external assignment.
   Every fixed point must obey these values. Enumerate every assignment of
   the remaining coordinates in deterministic chunks of 8192; test all
   remaining fixed-point equations. Every retained full target is checked
   against every compiled source row and every original parsed expression.
2. Independently repeat the declared orbit from zero regulated coordinates
   to identify the original selected target or cycling discovery outcome.
   Original target/residual/cost records are checked for exact agreement.
3. For each target, conjugate to zero deviations and compute exact physical
   support after input clamps. Restrict each row to this residual cube, then
   remove all functionally inessential parents by testing its complete truth
   table. These tests examine every context, not only singleton sensitivities.
4. Apply SCC localization to that exact functional graph. Its nominal
   slots now represent the exact surviving functional parents. The original
   syntactic-slot graph remains a valid overgraph but can merge distinct
   feedback components; it is not used for the mechanistic decomposition.
5. For every isolated SCC, retain its restricted error truth tables and
   enumerate all pin sets at the independently established minimum cost.
   A local subsystem can simplify further when its exterior errors are set
   to zero. Its original functional SCC label is retained, with any further
   local simplification visible in its exported isolated truth tables.
6. Combine local minimum witnesses. Explicitly inspect every compatible
   residual state to verify convergence to the target under this witness.
   The computed support reduction and transition checks establish convergence
   from every compatible initial state of the original full model.

No monotonicity of pinning success is assumed. SCC-local exact searches test
the actual synchronous transition map. The code does not infer phenotype
labels, asynchronous convergence, intervention feasibility in a laboratory,
or a random-ensemble Poisson law for these models.

Eligibility sensitivity
-----------------------
Every nonbaseline coordinate is eligible in the main full-state problem.
The exact source memory rule proliferation_next = STAT5_HIGH OR proliferation
is checked directly. A separate analysis prohibits pinning proliferation.
If that coordinate belongs to the target residual, every allowed pin subset
of its isolated SCC is checked, up to exhaustion when infeasible. This
correctly handles target proliferation=1 as well as target proliferation=0.
For a target with STAT5_HIGH=0, an unpinned isolated memory coordinate can
remain at the opposite target value, so target proliferation=1 is not an
automatic escape from the eligibility obstruction. If proliferation lies
outside the target residual, the original optimum and witness remain valid.

Raw named marker signatures are exported without assigning cell-lineage,
healthy/disease, or source-experiment phenotype labels.

Declared finite guards
----------------------
Source truth-table arity: at most 20 physical parents.
Fixed-point enumeration: at most 24 remaining free coordinates, chunk 8192.
SCC component optimization: size at most 18, pins at most 18, work at most
100000000 per complex component, with SCC's documented work accounting.
All-minimizer/restricted-eligibility enumeration: component size at most 12,
hence at most 4096 eligible subsets; each tested directly on all local states.
Whole-residual global-witness check: at most 20 residual coordinates.
Zero-orbit discovery: 2000 updates, matching the earlier update cap.

The frozen study fits every guard. If a changed input exceeds a guard or an
exact search remains unresolved, the command raises an error before claiming
a completed result; it does not silently omit that target. SCC component
records retain symbolic or exact-search reasons and lower-bound evidence.

Results of this bounded study
----------------------------
There are 46 CD4 fixed targets across six contexts, plus four macrophage
targets in its four fixed-point contexts. The macrophage BACTERIA condition
has zero fixed points after all 524288 remaining assignments are examined.
Across all eleven contexts, 805380 remaining assignments are examined.
The original ten fixed-target optima are preserved. When proliferation is
ineligible, seven CD4 targets remain feasible (all six ON targets and one IL2
target); the other 39 CD4 full-state target problems are infeasible.

Independent verification is maintained separately under verification/.
Its source parser, fixed-point enumeration, support, graph and optimization
oracle do not use production functions to generate the reference answers.
No full synthetic-ensemble replay or database-scale biological scan is part
of this package.
