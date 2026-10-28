INDEPENDENT Bernoulli NETWORK VERIFICATION

From the workflow root, run:

  python verification/validate_bernoulli_networks.py --workers 4

The default input is simulation/results/main. The default output is
results/independent_bernoulli_validation.json. Use --results and --output to validate a
separate fresh replay without changing delivered scientific records. NumPy is
required for the frozen PCG64 stream and independently evaluated state tables;
no production generator wrapper, support routine, SCC code or solver is imported.

The separate tiny reference-oracle check can be rerun with:

  python verification/independent_network_oracle.py

That check enumerates all 1,024 ordered two-node, two-slot, zero-target networks,
all 4,096 physical pin sets, and 9,216 compatible initial states. A literal
state-trajectory evaluator is compared against the indegree-peeling oracle.
The full validator also runs this check automatically.

After the full audit, run its one targeted completion:

  python verification/targeted_component_check.py

The full audit's 20-million-compatible-state guard stopped during cardinality
four for the 17-coordinate component of B, n=100, replicate=3449. This separate
command raises only that component's state guard to 100 million. It independently
finds local optimum four after 1,530 candidates and 22,544,384 compatible states;
the complete draw has optimum five. Its record is
results/independent_bernoulli_targeted.json. Combined, the full and targeted reports
verify all 249,208 exact draw costs and all 13,697 exact complex-component costs.
The 793 remaining component guards all exceed 18 coordinates, occur in 792
retained interval-valued draws, and are not asserted exact. All 250,000 final
feasible witnesses are verified in the full report.

Coverage and interpretation
---------------------------
Each full-run draw is reconstructed directly from Protocol_Baseline.json with the
specified PCG64 stream and call order. The complete generated arrays are hashed;
literal synchronous support iteration produces the independently checked residual.
Tarjan's algorithm reconstructs SCCs on the nominal slot graph. Both occurrences
of a repeated physical parent remain distinct for G and G-prime.

Acyclic singleton and marked-cycle components use exact symbolic certificates.
Complex components use separately constructed physical successors, nonzero
fixed-point rejection and indegree peeling, with all exterior coordinates zero.
The finite reference search guards are recorded in the output. Components beyond
these guards are not described as exhaustively optimized. Their valid bound
certificate and feasible witnesses are still checked. Every witness must pass;
an unverified witness causes validation failure rather than being silently omitted.
All non-G residuals of size at most 12 also receive direct whole-residual optimum
enumeration, independently checking the component cost combination.

The former double-row deletion/copy-cycle interval is reconstructed for every
originally unresolved draw. Every final bound must equal the intersection of
that earlier interval and the SCC interval. Every original exact optimum is
preserved. The 25 original and 25 refined raw-file hashes are checked against
their frozen manifests before network validation.

Exact coverage counters, any independent optimum guards, source hashes and all
raw/detail input hashes are recorded. Machine timings and scheduling are excluded
from the deterministic scientific JSON. The verifier reports only performed
checks; exported residual inspection is not confused with full source replay.

