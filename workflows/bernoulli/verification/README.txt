INDEPENDENT BERNOULLI VERIFICATION

The reference implementations do not import the production support, graph,
solver, random-generator wrapper or aggregation functions. NumPy PCG64 and
SeedSequence provide the prescribed random stream. The reference generator
implements the call sequence specified in the baseline protocol.

validate_bernoulli_networks.py regenerates the seeded networks, checks support
and nominal-slot SCC structure, and verifies cost bounds and feasible pin sets.
Its reports distinguish complete finite-state checks, copy-component checks,
preserved baseline bounds and unresolved cases. targeted_component_check.py
covers the recorded component outside the general verifier's state limit.

independent_statistics.py reads all 25 cells from the baseline, refined and SCC
records. It checks row identities, sample retention, interval nesting,
probability and moment exports, residual distributions, and sampling intervals.
The output includes the actual row coverage and numerical comparison counts.
All 25 complete SCC raw files are required for this check.

bernoulli_stats.py computes the complete-sample statistical exports.
Bernoulli_Analysis_Protocol.json documents their operational definitions and
pointwise uncertainty conventions. Numerical boundary tests independently
compare the probability calculations against SciPy.

Compressed records retain the original network identifiers and explicit
physical-coordinate mappings. Network SHA256 uses parents in little-endian
64-bit row-major bytes followed by the three rule entries as uint8 bytes.
