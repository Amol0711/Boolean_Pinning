General two-input joint-rule validation (JointRules)

This standalone package contains the fixed 18,000-network focused study and
its exact finite checks. It is separate from the 250,000-network Bernoulli
study and the biological computations.

Install requirements.txt with Python>=3.11. A C++17 compiler enables the original
accelerated certificate engine. From the extracted package root:

1. python validation/reproduce.py --workers 4
   Repeats all exact/main calculations in a new temporary directory and compares
   all23 scientific files bytewise, including all18 compressed main streams.
2. python verification/independent_exact_checks.py --output qa/replayed_exact_checks.json
   Rechecks exact rational coefficient/incidence/cycle identities and physical
   exhaustive small networks without overwriting the distributed validation report.
3. python verification/validate_exact_exports.py
   Independently verifies all228 production tiny outcomes and32support patterns.
4. python verification/validate_joint_simulations.py --workers 4
   Independently regenerates and verifies all18,000 networks, optima and pins.
5. python verification/independent_joint_statistics.py
   Independently recomputes every complete-sample summary and sampling interval.

Read validation/README.md for laws, seeds, call sequence, slot conventions,
search budgets, schemas and statistical methods. Read verification/README.txt
for independent oracle methods and finite verification guards. All18,000 costs
are exact, all106 complex SCC optima were independently verified, and no guard
was reached. Pointwise sampling intervals summarize the observed finite-size sample.

SHA256_MANIFEST.json covers every other distributed file, including nested manifests. Operational times and paths can differ on reproduction;
scientific payloads and seeds determine the byte comparisons.
