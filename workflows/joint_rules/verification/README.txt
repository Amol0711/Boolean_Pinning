JointRules independent verification
===========================
Run these commands from the extracted JointRules package root:

  python verification/independent_exact_checks.py --output qa/replayed_exact_checks.json
  python verification/validate_exact_exports.py
  python verification/validate_joint_simulations.py --workers 4
  python verification/independent_joint_statistics.py

Dependencies: Python 3, NumPy; the independent statistical comparator also uses
SciPy for independent distribution calculations. No reference answer is obtained
by importing the production solver, engine, generator or summarizer.

The exact gate enumerates all labeled matchings and all positive-probability
truth-law assignments for twelve fixtures of one to three vertices, then checks
coefficient, distinguished-incidence and marked-cycle identities using exact
integer/rational arithmetic. The finite incidence inequality is checked at each
available k using the zero-safe constant. These computations check finite identities and numerical fixtures.

Physical checks exhaust all 1,024 two-vertex ordered two-slot planted networks,
all 4,096 pin sets and all 9,216 compatible initial states. A separate literal
trajectory evaluator cross-checks the indegree-peeling success oracle. Another
300 seeded networks of three through six vertices check global optimum,
residual capture, SCC additivity and the G-prime cycle certificate. These counts
do not imply exhaustive coverage of every larger network.

The independent finite-state report is qa/independent_exact_checks.json. REPRESENTATION.json records its current hash. A deterministic rerun writes a separate report for comparison.

The exact-export comparator checks every one of the 228 production tiny
realizations, including its rational probability, support/certificate fields,
marked cycles, physical optimum and feasible witness, and recomputes all exact
summaries. The all-AND fixture explicitly verifies positive cost with no marked
singleton cycles.

The full replay independently implements the fixed PCG64 call sequence,
rebuilds all networks, checks their hashes, computes support by synchronous
iteration and SCCs by Tarjan's algorithm, and constructs isolated physical maps
with repeated parents sharing one physical bit. Pin success is tested by
functional-graph indegree peeling. Cyclic complex components have an independent
finite guard of 18 vertices, 100,000 pin candidates and 100,000,000 compatible
states. The JSON reports actual verified counts and any guards; it never labels
a guarded optimum as independently exhaustive. Feasible witnesses are checked
by physical support capture or full remaining-state enumeration. All non-G
residuals of size at most twelve are additionally globally optimized.

Outputs are deterministic scientific records under qa/. Timing is excluded.
The main replay verifies the fixed protocol, gate and production source hashes,
and all compressed stream hashes, and consumes every stream through its end.
The existing Bernoulli simulations, figures and biological calculations are not rerun
or modified by these checks.
