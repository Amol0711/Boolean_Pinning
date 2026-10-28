# SCC-localized permanent-pinning solver

This workflow supplies the deterministic solver, examples, and independent
finite checks. It uses the Python standard library. It does not replay the
250,000-network study.

From the repository root:

```sh
python tools/reproduce.py doctor scc-validation
python tools/reproduce.py run scc-validation --output runs/scc-validation-001
```

The runner executes in a new staged workspace and preserves the distributed
baseline. See [REPRODUCTION.md](../../docs/REPRODUCTION.md) for the commands and
outputs covered by the task. The [README.txt](README.txt) documents
the API, JSON schema, support conventions, and component-budget accounting;
[verification/README.txt](verification/README.txt) describes oracle coverage.

Inputs specify physical Boolean rules and a fixed target. Baseline coordinates
are clamped permanently; reported optimum costs count additional pins. Repeated
parent slots read the same physical coordinate during dynamics. Nominal-slot
support and exact physical support are separate available conventions.

Exact-search caps apply per complex component. If a cap prevents resolution,
the result retains lower/upper bounds, a feasible upper witness, and its stop
reason; `exact_cost` is null. A certificate flag in a general deterministic
input is not evidence that the input was sampled from a random ensemble.

Frozen evidence is in `results/solver_validation.json`,
`results/reduction_validation.json`, and `results/targeted_cap_check.json`.
These records concern finite enumerated/seeded cases, not arbitrary-size
exhaustive verification. See the [results index](../../docs/RESULTS_INDEX.md)
and [fresh validation record](../../docs/VALIDATION.md) for their separate
recorded and repository verification scopes. The root manifest governs current package
integrity; the computational inventory records the distributed solver and validation files.
