# Joint-rule simulations and validation

This workflow contains 18,000 networks across nine family/size cells and exact
finite-state validation fixtures. Each local truth table is sampled from a
specified eight-code probability law.

Run from the repository root:

```sh
python tools/reproduce.py doctor joint_rules-exact
python tools/reproduce.py run joint_rules-exact --output runs/joint_rules-exact-001
python tools/reproduce.py run joint_rules-statistics --output runs/joint_rules-statistics-001
```

`joint_rules-reproduce` regenerates the exact fixtures and all production draws,
then compares the scientific files. Full independent reconstruction of the
networks and pinning witnesses is available through
`verification/validate_joint_simulations.py`.

`validation/Protocol_JointRules.json` records the model laws, random streams,
network sizes, search limits and statistical definitions. Exact fixtures are in
`validation/results/exact/`; the network sample is in `validation/results/main/`.
All 18,000 recorded costs were resolved. Independent validation reports are in
`qa/`. Current per-file checksums are in `SHA256_MANIFEST.json`.

The separate J1 replication retains its own sample and does not replace these
observations. Fixed Poisson reference values are computed from the declared
probability laws. Moment intervals and Clopper–Pearson probability bounds describe
pointwise sampling uncertainty. Every bounded search retains its interval, stop
reason and feasible pinning witness.

See [validation/README.md](validation/README.md) for model definitions and output
schemas, and [verification/README.txt](verification/README.txt) for independent
checks and search limits. Python 3.11 or later, NumPy and SciPy are required; a
C++17 compiler enables the accelerated production backend.
