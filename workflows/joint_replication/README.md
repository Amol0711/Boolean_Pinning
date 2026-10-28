# Independent J1 replication

This workflow contains 20,000 J1 networks with 25,600 vertices and the separate
2,000-network comparison sample from the joint-rule study. The samples retain
their own random streams and summaries.

Run from the repository root:

```sh
python tools/reproduce.py doctor replication-preflight
python tools/reproduce.py run replication-preflight --output runs/replication-preflight-001
python tools/reproduce.py run replication-statistics --output runs/replication-statistics-001
```

`replication-reproduce` repeats all 20,000 production draws and compares the
scientific output files and ordered random-stream identities. The independent
full-network verifier is a separate computation described in [README.txt](README.txt).

`Protocol_J1_Replication.json` specifies the model, random calls, component search
limits and statistical conventions. `seed_manifest.csv` lists all 20,000 seeds.
`integrity_record.json` and `SHA256_MANIFEST.json` identify the current distributed
protocol, sources and validation reports. The production replay checks the
recorded NumPy 2.3.5 and SciPy 1.17.0 versions. It uses a C++17 compiler for the
accelerated backend.

The production entry point is `src/run_replication.py`. The imported
`src/joint_validation.py` supplies the shared generator, analyzer and summary
routines. The `original/` directory contains the comparison sample; its summary
and manifest describe the wider joint-rule computation. The replication results
are under `results/main/`.

All 20,000 replication costs were resolved. The sample contains no noncopy cyclic
component. Diagnostic rule frequencies do not change the fixed model probabilities
or reference values. Between-sample contrasts are exploratory because the J1
cell was selected after examining the first sample. Complete-sample moment
intervals and probability bounds describe pointwise sampling uncertainty.

See [README.txt](README.txt) for data schemas, commands and validation coverage.
