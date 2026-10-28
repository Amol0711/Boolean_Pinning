# Bernoulli-network simulations

This workflow contains 250,000 networks in 25 mixture/size cells, with 10,000
observations per cell. The SCC calculations produce 249,208 exact costs and
792 cost intervals. Every observation is retained in the statistical summaries.

| Root task | Computation |
|---|---|
| `bernoulli-quick` | Independent small-network checks and plot-data validation |
| `bernoulli-statistics` | Complete-sample statistical aggregation and independent verification |
| `bernoulli-boundary` | Statistical special-function boundary checks |
| `bernoulli-pilot` | Reproduce the retained 2,500-network pilot subset |
| `bernoulli-reproduce` | Reproduce all 250,000 networks and compare saved outputs |
| `bernoulli-figures` | Rebuild the native TikZ/PGFPlots figures |

Run these commands from the repository root.

```sh
python tools/reproduce.py doctor bernoulli-quick
python tools/reproduce.py run bernoulli-quick --output runs/bernoulli-quick-001
```

The root reproduction guide lists dependencies and commands. Python 3.11 or
newer is required; C++17 `g++` enables the optional numerical accelerator.

`simulation/Protocol_Bernoulli.json` specifies the seeds, generation procedure,
and component search budgets. `simulation/baseline/` contains the baseline
calculations on the same networks. The pilot is a subset of these observations.
The complete SCC records are in `simulation/results/main/`, and the statistical
exports are in `results/bernoulli_statistics/`.

Repeated nominal input slots remain distinct in graph calculations, while
repeated physical coordinates share one bit during state enumeration. Every
unresolved cost is represented by an interval and a feasible upper-bound pin
set. Statistical uncertainty is reported separately using the definitions in
`verification/Bernoulli_Analysis_Protocol.json`.

See `README.txt` for direct commands and `simulation/README.md` for record
schemas and computational conventions. `SHA256_MANIFEST.json` records the
current distributed file bytes. All compressed per-network data are retained.
