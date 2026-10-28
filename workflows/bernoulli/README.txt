BERNOULLI NETWORK COMPUTATIONS

Run the following commands from this workflow directory.

1. python simulation/reproduce.py --workers 4
   Reproduce the 250,000-network sample and compare the 50 compressed scientific
   records with the saved outputs. Use --out for a new destination or --pilot
   for the retained first 100 networks in each cell (2,500 total).
2. python verification/validate_bernoulli_networks.py --workers 4
   Independently regenerate the networks and check residual support, SCC
   certificates, cost bounds and pin witnesses within the declared state limits.
3. python verification/targeted_component_check.py
   Check the exact component beyond the general verifier's state limit.
4. python verification/bernoulli_stats.py
   Rebuild the complete-sample statistical exports.
5. python verification/bernoulli_stats.py --special-function-checks-only
   Check 10,001 binomial count cases and 200 chi-square boundary cases.
6. python verification/independent_statistics.py
   Independently check all raw aggregation and exported bounds against SciPy.
7. python figures/build_bernoulli_figures.py
   Rebuild three computational figures using TikZ/PGFPlots.
8. python figures/verify_bernoulli_plot_data.py
   Check every plotted record against the statistical exports.

Python 3.11 or newer is required. Install requirements.txt. The optional C++17
accelerator requires g++. Figure compilation requires pdflatex, TikZ/PGFPlots,
Latin Modern and pdftoppm.

The complete data are in simulation/results/main. Baseline calculations on the
same sampled networks are in simulation/baseline. Simulation and analysis
protocols record the sample design, search budgets and statistical definitions.
All unresolved costs remain intervals with feasible upper witnesses. Moment
envelopes are approximate and pointwise; binomial intervals are exact and
pointwise. Finite-size goodness-of-fit diagnostics use the specified reference
probabilities without fitting parameters to the observed costs.

SHA256_MANIFEST.json describes the current package. Runtimes, local paths and
invocations can differ on reproduction; the scientific records, seeds and
search budgets determine the data comparisons.
