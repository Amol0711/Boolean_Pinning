# Computational figures

This directory contains numerical plots, finite-state illustrations, their data,
and native LaTeX/TikZ/PGFPlots sources. Figure identifiers name repository assets.

Run from the repository root with Python requirements installed and assertions enabled.

```sh
python tools/verify_figures.py
python tools/build_figures.py --output runs/figures-001
```

The verifier checks an isolated copy. The build checks the data, regenerates the
numerical sources, and compiles all nine views. It requires `pdflatex`, `pdfinfo`,
`pdftotext`, `pdftoppm`, TikZ/PGFPlots, the standalone class and Latin Modern.
The output directory must be new. Both commands check the repository integrity
manifest and leave the distributed data unchanged.

| Native stem | Computational content |
|---|---|
| `fig1_cost_distribution` | Cost probability intervals and specified Poisson comparisons |
| `fig2_certificate_failure` | Residual and componentwise condition failure frequencies |
| `fig3_residual_and_moments` | Residual survival and normalized sample moments |
| `fig5_larger_network_support_grouped` | CD4 support bounds and grouped macrophage orbit |
| `fig6_cd4_target_mechanisms` | CD4 target counts and isolated component maps |
| `fig4_biological_certificates` | Synchronous trajectories and finite-state control examples |

Additional views are `fig1_cost_differences`, `fig5_larger_network_support`, and
`fig_condition_moment_overview`. Every view has a `.tex` source and a `.pdf`
output. `FIGURE_CONTENT.json` lists the views and the seven-coordinate update
map used in the exhaustive finite-state check. `FIGURE_INVENTORY.json` records
current file hashes. Statistical input hashes are in `bernoulli_plot_metadata.json`.

The probability checker verifies 45 cost bins, Poisson reference probabilities,
computational endpoints and pointwise Clopper–Pearson intervals. The final bin
includes costs of at least eight, and all 10,000 draws per mixture contribute.
The condition and moment checker verifies 50 failure counts, 25 moment records,
50 normalized ranges, and 563 residual-survival records. Moment envelopes use
the retained endpoint standard errors. The `1/n` curve is a slope guide.

The regulatory checks compare 46 fixed-target counts, nine isolated rules, five
copy arrows, support endpoints, and the 19-coordinate orbit with the retained
computational records. The four grouped profiles reconstruct all 95 orbit cells.
Small component maps are enumerated exhaustively. Full-model census and support
enumeration are separate workflow tasks described in `docs/REPRODUCTION.md`.

These checks validate the retained data and plotting transformations. They do
not regenerate the random-network samples or recompute moment standard errors.
The build compares each compiled view with its distributed PDF by extracted
text, dimensions and rendered pixels. Build logs and reports are written to the
requested output directory.
