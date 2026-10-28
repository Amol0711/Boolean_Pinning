Baseline figure package
=================

The three PDFs are vector graphics generated with editable TikZ/PGFPlots sources.
They use the complete, independently verified post-run refined analysis tables.
Every mixture/size cell retains all 10,000 draws. Refinement changes only valid
deterministic bounds for originally unresolved draws; it adds no random samples.

Files
  fig1_cost_distribution.pdf/.tex
  fig2_certificate_failure.pdf/.tex
  fig3_residual_and_moments.pdf/.tex
  data/*.csv                    Small plotting exports, regenerated from analysis
  plot_metadata.json            Input hashes, normalization, offsets and provenance

Rebuild all figures from this package's root
  python simulation/build_figures.py --analysis-dir results/refined

The script requires NumPy, SciPy, pdflatex, PGFPlots and pdftoppm. It compiles each
figure three times in an isolated temporary directory, copies only its final PDF
back here, and renders a diagnostic PNG under qa/figures. The shared legend in
Figure 3 requires multiple LaTeX passes. To recompile an edited standalone source,
run pdflatex three times from this figures directory, with data/ alongside it.

Interpretation
  Figure 1 uses exact bins 0,...,7 and one >=8 tail bin. A marker is a lower bound
  for the empirical bin probability. Its upward bar gives the pointwise upper
  bound. These intervals are computational identification ranges, not sampling
  confidence intervals. The upper endpoints across bins need not sum to one.

  Figure 2 uses actual nonzero failure proportions, with pointwise two-sided 95%
  Clopper--Pearson confidence intervals. No pseudo-count is used. The n^-1 curve
  has prefactor one and is an illustrative reference, not a fitted rate.

  Figure 3(a) uses the exact empirical residual survival P(R>=r). Panels (b,c)
  average powers of the per-draw lower and upper control-cost bounds. Their
  normalizers are mu and mu+mu^2, with mu=-log(1-lambda) for the corresponding
  mixture. Markers denote lower endpoints and bars denote upper endpoints.
  Horizontal plotting factors A=.88, B=.94, C=1, D=1/.94, E=1/.88 separate
  overlapping intervals; these are display offsets only. True sizes are retained
  in each plotting table's n column. Sampling uncertainty is not plotted in the
  cost distribution or cost moment panels.

The figure PDFs are single-page vector graphics.
