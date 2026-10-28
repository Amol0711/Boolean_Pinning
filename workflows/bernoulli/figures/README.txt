Bernoulli synthetic-study figure package
================================

Figures 1-3 are vector PDFs with editable TikZ/PGFPlots sources. They use the
complete frozen sample of 250000 networks after SCC-based reanalysis. Each
mixture/size cell retains all 10000 draws, including unresolved costs.

Files
  fig1_cost_distribution.tex/.pdf
  fig2_certificate_failure.tex/.pdf
  fig3_residual_and_moments.tex/.pdf
  build_bernoulli_figures.py              Data checks, plotting exports and builds
  verify_bernoulli_plot_data.py           Independent export-to-statistics correspondence checks
  data/bernoulli_*.csv                    Exact inputs used by the plotting sources
  bernoulli_plot_metadata.json            Input hashes, offsets and interval semantics
  ../qa/bernoulli_figures/                Rendered previews, logs and build checks

Rebuild from the workflow root
  python figures/build_bernoulli_figures.py

The full regeneration reads the independent statistics under
results/bernoulli_statistics/. It verifies all 25 cells, the denominator 10000 in each
cell, interval nesting and containment of G failure counts above G-prime failure
counts. It does not derive fresh sampling intervals or impute an unresolved
cost. The statistical definitions are in verification/Bernoulli_Analysis_Protocol.json.

Compile the supplied figure sources
  python figures/build_bernoulli_figures.py --compile-only

The build requires Python's standard library, pdflatex, PGFPlots, TikZ, Latin
Modern and pdftoppm. Three LaTeX passes run in isolated temporary directories
with a fixed SOURCE_DATE_EPOCH. No author metadata is embedded in the PDFs.

Figure 1
  The nine disjoint bins are costs 0 through 7 and one >=8 tail bin. The ordinate
  is empirical bin probability minus its specified Poisson probability. Every
  marker is a lower endpoint of that difference; its thick upward bar reaches the computational
  upper endpoint. Thin wider-cap bars show the pointwise exact binomial outer
  interval supplied by the independent statistics. The lower endpoint uses
  definite bin membership and the upper endpoint uses possible bin membership.
  These are pointwise intervals and do not give simultaneous coverage over bins.
  No unresolved network is discarded. The upper computational endpoints need
  not sum to one before subtraction. The dashed zero line is the Poisson
  reference using the prescribed degree-weighted parameter. The source CSVs
  retain both the original probability endpoints and the plotted differences.
  Within the A/B and C/D panels, points are offset by -0.18 and +0.18 bin units
  solely for visibility; the integer bin k is retained in each plotting table.

Figure 2
  Panels show G and G-prime failure probabilities separately. Five color/shape
  identities retain mixtures A-E. Failure bars are pointwise two-sided 95%
  Clopper-Pearson intervals. If a count is zero, its point estimate remains zero
  in the source data and no substitute positive point estimate is plotted.
  All 50 current certificate/cell combinations have positive failure counts.
  If future data contain zero counts, the generator places an open downward
  triangle at the CP upper endpoint and labels it explicitly; the lower endpoint
  remains zero and nonzero-series connections break at that cell. The dashed 1/n curve has unit prefactor
  and indicates an illustrative slope, not a fitted rate or universal constant.

Figure 3
  The left panel shows the exact empirical survival P(|R|>=r) at n=25600. All
  zero probabilities are retained in the survival column of the plotting CSVs;
  only the plot_survival column uses nan to omit zero on a logarithmic axis.
  The plotted curves therefore stop at the observed sample maxima.

  The moment panels divide the first raw moment by mu=-log(1-lambda_star) and
  the second raw moment by mu+mu^2. Markers and thick upward bars give the
  computational ranges. Thin wider-cap bars give approximate 95% sampling
  envelopes for the lower/upper moment endpoints, as defined by the independent
  statistical analysis. Their interpretation differs from the exact binomial
  intervals in Figures 1 and 2. The moment axes are linear and contain the full
  nonnegative outer envelopes. The horizontal reference 1 is the Poisson limit.

Horizontal display offsets
  For Figures 2 and 3, x=n times 0.86,0.93,1,1/0.93,1/0.86 for A-E respectively.
  Each plotting table retains the true n alongside the display coordinate x.
  All five mixtures use the same five sizes 100,400,1600,6400,25600.

The biological figures have their own workflow and data guides under
workflows/biological/; these numerical plots use only the Bernoulli inputs.
