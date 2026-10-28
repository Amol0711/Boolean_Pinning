Figure 5 — Support reduction and a verified synchronous cycle
===========================================================

Files
  fig5_larger_network_support.tex/.pdf   Editable vector figure
  data/fig5_cd4_apc.csv                  Panel A support cardinalities
  data/fig5_mac_bacteria_cycle.csv       Panel B complete varying-coordinate raster
  data/fig5_mac_bacteria_full_orbit.csv  All 321 coordinates at updates 0 through 15
  data/fig5_mac_bacteria_constant_coordinates.csv
                                       All 302 coordinates constant on the cycle
  data/fig5_mac_bacteria_status.csv      Retained original discovery-status record
  ../qa/figure5_build.json              Source, protocol, data, and output hashes
  ../qa/figure5_validation.json         Visual and data-preservation checks

Rebuild from the workflow root
  python scripts/build_larger_figure.py

The script requires Python, pdflatex, PGFPlots, TikZ, and Latin Modern. The native
BNET checker uses only the Python standard library. No author metadata is
embedded in the figure PDF. The PDF contains vector marks and embedded text.

Selection
  CD4_APC and MAC_BACTERIA were selected by their input identities before the
  seven focal-context computations. The eleven-context protocol distinguishes four endpoint contexts observed before focal-context selection from seven focal contexts selected before their computations.
  Panel B shows the actual recorded cycle;
  no successful context has been substituted for the bacterial condition.

Panel A
  Every ordinate is an exact coordinate-cube support cardinality. The curves
  bound coordinates that may differ from a fixed target; they are neither
  Hamming distances along one trajectory nor exact reachable-image supports.
  The initial support excludes all clamped coordinates. The selected additional
  four pins are FOXP3=0, GATA3=0, TBET=0, and proliferation=0. The latter is a
  modeled output coordinate, not automatically an experimentally available
  intervention. Environmental inputs are APC=MHC_II=1 and the other 32 inputs=0.
  The baseline support retains 17 coordinates. The selected controller's
  support is empty at round 13. The final stationary values are repeated once;
  CSV flags identify that repetition. Panel A's CSV contains the exact support sequence for that selected target.

  The four-pin controller attains the exact additional cost. Every smaller
  cardinality was excluded, but the enumeration of all optimal quadruples is
  incomplete. No complete optimal-controller list is inferred from this plot.

Panel B
  This is a binary state raster from one synchronous orbit, with no target
  subtraction and no support-envelope interpretation. The 302 regulated
  coordinates start at zero. Of the 19 input coordinates, Bacteria and
  External_Activator are 1 and the other 17 are 0. All inputs remain clamped.
  The first repeated full state is x(15)=x(11). The four distinct cycle states
  are at updates 11, 12, 13, and 14. The final displayed column repeats the first
  phase and is separated by a dashed line. Dark cells contain 1; white cells 0.

  All 19 varying coordinates are displayed separately and in source order.
  Only the common source-name prefix v_ is omitted from the labels. The other
  302 coordinates are constant across the four cycle phases: 283 regulated
  coordinates and all 19 inputs, with 171 zeros and 131 ones. The orbit CSV and
  constant-coordinate CSV preserve their full identities and values.

  Before plotting, the script checks the native model.bnet SHA-256, parses it
  using verification/independent_core.py (which does not import the production
  engine), and replays all 15 transitions from the declared initial state.
  The first cycle entry, minimal period, all full-state cycle vectors, and all
  four cycle edges must agree with the preserved certificate. The varying-
  coordinate raster is then an exact projection of that verified full orbit.

  A single cycling orbit does not prove that the context has no fixed points.
  It explains why the prescribed target-discovery procedure did not produce a
  fixed target. No fixed-target convergence or control-cost claim is assigned
  to this panel. External_Activator is a repository background input, and the
  aggregate Bacteria input is not a pure LPS experiment. These observations
  refer to the deterministic synchronous Boolean model used in this study.

Figure caption
  Support reduction and the recorded orbit in the two prespecified focal
  contexts. (A) CD4_APC clamps APC and MHC II at 1 and the other 32 inputs at 0.
  The dashed curve uses these input clamps alone; the solid curve also pins
  FOXP3, GATA3, TBET, and the modeled proliferation output at 0. Its empty
  coordinate-cube support at round 13 certifies convergence within 13 updates
  from every compatible state. Support cardinalities are bounds on possible
  deviations, not simulated Hamming distances. (B) MAC_BACTERIA clamps Bacteria
  and External_Activator at 1 and the other 17 inputs at 0. The raster displays
  all 19 coordinates varying on the four-cycle reached from zero regulated
  coordinates. The other 302 coordinates are constant on this cycle. Update
  15 repeats update 11; the source prefix v_ is omitted from row labels. This
  panel depicts one actual synchronous orbit, for which no fixed-target
  certificate is applied.
