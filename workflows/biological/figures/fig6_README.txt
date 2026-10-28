Figure 6 - Exact target costs and isolated component mechanisms
==============================================================

Files
  fig6_cd4_target_mechanisms.tex/.pdf   Editable TikZ/PGFPlots vector figure
  fig6_build.py                        Data checks, export and figure build
  data/fig6_cd4_cost_distribution.csv   Panel A exact target counts
  data/fig6_bar_labels.tex              Exact segment and total labels
  data/fig6_apc_isolated_rules.csv      Panel B isolated component truth tables
  data/fig6_source_records.json         Selected records and source hashes
  data/fig6_build_metadata.json         Figure source/data/output hashes

Compile the supplied figure sources from the workflow root
  python figures/fig6_build.py --compile-only

Regenerate plotting records from the complete biological workflow
  python figures/fig6_build.py

The figure uses Python's standard library, pdflatex, TikZ, PGFPlots, Latin
Modern and pdftoppm. The rebuild uses a fixed SOURCE_DATE_EPOCH and two LaTeX
passes in an isolated temporary directory. It embeds no author metadata.
Full data regeneration requires results/all_targets.json and the
preserved CD4 BNET source and independent grammar checker under larger_networks/.

Panel A
  This panel includes all 46 fixed targets enumerated in the six CD4 input
  contexts. Counts are exact model enumerations, not sampled frequencies.
  Every displayed cost is the exact minimum number of additional permanent
  target-valued pins. The 34 environmental input clamps are excluded from c_B.
  All regulated coordinates are eligible, including the modeled proliferation
  output. The plotted context labels identify the input assignments specified
  in larger_networks/protocol/contexts.json, not independently selected
  phenotype classes. The IL2 context sets CGC, IL2RB and IL2_e to 1; all other
  inputs are 0. The APC context sets APC and MHC_II to 1; the other inputs are 0.

  The bar labels are exact counts. The numbers at each bar's right endpoint
  give the total number of fixed targets for that context. OFF, ON, ECM and
  GALPHAQ each have two targets of cost 2 and four of cost 3. APC has three
  targets of cost 3, eight of cost 4 and four of cost 5. IL2 has one target of
  cost 2 and six of cost 3. The exporter recalculates every count from the
  complete per-target cost records and compares it with the context summary.

Panel B
  The selected APC target is CD4_APC_FP01, the original fixed point reached
  from zero regulated coordinates with the prescribed environmental clamps.
  The complete 188-coordinate target is preserved in fig6_source_records.json.
  Its baseline residual contains 17 coordinates. The SCC partition is computed
  from functional dependencies over the entire baseline residual coordinate
  cube. The five displayed boxes are its five cyclic SCCs; eight acyclic
  singleton components are omitted from the boxes and contribute zero cost.

  The rules inside a box are the isolated map F_C after all coordinates outside
  that component are fixed at this target. They therefore need not retain every
  dependency used to define the baseline SCC partition. The displayed variables
  are target-error coordinates. Each of the nine displayed coordinates has
  target value 0, so its target-error value equals its original Boolean value.
  Superscript + denotes one synchronous update. The variable p denotes the
  modeled proliferation coordinate.

  In the four-node FOXP3/SMAD3/TGFB/TGFBR component, the baseline residual rule
  TGFB^+ = FOXP3 AND proliferation becomes TGFB^+ = 0 upon isolation, because
  the external proliferation coordinate is fixed at target value 0. The
  persistent FOXP3 self-dependence still requires one pin. In the RORGT
  singleton, the baseline residual rule RORGT^+ = RORGT AND TGFBR becomes zero
  upon isolation, hence this cyclic SCC has cost zero. The component costs are
  1, 1, 1, 1 and 0, and their sum is the exact additional cost 4.

  The exporter checks every displayed symbolic rule against the production
  isolated truth table and independently evaluates all 78 relevant output
  values directly from the native BNET source using the separately preserved
  recursive-descent grammar checker. It does not infer biochemical edge signs,
  clinical relevance, or experimental accessibility from these equations.

Figure caption
  Target dependence and SCC localization in the CD4 model. (A) Exact counts of
  all fixed targets, grouped by their minimum additional pinning cost c_B in
  each input context. Environmental clamps are excluded from c_B. (B) The five
  cyclic SCCs of the 17-coordinate baseline residual for the original APC
  target, with their isolated target-error rules and exact component costs.
  Coordinates outside each component are fixed at the target when forming its
  isolated map; p denotes proliferation. The RORGT component has zero isolated
  cost despite its residual self-dependence. Eight acyclic singleton components
  also contribute zero, giving c_B = 4. These costs allow every regulated
  coordinate, including the modeled proliferation output, to be pinned.

Visual review
  The final PDF is one vector page, 372.147 by 241.336 points. It was rendered
  at 180 dpi and checked for all tick labels, exact bar counts, rule notation,
  component costs and footer placement. There are no clipped or overlapping
  labels. At five-inch display width, ordinary 7.3-point figure
  text prints at approximately 7.06 points. Mathematical scripts are smaller
  by the normal TeX rules. No overfull boxes or undefined references occur.
