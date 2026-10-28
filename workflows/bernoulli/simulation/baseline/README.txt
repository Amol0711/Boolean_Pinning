BASELINE NETWORK COMPUTATIONS

This directory retains the baseline calculations and deterministic cost-bound
refinement for the 250,000-network sample. The original sample contains five
mixtures, five sizes (100, 400, 1600, 6400 and 25600), and 10,000 replicates per
cell. Protocol_Baseline.json specifies the generator, seeds and search budgets.
Bound_Refinement_Protocol.json specifies the additional deterministic analysis
of unresolved records. The refinement draws no additional networks.

Raw baseline records are in results/raw and refined records are in
results/refined/raw. Each observation retains its identifiers, residual size,
certificate status, cost endpoints, search-work count and a feasible pin set.
Repeated nominal slots are retained; repeated references to a physical
coordinate share the same bit in state enumeration. Every observation appears
in the all-draw statistical summaries, including unresolved cost intervals.

Reproduction from this directory

  python simulation/reproduce.py --destination /path/to/new/output --workers 4

The destination must be new or empty. Add --validate-only for independent
small-network validation without generating the main sample. Add --figures to
compile the native TikZ/PGFPlots plots after a complete reproduction.

Direct validation commands

  python simulation/validate_engine.py
  python simulation/validate_refined_bounds.py
  python simulation/analyze_results.py --self-test

The baseline and refined output manifests record current source/protocol
checksums and unchanged raw-data checksums. SHA256_MANIFEST.json inventories the
current distributed files. No random draws are excluded or replaced.
