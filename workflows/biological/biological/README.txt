EXACT BOOLEAN-MODEL CERTIFICATES

This package contains two explicit coded Boolean models, a general local-cube
support engine, exact control certificates, and independent finite-state tests.
It requires Python 3.10 or later and the standard library only. No original
learning script is executed and no third-party dependency is needed.

REPRODUCTION

From this directory run

    python src/reproduce.py --output /absolute/path/to/a/new/results_directory

The output directory must not already exist. The command writes deterministic
JSON certificates, validation results, a CSV convergence curve, and a manifest
with SHA256 hashes of the model/engine/test inputs and generated outputs. To
compare a rerun with the delivered results, compare the contents of their files;
the output path is not embedded. The included results directory is a completed
reference run. It must not be used as the new output path.

MODEL AND CONTROL CONVENTIONS

Each network is a deterministic synchronous Boolean map. Pins are permanent,
take the specified full-state fixed-point target values, and are applied from
time zero. Global success means convergence from every initial state compatible
with those pins. All nodes are eligible and each pin has unit cost unless the
explicitly stated fixed-baseline convention is used.

Source/engine indices are zero based. The BNET files use one-based
neutral labels n01,...,nNN. No gene-label mapping or biological reconstruction
fidelity is asserted. The JSON files provide the exact mapping, complete target,
complete physical truth tables, and Boolean formulas. A truth-table index uses
bit j for parents_zero_based[j]. Parents are distinct physical variables.
Repeated occurrences of one variable always receive the same value. Relabeling
uses y_i=x_i XOR target_i before support reduction.

The original scripts also contain learning/individual-based mean-field routines
and automatic freezing of deterministic coordinates. Those routines are not
part of the certified dynamics. The extracted update assignments define the
Boolean map; only the chosen permanent pin coordinates are clamped.

29-NODE CODED IMPLEMENTATION

models/tlgl_coded_29.json and .bnet reproduce TRUE_RULES and TARGET from the
bundled provenance/source_files/tlgl_topology_rewiring_test.py. That
source explicitly describes its rule templates as inferred from the released
29-node QL script. The model JSON records local source hashes and provenance.
The templates and target were also checked against the frozen 29nodeQL.py
polynomial assignments and target at all 88 source-local Boolean corners. The
export contains 82 truth-table corners because source row n18 has two irrelevant
references and simplifies to NOT n17, reducing its eight source corners to two.
The source polynomial expressions are retained as neutral model data and checked
against the exported BNET rules during every reproduction run.

Physical identity rows n01 and n09 are mandatory if every initial state must
reach the specified target. With those pins the descending support sizes are
27,21,12,8,5. The remaining nodes are n15,n16,n17,n18,n19. In target-relative
coordinates their map is

    y15'=y16, y16'=y15, y17'=y15, y18'=y17, y19'=y18.

The only cycle is n15 <-> n16. Exactly one of its nodes must also be pinned.
Thus all optimal total sets are {n01,n09,n15} and {n01,n09,n16}; the total
optimum is three pins. If n01,n09 are imposed as free background inputs, the
additional optimum is one. Both optimal interventions bring every compatible
initial state to the target in at most four updates. An explicit trajectory
requiring four updates proves this bound is exact. A non-target fixed point
under only the two mandatory pins supplies an explicit failure witness.

The global statement follows from local-cube checks and the exact copy-cycle
argument, without enumerating 2^29 states. As an independent check, the reduced
five-node state graph is enumerated. These sets are not claimed to be new control
sets, and the supplied full-state target is not reinterpreted as a phenotype.

17-NODE CODED IMPLEMENTATION

models/drosophila_coded_17.json and .bnet contain the exact Boolean corners of the
frozen 17nodeQL.py assignments before its automatic-freezing loop. The source
hash, immutable commit URL, assignment expressions, and supplied target are in
the JSON. All 80 source-local Boolean corners are checked against the exported
BNET formulas during every reproduction run. This is the source-coded target,
without a claim that it is a particular published wild-type state.

Identity rows n01,n02,n03 are mandatory. With exactly those three pins, the
support sizes are 14,12,11,9,7 and remain seven. The residual nodes are
n04,n05,n10,n11,n13,n16,n17, and their map is nonlinear. The copy-cycle certificate
therefore does not apply. Exact enumeration of all 128 residual states shows
that they nonetheless converge to zero, with worst transient eight updates.
An independent enumeration of all 16,384 compatible states of the ORIGINAL
17-node map confirms global convergence and an exact worst transient of eleven
updates. The three mandatory pins therefore constitute the unique optimal total
set. The additional optimum is zero when they are free background inputs.

The support-capture bound (four) plus residual-transient bound (eight) gives the
valid analytic upper bound twelve. Eleven is the sharper exact full-state value.
The nonempty cube residual is conservative: a cube can contain state combinations
that never coexist in the exact reachable set. It does not establish a failure
of global convergence. No extra pin is introduced merely to empty that cube.

convergence_curves.csv compares iterated cube support with the union of deviation
coordinates over the EXACT successive reachable sets of all 16,384 compatible
initial states. The latter is a union across states, not maximum Hamming distance
within a single state. The CSV also gives the number of distinct reachable
states at every update from zero through eleven.

GENERAL ENGINE AND VALIDATION

src/boolean_certificates.py supports arbitrary-arity local truth tables, nonzero
fixed targets, repeated-reference compilation, exact cube supports, fixed
baseline pins, physical residual maps, and pure-copy residual certificates.
It does not assert a general noncopy-row approximation guarantee.

tests/independent_oracle.py enumerates physical compatible states and their
successor graph. It does not use support reduction or the copy-cycle criterion.
tests/validate_engine.py checks the following independently of the two examples.

  - All 256 pairs of a two-node target and target-compatible Boolean map.
  - 72 seeded three-, four-, and five-node maps with varying arity and targets.
  - A repeated-reference model, for 329 total test networks.
  - Exact support cubes against direct state evaluation, 2,052 cubes.
  - Global/reduced success equivalence for 7,713 baseline/additional-pin pairs.
  - 1,680 pure-copy certificates against exhaustive optimum and minimum sets.
  - Rejection of five invalid input cases.
  - The nonmonotone-feasibility example a'=a AND NOT b, b'=a: the empty pin set succeeds,
    but adding a permanent zero pin at b fails.
  - Every exported BNET local truth-table corner and every retained source
    polynomial local corner of both supplied models.

These finite tests validate the implementation on the stated cases. The seed
2026101401 is used only for the small validation suite. No new random biological
data are generated.

OPTIONAL SOURCE EXTRACTION

src/extract_models.py parses source text as AST data using a restricted integer
arithmetic interpreter. It never imports or runs the original source. To recreate
the model files with access to the original sources, use an existing EMPTY output
directory and the following command (replace all paths):

    python src/extract_models.py --source29-template provenance/source_files/tlgl_topology_rewiring_test.py --source17 provenance/source_files/17nodeQL.py --source29 provenance/source_files/29nodeQL.py --output /path/empty_model_directory

The optional --source29 argument enables the independent original-polynomial
comparison and was supplied when preparing this workflow. Normal certificate
reproduction uses the included explicit model files and requires no source ZIP.

This workflow contains model data, certificate code, and verification outputs.
