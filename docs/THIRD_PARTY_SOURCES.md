# Third-party sources and model provenance

Evidence checked on 17 October 2026. Paths below are relative to
`workflows/biological/`. Retain these scientific source credits when using the models. See [LICENSE_STATUS.md](../LICENSE_STATUS.md) for the
redistribution terms.

## Frozen small-model sources

`biological/provenance/source_files/{17nodeQL.py,29nodeQL.py,44nodeQL.py}`
come from [ziyuxuan67-maker/QL-and-RTS at commit
29b4f8c5e8e15c71735ab643813cfdb6de7d8d77](https://github.com/ziyuxuan67-maker/QL-and-RTS/tree/29b4f8c5e8e15c71735ab643813cfdb6de7d8d77).
The complete five-file tree and its README at that commit were inspected;
neither contains an explicit license grant. Public availability is not recorded
here as redistribution permission.

| Frozen source | SHA256 |
|---|---|
| `17nodeQL.py` | `be2aca286b7989cccd390809bf4880c301f983f4ccd7ac50917bb7ba888cd3f4` |
| `29nodeQL.py` | `22d2e8cc8d058d6e2c391838d521331bad218a777bc65b1c1670836be6dc7b25` |
| `44nodeQL.py` | `d161697c5b440a7f201eb232e2cf92d12d0a5e2aabd334c1d1472f5dfd2c7898` |

The 17- and 29-node JSON/BNET models under `biological/models/` encode the
frozen updates and supplied targets. They are particular coded implementations,
not assertions of exact agreement with every published biological model.
The 44-node source is source provenance with unresolved model-identity and
rule issues; it is not a certified biological benchmark here.

The T-LGL template copy
`biological/provenance/source_files/tlgl_topology_rewiring_test.py` comes from
the supplied source collection, with SHA256
`1f6ae48680cff8ebab3953e3fb9c5c64cf337b7fc90154f7375b9988520629f5`.
Its source description identifies templates inferred from the released QL
script. `biological/provenance/tlgl_archive_to_frozen_code.json` records the
original member and exhaustive comparison. No separate redistribution grant
was found in the supplied copy.

The named T-LGL comparison is Zhu, Cao, Lin, Lam, and Azuma, *Toward Stabilizable
Large-Scale Boolean Networks by Controlling the Minimal Set of Nodes*, IEEE
Transactions on Automatic Control 69 (2024), 174–188,
[doi:10.1109/TAC.2023.3269321](https://doi.org/10.1109/TAC.2023.3269321),
Example 6.2, equation (17). The coded Fas rule differs from that comparison;
the source-coded target has Caspase and Apoptosis equal to zero. The model must
not be relabeled as an apoptotic target or silently replaced by the comparison
rule. Exact differences and node mappings are retained in
`biological/provenance/node_maps_and_conflicts.json`.

## Drosophila comparator

`biological/provenance/source_files/drosophila_single_cell.txt` credits Albert
and Othmer (2003), *The topology of the regulatory interactions predicts the
expression pattern of the segment polarity genes in Drosophila melanogaster*,
Journal of Theoretical Biology 223, 1–18. Its SHA256 is
`bc5bba3c57b0ca749da041f89ca85235ecbc727c0a5444189ad807c86686e642`.

The frozen bytes have Git blob SHA1 `0dc263e43b88c3c1f07d495931b8f2e2d5af2853`,
matching the file in [tjparmer/influence_maximization_BNs at commit
7e8af0b63d8c4e59a57e1b21cb74a82a0a2f0dcf](https://github.com/tjparmer/influence_maximization_BNs/blob/7e8af0b63d8c4e59a57e1b21cb74a82a0a2f0dcf/drosophila_single_cell.txt).
This matching commit was inspected in the source check; it is not asserted to be
the download commit. Credit the accompanying work by Parmer, Rocha,
and Radicchi, *Influence maximization in Boolean networks*, Nature
Communications 13 (2022), 3457,
[doi:10.1038/s41467-022-31066-0](https://doi.org/10.1038/s41467-022-31066-0).

The inspected repository tree and README contain no explicit license. Its
README links [Zenodo record 6581810](https://zenodo.org/records/6581810), but the
record's particular license could not be verified in the source check; no license
is inferred from Zenodo's defaults. The comparator substitutes the SMO formula
into CIA/CIR, whereas the certified coded implementation updates dynamic SMO.
The retained comparison records both differences without replacing the model.

## Larger biological models

The two `larger_networks/data/*/model.bnet` files and accompanying metadata came
through [Biodivine Boolean Models at commit
314b2fed94689c38d3aaaf8101ed622b3501c7f7](https://github.com/sybila/biodivine-boolean-models/tree/314b2fed94689c38d3aaaf8101ed622b3501c7f7).
Credit BBM's source collection and the original model publications. Its README
identifies the associated report, Pastva, Safranek, and colleagues, *Repository
of logically consistent real-world Boolean network models*,
[doi:10.1101/2023.06.12.544361](https://doi.org/10.1101/2023.06.12.544361).

| Model | Frozen BBM identity | BNET SHA256 |
|---|---|---|
| CD4 | id-050; 154 dynamic nodes and 34 inputs | `24722017d5be15ab31cf9a8c3dfc31139a80313c7a69196f8ab77b0945444209` |
| Macrophage | id-001; 302 dynamic nodes and 19 inputs | `ea60f90ecab55c939f44648dc24608d5779916f7d1df54200f1e85681a4951d4` |

Each model's `README.md` and `metadata.json` retain the complete original
bibliographic citation and Cell Collective link. The CD4
`download_manifest.json` and macrophage `immutable_download_manifest.json`
contain exact pinned file URLs and hashes. Keep those records with the models.

For CD4, credit Conroy and colleagues, *Design, Assessment, and in vivo
Evaluation of a Computational Model Illustrating the Role of CAV1 in CD4+
T-lymphocytes*, Frontiers in Immunology 5 (2014), 599,
[doi:10.3389/fimmu.2014.00599](https://doi.org/10.3389/fimmu.2014.00599).
The publisher's article notice specifies CC BY with author and publication
attribution. For macrophage, credit Raza and colleagues, *A logic-based diagram
of signalling pathways central to macrophage activation*, BMC Systems Biology
2 (2008), 36, [doi:10.1186/1752-0509-2-36](https://doi.org/10.1186/1752-0509-2-36).
The publisher specifies CC BY 2.0 for that article. The tested macrophage
network is the Cell Collective/BBM Boolean realization associated with the
pathway map; it is not described as a direct transcription of Boolean rules
from the original article.

The pinned [BBM README](https://github.com/sybila/biodivine-boolean-models/blob/314b2fed94689c38d3aaaf8101ed622b3501c7f7/README.md)
reserves model copyright to the respective authors and/or publishers, while
allowing reproduction of other information. These article and repository
notices do not settle the licensing scope of every later model export. That
specific redistribution question remains open.

## Derived records and software

Model expressions and truth tables also occur in the small-model JSON exports,
`biological/provenance/source_rule_validation.json`, and residual rows in
`results/all_targets.json` and component certificates. Retain their provenance;
removing original source copies alone would not resolve the scope question for
these representations.

The folders called `vendor` contain shared project solver/engine
modules reused across workflows. No third-party Python library distribution is
bundled. External Python and TeX dependencies are obtained separately. The
project-written software and generated results still require an explicit
rights-holder license decision; this provenance record supplies none.
