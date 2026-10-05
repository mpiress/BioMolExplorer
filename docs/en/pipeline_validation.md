# Validation of pipeline stages and connections

[Documentation](../README.md) · English · [Português](../pipeline_validation.md)

The review covers the 14-operation catalog, file contracts, block relationships, input selection and materialization, individual/merge processing, branch execution, persistence and reuse. The [user manual](user_manual.md) teaches the operational workflow; this page records how to verify it and the limits of the evidence.

## Contract matrix

| Operation | Inputs and conditions | Output for other stages |
| --- | --- | --- |
| Import files | Authorized assets validated according to their selected type | Explicitly published type |
| Retrieve compounds | Target/filters; ChEMBL query and optional PubChem expansion | Compounds and original ChEMBL context |
| Expand similar compounds | Retrieval with original ChEMBL downloads | Compounds |
| Retrieve structures | PDB criteria and external query | Raw structures and metadata |
| Retrieve ZINC | Address list of type `other` | Compounds |
| Prepare structures | Raw PDBs and reference records | Prepared receptors/ligands, centers and metadata |
| ADMET | Compounds | Evaluated/filtered compounds and visualizations |
| Fingerprints | Compounds | Fingerprints with a determined type and width |
| Similarity | Compatible fingerprints | `source,target,value` relationships and molecular lineage |
| Graphs | Relationships from similarity blocks; external data through the form | MCC compounds and visualizations |
| Redocking | Raw if `prepare_complex=true`; prepared if `false` | Structures, evaluation metadata and prepared files |
| Vina | Prepared structures and compounds | Poses identified by full reference and compound |
| DOCK6 | Prepared structures, candidates and corresponding Vina poses | MOL2 scores and other protocol results |
| Consensus | Vina poses and DOCK6 results with matching identifiers | Score table and comparison |

`tests/test_connection_matrix.py` iterates over every native producer, ready-made output and import type against every port, including both redocking configurations. Permitted connections are applied and submitted to pipeline ordering. Rejected connections leave configuration unchanged. Disabled sources, self-links and the internal consensus alias are also checked. Existing flow tests verify cycles, deletion, duplication and normalization of older projects.

Type compatibility does not replace file contracts: execution checks headers, content, paths and dependencies. Graphs only receive a canvas connection from **Calculate similarity**, preserving lineage; external CSVs can be selected directly in the form. ChEMBL expansion requires original downloads; generic ready-made results do not publish this context.

## Corrections and covered regressions

| Case | Verified behavior |
| --- | --- |
| Repeated selection | Stage configuration receives popup choices; confirmation preserves mode and references |
| Individual versus merge | Real calculations retain separate datasets or create a combined set according to the choice |
| Several Vina sources in consensus | The alias follows the entire chosen group without losing earlier sources |
| Repeated result name | An ambiguous selector fails explicitly; paths distinguish files |
| Preparation → docking | Metadata, native/legacy centers, ligands and MOL2/noH auxiliaries accompany the correct receptor |
| Prepared redocking | Four-field records are normalized; the port accepts prepared data only with preparation disabled |
| Vina with multiple references | PDB/ligand/residue/chain distinguish poses; one reference's results do not make another get skipped |
| Individual DOCK6 | Only combinations with corresponding receptors, candidates and poses execute |
| DOCK6 compound selection | The table restricts conformers passed to the engine |
| Individual consensus | Corresponding results are paired instead of crossing every file |
| Consensus scores | Integers and scientific notation are read; missing/nonfinite values are rejected |
| One compound or constant scores | Normalizations are finite and zero; native output can be imported again |
| Restart/reuse | Compatible inputs/configurations reuse results; changed data or tampered hashes prevent adoption |
| Graphs/fragments | Separated components/vertices; valid SMILES extracted from the reference, with cut stereochemistry verified |
| Language | Independent sessions, login switching, translations and preservation of parameters/user data |

Docking cases are in `test_docking_handoffs.py`; engines are replaced at the command boundary, but copies, metadata, matching and score parsing/normalization are real. `test_input_processing.py` runs fingerprints → similarity → graphs in real workers with two datasets and both modes. Other tests cover real ADMET, cancellation, timeout, isolation, permissions, import/export, history and caching.

## Reproduce verification

Run from the repository root using the scientific environment's Python with the interface installed:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/biomol-mpl PYTHONPATH=src:tests \
  python -m unittest discover -s tests -q
python scripts/build_docs.py
```

To review only the added contracts:

```bash
PYTHONPATH=src:tests python -m unittest \
  test_connection_matrix test_docking_handoffs test_localization -v
```

The connection test is exhaustive for the declared catalog. This does not execute every possible scientific parameter combination or establish provider availability and compatibility with every external engine version.

## Validate the protocol in the deployment environment

1. Confirm the worker's Python and installations of the executables used.
2. Choose a reference complex and a small set of known candidates.
3. Run retrieval/import, preparation and redocking; check ligand, chain, centers, RMSD and logs.
4. Run Vina and DOCK6 with the same receptors/identifiers; check that candidates and poses correspond to selected files.
5. Generate consensus and compare raw scores with engine files; inspect repulsion weight and normalization effects.
6. Repeat with two datasets in individual and merge modes; verify composition at each stage.
7. Restart, choose reuse and compare with a fresh run when validating environment changes.

External queries depend on network, limits and availability of ChEMBL, PubChem, PDB and ZINC. Deterministic tests do not perform a campaign of live queries. Preparation/docking need validation with real engines and data in the destination installation. Reuse after environment changes is an explicit researcher choice; to recalculate with the current installation, choose **No, run again**. ADMET rules and docking scores are computational results, not experimental proof.
