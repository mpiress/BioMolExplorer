# Validation of pipeline stages and connections

[Documentation](../README.md) · English · [Português](../pipeline_validation.md)

The review covers the operation catalog, file contracts, block relationships, input selection and materialization, individual/merge processing, branch execution, persistence and reuse. The [user manual](user_manual.md) teaches the operational workflow; this page records how to verify it and the limits of the evidence.

## Contract matrix

| Operation | Inputs and conditions | Output for other stages |
| --- | --- | --- |
| Import files | Authorized assets validated according to each file's type | One or more types published by the selected set |
| Retrieve compounds | Target/filters; ChEMBL query and optional PubChem expansion | Compounds and original ChEMBL context |
| Retrieve PubChem | Manual reference, file or connected compounds | New PubChem similar compounds with identifiers and SMILES |
| Expand similar compounds | Retrieval with original ChEMBL downloads | Compounds |
| Retrieve structures | PDB criteria and external query | Raw structures and metadata |
| Retrieve ZINC | Tranche list of type `zinc_urls` or `other` | `compounds.csv`, MOL2 conformations and report |
| Prepare for docking | Raw or redocking PDB receptor; one or more compound sources | Receptors, centers, metadata and prepared candidates for Vina, DOCK6 or both |
| ADMET | Compounds | Evaluated/filtered compounds and visualizations |
| Fingerprints | Compounds | Fingerprints with a determined type and width |
| Similarity | Compatible fingerprints | `source,target,value` relationships and molecular lineage |
| Graphs | Relationships from similarity blocks; external data through the form | MCC compounds and visualizations |
| Redocking | Raw if `prepare_complex=true`; prepared if `false` | Structures, evaluation metadata and prepared files |
| Vina | Prepared structures and compounds | Poses identified by full reference and compound |
| DOCK6 | Prepared structures and candidates from any molecular source; optional Vina poses | MOL2 scores and other protocol results |
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
| Preparation → docking | Metadata, native/legacy centers and MOL2/noH auxiliaries accompany the correct receptor; candidates arrive separately without reference ligands |
| Prepared redocking | Four-field records are normalized; the port accepts prepared data only with preparation disabled |
| Vina with multiple references | PDB/ligand/residue/chain distinguish poses; one reference's results do not make another get skipped |
| Individual DOCK6 | Only combinations with corresponding receptors, candidates and poses execute |
| DOCK6 compound selection | The table restricts conformers passed to the engine |
| Individual consensus | All selected results are merged before intersection by receptor and compound |
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

## Redocking and language regressions

The review covers pandas structured records (which cannot be sliced like lists), configured selection without another prompt, missing pairs, invalid residues, missing tools, paths with spaces, cofactors and prepared centers. Preparation errors propagate their cause instead of returning an empty set.

Tests traverse forms for every operation in English, verify nested worker messages and preserve scientific values and custom names. Template labels, default stage names and selectors follow the session language. Messages originally written in English are translated for Portuguese sessions.

Automated integration tests replace external tools at the execution boundary. Following the October 6, 2026 log review, the real 4M0E / 1YL / 604 / A case was also executed using Chimera, Open Babel, Vina and PyMOL: with the original options, including receptor and ligand minimization, redocking completed with an RMSD of approximately 0.149 Å. A complementary run without minimization produced approximately 0.215 Å. Both runs used temporary copies of the input; these results verify this case's execution flow and do not establish scientific tolerances for other complexes.

## Results and 3D visualization

`tests/test_redocking_results.py` covers simulation grouping, imported collections with identical identifiers, duplicate artifacts, invalid RMSD values, repeated reference/pose filenames in ZIP archives, permissions, revoked sessions and dialog actions in Portuguese and English. `tests/test_pdb_viewer.py` validates formats and viewer access.

Browser checks also covered a MOL2 ligand (35 atoms), a PDBQT receptor (4,142 atoms) and PDBQT poses (21 atoms), including visible geometry, rotation, zoom, four representations and PNG export, without JavaScript errors or external requests. These checks complement automated tests and do not replace scientific evaluation of the structures. To repeat with your own file, run:

```bash
python scripts/validate_pdb_viewer.py --structure /path/to/structure.mol2
```

See [Logs and diagnostics](logging.md) for the common format, execution context, failure codes, job summary and `python -m biomolexplorer.log_report` command.


## Interoperable docking: October 7, 2026 validation

The real `scripts/validate_docking.py` check used the 1ABE receptor and charged reference ligand from DOCK6's distributed `ligand_sampling_demo`, plus ethanol supplied as SMILES. Vina and DOCK6 ran independently, followed by Vina → DOCK6 and DOCK6 → Vina. All four runs retained `REF` and `ETHANOL`, produced finite scores and exported poses. All three consensus cases (independent and both handoffs) intersected the two compounds against `1ABE_A`; the empty-intersection case returned `skipped_reason` without a score table. Values are preserved in the [JSON report](../validation/docking_2026-10-07.json).

| Compound | Independent Vina | Independent DOCK6 | Vina → DOCK6 | DOCK6 → Vina |
| --- | --- | --- | --- | --- |
| REF | -6.511 | -12.755546 | -23.418610 | -6.509 |
| ETHANOL | -2.673 | -13.771717 | -11.354554 | -2.673 |

The check uses Vina exhaustiveness 1, two poses and rigid DOCK6 search. It validates execution and software contracts for this case without defining scientific quality criteria or comparing scores between different engines. Consensus adjusts the DOCK6 score using repulsion weighting, so it may differ from the raw scores above.

To repeat, activate the scientific environment and choose a new output directory:

```bash
PYTHONPATH=src python scripts/validate_docking.py --dock6-root /path/to/dock6 --output /tmp/docking-check
PYTHONPATH=src:tests python -m unittest test_docking_interoperability test_docking_handoffs test_connection_matrix test_results_display test_localization
```

Regressions cover isolated inputs and CSV aliases, poses in both directions, identity of individually selected consensus poses, multiple batches, different receptors, matching IDs with different structures, audited removal and authorized pose access. Empty curated summaries remain valid consensus inputs; auxiliary tables cannot resurrect their compounds. DOCK6 calculations use copies in short temporary paths and export artifacts back to the project, avoiding legacy helper truncation; subprocess failures retain their diagnostics across the parallel executor.

Consensus PDBQT and MOL2 poses were also checked in the browser: rotation, zoom, four representations, PNG export and no external requests or JavaScript errors. PDBQT loading removes torsion records that would otherwise be mistaken for model boundaries, preserves atoms in each pose and identifies files with Vina scores as ligands even after consensus renaming.

The final general suite ran 442 tests with no failures or errors; two localhost-server checks were skipped due to sandbox restrictions. Browser tests ran with authorized local access and verified all 14 atoms of the Vina pose and 20 atoms of the DOCK6 pose. A complementary **flexible** DOCK6 calculation reused grids and preparation from the independent case and produced scores -38.876862 for REF and -15.546004 for ETHANOL.

To check local links and anchors after documentation changes, run `python scripts/build_docs.py` and `python scripts/validate_docs.py`. Validation covered 27 HTML pages: `index.html`, 24 guides and two language portals.

## Prepare for docking: validation on October 8, 2026

Regressions in `test_preparation_settings` and `test_docking_preparation` cover raw and reused receptors, selection without reference-ligand files, independent settings, combined ChEMBL/PubChem/ZINC sources, each engine's exports from one preparation and prepared-file reuse without further conversion. Invalid centers, missing formats and mixed raw/prepared receptors are rejected. Connections respect selected engines, and rejected connections preserve receptor settings.

Activate the scientific environment to repeat focused tests and validation with real tools:

```bash
PYTHONPATH=src:tests python -m unittest test_preparation_settings test_docking_preparation test_docking_handoffs test_connection_matrix test_docking_interoperability
PYTHONPATH=src python scripts/validate_docking.py --dock6-root /path/to/dock6 --prepare-inputs --output /tmp/preparation-check
```

The second command reuses the tutorial receptor, prepares candidates with Chimera and Open Babel and runs Vina, DOCK6 and both engine sequences. Selection and multiple-source tests use local data; they do not query live providers.

All 48 focused tests passed. The real `--prepare-inputs` run prepared REF and ETHANOL, preserved the receptor byte for byte and validated independent Vina and DOCK6 runs, both engine sequences and three consensuses. See the [JSON report](../validation/preparation_2026-10-08.json). The case with no common compounds returned a skip reason without scores. This run verifies this case's software contracts; it does not establish scientific quality for other complexes.

The final general suite ran 465 tests without failures or errors, with two localhost-server checks skipped due to sandbox restrictions. Documentation review checked local links in 26 Markdown documents and files, IDs and anchors across 27 HTML pages.

## Import and ZINC tranches: review on October 8, 2026

All 75 focused tests passed. `test_import_files` covers disk/project selection, the table, removal without deleting originals, per-file types, type-change validation, receptor/compound compatibility and read-only mode. `test_zinc_tranches` covers URI lists and scripts, HTTPS, compressed formats, headerless SMI, column order, individual MOL2 files, duplicates, conflicts, redirects and multiple-list materialization. It also checks binding-site placement of library conformations for DOCK6 and original-file preservation.

A real download of the public `AAIA.smi` tranche produced 11 compounds with ZINC identifiers and valid SMILES. The [JSON report](../validation/zinc_2026-10-08.json) records the URL and hash. Compressed MOL2 reading and conformation contracts were checked with deterministic responses; this live test did not download an entire 3D library.

Repeat focused tests with:

```bash
PYTHONPATH=src:tests python -m unittest test_import_files test_zinc_tranches test_preparation_settings test_localization test_docking_preparation test_docking_handoffs test_guided_ui test_connection_matrix
```

The general suite for this review ran 482 tests without failures or errors; two localhost-server checks were skipped due to sandbox restrictions. Local links in 26 Markdown documents and files, IDs and anchors in 27 HTML pages were also validated.

## Help and docking configuration: review of October 8, 2026

`test_pipeline_guidance` covers rejection of ChEMBL/PubChem → similarity through the canvas action, preserving connections and history and creating a valid connection after rejection. It also verifies guidance for all blocks for viewers and translation of explanations.

Compound selection tests cover both engines: identifier filtering, the all-compounds default, per-file selection, file changes, popup state preservation, copying prepared molecular files, separate materializations for different selections and a missing-code message. In real pipeline service execution, the scientific supervisor is replaced by a deterministic response: explicit inputs execute without another popup; automatic connections still wait for selection.

```bash
PYTHONPATH=src:tests python -m unittest test_pipeline_guidance test_connection_matrix test_visual_flow test_file_selection test_pipeline_execution test_localization
```

The full suite for this review ran 495 tests with no failures or errors; two local server checks were skipped due to sandbox restrictions. All eight new tests passed. The 27 HTML pages were checked for local files, IDs and anchors.
