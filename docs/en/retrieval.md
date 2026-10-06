# Flexible structure and compound retrieval

[Documentation](../README.md) · English · [Português](../retrieval.md)

Retrieval stages can start independent branches. Choose structures, activity-associated compounds, direct compound searches or similarity expansion. Review retrieved files before passing them to the next stage.

## PDB structures

EC numbers and enzyme names are optional. `target` names the output collection under `PDB/<collection>/` and defaults to `Estruturas`. Supply at least one actual search criterion:

| Criterion | Example | Meaning |
| --- | --- | --- |
| `pdb_query` | acetylcholinesterase, EGFR | RCSB textual annotations |
| `pdb_ids` | 1CRN, 4HHB | Explicit structures |
| `uniprot_ids` | P00533 | Associated reference protein |
| `ligand_ids` | ATP, ADP | Structures containing those chemical components |
| `pdb_ec` | 3.1.1.7 | Enzyme classification, when known |
| Organism, polymer, method, resolution | Homo sapiens, Protein, 2.5 Å | Attribute search/refinement |

Different criteria use AND; values within a list are alternatives (OR). Clear fields to widen the search. RCSB free text may match alternative words; quotes request a phrase. See the [official RCSB search specification](https://search.rcsb.org/). For legacy calls without other criteria, a nondefault collection name also serves as search text.

`must_have_ligand` remains true by default. Disable it to retrieve structures without ligands; these require an appropriate ligand reference for the implemented redocking workflow. Detected nonwater components can include ions/additives, so inspect ligand identity.

The default limit is 100 structures. Search pagination and downloads have finite timeouts and transient retries. Up to four downloads run concurrently. Coordinate responses are parsed before atomic publication of `.pdb` files. LINK/SSBOND removal preserves the previous preparation contract.

Outputs include successful `<PDB>.pdb` files, `pdb_codes.csv` with `PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION`, and `retrieval_report.json` recording criteria/query, available total count, limits, selected IDs and individual download outcomes. Ligand-free structures generate no complex metadata rows. Individual failures retain successful downloads; no matches or no successful downloads produce explicit errors. mmCIF-only entries are reported as unavailable in PDB format; automatic conversion is not introduced into the PDB pipeline.

## ChEMBL search modes

Choose **Search ChEMBL by** first:

| `search_mode` | Query | Meaning |
| --- | --- | --- |
| `target` | Name, CHEMBL220 or P00533 | Automatic target detection, retaining the previous workflow |
| `target_name` | Partial name | Target preferred-name matching |
| `target_id` | CHEMBL220, CHEMBL240 | Explicit target IDs |
| `uniprot` | P00533, P22303 | UniProt-associated targets |
| `target_text` | Free text | Target search endpoint |
| `molecule_id` | CHEMBL25, CHEMBL50 | Explicit compounds, without requiring activity evidence |
| `molecule_name` | aspirin | Partial preferred-name compound matching |
| `similarity` | Compound ID or SMILES | Structurally similar compounds |
| `substructure` | c1ccccc1 | Compounds containing a SMILES fragment |

Target and compound ChEMBL IDs identify different resources. Automatic UniProt detection now validates accession syntax. The modes follow the [official ChEMBL API resources](https://www.ebi.ac.uk/chembl/api/data/docs).

Target searches run targets → activities → molecules, with optional ChEMBL expansion. Direct compound searches bypass activities. Structural similarity/substructure does not establish activity against a target; the report records `activity_evidence: false` for direct searches. The interface shows controls/filters relevant to the mode.

`expand_chembl` toggles similarity expansion after target retrieval. `include_pubchem` independently enables PubChem expansion, including for direct compound searches. Direct `similarity` defaults to 70%; target-associated ChEMBL expansion retains its own filter threshold, initially 65%.

## Filters, limits and existing configurations

Updated defaults no longer silently require human targets or natural/non-natural products. Organism/type can be empty and natural product supports **Any**. Other optional filter fields can be cleared. Saved overrides retain their values; configurations without overrides use the updated resources displayed for review.

Activity defaults remain Ki/IC50, nM, assay B, nonnull pChEMBL and maximum value 5000. Clear the maximum for unrestricted values; also clear units to search without a fixed unit. A numeric threshold requires explicit standard units. Missing numeric values may remain when no numeric threshold applies. Exported relations (`=`, `<`, `>`) must still be considered when interpreting activities. Molecule type belongs in the molecule filter group rather than target/activity filters.

Backend `chembl_filters` keys `target`, `bioactivity`, `molecules` and `similars` replace their respective default groups. `{}` clears that group's restrictions. The guided interface exposes the available filters without requiring JSON editing.

Initial limits are 25 targets, 1000 activities per target and 1000 records per direct search or ChEMBL expansion reference. These bound returned records before some local filtering; a smaller final dataset is not proof of an exhaustive search. Raise limits or refine criteria as needed. Multiple expansion references can generate many requests; both expansion providers can be disabled.

Targets/activities are queried afresh; molecule-record and PubChem caches remain. Application jobs have isolated output folders. When calling wrappers directly, use a fresh folder for each query/filter combination because older molecular exports can remain in reused directories.

The consolidated `compounds/<collection>/compounds.csv` retains `molecule_chembl_id,canonical_smiles,molecule_properties,source`, with existing structure validation/deduplication and compatibility with ADMET, fingerprints and graphs. Structural query syntax uses a deterministic `consulta_<hash>` folder, avoiding path separators in SMILES; the original query remains in `retrieval_report.json`.

## PubChem, ZINC and local files

PubChem expansion accepts downloaded ChEMBL collections or a selected curated compound table. Thresholds/limits per reference, cache, rate control, retries and relation exports remain available. Additional structures pass through consolidated validation/deduplication; similarity does not establish target activity.

ZINC retrieval still accepts a URL list. Its table parser now accepts repeated spaces/tabs and verifies exactly two SMILES/identifier columns, instead of silently misaligning fields. Retries and one-at-a-time URL processing remain. Local tables can be imported without mandatory retrieval.

## Backend examples

Use parameter JSON with `biomolexplorer <operation> --parameters file.json --output folder`.

For `retrieve_structures`, no collection name or EC is needed:

```json
{
  "pdb_ids": ["1CRN", "4HHB"],
  "must_have_ligand": false,
  "max_records": 10
}
```

For compounds without a target:

```json
{
  "search_term": "CHEMBL25, CHEMBL50",
  "search_mode": "molecule_id",
  "max_records": 100,
  "include_pubchem": false
}
```

## Verification

Offline regression tests cover criteria validation, modes/IDs, bounded pagination, timeouts, structural URL encoding, compatible output, partial failures and guided controls. Small live queries confirmed RCSB ID/UniProt/ligand attributes and ChEMBL preferred-name, free-text, UniProt, similarity and substructure endpoints. Updated wrappers downloaded/parsed 1CRN and produced consolidated compound tables by ID, name, similarity, substructure and target.

These checks establish contracts and representative behavior. Provider availability/coverage varies, and bounded searches are not exhaustive database reviews.

```bash
PYTHONPATH=src:tests python -m unittest test_flexible_retrieval test_chembl_rest test_compound_retrieval test_guided_ui
```

## Review structures and ligands in the interface

Results list the PDB structures and hide `pdb_codes.csv` and `retrieval_report.json`; both files remain stored and retain their pipeline contracts. Each structure provides three actions:

- **Ligands** lists detected records. Remove unwanted entries or add a ligand code, residue number and chain present in the downloaded PDB. Click **Save ligands** to confirm.
- **View 3D structure** opens the WebGL explorer in a browser tab, with mouse rotation, wheel zoom and right-drag panning. It offers ribbons, sticks, spheres and lines, chain/element/sequence coloring, chain/model selection, ligand/water toggles, fullscreen and PNG export. The complete file remains available for download.
- **Open in RCSB PDB** opens the structure’s official page by PDB identifier.

These actions are also available when choosing PDBs before the next stage. Viewers can inspect records; editors/owners can change them when calculations are inactive or paused for file selection. Changes enter project history, preserve other structures and update artifact integrity. Changed inputs invalidate reuse of affected consumers. Concurrent edits require reopening the list. Removing every ligand leaves the structure without a complex reference; add a valid reference or exclude that PDB before preparation/redocking. Running retrieval again can detect candidates afresh.


## ChEMBL selectors

**Organism** suggests scientific names, including `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`, `Escherichia coli`, `Saccharomyces cerevisiae` and `Danio rerio`. Search and free typing support other names; organisms are not a closed vocabulary. Leave it blank or choose **Any** to remove this filter. The example and help text identify it as the target organism.

**Assay type** displays official codes with descriptions: B — Binding, F — Functional, A — ADMET, T — Toxicity, P — Physicochemical and U — Unassigned, plus **Any**. API codes remain unchanged across interface languages. [Official ChEMBL vocabulary](https://chembl.gitbook.io/chembl-data-deposition-guide/file-structure/field-names-and-data-types-minimal-data-submission/assay.tsv).

**Activity measures** contains 6,433 distinct `standard_type` values found in a public ChEMBL query on October 5, 2026. This snapshot ships with the application; opening the form does not require another online lookup. Options use up to three columns, 24 items per page and name search. Selection survives filtering and pagination; **Clear selection** accepts any measure. **Other measure** accepts exact names of newer or specialized types missing from this snapshot. Scientific values, including case and spaces, are preserved. [Public catalog source](https://www.ebi.ac.uk/chembl/elk/es/chembl_activity/_search).

Activity vocabulary grows with the database. The snapshot is recorded in `resources/crawlers/activity_types.json`, rather than an immutable enum. Inhibition, percentage activity and kinetic parameters may require removing the nM filter, adjusting numeric limits and permitting missing pChEMBL values. Choosing a measure does not automatically change other filters. API and pipeline record limits continue to apply.

Names containing literal commas, such as `K(p,uu,brain)`, are queried using separate exact filters while preserving the total activity limit per target.

The ligand dialog also offers **Ligands present in the structure**, suggesting nonwater residues from the file and filling code, residue and chain when added. Manual entry remains available.

## Interactive PDB viewer

Click **View 3D structure** beside a PDB in results or file selection. The explorer opens in a browser tab on both desktop and web. This also supports Linux, where Flet’s native WebView is unavailable.

- Left-drag to rotate; use the mouse wheel to zoom.
- Right-drag to pan. **Recenter** fits the camera to the selection.
- Choose ribbons with ligands, sticks, spheres or lines; color by chain, element or sequence.
- Select a chain or another model in multimodel PDBs. The first model is shown initially.
- Toggle ligands and water. Click an atom to inspect its residue, chain, atom name and element.
- Use **Fullscreen** and **Save image** to export the current view as PNG.

[3Dmol.js 2.5.5](https://github.com/3dmol/3Dmol.js/releases/tag/2.5.5) ships with the application, including license and integrity hashes. Viewing uses browser WebGL and perspective, without CDN requests, RCSB downloads or uploading structures to outside providers. A WebGL-capable browser is required; initialization failures show guidance.

Web viewer routes share the Flet origin and port. Desktop uses an application-managed loopback-only server on `127.0.0.1`, stopped at shutdown. Opaque temporary links recheck session and project permissions when serving the page and PDB. Logout, access revocation and expiration deny subsequent requests; already loaded tabs retain their local copy, as with a download. Pages and structures use `Cache-Control: no-store`; project folders are not published as assets. The opening limit remains 32 MB. Viewing does not change coordinates, ligand curation or pipeline results.

The bundled JavaScript library requires no separate application or Python package. Web hosting uses FastAPI/Uvicorn already supplied with the Flet web `ui` extra.

### Verify the viewer

Authorization and integration tests are in `tests/test_pdb_viewer.py`. To verify WebGL, gestures and export with an installed Chrome/Chromium, run in the UI environment:

```bash
PYTHONPATH=src python scripts/validate_pdb_viewer.py --chrome /usr/bin/google-chrome
```

The script uses only a temporary loopback server and browser profile, verifies rotation, zoom, representations and PNG export and rejects external service requests. Its default screenshot is `/tmp/biomol-pdb-modern.png`; use `--screenshot` for another destination.
