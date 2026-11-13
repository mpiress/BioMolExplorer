<p align="center"><img src="imgs/new_BioMolLogo.png" alt="BioMolExplorer" width="480" style="max-width:100%;height:auto"/></p>

# BioMolExplorer

BioMolExplorer is research software for molecular retrieval, chemical analysis and structure-based screening. A Flet workspace connects Python scientific services in a configurable directed acyclic graph (DAG); projects retain inputs, parameters, results and provenance.

[Platform overview](index.html) · [English manual](docs/en/user_manual.md) · [Manual em português](docs/user_manual.md) · [GitHub](https://github.com/mpiress/BioMolExplorer)

## Research workflows

| Area | Available capabilities |
| --- | --- |
| Data and inputs | ChEMBL retrieval and expansion, PubChem reference searches, PDB queries, ZINC tranche downloads, validated local files and supported prepared results |
| Molecular analysis | RDKit descriptors, BBB/HIA heuristics, physicochemical rules, chemical review alerts; Morgan/MACCS/pharmacophore fingerprints, similarity, graphs, components and scaffolds |
| Structure-based screening | Shared Chimera receptor/ligand preparation; redocking with Vina or DOCK6 6.11 and RMSD; candidate docking, compound-level multi-target ECR consensus and available DOCK6 footprints |
| Native Amber Score | Selected receptor–pose pairs, implicit-solvent refinement/energy evaluation through DOCK6 6.11; no conventional MD or trajectory-based RMSF |
| Results and workspace | Tables, plots, molecular 2D/3D, typed connections, lazy input selection, cooperative pause/cancellation, verified incremental reuse, history, sharing and portable backups |

The workspace provides **16 block types**: 15 scientific operations plus import. The CLI exposes the 15 scientific operations. A typical branch is `PubChem → fingerprints → similarity → graphs → ADMET`; a structural branch is `PDB → redocking → prepared receptors + candidates → Vina/DOCK6 → consensus`. Native Amber Score can consume selected receptors and three-dimensional ligand poses. These are configurable examples, not validated universal protocols.

## Environment and first run

Use Python 3.12 and the scientific dependencies declared in [requirements.yml](requirements.yml). Active preparation protocols require classic UCSF Chimera 1.17 and Open Babel; Vina protocols need AutoDock Vina, and DOCK6 protocols need a complete DOCK6 6.11 installation with the relevant accessory executables and parameter files. The native SES surface implementation does not require a standalone DMS executable.

Follow the [step-by-step installation guide](docs/en/installation.md). From the repository root:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
python -m pip check
biomolexplorer-ui --help
biomolexplorer-ui --web --language en
```

For a local web session with explicit worker/tool paths:

```bash
PYTHONPATH=src python -m biomolexplorer.ui.app --web --host 127.0.0.1 --port 8550 \
  --worker-python /path/to/BioMolExplorer/bin/python --dock6-path /path/to/dock6
```

These paths are placeholders. Configure the actual installations; the Python package does not install Chimera or DOCK6. Optional trained ADMET models have separate environment requirements described in the installation guide.

## Scientific CLI example

Create `pubchem.json` with the following reference query, then execute it in a separate output directory:

```json
{"target":"Aspirin","reference_source":"manual","reference_type":"cid","reference":"2244","threshold":90,"max_records":10}
```

```bash
PYTHONPATH=src python -m biomolexplorer retrieve_pubchem --parameters pubchem.json --output results/aspirin
```

This example requires network access and is not executed by the documentation build. Consult the [API parameter reference](docs/en/backend_usage.md#exact-active-parameter-reference) for required/optional keys and the manual for a short example of every block.

## Persistence and architecture

Application services live in `src/biomolexplorer/`; scientific entry points are under `src/wrappers/`, with implementations in `src/caad/`, `src/crawlers/` and `src/kernel/`. Native commands run through worker services. The DAG enforces compatible dependencies and saves execution priority.

Reuse verifies inputs, settings, processing/tool versions and output integrity; it is not based on file existence alone. Destructive edits require confirmation. **Remove All** removes computed results while preserving configuration; **Clear All** also resets selections/scientific defaults while retaining blocks and connections. Public artifact aliases use `TG[TARGET]-SRC[SOURCE]-OP[OPERATION]-TYPE-ID[hash].ext`; native paths and complete lineage remain in metadata. See [project contracts](docs/en/projects.md#execution-lazy-selection-and-result-preservation).

## Documentation

- [Browser portal](docs/index.html) · [Portal em português](docs/pt/index.html) · [Source index](docs/README.md)
- [User manual](docs/en/user_manual.md) · [Installation](docs/en/installation.md) · [Backend and parameters](docs/en/backend_usage.md)
- [Architecture](docs/en/architecture.md) · [Project history and migration](docs/en/projects.md) · [Preparation](docs/en/molecular_preparation.md)
- [Documentation review and validation limits](docs/validation/documentation_review_2026-10-31/README.md)

Markdown is the source for the static guides. Existing documentation tools use the Python standard library:

```bash
python scripts/build_docs.py
python scripts/validate_docs.py
```

The portal uses its existing template and registry. This audit changes documentation only; source scripts and their registry remain untouched. Documentation checks do not certify scientific correctness, environment solvability or release readiness. Docking scores are not experimental affinities; local ADMET heuristics and alerts require hypothesis-specific interpretation. Native Amber Score is not a conventional molecular dynamics engine.

## License and related research

Application code is distributed under the [MIT license](LICENSE); bundled Forcefield_PTM parameter data have separate **CC BY-NC 4.0** terms recorded in LICENSE and their NOTICE/manifest. The following is an associated Alzheimer’s disease dataset, **not a DOI assigned to the software release**. An earlier documentation audit recorded HTTP 404 for the supplied DOI resolver; this source-code audit did not revalidate that external record. Confirm bibliographic details before citation or reuse. The identifier is retained as text rather than a broken navigation link.

**Project-supplied reference (external record not revalidated by this source-code audit):**

Pires da Silva, M., Alves de Oliveira, T., Habib Bechelane Maia, E., Oliveira Mendes, G., Cristina Moreira Damázio, L., Brito Barbosa, D., Andrade Leite, F. H., Falkoski, L., Flores de Souza Marra, I., Siqueira Valle, M., Silva Matos Andrade, L., Čmelo, I., Fayne, D., Batista de Carvalho, P., Marques da Silva, A., & Gutterres Taranto, A. (2025). *Data-driven chemical domain for polypharmacology agents: Focus on Alzheimer’s disease*. Mendeley Data. `10.17632/5njg46dfj4.3`


## Authors

| [<img loading="lazy" src="imgs/michel.jpg" width=150><br><sub> Michel Pires da Silva</sub>](http://lattes.cnpq.br/1449902596670082) |  [<img loading="lazy" src="imgs/alisson.png" width=150><br><sub> Alisson Marques da Silva</sub>](http://lattes.cnpq.br/3856358583630209) |  [<img loading="lazy" src="imgs/alex.png" width=150><br><sub> Alex Gutterres Taranto</sub>](http://lattes.cnpq.br/4759006674013596) |
| :---: | :---: | :---: |

## Current implementation reference

[Source-derived command/parameter inventory](docs/en/current_reference.md) · [Compound-level consensus](docs/en/docking_consensus.md) · [Native AMBER structural support](docs/en/amber_score_structural_support.md) · [Documentation audit](docs/documentation_audit_2026-11-13.md).

Native Amber Score includes coherent PTM profiles and a cached execution-local compatibility build when needed; this does not alter the licensed DOCK6 installation. Unknown covalent chemistry still requires validated parameters. The update button installs a validated, commit-pinned package into the host Python environment; restart remains required, and scientific dependencies are not upgraded automatically.
