<p align="center"><img src="imgs/new_BioMolLogo.png" alt="BioMolExplorer" width="480" style="max-width:100%;height:auto"/></p>

# BioMolExplorer

BioMolExplorer combines molecular retrieval, local analysis and docking in a Flet workspace. Projects hold pipeline configuration, uploaded files, results, logs and version history. Use the [platform overview](index.html), [English manual](docs/en/user_manual.md) or [manual em português](docs/user_manual.md).

## Available functionality

| Area | Active operations |
| --- | --- |
| Retrieval | ChEMBL compounds and similar expansion, PubChem reference searches, PDB structures, ZINC tranche downloads |
| Local inputs | Import compound tables, structures and supported prepared results |
| Molecular analysis | Morgan, MACCS and pharmacophore fingerprints; configurable similarity; molecular graphs, components and common fragments |
| ADMET | RDKit descriptors, BBB/HIA heuristics, physicochemical rules, structural alerts and optional ADMET-AI predictions |
| Preparation | Chimera receptor/ligand preparation, gas or am1 charges, Open Babel format conversion |
| Docking | Vina redocking with RMSD; independent Vina and DOCK6 docking; consensus of matching compound/receptor results |
| Inspection | Compound tables, ADMET plots, graph exploration, molecular 2D/3D, pose/receptor comparisons and DOCK6 footprints |
| Workspace | Accounts, project sharing, history/restore, execution monitoring, cancellation, compatible result reuse, downloads and project migration |

The [operation registry](src/biomolexplorer/operations.py) defines the scientific API; the [catalog](src/biomolexplorer/catalog.py) adds the import block and UI metadata. The manual includes a short example for every active block. Geometric contacts, local ADMET heuristics and docking scores require interpretation within a research protocol.

## Installation and startup

From the repository root, create the scientific environment and install the UI:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
python -m biomolexplorer.ui.app --language pt
```

Install UCSF Chimera 1.17, Open Babel and Vina for the protocols that use them. DOCK6 docking also requires DOCK6 6.11, its accessory executables and parameter files. See [installation](docs/en/installation.md) ([Português](docs/installation.md)) for executable checks. Installing the Python package does not install Chimera or DOCK6. Molecular surfaces use the local Python SES implementation; a standalone DMS executable is unnecessary.

Web mode and a configured DOCK6 installation root:

```bash
python -m biomolexplorer.ui.app --web --host 127.0.0.1 --port 8550 \
  --worker-python /path/to/BioMolExplorer/bin/python \
  --dock6-path /path/to/dock6
```

Optional trained ADMET models use an isolated environment created with `python scripts/setup_admet_ai.py`; local RDKit analysis remains available without it.

## Scientific CLI

The same registry is exposed through `python -m biomolexplorer`. Supply a JSON parameter file and a separate output directory:

```bash
cat > pubchem.json <<'JSON'
{"target":"Aspirin","reference_source":"manual","reference_type":"cid","reference":"2244","threshold":90,"max_records":10}
JSON
python -m biomolexplorer retrieve_pubchem --parameters pubchem.json --output results/aspirin
```

This example queries PubChem and requires network access. Consult [backend usage](docs/en/backend_usage.md) for other operations and their input contracts. Legacy scripts under `workflow/` remain in the repository; the registry and workspace are the supported starting points described here.

## Source layout

```text
src/biomolexplorer/  Application services, pipeline, results and Flet UI
src/wrappers/        Scientific operation entry points
src/caad/            Docking and molecular analysis implementations
src/crawlers/        Retrieval implementations
src/kernel/          Shared scientific utilities
workflow/            Standalone workflow scripts
scripts/             Setup, documentation and validation tools
docs/                Bilingual Markdown guides and generated HTML
tests/               Software tests and scientific reference fixtures
```

## Documentation and validation

- [User manual](docs/en/user_manual.md) · [Manual do usuário](docs/user_manual.md)
- [Documentation portal](docs/index.html) · [Portal em português](docs/pt/index.html)
- [Architecture](docs/en/architecture.md) · [Pipeline contracts](docs/en/pipeline_validation.md)
- [Preparation with Chimera](docs/en/chimera_migration.md) · [Native molecular surfaces](docs/en/dms_migration.md)

Update both Markdown languages before regenerating HTML. The documentation tools use the Python standard library:

```bash
python scripts/build_docs.py
python scripts/validate_docs.py
```

Scientific execution and reference-case validation have separate requirements from documentation checks. See the scope stated in each validation report.

## License

[MIT](LICENSE).

# 📄 How to Cite Our Work

If you use the data from this study in your research, please cite the dataset as follows:

**Full Reference:**

Pires da Silva, M., Alves de Oliveira, T., Habib Bechelane Maia, E., Oliveira Mendes, G., Cristina Moreira Damázio, L., Brito Barbosa, D., Andrade Leite, F. H., Falkoski, L., Flores de Souza Marra, I., Siqueira Valle, M., Silva Matos Andrade, L., Čmelo, I., Fayne, D., Batista de Carvalho, P., Marques da Silva, A., & Gutterres Taranto, A. (2025). *Data-driven chemical domain for polypharmacology agents: Focus on Alzheimer’s disease*. Mendeley Data. https://doi.org/10.17632/5njg46dfj4.3

**Example of in-text citation:**

> (Pires da Silva et al., 2025)

<p> </p>
<p> </p>

## Authors

| [<img loading="lazy" src="imgs/michel.jpg" width=150><br><sub> Michel Pires da Silva</sub>](http://lattes.cnpq.br/1449902596670082) |  [<img loading="lazy" src="imgs/alisson.png" width=150><br><sub> Alisson Marques da Silva</sub>](http://lattes.cnpq.br/3856358583630209) |  [<img loading="lazy" src="imgs/alex.png" width=150><br><sub> Alex Gutterres Taranto</sub>](http://lattes.cnpq.br/4759006674013596) |
| :---: | :---: | :---: |
