# Visual examples / Exemplos visuais

These assets illustrate actual BioMolExplorer outputs. They contain tutorial or test data, not private research-project results. They are presentation examples rather than experimental evidence of affinity, permeability or safety.

As imagens ilustram saídas do BioMolExplorer com dados de tutorial ou teste, sem resultados privados de projetos. Não constituem evidência experimental de afinidade, permeabilidade ou segurança.

| Asset | Origin and interpretation / Origem e interpretação |
| --- | --- |
| `workspace-pipeline.png` | Actual 1700 × 1100 Chrome capture of the current Flet workspace, using an isolated demonstration project. PubChem → fingerprints → similarity → graphs → ADMET connections are configured; no calculations were run. / Captura real do workspace com projeto temporário, sem cálculos científicos. |
| `workspace-configuration.png` | Actual guided configuration dialog for the same project’s fingerprint block. / Formulário real de configuração do bloco de fingerprints do mesmo projeto. |
| `admet-egg.png` | Export from `biomolexplorer.egg_plot.egg_png`, using the current RDKit descriptor/local BBB classifier on the 11 structures in `admet-example.json`. Orange circles = BBB+, blue squares = BBB−; regions use the application's local heuristic. / Exportação do módulo atual, com propriedades e classes locais dos 11 compostos do JSON. |
| `interactions-3d.png` | Actual WebGL screenshot from `scripts/validate_docking_viewer.py`, using the benzene/PHE geometric test scene. The displayed score belongs to a test fixture; it is not a native docking result. Colored traces demonstrate overlapping π–π, hydrophobic and geometric van der Waals contacts. / Captura real do visualizador com cena geométrica de teste; o score do fixture não representa cálculo de docking. |
| `dock6-footprint.pdf` | Native plotting output after real Vina → DOCK6 execution on the tutorial 1ABE receptor/reference ligand (`REF`). Reference/pose van der Waals and electrostatic contributions are plotted separately. / Saída do fluxo real no receptor e ligante do tutorial 1ABE. |
| `dock6-footprint.png` | Raster export of the original PDF using `pdftoppm -png -singlefile -scale-to 1800`. The original curves and values are preserved. / Exportação raster do PDF original, preservando curvas e valores. |

The graph example `../Fig1.png` remains the previously published illustration.

## Reproduce / Reproduzir

Use the scientific environment described in the [installation guide](../../docs/en/installation.md). The EGG data include submitted SMILES, RDKit descriptors and local BBB/HIA classes. Recalculate them with `caad.admet.MoleculeEvaluator.calculate_properties`, `classify_bbb` and `classify_hia`; build plot nodes with `x=TPSA`, `y=WLOGP`, `properties=<calculated properties>` and call `egg_png` with a title and nodes.

For the native example, run:

```bash
PYTHONPATH=src python scripts/validate_docking.py --dock6-root /path/to/dock6 --output /tmp/biomol-validation
pdftoppm -png -singlefile -scale-to 1800 /tmp/biomol-validation/vina_to_dock6/Validation/1ABE_A/Dock6/footprint/plots/REF.pdf /tmp/dock6-footprint
```

Exact scores can depend on installed engine versions and settings. The [review evidence](../../docs/validation/dock6_review_2026-10-10.json) records the original functional run. The structural scene is included in `interaction-fixture/` (synthetic geometry and fixture score). Reproduce its viewer capture with:

```bash
PYTHONPATH=src python scripts/validate_docking_viewer.py --docking-root imgs/examples/interaction-fixture --screenshot /tmp/interactions-3d.png
```

This uses an installed Chrome/Chromium and a temporary local workspace; it does not require a native docking calculation.

DOCK6 tutorial source: [ligand sampling tutorial](https://dock.compbio.ucsf.edu/DOCK_6/tutorials/ligand_sampling_dock/ligand_sampling_dock.html). See the [ADMET guide](../../docs/en/admet_research.md) and [viewer guide](../../docs/en/molecular_viewer.md) for interpretation limits.

## Workspace captures / Capturas do workspace

Open an isolated example project in the current Flet application, in English.
Add Retrieve PubChem, Generate fingerprints, Calculate similarity, Filter by graphs
and Evaluate ADMET, then connect each output to the next compatible input.
Use PubChem CID 2244, a 90% retrieval cutoff and at most 100 records; Morgan
fingerprints with radius 2 and 2,048 bits; Tanimoto relationships at 70%; an
unlimited common-fragment search; and the oral ADMET profile with local descriptors.
Do not run the pipeline. Capture the canvas in a 1700 × 1100 browser viewport;
select Generate fingerprints and open Configure selection for the form capture.
The visible demo identity is synthetic and the application uses temporary storage.

Abra um projeto de demonstração isolado no Flet, conecte os cinco blocos e configure
os parâmetros acima, sem executar o pipeline. Capture o canvas e o formulário de
fingerprints. As imagens mostram controles reais, não resultados científicos.
