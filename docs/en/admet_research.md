# Hypothesis-driven ADMET exploration

[Documentation](../README.md) · English · [Português](../admet_research.md)

Implementation review: October 10, 2026. This stage supports reproducible compound exploration. Physicochemical rules, structural alerts and model predictions do not establish safety or efficacy.

## Earlier behavior and compatibility

The earlier module calculated six RDKit descriptors, classified BBB/HIA with local ellipses and penalties, and estimated P-gp from a feature count. PAINS and small custom SMARTS sets triggered automatic exclusions described as toxicity. It did not provide experimental metabolism, excretion or toxicity models, research profiles, or a complete explanation of each rejected molecule.

The current implementation separates calculated descriptors, QSPR estimates, selection rules, structural alerts and optional trained models. Existing results are not automatically recalculated. Configurations without a profile retain `legacy` behavior; new blocks start with broad exploration and descriptive output.

## Methods and interpretation

| Method | Use | Implementation and limits |
| --- | --- | --- |
| Wildman–Crippen and TPSA | Lipophilicity and polarity | RDKit calculates WLOGP/TPSA from the submitted structure. WLOGP is not logD at a specified pH. |
| Lipinski | Small molecules under an oral hypothesis | Counts MW > 500, WLOGP > 5, HBD > 5 and HBA > 10. These are not universal restrictions for natural products, macrocycles or other routes. |
| Veber | Flexibility and polarity under an oral hypothesis | Counts RB > 10 and TPSA > 140 Å². The original evidence concerns rat bioavailability, not measured human absorption. [Veber et al., 2002](https://pubmed.ncbi.nlm.nih.gov/12036371/). |
| Extended rule of three | Fragment libraries | MW ≤ 300, WLOGP ≤ 3, HBD/HBA ≤ 3, plus RB ≤ 3 and TPSA ≤ 60. The extension is explicit and custom ranges can replace it. [Congreve et al., 2003](https://pubmed.ncbi.nlm.nih.gov/14554012/). |
| QED | Continuous characterization of drug-like chemical space | Range 0–1, without a mandatory threshold. RDKit QED can differ from the publication's implementation. It does not measure potency, exposure or toxicity. [Bickerton et al., 2012](https://www.nature.com/articles/nchem.1243), [RDKit documentation](https://www.rdkit.org/docs/source/rdkit.Chem.QED.html). |
| RDKit-adapted ESOL | Estimated aqueous solubility | `ESOL_logS`, log10 mol/L, uses the RDKit refit documented by Pat Walters: MW, MolLogP, RB and aromatic proportion. It is neither the original CLOGP equation nor SwissADME ESOL. Crystal form, pH and formulation are not modeled. [Implementation and Delaney reference](https://github.com/PatWalters/solubility). |
| Adapted Potts–Guy | Skin-permeation exploration | `logKp = −6.3 + 0.71 × WLOGP − 0.0061 × MW`, log10 cm/s, explicitly substitutes RDKit logP. Vehicle, skin condition and ionization are not modeled. [Potts & Guy, 1992](https://pubmed.ncbi.nlm.nih.gov/1608900/). |
| PAINS/BRENK/NIH | Review of interference, reactivity and structural motifs | Matching pattern names appear in the CSV. Annotation is the default; exclusion is an explicit choice. A match does not prove toxicity. [RDKit catalogs and references](https://www.rdkit.org/docs/GettingStartedInPython.html#filtering-unwanted-substructures). |

Additional fields include heavy atoms, rings, aromatic rings, fraction of sp³ carbons, formal charge, disconnected fragments and structure warnings. Salts, mixtures and ionized structures do not automatically represent the predominant species in an assay.

### BBB, HIA and the EGG

The published BOILED-Egg describes passive intestinal and brain permeation in WLOGP–TPSA space. It does not establish unbound brain exposure, active transport or complete oral bioavailability. [Daina & Zoete, 2016](https://pubmed.ncbi.nlm.nih.gov/27218427/).

BioMolExplorer retains its historical local geometry: HIA center (75, 2), radii (75, 3); BBB center (42, 2.3), radii (47, 2.2). The BBB classifier adds penalties, so a point inside the yolk can still be BBB−. The legend distinguishes regions from classifications. These formulas have not been validated here as a reproduction of the published model.

SwissADME colors indicate predicted P-gp substrate status. Our orange circles and blue squares indicate local BBB+/BBB−. New analyses do not infer P-gp from descriptors: `PGP=not_evaluated`. Legacy mode retains its identified historical heuristic. [SwissADME, Daina et al., 2017](https://www.nature.com/articles/srep42717).

## Libraries evaluated

| Resource | Role | Integration decision |
| --- | --- | --- |
| RDKit | Deterministic local descriptors and interpretable catalogs | Required foundation; not a universal ADMET predictor. |
| ADMET-AI | Local Chemprop models trained on TDC datasets | Optional adapter for version 2.x, pinned to 2.0.1; no remote inference or DrugBank percentiles. [Authors' repository](https://github.com/swansonk14/admet_ai), [2.0.1 distribution](https://pypi.org/project/admet-ai/2.0.1/). |
| Therapeutics Data Commons | Endpoint definitions, datasets and evaluation | Reference for endpoint meaning; installing TDC does not supply a validated model by itself. [TDC overview](https://tdcommons.ai/single_pred_tasks/overview/). |
| DeepChem | Training and evaluating task-specific models | Future integration requires datasets, splits, architecture and external validation. [Official documentation](https://deepchem.readthedocs.io/en/stable/index.html). |
| OpenADMET | Developing property models using experimental datasets | Candidate for endpoint-specific integration and validation, not a current backend option. [Project](https://www.openadmet.org/), [inference example](https://demos.openadmet.org/en/latest/demos/05_Ensemble_Model_Inference/05_Model_Ensemble_Inference.html). |
| SwissADME | Scientific reference and web tool | Used for interpretation; no scraping dependency or numerical-equivalence claim. |

ADMET-AI v1 and v2 were trained separately and can return different values. Reports record the version. Endpoint names, task types, species and units come from installed model metadata. There is no universal 0.5 cutoff, and neural outputs do not replace local BBB classes.

### Optional trained-model coverage

The adapter adds available outputs with the `ML_` prefix. Depending on the installed version, they cover absorption/permeability and solubility; distribution, protein binding and distribution volume; CYP metabolism; half-life and clearance; and endpoints such as hERG, AMES and DILI. Availability does not establish reliability for your chemical series. Check each endpoint's scale and meaning. [TDC ADME](https://tdcommons.ai/single_pred_tasks/adme/), [TDC toxicity](https://tdcommons.ai/single_pred_tasks/tox/).

`ML_Pgp_Broccatelli` concerns **P-gp inhibition**, not substrate classification, and does not color the EGG or imply efflux. `ML_BBB_Martins` remains separate from local `BBB`. Custom ranges can select any returned `ML_` output; a missing value fails a required range.

Install the optional predictor separately from the main scientific environment:

```bash
python scripts/setup_admet_ai.py
```

The script creates `.admet-env` using the Python that invokes it, installs CPU PyTorch and ADMET-AI 2.0.1, checks an inference and only then writes `.ready`. Use a compatible Python; the verified environment uses Python 3.12. Workers discover the repository-local environment automatically. For a different directory or deployment, set `BIOMOL_ADMET_PYTHON` to its Python executable. Alternatively install `.[admet-ml]` in a dedicated worker. The configuration panel shows availability before execution. See [installation](installation.md#optional-admet-ai-predictor) for troubleshooting.

The adapter predicts unique SMILES in bounded batches, checks identity/order and finite values, and reuses predictions for duplicate structures. Real inference for ethanol and acetic acid returned 41 endpoints in the tested isolated environment. This is functional verification; applicability domain and uncertainty remain unevaluated for a research library.

## Configure a research hypothesis

| Profile | Starting criteria | Interpretation |
| --- | --- | --- |
| Broad exploration | All valid compounds, annotated alerts | Preserve diversity while studying the series. |
| Oral candidates | HIA+, at most one Lipinski violation, no Veber violations | Passive absorption does not include first-pass metabolism. |
| CNS access | Local BBB+ | Investigate efflux, pKa/logD, binding and unbound brain exposure. |
| Peripheral restriction | Local BBB− | Does not prove absence of brain exposure. |
| Fragments | Extended rule of three, no violations | Adapt to fragment binding and growth strategy. |
| Natural products / beyond rule of five | No mandatory oral cutoff | Macrocycles, transport and alternative routes may require specific models. |
| Custom | Freely combined criteria | Record rationale and units for each range. |
| Legacy | Historical pattern exclusions | Compatibility only; not validated toxicity assessment. |

Profiles supply editable starting values. BBB and HIA accept unrestricted, positive or negative. Rules, violation allowances, catalogs and alert handling are configurable. Ranges are inclusive; criteria combine with **AND**. An empty minimum or maximum leaves that side unrestricted.

**Annotate all compounds** forwards every evaluated compound and records hypothesis selection. **Filter by hypothesis** filters the main CSV. Both preserve `_evaluated.csv`, `ADMET_selected` and `ADMET_reasons`; new analyses also export `_selected.csv`, `_rejected.csv`, BBB± and HIA± subsets. `_admet_report.json` records configuration, versions, methods, counts and limits. Invalid molecular records remain in the pipeline's quality report.

For example, peripheral oral exploration can combine BBB−, HIA+, oral rules and `ESOL_logS ≥ −4`. That threshold is a research choice. Natural-product exploration can retain all compounds and annotate QED, solubility and alerts without automatic exclusion.

## Visualization and interpretation

The interactive EGG measures its panel and redraws axes, regions and points at the available width. The chart uses a 560 px panel; details and downloads move below it on narrower screens. Coordinates are preserved without jitter, and limits include outliers. Hover identifies compound, TPSA, WLOGP and classes; clicking opens structure/properties. The selector provides access to coincident points.

Orange circles and blue squares reduce dependence on color perception. Subtle grids, units and explicit legends help interpretation. PNG export uses the same geometry and encoding. These interface decisions do not constitute an experimental usability evaluation. [Crameri et al., 2020](https://www.nature.com/articles/s41467-020-19160-7), [Nature Methods, Points of View, 2026](https://www.nature.com/articles/s41592-026-03143-5).

[View an exported EGG example](../../imgs/examples/admet-egg.png) and its [input structures and calculated properties](../../imgs/examples/admet-example.json). This is illustrative local screening, not a validated SwissADME prediction.

## Current limits and scientific development

Rules, alerts and neural outputs are not experimental measurements. Validated compound-level applicability domains, prediction intervals, pKa/logD, pH-dependent species permeation, experimental aggregation, unbound brain exposure and assay-specific data integration are not implemented. Reports identify unevaluated uncertainty instead of supplying fictitious confidence.

For a specific hypothesis, gather data relevant to species, tissue, route and assay; use splits that prevent series leakage; compare models with baselines; evaluate calibration, extrapolation and external performance; and confirm selections in suitable experiments.

## Configuration layout and result management

Fields are grouped by objective/output and selection criteria, using two columns when space permits and one on narrow screens. Alerts and custom ranges are collapsible; longer explanations appear in contextual help. This applies proximity, consistent spacing and progressive disclosure. [Form grouping](https://www.nngroup.com/articles/form-design-white-space/), [progressive disclosure](https://www.nngroup.com/articles/progressive-disclosure/).

Execution honors explicit file selectors and saved filenames. Input curation is requested when a connected input still uses automatic selection. **Remove all results** requires confirmation, includes auxiliary files and invalidates reused references; run records and logs remain, and history can restore files. Active runs block removal. See the [workspace guide](frontend.md#navigation-canvas-and-completed-runs).

Schema 4 exports declare each output in the report: BBB± and HIA± independently partition all evaluated compounds; `rules-<rule>.csv` applies only the chosen rules; `alerts.csv` contains Boolean columns per catalog. The main table retains the global hypothesis selection. Legacy batch-qualified selections migrate only when the producer and output role match unambiguously.
