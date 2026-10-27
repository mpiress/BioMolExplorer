# Molecular preparation for docking and redocking

UCSF Chimera is the only preparation tool for Prepare for docking, Redocking,
Vina docking and DOCK6 docking. Historical `auto` and `python` configurations
are normalized to `legacy_chimera`, including saved projects and resumed runs.
AmberTools and Meeko are not selectable preparation backends for these blocks.
ChimeraX does not interpret the Chimera scripts.

One shared receptor section and one shared ligand section control solvent and
hydrogen removal, hydrogen addition, minimization and charges. Pair identities,
chains and cofactors remain associated with the selected structures.

- `gas`: Chimera Gasteiger charges.
- `am1`: Chimera AM1-BCC charges.
- Minimization uses Chimera. Standard protein residues retain Chimera's protein
  charge templates; gas/am1 controls nonstandard components and ligands.

Chimera may invoke utilities bundled with its own distribution for AM1-BCC.
The application does not route docking/redocking to the separate
BioMolExplorer-AmberTools environment. That environment remains available to
molecular-dynamics protocols requiring it.

ChEMBL, PubChem, ZINC and uploaded compounds follow the same preparation route.
RDKit may generate an initial SMILES conformer; Open Babel converts formats.
Chimera assigns final charges and performs minimization. MOL2/PDBQT exports
retain the chosen charge method. Vina and DOCK6 remain the docking engines.

Install Chimera, Open Babel and the selected docking engine in the worker's
executable search path. The worker interpreter's directory is included in PATH.
Explicitly supplied prepared receptors may be reused. Cache signatures include
chemistry options and tool identities; native historical preparations are not
silently relabeled as Chimera results. Pair failures remain isolated and logged.

Chimera runs structural Dock Prep before AddH/AddCharge, resolving alternate
locations for receptors and ligands and repairing incomplete receptor side
chains. Minimization then uses `prep false` to preserve the chosen charges.

Available local CCD definitions guide carbon hybridization and oxime/nitrone
atom types before protonation. Failed charge assignment may retrieve a
component definition and retry in Chimera with the same charge method.
The application does not substitute gas charges for AM1 failures.

MOL2 export may repair heavy-bond orders from the component definition while
preserving prepared coordinates and partial charges. Missing atoms are never
invented. An internal `.export.json` records missing atoms and the scope when
the reference represents only an observed fragment.
