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
No block uses a standalone AmberTools environment. Native Amber Score uses the utilities bundled with DOCK6 6.11.

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

## Concurrency and resource use

Vina, DOCK6 and Prepare for docking share candidate preparation. Independent molecules use a bounded queue and private temporary Chimera directories. Charge methods, pH, minimization and chemical validation are unchanged. When both formats are requested, they are exported from the same prepared conformation without repeating charge calculation.

Admission respects the worker CPU budget (normally two in the workspace), half the effective CPUs, available-memory estimates and a ceiling of four. CPU affinity/quotas and cgroup-v2 limits are considered. The estimate reserves at least 1 GiB and, when larger, 25% of available memory, with a 1 GiB allowance per process. Unknown or small available memory leads to serial preparation. This is a conservative admission policy, not a hard memory limit or a guarantee against OOM for every molecule.

Linux file leases coordinate candidate-preparation Chimera processes across workers of the same OS user. Process-group cancellation releases leases automatically. Numerical-library threads inside Chimera are limited to one. Manually started Chimera sessions and other OS users are outside this coordination scope.

Set `BIOMOL_PREPARATION_WORKERS=1` in the application's environment to reduce concurrency without changing scientific settings. Higher values never exceed the worker budget or safety limits. The workspace supplies `BIOMOL_CPU_WORKERS` from its configuration; this variable sets the reference budget for CLI use. No scientific numeric option was added to forms.

Per-molecule native files and verified cache manifests remain durable. The growing `compounds.csv` index is coalesced with a five-second minimum interval between periodic writes and forced on completion/handled fatal error, reducing repeated I/O. This does not promise a CSV checkpoint for every molecule; valid rows can be rebuilt from verified cache after interruption. Explicitly disabled reuse remains disabled. Exact version equivalences were audited so scheduling changes do not invalidate existing scientific caches.

Speedup depends on molecules, hardware and cache state. Parallelism does not repair chemical typing or convergence failures. See the [benchmark scope and limitations](../validation/chimera_parallel_2026-11-02/README.md).
