# Molecular preparation with Chimera

> Current preparation guide; its historical URL has been preserved.

## Active path

Receptor and ligand preparation uses UCSF Chimera. Controls expose solvent and
hydrogen removal, hydrogen addition, `gas` or `am1` charges, minimization and
cofactors. Receptor and ligand settings are independent. Open Babel applies pH
during conversion; pH is not a Chimera template argument. This guide describes
the current code and provides no alternative backend setup.

## Entry points

| Operation | Implementation | Purpose |
| --- | --- | --- |
| Standalone preparation | `wrappers.redocking.prepare_structures` | Receptors and candidates for Vina, DOCK6 or both |
| Redocking | `wrappers.redocking.perform_redocking` | Explicit pairs, preparation and RMSD with Vina or DOCK6 6.11 |
| Vina docking | `wrappers.docking.perform_docking_vina` | Candidates and receptor, reusing compatible prepared structures |
| DOCK6 docking | `wrappers.docking.perform_docking_dock6` | Candidates, receptor, spheres, grids, docking and footprint |
| Consensus | `wrappers.docking.generate_consensus` | Combine scores for compound/receptor pairs present in both engines |

`biomolexplorer.redocking_config` validates pairs and applies template options.
`biomolexplorer.redocking_preparation` isolates pair failures and reuses receptors.
`biomolexplorer.docking_preparation` prepares candidates and validates receptor
formats. The docking backend is implemented in `src/caad/docking.py`.

## Resources and results

Preparation templates live in `src/biomolexplorer/resources/chimera/`:

| Template | Responsibility |
| --- | --- |
| `prepare_complex.template` | Select chain and cofactors for the working complex |
| `prepare_receptor.template` | Prepare the receptor and export MOL2 and PDB without H |
| `prepare_ligand.template` | Separate and prepare the reference ligand as PDB and MOL2 |
| `prepare_better_conform.template` | Prepare a conformation for the DOCK6 protocol |

Open Babel produces formats such as PDBQT. Centers and metadata identify the
binding site. Reused receptors must contain the formats required by the engine:
Vina requires PDBQT; DOCK6 also uses MOL2 and PDB without H. The SES surface is
computed by `biomolexplorer.molecular_surface` using NumPy and SciPy and written
in the format consumed by `sphgen`. No standalone DMS executable is needed.

## Execute and diagnose

```text
Connect receptor and compounds to Prepare for docking.
Choose the engine and review receptor and ligand settings.
Run and inspect prepared files, centers and metadata.
Connect outputs to the corresponding engine and inspect its results.
```

Chimera and protocol executables must be accessible to the worker. Command failures
retain diagnostic information; inspect stage logs before changing parameters.
Cancellation stops scientific processes. Files alone do not establish chemical
validity: assess structures, RMSD and scores against suitable study references.

See [installation](installation.md), [redocking](redocking_configuration.md),
[user manual](user_manual.md) and [molecular surfaces](dms_migration.md).
