# Chimera-free main preparation path

> Native migration history. Since 2026-10-26, docking and redocking use Chimera-only preparation again. See the [current configuration](molecular_preparation.md). The descriptions below document the earlier implementation.

Updated October 23, 2026.

## Current status

In-process Python/Gemmi structure selection, trusted-graph crystal ligand extraction, ff14SB receptor preparation, direct ligand charging and sander minimization replace Chimera on the standard path. Meeko/AutoDockAdapter exports Vina inputs; DOCK6 receives charged SYBYL MOL2.

legacy_chimera remains an explicit option for validated special systems, historical scientific protocols and custom Chimera templates. Missing native tools or trusted chemistry produce actionable errors, not silent legacy fallback. The historical caad.redocking preparation API delegates to the shared caad.docking adapter.

## Scientific interpretation

- STRUCTURE_SELECTION: operational equivalence; selected heavy atoms, identifiers and coordinates are preserved.
- CRYSTAL_LIGAND_PREPARATION: CCD or named SDF/MOL2 topology is matched to crystal coordinates; bond orders are never invented from distances.
- RECEPTOR_PREPARATION: ff14SB is retained; protonation decisions, H placement and LEaP-built missing atoms are documented.
- CHARGE_ASSIGNMENT: direct Antechamber Gasteiger/AM1-BCC; RDKit is allowed for standard non-minimized Vina. Amber minimization preserves the selected charge method.
- MINIMIZATION: controlled MMTK-to-sander change, vacuum ff14SB/GAFF, 100 SD + 10 CG cycles. Numerical identity is not required.
- FORMAT_CONVERSION: in-process adapters, no external process just to select atoms.
- LEGACY_FALLBACK: explicit compatibility only.

Meeko poses can recover their explicit chemical graph/index mapping. Unannotated PDBQT needs the corresponding chemical source. The legacy reference-ligand template was also corrected: reopening PDB before MOL2 export lost assigned charges. MOL2 now comes directly from the charged molecule; old accidental zero-charge results are not a scientifically valid reference.

## Validation and limitations

Real 1LGU/AAH checks cover paired structural selection, receptor preparation, Vina redocking, continuous DOCK6 scoring, native SES/spheres/box/grid and flexible DOCK6 docking. Short minimization checks energy trends, RMSD, clashes and bond lengths. Cold/warm reuse and absence of Chimera calls are measured.

See [real results](../validation/native_structures_2026-10-23.json) and [configuration/scientific protocol](molecular_preparation.md). [Earlier reports](../validation/chimera_remaining_2026-10-23.md) are historical, not the current dependency status.

Metals/coordination, covalent ligands, radicals, unsupported charge states and residues without parameters require dedicated handling. Model-compound-pKa protonation is not site-specific prediction. Separate vacuum minimization can worsen a fixed-pose interaction score even while lowering intramolecular energy; docking pose quality must be assessed separately. One protein target does not establish systematic non-regression or universal chemical coverage.

Sources: [wwPDB CCD](https://www.wwpdb.org/data/ccd), [Amber sander tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php). Dynamics keeps its separate DOCK6/AmberTools protocol and remains Chimera-free.
