> Registro histórico. As instruções de ambiente/protocolos desta validação foram substituídas pela [revisão de Amber Score nativo de 2026-10-29](native_amber_score_2026-10-29/README.md). Não representam as opções atuais da aplicação.

# AmberTools isolated-environment validation — 2026-10-24

- Main runtime: `BioMolExplorer`, Python 3.12.3. Imports verified: RDKit, Open Babel, Meeko, Vina, Gemmi, NumPy and Pandas.
- Amber runtime: `BioMolExplorer-AmberTools`, AmberTools 26.0.
- Found in the Amber environment: `antechamber`, `parmchk2`, `tleap`, `sander`, `cpptraj`, `ambpdb`, `sqm`.
- Optional and not installed: `pmemd`, `reduce`. Neither is called by the current application.
- Invocation verified as `conda run --no-capture-output -n BioMolExplorer-AmberTools /absolute/environment/bin/<tool> ...`; the main Python process remained in `BioMolExplorer`. The absolute executable prevents DOCK6's bundled `tleap` from shadowing AmberTools when its directory is first on `PATH`.

## AM1-BCC, tleap and minimization

Ethanol was prepared through the application API. Antechamber produced nine finite AM1-BCC charges with total `0.000001`; parmchk2/tleap produced MOL2, frcmod, prmtop and inpcrd artifacts. A short sander minimization completed with final potential energy `-1.2374 kcal/mol`.

## Redocking, DOCK6 and Vina

The repository's real 1LGU/AAH validation completed receptor and ligand preparation, fixed-pose scoring, Vina docking, and flexible DOCK6 redocking (SES, spheres, box, grid and scored poses). The native receptor retained all 1,299 heavy atoms and their identifiers/coordinates before minimization.

- Native Vina: score `-5.588 kcal/mol`; heavy-atom pose RMSD `0.3782 Å`.
- Native fixed-pose DOCK6 continuous score: `6.693309`.
- Flexible DOCK6: scored poses produced for legacy, native and native-minimized preparations.

Three receptors from the failing Alzheimer project were parameterized again with the corrected PDB-to-LEaP residue mapping. `tleap` completed for 4M0E (532 residues, three disulfides), 2QP8 (390 residues, three disulfides), and 6UWV (391 residues, three disulfides), with finite integral net charges. The generated bond records retain the original chain/residue identity and use a separate dense address in the deliberately renumbered LEaP input.

The cached 4M0E receptor was also resolved using the existing downstream docking reference. Its canonical locator and artifact ID remained unchanged after materialization below `cache/redocking-0`, eliminating the observed `A identidade do arquivo selecionado mudou` failure.

1P0I/BUA and 7MYI/ZQS contain inter-residue bonds to the selected ligand. They remain explicitly excluded from the ordinary noncovalent workflow because accepting them would require a dedicated covalent/coordination parameterization protocol.
When a partial redocking report explicitly identifies such a failed receptor, downstream Vina/DOCK6 variants now omit only that exact receptor and continue with the successfully prepared selections. Unreported missing or identity-mismatched artifacts remain hard errors; no stale receptor is substituted.

## Short MD and trajectory analysis

A 20-step sander trajectory was run through the isolated environment. It produced restart, NetCDF trajectory and output files. cpptraj read four frames and produced finite heavy-atom RMSD values: `0.0001`, `0.0308`, `0.0487`, `0.0425 Å`.

## Automated tests

The final focused AmberTools/artifact/preparation/redocking/docking suite passed 91 tests. Separate real-tool checks covered AM1-BCC, ff14SB/GAFF topology generation, minimization, short MD, cpptraj, Vina and DOCK6. The complete repository suite previously ran 991 tests; unrelated failures remain in pre-existing graph/UI/cache fixtures and curated inputs.
