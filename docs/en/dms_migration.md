# Native DMS algorithm port

[Documentation](../README.md) · English · [Português](../dms_migration.md)

BioMolExplorer now computes the molecular surface in Python. `Dock6.prepare_surface()` directly calls `biomolexplorer.molecular_surface.generate_surface()` without installing or executing `dms`/`dmsd`. The `{receptor}.dms` artifact, density/probe parameters and `sphgen → sphere_selector → showbox → grid → DOCK6` sequence are retained.

## Original source and license

The [official UCSF distribution](https://www.cgl.ucsf.edu/Overview/ftp/dms.zip), linked from the [UCSF software page](https://www.cgl.ucsf.edu/Overview/software.html), was downloaded and compiled. ZIP SHA-256: `0699283fc4902d4073fba7ccb393e830c9d8b9e1a8b2c4a883b9f0dd102dc011`.

| Original source | Ported responsibility |
| --- | --- |
| `ms.c` | Parameters and accepted ranges for the application path |
| `input.c`, `radii.proto` | PDB eligibility, identifiers, atom-name prefix radii |
| `compute.c` | Ordering, scheduling, neighbor removal, near-coincident probe collapse |
| `dmsd/server.c` | Neighbors, probe centers, occlusion, sphere sampling, contact/toroidal/concave patches |
| `dmsd/viewat.c`, `lookat.c` | Geometry coordinate transforms |
| `output.c` | Atom/point formatting, patch types, areas and normals |

The original permits redistribution and modification with attribution. Its full license is retained in `src/biomolexplorer/resources/dms/LICENSE` and distributed with the Python package. Originally developed by the UCSF Computer Graphics Laboratory with NIH National Center for Research Resources support, grant P41-RR01081. This derived component retains that license; the project's general MIT license does not replace its terms.

## Implemented behavior

This is a rolling-probe solvent-excluded surface, rather than a SASA approximation. It includes convex contact patches (`SC0`), pairwise toroidal patches including spindle tori (`SS0`), and concave spherical triple patches (`SR0`), with occlusion and near-coincident probe collapse. Atom association, per-point areas and oriented normals follow the original rules.

The UCSF radius table is bundled. An existing `./radii` file or explicit `radii_path` overrides it. PDB eligibility, insertion codes and chain identifiers are preserved. Optional inclusion of HETATM also accepts nonstandard ATOM residues, as with `dms -a`. Output is written atomically after successful computation; point counts and area are available for logging.

Historical `PI=3.141592`, layered sampling, internal `2.75 × density`, rounding and six-decimal C protocol precision are intentional. NumPy and SciPy's `cKDTree` implement spatial searches and filtering. The port covers the entire application call: `dms receptor.noH.pdb -d <density> -n -w <radius> -v -o receptor.dms`. Distributed DMS server administration, its complete CLI and unused residue selectors are outside this application path. Computation runs locally in the Python worker.

## Usage and environment

NumPy and SciPy are already scientific environment dependencies; no new library is needed. `requirements.yml` is the available Conda manifest. `environment.yml` had already been removed before this port and was not recreated. The installer no longer downloads, compiles or installs DMS.

```python
from biomolexplorer.molecular_surface import generate_surface

summary = generate_surface(
    "receptor.noH.pdb", "receptor.dms",
    density=0.5, probe_radius=1.4,
)
print(summary)
```

Accepted `density`: 0.1–10; `probe_radius`: 1–201 Å, matching the original. Defaults remain 0.5 and 1.4 Å. `normals=False` omits normals; `include_hetero=True` accepts all residues and HETATM; `radii_path` selects a custom table. DOCK6 integration retains normals and excludes HETATM, matching its previous invocation.

## Equivalence evidence

The C distribution was built in a temporary directory. Only an unused obsolete `sbrk` declaration conflicting with modern headers was removed; geometry was unchanged.

`tests/fixtures/dms` retains 31 C-generated reference surfaces, with parameters and checksums in `manifest.json`: an isolated sphere, pair, triangle, tetrahedron, spindle torus, buried atoms, 20 mixed-radius coordinate sets, four symmetric geometries (square, pentagon, cube and octahedron) and the real 1CRN protein (331 atoms). Tests compare every field and record multiplicity. Only record order and rounded negative zero are normalized; C server completion order is unspecified.

At the defaults, both implementations produced 2,479 contact, 2,384 saddle and 1,390 concave points for 1CRN (6,253 total), with identical formatted fields. Additional protein comparisons at `(density=0.1, probe_radius=1.0)` and `(1.0, 2.0)` also matched. The default native calculation took about 8.6 seconds on the review environment; this is a single observation, not a performance guarantee for large receptors.

```bash
PYTHONPATH=src:tests python -m unittest test_molecular_surface test_docking_handoffs
PYTHONPATH=src python scripts/validate_native_surface.py --dms-executable /path/to/dms
```

The first command needs no DMS installation. The second audits an independently built C reference. Fixture provenance/build details are in `tests/fixtures/dms/README.md`. The integration test runs the real generator and verifies the remaining sphgen/sphere_selector command and file handoffs.

## Validation limits

Collinear/coincident atoms can cause undefined C behavior and `NaN` output. The native port explicitly handles blocked probe circles and never writes non-finite geometry; separate tests cover degeneracies. Equivalence to invalid original `NaN` records is not intended.

Real sphgen/DOCK6 was not run during this review because those executables are unavailable here. Evidence establishes surface data equivalence for the tested inputs and preserved integration contracts; it does not establish identical scores/poses or performance for every receptor. The `.dms` generation step is replaced; Chimera, DOCK6 and its accessories remain dependencies of the other existing steps.
