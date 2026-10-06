# UCSF DMS reference surfaces

These compressed outputs were generated with the official C distribution from
https://www.cgl.ucsf.edu/Overview/ftp/dms.zip (downloaded 2026-10-05).
Archive SHA-256: `0699283fc4902d4073fba7ccb393e830c9d8b9e1a8b2c4a883b9f0dd102dc011`.
Copyright and redistribution conditions: [UCSF license](../../../src/biomolexplorer/resources/dms/LICENSE).
Originally developed by the UCSF Computer Graphics Laboratory under NIH National
Center for Research Resources grant P41-RR01081.

The only source change required to compile with a modern GNU toolchain was
removing the unused declaration `extern char *sbrk(unsigned long);` in
`dmsd/server.c` (incompatible with unistd.h). No geometry was changed.
Build flags: `-O2 -std=gnu89 -include stdlib.h -include string.h -include unistd.h
-Wno-implicit-function-declaration`. `LIBDIR` points to a temporary directory
containing `dms/dmsd` and `dms/radii` (copied from `radii.proto`).

Each output uses `dms <name>.pdb -d <density> -w <probe_radius> -n -o output.dms`,
with the parameters and uncompressed output checksums in `manifest.json`.
Random cases use NumPy default_rng(5827), 8 normal(0, 2.2) coordinate triplets,
and the atom names C, O, N, S, CA, FE, 1H, P; they exercise mixed radii,
partial/buried surfaces and spindle tori. The analytical fixtures cover a sphere,
pair, triangle, tetrahedron, spindle torus, buried atoms and symmetric square,
pentagon, cube and octahedron (coincident-probe/angle degeneracies). 1CRN contains the
327 ATOM records of https://files.rcsb.org/download/1CRN.pdb; it is a real protein.

Run `python -m unittest discover -s tests -p test_molecular_surface.py` with
`PYTHONPATH=src` in the scientific environment. The executable is unnecessary
for these tests. To compare an independently built reference again, run
`scripts/validate_native_surface.py --dms-executable /path/to/dms`.
Record ordering and the sign of printed zero are ignored; all other characters,
including atom identifiers, coordinates, patch types, areas and normals, must
match. Each record's multiplicity is checked.
