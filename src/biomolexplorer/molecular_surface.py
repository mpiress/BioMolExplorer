"""Native Python port of UCSF DMS's rolling-probe molecular surface.

Derived from compute.c, input.c, output.c and dmsd/server.c/viewat.c/lookat.c
in https://www.cgl.ucsf.edu/Overview/ftp/dms.zip. Copyright (c) 2002
The Regents of the University of California. Full redistribution conditions
and disclaimer are retained in resources/dms/LICENSE, distributed with this
module. Originally developed by the UCSF Computer Graphics Laboratory under
NIH National Center for Research Resources grant P41-RR01081.

This implements the PDB -> DMS surface path used by DOCK6, with no executable,
network server, or approximation of the surface by solvent-accessible dots.
The historical sampling constant and six-decimal protocol rounding are
intentional: changing them changes the input consumed by sphgen.
"""
from dataclasses import dataclass
from itertools import combinations
from importlib.resources import files
from math import acos, atan2, cos, sin, sqrt, isfinite
from pathlib import Path
import os
import tempfile

import numpy as np
from scipy.spatial import cKDTree

_PI = 3.141592
_RESIDUES = set('ALA ARG ASN ASP CPR CYS CYX CYZ GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()) | {'  A', '  C', '  G', '  T', '  U'}


@dataclass(frozen=True)
class SurfaceSummary:
    atoms: int
    contact_points: int
    saddle_points: int
    concave_points: int
    area: float


@dataclass
class _Atom:
    residue: str
    sequence: str
    name: str
    coord: np.ndarray
    radius: float
    batches: list


def _quantize(value):
    return float(f'{value:.6f}')


def _dist2(vector):
    return vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]


def _unit(vector):
    length = sqrt(_dist2(vector))
    return vector / length if length else None


def _view(origin, target, up):
    z = _unit(target-origin)
    x = _unit(np.cross(up-origin, target-origin))
    if z is None or x is None:
        return None
    y = _unit(np.cross(z, x))
    return np.column_stack((x, y, z))


def _local(points, origin, matrix):
    # Preserve the original homogeneous transform's arithmetic order.
    offset = -origin[0]*matrix[0] - origin[1]*matrix[1] - origin[2]*matrix[2]
    return points[..., 0, None]*matrix[0] + points[..., 1, None]*matrix[1] + points[..., 2, None]*matrix[2] + offset


def _local_look(points, origin, matrix):
    # lookat/probe_angle accumulate translation before the rotation terms.
    offset = -(origin[0]*matrix[0] + origin[1]*matrix[1] + origin[2]*matrix[2])
    return ((offset + points[..., 0, None]*matrix[0]) + points[..., 1, None]*matrix[1]) + points[..., 2, None]*matrix[2]


def _world(points, origin, matrix):
    return ((origin + points[..., 0, None]*matrix[:, 0]) + points[..., 1, None]*matrix[:, 1]) + points[..., 2, None]*matrix[:, 2]


def _rotate(points, matrix):
    return (points[..., 0, None]*matrix[:, 0] + points[..., 1, None]*matrix[:, 1]) + points[..., 2, None]*matrix[:, 2]


def _radii(path):
    text = Path(path).read_text() if path else files('biomolexplorer').joinpath('resources/dms/radii').read_text()
    result = []
    for line in text.splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        name, value = line.split()
        value = float(value)
        if not isfinite(value) or value <= 0:
            raise ValueError('Atomic radii must be positive and finite')
        result.insert(0, (name, value))
    return result


def _read_atoms(path, radii, include_hetero):
    atoms = []
    default = next((r for n, r in radii if n == 'default'), None)
    with Path(path).open() as stream:
        for line_number, line in enumerate(stream, 1):
            record = line[:6].strip()
            if record == 'END':
                break
            if record not in ('ATOM', 'HETATM'):
                continue
            hetero = record == 'HETATM'
            residue = line[17:20]
            if (hetero and not include_hetero) or (not include_hetero and residue not in _RESIDUES):
                continue
            try:
                name = ''.join(line[12:16].split())
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                sequence = str(int(line[22:26]))
                sequence += ''.join(c for c in (line[26:27]+line[21:22]) if c.isalnum())
                if hetero:
                    sequence += '*'
                radius = next((r for n, r in radii if n != 'default' and name.startswith(n)), default)
                if not name or radius is None or not np.isfinite(coord).all():
                    raise ValueError('invalid name, coordinates or missing radius')
                atoms.append(_Atom(residue.strip(), sequence, name, coord, radius, []))
            except (ValueError, IndexError) as error:
                raise ValueError(f'Invalid PDB atom on line {line_number}: {error}') from error
    if not atoms:
        raise ValueError('PDB contains no eligible DMS atoms')
    return atoms


class _Surface:
    def __init__(self, atoms, probe, density):
        self.atoms = sorted(atoms, key=lambda a: tuple(a.coord))
        self.xyz = np.array([a.coord for a in self.atoms])
        self.radii = np.array([a.radius for a in self.atoms])
        self.probe = probe
        self.arc = 1/sqrt(sqrt(3)*density*2.75)
        tree = cKDTree(self.xyz)
        self.nb = []
        maximum = max(self.radii)
        for i, atom in enumerate(self.atoms):
            candidates = tree.query_ball_point(atom.coord, atom.radius+maximum+2*probe)
            neighbors = []
            for j in sorted(candidates):
                if i == j:
                    continue
                v = atom.coord-self.xyz[j]
                distance = v[2]*v[2]+v[1]*v[1]+v[0]*v[0]
                if distance < (atom.radius+self.radii[j]+2*probe)**2:
                    neighbors.append((distance, j))
            self.nb.append([j for _, j in sorted(neighbors)])
        self.probes = []
        self.pair_probes = {}
        self.spheres = {}

    def sphere(self, radius):
        if radius in self.spheres:
            return self.spheres[radius]
        layers = int(_PI/(self.arc/radius)+0.5)+1
        step = _PI/layers
        phi = 0.0
        points = []
        for layer in range(layers):
            rsin = radius*sin(phi)
            z = radius*cos(phi)
            dtheta = 2*_PI if rsin == 0 else self.arc/rsin
            count = max(1, int(2*_PI/dtheta+0.5))
            dtheta = 2*_PI/count
            theta = 0.0 if layer % 2 else _PI
            for _ in range(count):
                points.append((rsin*cos(theta), rsin*sin(theta), z))
                theta += dtheta
                if theta > 2*_PI:
                    theta -= 2*_PI
            phi += step
        result = np.array(points), 4*_PI*radius*radius/len(points)
        self.spheres[radius] = result
        return result

    def add(self, index, kind, points, normals, area):
        if len(points):
            if not np.isfinite(points).all() or not np.isfinite(normals).all() or not isfinite(area):
                raise ValueError('Non-finite molecular surface geometry')
            self.atoms[index].batches.append((kind, points, normals, _quantize(area)))

    def positions(self, i, j, k):
        a, b, c = self.xyz[[i, j, k]]
        matrix = _view(a, b, c)
        if matrix is None:
            # Collinear atoms do not define two discrete probe centers. Test
            # whether the third inflated sphere hides the entire pair circle.
            distance2 = _dist2(b-a)
            if not distance2:
                return 0, ()
            distance = sqrt(distance2)
            r0, r1, r2 = self.radii[[i, j, k]]+self.probe
            length = (r0*r0-r1*r1+distance2)/(2*distance)
            circle2 = r0*r0-length*length
            if circle2 > 0:
                center = a+(b-a)*length/distance
                if _dist2(c-center)+circle2 <= r2*r2:
                    return -1, ()
            return 0, ()
        z1 = _local(b, a, matrix)[2]
        y2, z2 = _local(c, a, matrix)[1:]
        r0, r1, r2 = self.radii[[i, j, k]]+self.probe
        cz = (r0*r0-r1*r1+z1*z1)/(2*z1)
        dz1, dz2 = cz-z1, cz-z2
        cy = (r1*r1-r2*r2-dz1*dz1+dz2*dz2+y2*y2)/(2*y2)
        t = r0*r0-cy*cy-cz*cz
        if t <= 0:
            remaining = r0*r0-cz*cz
            if remaining >= 0:
                cy = -sqrt(remaining)
                if (y2-cy)**2+(z2-cz)**2 <= r2*r2:
                    return -1, ()
            return 0, ()
        if k <= j:
            return 0, ()
        cx = sqrt(t)
        return 1, _world(np.array([[cx, cy, cz], [-cx, cy, cz]]), a, matrix)

    def generate_probes(self):
        bad = []
        for i, neighbors in enumerate(self.nb):
            for j in neighbors:
                if j <= i:
                    continue
                pair = []
                others = set(self.nb[j])
                for k in neighbors:
                    if k not in others:
                        continue
                    status, positions = self.positions(i, j, k)
                    if status == -1:
                        bad.append((i, j))
                        pair = []
                        break
                    if not status:
                        continue
                    occluders = sorted((set(neighbors) | others | set(self.nb[k])) - {i, j, k})
                    for position in positions:
                        delta = self.xyz[occluders]-position
                        distances = 1e-12+delta[:, 0]**2+delta[:, 1]**2+delta[:, 2]**2
                        if np.any(distances <= (self.radii[occluders]+self.probe)**2):
                            continue
                        pair.append(((i, j, k), position))
                self.probes.extend(pair)
        # Both server and controller prepend lists: reverse chronological order.
        self.probes.reverse()
        for i, j in bad:
            self.nb[i].remove(j)
            self.nb[j].remove(i)
        for number, (triplet, _) in enumerate(self.probes):
            for pair in combinations(triplet, 2):
                self.pair_probes.setdefault(tuple(sorted(pair)), []).append(number)

    def contacts(self):
        for i, atom in enumerate(self.atoms):
            points, area = self.sphere(atom.radius)
            points = points+atom.coord
            keep = np.ones(len(points), dtype=bool)
            for j in self.nb[i]:
                distance2 = _dist2(atom.coord-self.xyz[j])
                if distance2 == 0:
                    # Larger coincident atom hides the smaller one.
                    if self.radii[j] > atom.radius:
                        keep[:] = False
                    continue
                distance = sqrt(distance2)
                r0, r1 = atom.radius+self.probe, self.radii[j]+self.probe
                length = (r0*r0-r1*r1+distance2)/(2*distance)
                clip = atom.radius**2+distance2-2*distance*(length*atom.radius/r0)
                delta = self.xyz[j]-points
                keep &= delta[:, 0]**2+delta[:, 1]**2+delta[:, 2]**2 >= clip
            prototype, _ = self.sphere(atom.radius)
            self.add(i, 'SC0', points[keep], prototype[keep]/atom.radius, area)

    def tsection(self, i, j, center, matrix, tradius, clip, start, end, check=True):
        p = self.probe
        if check and tradius < p:
            xy = sqrt(p*p-tradius*tradius)
            # DMS changes this shared clip in place, including across sections.
            # Retain that behavior for downstream reproducibility.
            clip[1] = -clip[1]
            if clip[1] < -xy < clip[0]:
                self.tsection(i, j, center, matrix, tradius, [-xy, -clip[1]], start, end, False)
            if clip[1] < xy < clip[0]:
                self.tsection(i, j, center, matrix, tradius, [clip[0], -xy], start, end, False)
            return
        maxa, mina = -acos(max(-1., min(1., clip[0]/p))), -acos(max(-1., min(1., -clip[1]/p)))
        arange = maxa-mina
        if end < start:
            end += 2*_PI
        erange = end-start
        area = erange*tradius*arange*p-erange*p*p*(sin(maxa+_PI/2)+sin(-_PI/2-mina))
        ne = int(erange*tradius/self.arc+0.5)
        na = int(arange*p/(self.arc*sqrt(3)/2)+0.5)
        if ne <= 0 or na <= 0:
            return
        aincr = arange/na
        a = mina+aincr/2
        mida = (mina+maxa)/2
        startm, offset = 0, 0
        points, normals = [], []
        for _ in range(na):
            z = p*cos(a)
            xy = p*sin(a)+tradius
            ne = int(erange*xy/self.arc+0.5)
            if ne > 0:
                eincr = erange/ne
                e = start+eincr/2*(offset+1)
                offset = 1-offset
                for _ in range(ne):
                    ce, se = cos(e), sin(e)
                    x, y = xy*ce, xy*se
                    points.append((x, y, z))
                    normals.append(((tradius*ce-x)/p, (tradius*se-y)/p, -z/p))
                    e += eincr
            a += aincr
            if a > mida and startm == 0:
                startm = len(points)
        if points:
            points = _world(np.array(points), center, matrix)
            normals = _rotate(np.array(normals), matrix)
            area /= len(points)
            self.add(j, 'SS0', points[:startm], normals[:startm], area)
            self.add(i, 'SS0', points[startm:], normals[startm:], area)

    def tori(self):
        for i, neighbors in enumerate(self.nb):
            for j in neighbors:
                if j <= i:
                    continue
                a, b = self.xyz[[i, j]]
                d2 = _dist2(a-b)
                if not d2:
                    continue
                d = sqrt(d2)
                r0, r1 = self.radii[[i, j]]+self.probe
                length = (r0*r0-r1*r1+d2)/(2*d)
                lsq = r0*r0-length*length
                if lsq <= 0:
                    continue
                tradius = sqrt(lsq)
                center = a+(b-a)*length/d
                clip = [length*(1-self.radii[i]/r0), (d-length)*(1-self.radii[j]/r1)]
                target = a if clip[0] > 0 else center+(center-b)
                av, bv, cv = target-center
                l = sqrt(av*av+cv*cv)
                dd = sqrt(l*l+bv*bv)
                matrix = np.array([[cv/l if l else 1, -av*bv/(l*dd) if l else 0, av/dd], [0, l/dd, bv/dd], [-av/l if l else 0, -bv*cv/(l*dd) if l else 1, cv/dd]])
                angles = []
                for number in self.pair_probes.get((i, j), []):
                    triplet, position = self.probes[number]
                    local = _local_look(position, center, matrix)
                    angle = atan2(local[1], local[0])
                    third = next(k for k in triplet if k not in (i, j))
                    local = _local_look(self.xyz[third], center, matrix)
                    diff = atan2(local[1], local[0])-angle
                    while diff < 0:
                        diff += 2*_PI
                    angles.append((angle if angle >= 0 else angle+2*_PI, diff >= _PI))
                angles.sort(key=lambda item: (item[0], not item[1]))
                if not angles:
                    shared = set(neighbors).intersection(self.nb[j])
                    if not shared or all(_dist2(np.cross(self.xyz[k]-a, b-a)) == 0 for k in shared):
                        self.tsection(i, j, center, matrix, tradius, clip, 0, 2*_PI)
                else:
                    for index, (start, active) in enumerate(angles):
                        if active:
                            self.tsection(i, j, center, matrix, tradius, clip, start, angles[(index+1) % len(angles)][0])

    def concave(self):
        real = np.ones(len(self.probes), dtype=bool)
        if not len(real):
            return
        tree = cKDTree([pos for _, pos in self.probes])
        for index, (triplet, center) in enumerate(self.probes):
            if not real[index]:
                continue
            for other in sorted(tree.query_ball_point(center, 1e-4)):
                if other <= index:
                    continue
                other_triplet, other_center = self.probes[other]
                if _dist2(center-other_center) > 1e-8:
                    continue
                shared = [k for k in triplet if k in other_triplet]
                if len(shared) != 2:
                    continue
                plane = np.cross(self.xyz[shared[0]]-center, self.xyz[shared[1]]-center)
                constant = -float(np.dot(plane, center))
                first = next(k for k in triplet if k not in shared)
                second = next(k for k in other_triplet if k not in shared)
                s1 = constant+float(np.dot(plane, self.xyz[first]))
                s2 = constant+float(np.dot(plane, self.xyz[second]))
                if (s1 < 0 and s2 < 0) or (s1 > 0 and s2 > 0):
                    real[other] = False
        prototype, area = self.sphere(self.probe)
        for index, (triplet, center) in enumerate(self.probes):
            if not real[index]:
                continue
            points = prototype+center
            normals = -prototype/self.probe
            keep = np.ones(len(points), dtype=bool)
            for at, up in combinations(range(3), 2):
                matrix = _view(center, self.xyz[triplet[at]], self.xyz[triplet[up]])
                if matrix is None:
                    keep[:] = False
                    break
                side = _local(self.xyz[triplet[3-at-up]], center, matrix)[0]
                x = _local(points, center, matrix)[:, 0]
                keep &= x <= 0 if side < 0 else x >= 0
            delta = points[:, None, :]-self.xyz[list(triplet)]
            distance = delta[:, :, 0]**2+delta[:, :, 1]**2+delta[:, :, 2]**2
            nearest = distance.argmin(axis=1)
            for at, atom in enumerate(triplet):
                selected = keep & (nearest == at)
                self.add(atom, 'SR0', points[selected], normals[selected], area)


def generate_surface(pdb_path, output_path, *, density=0.5, probe_radius=1.4,
                     normals=True, include_hetero=False, radii_path=None):
    """Write a DMS surface atomically and return counts and unrounded area.

    Defaults correspond to ``dms input.pdb -d .5 -n -w 1.4``. As in DMS,
    ATOM amino/nucleic residues are eligible, HETATM requires explicit opt-in,
    radii use atom-name prefixes, and an existing ``./radii`` overrides the
    bundled UCSF table unless an explicit radii_path is provided.
    """
    density, probe_radius = float(density), float(probe_radius)
    if not isfinite(density) or not 0.1 <= density <= 10:
        raise ValueError('DMS density must be between 0.1 and 10')
    if not isfinite(probe_radius) or not 1 <= probe_radius <= 201:
        raise ValueError('DMS probe radius must be between 1 and 201 angstrom')
    if radii_path is None and Path('radii').is_file():
        radii_path = Path('radii')
    atoms = _read_atoms(pdb_path, _radii(radii_path), include_hetero)
    surface = _Surface(atoms, _quantize(probe_radius), _quantize(density))
    surface.generate_probes()
    surface.contacts()
    surface.tori()
    surface.concave()
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = dict.fromkeys(('SC0', 'SS0', 'SR0'), 0)
    area = 0.0
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, prefix=f'.{output.name}.', suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            for atom in atoms:
                prefix = f'{atom.residue:>3} {atom.sequence:>4} {atom.name[:4]:>4}'
                x, y, z = atom.coord
                stream.write(f'{prefix}{x:8.3f} {y:8.3f} {z:8.3f} A\n')
                for kind, points, vectors, point_area in reversed(atom.batches):
                    counts[kind] += len(points)
                    area += len(points)*point_area
                    for point, vector in zip(points, vectors):
                        x, y, z = point
                        row = f'{prefix}{x:8.3f} {y:8.3f} {z:8.3f} {kind:>3} {point_area:6.3f}'
                        if normals:
                            nx, ny, nz = vector
                            row += f' {nx:6.3f} {ny:6.3f} {nz:6.3f}'
                        stream.write(row+'\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return SurfaceSummary(len(atoms), counts['SC0'], counts['SS0'], counts['SR0'], area)
