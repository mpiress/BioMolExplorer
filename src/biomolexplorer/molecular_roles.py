"""Content-level molecular roles used by MD's typed file selectors.

Polymer residue evidence is a role check, not a topology/force-field validation.
Scientific preparation and per-pair validation remain mandatory at execution.
"""
from functools import lru_cache
from pathlib import Path

POLYMER_RESIDUES=frozenset(('ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY HIS HID HIE HIP '
    'ILE LEU LYS LYN MET MSE PHE PRO SER THR TRP TYR VAL SEC PYL '
    'A C G U T I DA DC DG DT DU DI ADE CYT GUA THY URA').split())


def _role(path,kind):
    path=Path(path)
    try:
        stat=path.stat()
        return _content_role(str(path.resolve()),kind,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
    except OSError:return False


@lru_cache(maxsize=256)
def _content_role(path,kind,*signature):
    with Path(path).open(errors='replace') as stream:
        if kind=='receptor':
            return any(line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() in POLYMER_RESIDUES for line in stream)
        lines=iter(stream)
        for line in lines:
            if line.strip()=='@<TRIPOS>MOLECULE':
                next(lines,'');next(lines,'')
                return next(lines,'').strip().upper()=='SMALL'
    return False


def macromolecule_pdb(path):
    return Path(path).suffix.lower()=='.pdb' and _role(path,'receptor')


def small_molecule_file(path):
    path=Path(path)
    if '.dockprep.' in path.name or '.noH.' in path.name:return False
    if path.suffix.lower()=='.mol2':return _role(path,'ligand')
    return path.suffix.lower() in ('.csv','.sdf','.pdbqt')
