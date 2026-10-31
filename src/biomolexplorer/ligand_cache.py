"""Verified, per-compound preparation checkpoints, scoped to one project."""
from contextlib import contextmanager
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import tempfile

from .stage_cache import file_digest
from .storage import write_json


def validate_prepared_ligand(path, format):
    """Reject non-finite or structurally empty native-engine inputs.

    Native Vina aborts the process, rather than raising a Python exception,
    for some malformed PDBQT geometries.  This validation therefore belongs
    at the preparation/cache boundary, before native code sees the file.
    """
    import math
    path=Path(path)
    if not path.is_file() or not path.stat().st_size:
        raise ValueError('Arquivo preparado ausente ou vazio: '+path.name)
    if format=='pdbqt':
        atoms=[]
        for line in path.read_text(errors='replace').splitlines():
            if not line.startswith(('ATOM  ','HETATM')):continue
            try:
                coordinates=[float(line[start:end]) for start,end in ((30,38),(38,46),(46,54))]
                fields=line.split();charge=float(fields[-2]);atom_type=fields[-1]
                if not all(math.isfinite(value) for value in coordinates+[charge]) or not atom_type:raise ValueError()
            except (ValueError,IndexError):
                raise ValueError('PDBQT do ligante contém coordenadas, cargas ou tipos atômicos inválidos: '+path.name) from None
            atoms.append(line)
        if not atoms:raise ValueError('PDBQT do ligante não contém átomos: '+path.name)
    elif format=='mol2':
        from .dock6_validation import validate_mol2
        validate_mol2(path)
    else:raise ValueError('Formato de ligante preparado não reconhecido: '+str(format))
    return path


def preparation_identity():
    """Fingerprint the adapter and the installed preparation tools once per batch."""
    files = {}
    root = Path(__file__).parent
    for name in ('ligand_cache.py', 'docking_preparation.py', 'docking_data.py',
                 'visualizations.py', 'molecular_parsing.py', 'redocking_config.py', 'chimera_mol2.py', 'chimera_charge_repair.py'):
        files[name] = file_digest(root / name)
    for name in ('chimera', 'obabel'):
        executable = shutil.which(name)
        if not executable:
            files[name] = None
            continue
        executable = Path(executable).resolve()
        files[str(executable)] = file_digest(executable)
        if name == 'chimera':
            installation = executable.parent.parent
            # Chimera's launcher delegates the scientific work to these modules.
            for folder in ('share/DockPrep', 'share/Rotamers', 'share/AddCharge', 'share/AddH', 'share/MMMD',
                           'share/WriteMol2', 'share/ReadMol2', 'share/AmberInfo',
                           'share/chimera', 'lib'):
                for path in sorted((installation / folder).rglob('*')):
                    if path.is_file() and path.suffix in ('.py', '.so', '.dat'):
                        files[str(path)] = file_digest(path)
            for amber in sorted((installation / 'bin').glob('amber*')):
                for tool in ('antechamber', 'sqm', 'bondtype', 'atomtype', 'parmchk2'):
                    path = amber / 'bin' / tool
                    if path.is_file():
                        files[str(path)] = file_digest(path)
                for path in sorted((amber / 'dat').rglob('*')):
                    if path.is_file():
                        files[str(path)] = file_digest(path)
        else:
            for path in sorted((executable.parent.parent / 'lib').glob('libopenbabel*.so*')):
                if path.is_file():
                    files[str(path)] = file_digest(path)
            for path in sorted((executable.parent.parent / 'share' / 'openbabel').rglob('*')):
                if path.is_file():
                    files[str(path)] = file_digest(path)
    versions = {}
    for name in ('rdkit', 'openbabel', 'openbabel-wheel'):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return dict(files=files, versions=versions)


class LigandCache:
    def __init__(self, destination,backend='legacy_chimera'):
        self.root = Path(os.environ.get('BIOMOL_PREPARATION_CACHE') or
                         Path(destination) / 'cache' / 'ligand-preparation')
        self.root.mkdir(parents=True, exist_ok=True)
        self.reuse = os.environ.get('BIOMOL_REUSE_PREPARATIONS', '1') != '0'
        self.identity=hashlib.sha256(json.dumps(preparation_identity(),sort_keys=True).encode()).hexdigest()

    def key(self, row, dataset, options, ph, engines):
        pose = row.get('conformer_file') or row.get('pose_file')
        structure = dict(smiles=row.get('canonical_smiles', ''))
        if pose:
            source = Path(pose)
            if not source.is_absolute():
                source = Path(dataset).parent / source
            structure = dict(pose_sha256=file_digest(source), format=source.suffix.lower(),
                             index=int(row.get('conformer_index') or 1))
        request = dict(version=1, code=row['molecule_chembl_id'], structure=structure,
                       options=options, ph=float(ph), engines=engines, software=self.identity)
        return hashlib.sha256(json.dumps(request, sort_keys=True, allow_nan=False).encode()).hexdigest()

    @contextmanager
    def locked(self, key):
        # An OS lock is released even when the worker or its process group is cancelled.
        import fcntl
        with (self.root / (key + '.lock')).open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            yield

    def restore(self, key, row, destination, formats):
        if not self.reuse:
            return False
        entry = self.root / key
        try:
            manifest = json.loads((entry / 'manifest.json').read_text())
            if manifest['key'] != key:
                return False
            sources = {fmt: entry / ('ligand.' + fmt) for fmt in formats}
            if not all(source.stat().st_size and file_digest(source) == manifest['files'][fmt]
                       for fmt, source in sources.items()):
                return False
            for fmt,source in sources.items():validate_prepared_ligand(source,fmt)
        except (OSError, ValueError, KeyError, TypeError):
            return False
        for fmt, source in sources.items():
            output = Path(destination) / (row['molecule_chembl_id'] + '.lig.' + fmt)
            shutil.copy2(source, output)
            row['prepared_' + fmt] = str(output.resolve())
        if manifest.get('preparation_metadata'):row['preparation_metadata']=manifest['preparation_metadata']
        return True

    def save(self, key, row, formats,molecular_details=None):
        entry = self.root / key
        entry.mkdir(exist_ok=True)
        # Publish the manifest last; incomplete or corrupt entries are never reused.
        hashes = {}
        for fmt in formats:
            source = Path(row['prepared_' + fmt])
            validate_prepared_ligand(source,fmt)
            with tempfile.NamedTemporaryFile(dir=entry, delete=False) as stream:
                temporary = Path(stream.name)
            try:
                shutil.copy2(source, temporary)
                hashes[fmt] = file_digest(temporary)
                temporary.replace(entry / ('ligand.' + fmt))
            finally:
                temporary.unlink(missing_ok=True)
        write_json(dict(version=1, key=key, software_sha256=self.identity, files=hashes,
            preparation_metadata=row.get('preparation_metadata'),molecular_details=molecular_details), entry / 'manifest.json')
