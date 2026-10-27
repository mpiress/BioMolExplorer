"""Immutable, verified per-pair redocking snapshots, published by completed runs."""
import csv
import hashlib
import json
import math
import re
from pathlib import Path
import shutil

from .paths import resolve_path
from .redocking_config import DEFAULTS, pair_key
from .stage_cache import file_digest, implementation_digest
from .storage import write_text
from .calculation_cache import copy_artifact

NAME = 'redocking_cache.json'
PARAMETERS = dict(pH=7.4, sizeof_box=[24, 24, 24], exhaustiveness=20, num_modes=10,
                  prepare_complex=True, charge_type='gas')


def raw_structure(folder,code):
    from .structure_selection import source_structure
    try:return source_structure(folder,code)
    except ValueError:return Path(folder)/(code+'.pdb')


def raw_name(name,code):
    path=Path(name)
    return path.stem==code and path.suffix.lower() in ('.pdb','.cif','.mmcif')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def legacy_indexes(batch):
    """Adopt intact native results predating the pair index, without editing them."""
    configuration = batch.get('configuration', {})
    if configuration.get('provided_results'):
        return
    parameters = configuration.get('parameters', {})
    for name in batch.get('artifacts', []):
        metadata = Path(name)
        if metadata.name != 'pdb_codes.csv' or metadata.parent.parent.name != 'structures':
            continue
        source = metadata.parent; root = source.parent.parent
        with metadata.open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            row.setdefault('STATUS', 'succeeded')
        rows = [row for row in rows if all(re.fullmatch(r'[A-Za-z0-9_+-]{1,40}', row.get(k, ''))
                for k in ('PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN'))]
        legacy_settings={key:dict(value) for key,value in parameters.get('preparation_pairs',{}).items()}
        # Never relabel native/unknown historical preparation as Chimera.
        if any(legacy_settings.get('|'.join(str(row[name]) for name in ('PDB_CODE','LIGAND','RESNUM','CHAIN')),{}).get('preparation_backend')!='legacy_chimera' for row in rows):continue
        for row in rows:
            key='|'.join(str(row[name]) for name in ('PDB_CODE','LIGAND','RESNUM','CHAIN'))
            legacy_settings.setdefault(key,{}).setdefault('preparation_backend','legacy_chimera')
        cache = PairCache(source, root / source.name, [], legacy_settings,
                          {key: parameters.get(key, value) for key, value in PARAMETERS.items()},
                          template_overrides=configuration.get('templates', {}))
        yield root / NAME, cache.publish(rows, write=False)


def stage_sources(store, project_id, stage_id, destination, source=None, configuration=None,
                  only_codes=None, exclude_codes=()):
    """Copy only manifest-verified files from committed runs of this block."""
    with store.connect() as db:
        runs = db.execute("SELECT stages FROM runs WHERE project_id=? AND status='succeeded' ORDER BY created DESC",
                          (project_id,)).fetchall()
    count = 0
    seen = set()
    source = Path(source or destination)
    original_hashes = {}
    def original_matches(name, code, hashes):
        if not raw_name(name,code):return False
        path = store.scoped_path(project_id, name)
        if not path.is_file():return False
        if path not in original_hashes:original_hashes[path] = file_digest(path)
        return original_hashes[path] in hashes
    for run in runs:
        for stage in json.loads(run['stages']):
            if stage['id'] != stage_id or stage['status'] != 'succeeded' or stage.get('outputs_removed'):
                continue
            if stage.get('operation', stage.get('configuration', {}).get('operation')) != 'redocking' or stage.get('provided'):
                continue
            candidates = list(stage.get('batches') or [stage])
            manifest = stage.get('artifact_manifest', {})
            for batch in candidates:
                previous = batch.get('configuration', stage.get('configuration', {}))
                if configuration is not None:
                    from .bindings import sources
                    identities = lambda config: {(r.get('stage'), r.get('asset'))
                        for r in sources(config.get('bindings', {}).get('base_input_path') or {})
                        if r.get('stage') or r.get('asset')}
                    if identities(configuration) != identities(previous):
                        continue
                    if not identities(configuration) and configuration.get('parameters', {}).get('base_input_path') != previous.get('parameters', {}).get('base_input_path'):
                        continue
                indexes = [Path(p) for p in batch.get('artifacts', []) if Path(p).name == NAME]
                try:
                    if indexes:
                        snapshots = []
                        for index in indexes:
                            index = store.scoped_path(project_id, index)
                            if file_digest(index) == manifest.get(str(index.resolve())):
                                snapshots.append((index, json.loads(index.read_text())))
                    else:
                        # Legacy adoption requires a complete original manifest,
                        # including raw structures, metadata, prepared files and poses.
                        if not manifest or any(file_digest(store.scoped_path(project_id, p)) != sha
                                               for p, sha in manifest.items()):
                            continue
                        snapshots = list(legacy_indexes(batch))
                except (OSError, ValueError, KeyError, TypeError):
                    continue
                for index, data in snapshots:
                    try:
                        if data.get('version') != 1:
                            continue
                        if data.get('target') != Path(destination).name:
                            continue
                        entries = []
                        for entry in data['entries']:
                            if entry['pair_key'] in seen:
                                continue
                            code = entry['record'][0]
                            if code in exclude_codes or only_codes is not None and code not in only_codes:
                                continue
                            if not raw_structure(source,code).is_file():
                                # An unselected earlier receptor may be retained,
                                # but never resurrect a source deleted by curation.
                                originals = batch.get('input_files', stage.get('input_files', {})).get('base_input_path', [])
                                parameters = previous.get('parameters', {})
                                if not originals and parameters.get('base_input_path'):
                                    originals = [str(Path(parameters['base_input_path']) / Path(destination).name / f'{code}.pdb')]
                                raw_hashes = {sha for name, sha in entry['files'].items() if raw_name(name,code)}
                                if not any(original_matches(p, code, raw_hashes) for p in originals):
                                    continue
                            entries.append(entry)
                        data['entries'] = entries
                        if not data['entries']:
                            continue
                        files = {f for entry in data['entries'] for f in entry['files']}
                        paths = {}
                        for relative in files:
                            path = store.scoped_path(project_id, index.parent / relative)
                            path.relative_to(index.parent)
                            if file_digest(path) != manifest.get(str(path.resolve())):
                                raise ValueError(relative)
                            paths[relative] = path
                        seed = Path(destination) / 'cache' / f'redocking-{count}'
                        for relative, path in paths.items():
                            target = seed / relative
                            target.parent.mkdir(parents=True, exist_ok=True)
                            copy_artifact(path, target)
                        seed.mkdir(parents=True, exist_ok=True)
                        write_text(json.dumps(data), seed / NAME)
                        seen.update(entry['pair_key'] for entry in data['entries'])
                        count += 1
                    except (OSError, ValueError, KeyError, TypeError):
                        continue
    return count


class PairCache:
    def __init__(self, source, output, records, settings, parameters, template_overrides=None):
        self.source, self.output = Path(source), Path(output)
        self.root = self.output.parent
        self.prepared = self.source / 'Prepared'
        self.settings, self.parameters = settings, parameters
        self.code_hashes = {}
        self.native_identities={}
        self._centers_stamp = None
        self._centers_rows = []
        templates = {}
        from .templates import quiet_template, RESOURCE_ROOT
        for name in ('chimera/prepare_complex.template', 'chimera/prepare_receptor.template',
                     'chimera/prepare_ligand.template', 'vina/config.template'):
            if template_overrides is None:
                text = resolve_path('src/scripts/' + name).read_text()
            else:
                text = template_overrides.get(name, (RESOURCE_ROOT / name).read_text())
            templates[name] = hashlib.sha256(quiet_template(name, text).encode()).hexdigest()
        self.science = dict(code=implementation_digest(), templates=templates)
        self.records = records
        self.restored = {}
        from .calculation_cache import CalculationCache
        self.checkpoints = CalculationCache('redocking-pair', self.source)

    def keys(self, record, config=None):
        from .redocking_config import chimera_preparation
        config = chimera_preparation(self.settings.get(pair_key(record), {}) if config is None else config)
        code, _, _, chain = record[:4]
        if code not in self.code_hashes:
            raw = raw_structure(self.source,code)
            self.code_hashes[code] = file_digest(raw) if raw.is_file() else None
        receptor = dict(source=self.code_hashes[code], chain=chain,
                        options=dict(DEFAULTS, **config.get('receptor', {})),
                        cofactors=sorted(set(config.get('cofactors', []))),
                        pH=self.parameters['pH'], charge_type=self.parameters['charge_type'],
                        prepare_complex=self.parameters['prepare_complex'], science=self.science,
                        backend=config.get('preparation_backend','auto'),structural_selection=config.get('structure',{}))
        if self.parameters['prepare_complex'] and config.get('preparation_backend','auto')!='legacy_chimera':
            from .amber_preparation import parameterization_identity
            from .structure_selection import component_path
            prefix=config.get('ambertools_path','')
            environment=config.get('ambertools_env','')
            minimize=dict(DEFAULTS,**config.get('receptor',{}))['minimize']
            key=(prefix,minimize,environment)
            if key not in self.native_identities:self.native_identities[key]=parameterization_identity(prefix,minimize,environment)
            receptor['native_protocol']=self.native_identities[key]
            receptor['cofactor_topologies']={name:file_digest(component_path(name,self.source,config)) for name in config.get('cofactors',[])}
        # Supplied prepared structures are also scientific inputs.
        if not self.parameters['prepare_complex']:
            receptor['prepared'] = file_digest(self.prepared / f'{code}_{chain}.dockprep.pdbqt')
            if self.parameters.get('redocking_engine')=='dock6':
                receptor['prepared_mol2']=file_digest(self.prepared/f'{code}_{chain}.dockprep.mol2')
        rkey = digest(receptor)
        pair = dict(receptor=rkey, identity=pair_key(record), parameters=self.parameters,
                    options=dict(DEFAULTS, **config.get('ligand', {})),
                    ligand_chain=config.get('ligand_chain',chain),insertion_code=config.get('ligand_insertion_code',''))
        if config.get('preparation_backend')=='legacy_chimera':
            from .chimera_charge_repair import local_definition
            pair['chemical_topology']=local_definition(self.source,record[1],config)
        if self.parameters['prepare_complex'] and config.get('preparation_backend','auto')!='legacy_chimera':
            pair['chemical_topology']=file_digest(component_path(record[1],self.source,config))
            if dict(DEFAULTS,**config.get('ligand',{}))['minimize']:
                key=(prefix,True,environment)
                if key not in self.native_identities:self.native_identities[key]=parameterization_identity(prefix,True,environment)
                pair['minimization']=self.native_identities[key]
        if not self.parameters['prepare_complex']:
            pair['prepared'] = file_digest(self.prepared / (self.prefix(record) + '.lig.pdbqt'))
            if self.parameters.get('redocking_engine')=='dock6':
                pair['prepared_mol2']=file_digest(self.prepared/(self.prefix(record)+'.lig.mol2'))
            pair['centers'] = self.centers(record)
        return rkey, digest(pair)

    @staticmethod
    def prefix(record):
        return f'{record[0]}_{record[1]}_{record[2]}{record[3]}'

    def centers(self, record):
        path = self.prepared / 'centers.csv'
        stat = path.stat();stamp = (stat.st_ino, stat.st_mtime_ns, stat.st_size)
        if stamp != self._centers_stamp:
            with path.open(newline='') as stream:
                self._centers_rows = list(csv.DictReader(stream))
            self._centers_stamp = stamp
        rows = self._centers_rows
        key = self.prefix(record)
        if rows and key not in rows[0]:
            key = f'{record[0]}_{record[1]}_{record[2]}_{record[3]}'
        values = [float(row[key]) for row in rows]
        if len(values) != 3 or not all(math.isfinite(v) for v in values):
            raise ValueError(key)
        return values

    def restore(self):
        """Restore compatible pairs and shared receptors; never consume failures."""
        import pandas as pd
        selected = {pair_key(r): r for r in self.records}
        receptor_keys = {(str(r[0]), str(r[3])): self.keys(r)[0] for r in self.records}
        ready, centers, receptor_files, hashes = set(), {}, {}, {}
        # Per-pair publication survives failure/cancellation of the containing
        # run; only fully validated pairs are entered in this shared index.
        for record in self.records:
            try:
                _,pkey=self.keys(record)
                entry=self.checkpoints.metadata(pkey)
                if not entry:continue
                targets={name:(self.source/Path(name).name if raw_name(name,record[0]) else
                               self.prepared/Path(name).name if 'Prepared' in Path(name).parts else
                               self.output/Path(name).name) for name in entry['files']}
                if self.checkpoints.restore(pkey,targets) is None:continue
                self.restored[pair_key(record)]=entry
                ready.add((record[0],record[3]))
                centers[self.prefix(record)]=entry['center']
            except (OSError,ValueError,KeyError,TypeError):continue
        for index in sorted((self.source / 'cache').glob('redocking-*/' + NAME)):
            try:
                entries = json.loads(index.read_text())['entries']
            except (OSError, ValueError, KeyError):
                continue
            for entry in entries:
                try:
                    record = entry['record']; key = pair_key(record)
                    if (record[0],record[3]) not in receptor_keys:continue
                    paths = {name: (index.parent / name).resolve() for name in entry['files']}
                    for name, path in paths.items():
                        path.relative_to(index.parent.resolve())
                        if path not in hashes:hashes[path] = file_digest(path)
                        if hashes[path] != entry['files'][name]:
                            raise ValueError(name)
                    receptor_prefix = f'{record[0]}_{record[3]}.'
                    receptor = (record[0], record[3])
                    raw = raw_structure(self.source,record[0])
                    old_raw = next((p for p in paths.values() if p.name == raw.name), None)
                    if not raw.exists() and old_raw is not None:
                        self.code_hashes[record[0]] = hashes[old_raw]
                    config = self.settings.get(key, {}) if key in selected else entry['config']
                    rkey, pkey = self.keys(record, config)
                    expected = receptor_keys.get(receptor, rkey)
                    if rkey != entry['receptor_key'] or rkey != expected:
                        continue
                    # Unselected receptors enter the cumulative table only when
                    # their previous poses match this run's scientific settings.
                    if receptor not in receptor_keys and pkey != entry['pair_key']:
                        continue
                    receptor_digest = next(entry['files'][name] for name, path in paths.items()
                                           if path.name == receptor_prefix + 'dockprep.pdbqt')
                    if receptor in receptor_files and receptor_files[receptor] != receptor_digest:
                        continue
                    receptor_files[receptor] = receptor_digest
                    for name, path in paths.items():
                        if receptor not in ready and path.name.startswith(receptor_prefix):
                            self.prepared.mkdir(parents=True, exist_ok=True)
                            copy_artifact(path, self.prepared / path.name)
                    ready.add((record[0], record[3]))
                    if key not in selected:continue
                    if pkey != entry['pair_key'] or key in self.restored:
                        continue
                    for name, path in paths.items():
                        if path.name.startswith(receptor_prefix):continue
                        destination = self.source if path.name == raw.name else self.output if Path(name).parent.name == self.output.name else self.prepared
                        destination.mkdir(parents=True, exist_ok=True)
                        if destination == self.source and raw.exists():continue
                        copy_artifact(path, destination / path.name)
                    centers[self.prefix(record)] = entry['center']
                    self.restored[key] = entry
                except (OSError, ValueError, KeyError, TypeError, StopIteration):
                    continue
        if centers:
            path = self.prepared / 'centers.csv'
            previous = pd.read_csv(path).to_dict('list') if path.exists() else {}
            previous.update(centers)
            pd.DataFrame(previous).to_csv(path, index=False)
        return ready

    def publish(self, rows, write=True):
        entries = []
        hashes = {}
        for row in rows:
            if row['STATUS'] != 'succeeded':
                continue
            record = [row[k] for k in ('PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION')]
            key = pair_key(record)
            config = self.settings.get(key, self.restored.get(key, {}).get('config', {}))
            try:
                rmsd = float(row['RMSD'])
                if not math.isfinite(rmsd) or rmsd < 0:
                    continue
                receptor_key, pkey = self.keys(record, config)
                prefix = self.prefix(record)
                required = [self.prepared / f'{record[0]}_{record[3]}.dockprep.pdbqt',
                            self.prepared / (prefix + '.lig.pdbqt'), self.output / (prefix + '.lig.pdbqt')]
                if not all(p.is_file() and p.stat().st_size for p in required):
                    continue
                files = [p for p in self.prepared.iterdir() if p.is_file() and
                         p.suffix in ('.pdb', '.pdbqt', '.mol2') and
                         p.name.startswith((f'{record[0]}_{record[3]}.', prefix + '.'))]
                files += [p for p in self.output.iterdir() if p.is_file() and p.name.startswith(prefix + '.')
                          and p.suffix in ('.pdb', '.pdbqt', '.mol2', '.csv')]
                raw = raw_structure(self.source,record[0])
                if raw.is_file():files.append(raw)
                for path in files:
                    if path not in hashes:hashes[path] = file_digest(path)
                def relative(path):
                    if path.is_relative_to(self.root):return path.relative_to(self.root).as_posix()
                    return (Path('structures')/self.source.name/path.relative_to(self.source)).as_posix()
                entry=dict(record=record, config=config, receptor_key=receptor_key,
                    pair_key=pkey, row=row, center=self.centers(record),
                    files={relative(p): hashes[p] for p in files})
                entries.append(entry)
                if write:
                    # Missing experimental metadata may be NaN in pandas; it is
                    # not a scientific cache parameter and is stored as null.
                    safe=json.loads(json.dumps(entry,default=str))
                    def finite(value):
                        if isinstance(value,dict):return {k:finite(v) for k,v in value.items()}
                        if isinstance(value,list):return [finite(v) for v in value]
                        return None if isinstance(value,float) and not math.isfinite(value) else value
                    self.checkpoints.remember(pkey,{relative(p):p for p in files},finite(safe))
            except (OSError, ValueError, KeyError, TypeError):
                continue
        data = dict(version=1, target=self.source.name, entries=entries)
        if write:
            self.root.mkdir(parents=True, exist_ok=True)
            write_text(json.dumps(data, default=str), self.root / NAME)
        return data
