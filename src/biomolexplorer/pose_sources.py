"""Pose discovery for validation, by score table or by reference identifier.

Docking and consensus publish score tables whose pose paths are relative to the
table.  A payload that only carries identifiers (receptor_id, ligand_id,
pose_id, compound_id) is still resolvable: the pose and receptor files of the
upstream results are matched by identifier instead of raising an error.

Only the standard library is used here, so the UI, the catalog and the executor
share this resolution without loading the scientific stack.
"""
import csv
import re
from pathlib import Path

POSE_COLUMNS = ('pose_file','conformer_file','dock6_pose','vina_pose','prepared_mol2','prepared_pdbqt')
RECEPTOR_COLUMNS = ('receptor_file','dock6_receptor_file','vina_receptor_file')
REFERENCE_COLUMNS = ('reference_file','vina_reference_file','dock6_reference_file')
LIGAND_COLUMNS = ('ligand_id','LIGAND','molecule_chembl_id','compound_id','zinc_id','PubChem_CID')
POSE_SUFFIXES = ('.mol2','.sdf','.pdbqt','.pdb')
RECEPTOR_MARKERS = ('.dockprep.','.noh.','.prepared.','.amber.','receptor')
POSE_MARKERS = ('.lig.','_scored','pose','ligand')
# A four-character PDB code followed by the ligand code and residue/chain.
POSE_NAME = re.compile(r'^(?P<receptor>[A-Za-z0-9]{4})_(?P<ligand>[A-Za-z0-9]{1,8})(?:_(?P<residue>-?\d+)_?(?P<chain>[A-Za-z0-9])?)?')


def _text(value):
    return '' if value is None else str(value).strip()


def receptor_like(path):
    name = Path(path).name.lower()
    return any(marker in name for marker in RECEPTOR_MARKERS) and '.lig.' not in name


def pose_like(path):
    path = Path(path)
    if path.suffix.lower() not in POSE_SUFFIXES:
        return False
    name = path.name.lower()
    if receptor_like(path):
        return False
    return path.suffix.lower() in ('.mol2','.sdf') or any(marker in name for marker in POSE_MARKERS) or path.suffix.lower() == '.pdbqt'


def identifiers_from_name(path):
    """Receptor, ligand and pose identifiers encoded in a pose filename."""
    stem = Path(path).name.split('.',1)[0]
    match = POSE_NAME.match(stem)
    if not match:
        return {'pose_id':stem}
    receptor = match.group('receptor').upper()
    chain = match.group('chain') or ''
    return {'pose_id':stem,'receptor_id':receptor+('_'+chain if chain else ''),
            'ligand_id':match.group('ligand')}


def _resolve(base, value):
    path = Path(value)
    return path if path.is_absolute() else (Path(base)/path)


def table_records(path):
    """Pose records of one score table, with paths relative to the table."""
    path = Path(path); base = path.parent; records = []
    try:
        with path.open(encoding='utf-8-sig',newline='') as stream:
            rows = list(csv.DictReader(stream))
    except (OSError,UnicodeError,csv.Error):
        return records
    for index,row in enumerate(rows):
        row = {key:value for key,value in row.items() if key}
        shared = {'receptor_id':_text(row.get('receptor_id')) or _text(row.get('PDB_CODE')),
                  'ligand_id':next((_text(row[key]) for key in LIGAND_COLUMNS if _text(row.get(key))),''),
                  'compound_id':_text(row.get('molecule_chembl_id')) or _text(row.get('compound_id')),
                  'canonical_smiles':_text(row.get('canonical_smiles')),
                  'table':str(path),'row':index}
        # Keep preparation lineage separate from the docked pose.  These
        # columns are published by the docking preparation dataset.
        for field in ('prepared_mol2','prepared_pdbqt'):
            if _text(row.get(field)):
                shared[field] = str(_resolve(base,row[field]))
        if _text(row.get('conformer_index')).lstrip('-').isdigit():
            shared['model'] = int(row['conformer_index'])
        receptor = next((_resolve(base,row[key]) for key in RECEPTOR_COLUMNS if _text(row.get(key))),None)
        reference = next((_resolve(base,row[key]) for key in REFERENCE_COLUMNS if _text(row.get(key))),None)
        poses = [(key,_resolve(base,row[key])) for key in POSE_COLUMNS if _text(row.get(key))]
        engines = {'vina_pose':'vina','dock6_pose':'dock6'}
        chosen = [entry for entry in poses if entry[0] in engines] or poses[:1]
        if not chosen:
            records.append(dict(shared,pose=None,receptor=str(receptor) if receptor else None,
                                reference=str(reference) if reference else None,
                                engine=_text(row.get('engine')),score=_text(row.get('score')),
                                pose_id=_text(row.get('pose_id'))))
            continue
        for column,pose in chosen:
            engine = engines.get(column) or _text(row.get('engine'))
            score = _text(row.get(engine)) if engine and _text(row.get(engine)) else _text(row.get('score'))
            specific = next((_resolve(base,row[key]) for key in RECEPTOR_COLUMNS
                             if key.startswith(engine or '~') and _text(row.get(key))),None)
            records.append(dict(shared,pose=str(pose),
                                receptor=str(specific or receptor) if (specific or receptor) else None,
                                reference=str(reference) if reference else None,
                                engine=engine,score=score,
                                pose_id=_text(row.get('pose_id')) or Path(pose).name.split('.',1)[0]))
    return records


def score_tables(paths):
    """Tables that describe poses, newest/most specific layouts included."""
    tables = []
    for path in paths:
        path = Path(path)
        if path.suffix.lower() != '.csv' or not path.is_file():
            continue
        try:
            with path.open(encoding='utf-8-sig',newline='') as stream:
                header = next(csv.reader(stream),[])
        except (OSError,UnicodeError,csv.Error):
            continue
        fields = {name for name in header if name}
        if fields & set(POSE_COLUMNS) or ({'receptor_id','score'} <= fields) or ({'molecule_chembl_id'} <= fields and fields & {'vina','dock6'}):
            tables.append(path)
    return tables


def _candidates(paths):
    poses = []; receptors = []
    for path in paths:
        path = Path(path)
        if path.suffix.lower() not in POSE_SUFFIXES or not path.is_file():
            continue
        if receptor_like(path) or path.suffix.lower() == '.pdb' and not pose_like(path):
            receptors.append(path)
        elif pose_like(path):
            poses.append(path)
    return sorted(poses),sorted(receptors)


def _stem(path):
    """Name without any suffix: ``.mol2`` must not match the token ``MOL2``."""
    return Path(path).name.split('.',1)[0].lower()


def _match(candidates, *tokens):
    """The most specific candidate whose name carries all given identifiers."""
    tokens = [token.lower() for token in tokens if token]
    if not tokens:
        return None
    matches = [path for path in candidates if all(token in _stem(path) for token in tokens)]
    matches.sort(key=lambda path:(len(_stem(path)),path.as_posix()))
    return matches[0] if matches else None


def collect(paths, receptor_paths=()):
    """Pose records from score tables, completed by identifier when needed."""
    paths = [Path(p) for p in dict.fromkeys(map(str,paths))]
    receptor_paths = [Path(p) for p in dict.fromkeys(map(str,receptor_paths))]
    poses,receptors = _candidates(paths)
    extra_poses,extra_receptors = _candidates(receptor_paths)
    receptors = list(dict.fromkeys(receptors+extra_receptors+extra_poses))
    records = []
    for table in score_tables(paths):
        records.extend(table_records(table))
    seen = set()
    for record in records:
        if not record.get('pose') or not Path(record['pose']).is_file():
            # Identifier-based recovery: consensus payloads may carry only IDs.
            receptor_id = record.get('receptor_id') or ''
            # A receptor identifier is PDB_CHAIN; the pose filename may carry
            # only the PDB code, so both spellings are tried.
            code = receptor_id.split('_')[0]
            found = (_match(poses,record.get('pose_id')) or
                     _match(poses,receptor_id,record.get('ligand_id')) or
                     _match(poses,receptor_id,record.get('compound_id')) or
                     _match(poses,code,record.get('ligand_id')) or
                     _match(poses,code,record.get('compound_id')) or
                     _match(poses,record.get('compound_id')) or _match(poses,record.get('ligand_id')))
            record['pose'] = str(found) if found else None
            record['resolved_by'] = 'identifier' if found else 'unresolved'
        if not record.get('receptor') or not Path(record['receptor']).is_file():
            found = _match(receptors,record.get('receptor_id')) if record.get('receptor_id') else None
            if found is None and len(receptors) == 1:
                found = receptors[0]
            record['receptor'] = str(found) if found else None
        if record.get('pose'):
            record.setdefault('pose_id',Path(record['pose']).name.split('.',1)[0])
            for key,value in identifiers_from_name(record['pose']).items():
                if not record.get(key):
                    record[key] = value
        seen.add(record.get('pose'))
    for pose in poses:
        # Poses published without a score table are validated all the same.
        if str(pose) in seen:
            continue
        record = dict(identifiers_from_name(pose),pose=str(pose),engine='',score='',
                      compound_id='',canonical_smiles='',table='',row=None,resolved_by='scan')
        receptor = _match(receptors,record.get('receptor_id')) if record.get('receptor_id') else None
        if receptor is None and len(receptors) == 1:
            receptor = receptors[0]
        record['receptor'] = str(receptor) if receptor else None
        record['compound_id'] = record.get('ligand_id','')
        records.append(record)
    unique = {}
    for record in records:
        if not record.get('pose'):
            continue
        key = (record['pose'],record.get('engine',''),record.get('receptor') or '')
        unique.setdefault(key,record)
    return list(unique.values())


def from_directory(directory, receptor_directory=None):
    directory = Path(directory)
    paths = [p for p in sorted(directory.rglob('*')) if p.is_file()] if directory.is_dir() else [directory]
    receptors = []
    if receptor_directory:
        receptor_directory = Path(receptor_directory)
        receptors = ([p for p in sorted(receptor_directory.rglob('*')) if p.is_file()]
                     if receptor_directory.is_dir() else [receptor_directory])
    return collect(paths,receptors)


def label(record):
    parts = [record.get('receptor_id') or '—',record.get('compound_id') or record.get('ligand_id') or
             Path(record.get('pose') or '').name]
    if record.get('engine'):
        parts.append(record['engine'])
    return ' · '.join(str(part) for part in parts if part)
