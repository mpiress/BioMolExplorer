"""Targets and compounds offered by the blocks upstream of a pose validation.

The catalog is built from the artifact registry and from the consolidated score
tables the docking blocks already published: the editor never loads execution
JSON and never walks the pose tree, so opening the block stays cheap on
projects with hundreds of thousands of artifacts.

Targets are listed by PDB code, as the prepared receptors of redocking and
preparation blocks; compounds are the molecules that actually produced the
docking results, read from those tables.
"""
import csv
import re
from pathlib import Path

RECEPTOR_TYPES = ('prepared_structures','receptor_pdb','structures')
# Semantic ports that publish receptors, plus the per-file ports of the oldest
# executions, whose prepared receptors are only recognizable by their name.
RECEPTOR_PORTS = ('structures','receptors','prepared_receptors','receptor_structures')
PREPARED_PATTERNS = ('%.noH.pdb','%.dockprep.pdbqt','%.prepared.pdb','%.amber.pdb')
POSE_TYPES = ('vina','dock6','molecular_poses','scores','other')
# Prefer a receptor file that is already free of the crystal ligand, solvent
# and hydrogens; a raw complex PDB is the last resort and needs exclusions.
RECEPTOR_PREFERENCE = ('.noh.pdb','.dockprep.pdbqt','.prepared.pdb','.amber.pdb','.pdb')
RECEPTOR_NAME = re.compile(r'^(?P<code>[A-Za-z0-9]{4})(?:_(?P<chain>[A-Za-z0-9]{1,4}))?$')
METADATA_NAMES = ('pdb_codes.csv',)
MAX_TABLES = 8


def _rows(store, project_id, stage_id, *, types=None, ports=None, patterns=None, names=None, newest_only=False):
    """Registry rows of one stage, filtered in SQLite on the stage index.

    Only the columns the catalog needs are read: a docking stage can own
    hundreds of thousands of rows whose provenance JSON is far larger.
    """
    clauses = []; values = [project_id,stage_id]
    for column,items in (('artifact_type',tuple(types or ())),('output_port',tuple(ports or ()))):
        if items:
            clauses.append(column+' IN ('+','.join('?' for _ in items)+')')
            values.extend(items)
    for pattern in tuple(patterns or ())+tuple('%/'+name for name in (names or ())):
        clauses.append('storage_path LIKE ?')
        values.append(pattern)
    query = ('SELECT output_port,artifact_type,storage_path,state,indexed_run_id,created '
             'FROM artifact_instances WHERE project_id=? AND stage_id=?'
             +(' AND ('+' OR '.join(clauses)+')' if clauses else ''))
    with store.connect() as db:
        rows = [dict(row) for row in db.execute(query,values)]
    rows = [row for row in rows if row['state'] == 'available']
    if newest_only and rows:
        newest = max(rows,key=lambda row:row['created'])['indexed_run_id']
        rows = [row for row in rows if row['indexed_run_id'] == newest]
    return rows


def _receptor_identity(path):
    """``4M0E_A`` from ``4M0E_A.dockprep.pdbqt`` or ``4M0E_A.noH.pdb``."""
    stem = Path(path).name.split('.',1)[0]
    match = RECEPTOR_NAME.match(stem)
    if not match:
        return None,None,None
    return stem,match.group('code').upper(),(match.group('chain') or '')


def _preference(path):
    name = Path(path).name.lower()
    for index,suffix in enumerate(RECEPTOR_PREFERENCE):
        if name.endswith(suffix):
            return index
    return len(RECEPTOR_PREFERENCE)


def _metadata(paths):
    """PDB_CODE/LIGAND/RESNUM/CHAIN/RESOLUTION/RMSD of the upstream complexes."""
    records = {}
    for path in paths:
        try:
            with Path(path).open(encoding='utf-8-sig',newline='') as stream:
                for row in csv.DictReader(stream):
                    code = str(row.get('PDB_CODE') or '').strip().upper()
                    chain = str(row.get('CHAIN') or '').strip()
                    if not code:
                        continue
                    entry = records.setdefault(code+('_'+chain if chain else ''),{})
                    for key in ('LIGAND','RESNUM','RESOLUTION','RMSD','STATUS'):
                        if row.get(key) not in (None,'') and not entry.get(key):
                            entry[key] = str(row[key]).strip()
        except (OSError,UnicodeError,csv.Error):
            continue
    return records


def targets(store, project_id, stage_ids):
    """Prepared receptors available upstream, one entry per PDB code and chain."""
    found = {}; metadata_paths = []
    for stage_id in dict.fromkeys(stage_ids):
        rows = _rows(store,project_id,stage_id,types=RECEPTOR_TYPES,ports=RECEPTOR_PORTS,
                     patterns=PREPARED_PATTERNS,names=METADATA_NAMES)
        metadata_paths.extend(row['storage_path'] for row in rows
                              if Path(row['storage_path']).name in METADATA_NAMES)
        for row in rows:
            path = Path(row['storage_path'])
            if path.suffix.lower() not in ('.pdb','.pdbqt'):
                continue
            if '.lig.' in path.name.lower() or '.complex.' in path.name.lower():
                continue
            identity,code,chain = _receptor_identity(path)
            if not identity:
                continue
            entry = found.setdefault(identity,{'receptor_id':identity,'pdb_code':code,'chain':chain,
                                               'stage':stage_id,'files':[],'created':{}})
            entry['files'].append(str(path))
            # Several executions publish the same receptor; the newest wins.
            entry['created'][str(path)] = max(entry['created'].get(str(path),0),row['created'] or 0)
    records = _metadata(metadata_paths)
    # A raw complex is published without a chain in its name; keep it as a
    # lower-preference file of the chain it belongs to instead of a second
    # target for the same PDB code.
    chainless = {entry['pdb_code']:entry for entry in found.values() if not entry['chain']}
    chained = {entry['pdb_code'] for entry in found.values() if entry['chain']}
    for code,entry in chainless.items():
        if code in chained:
            for other in found.values():
                if other['pdb_code'] == code and other['chain']:
                    other['files'].extend(entry['files'])
            found.pop(entry['receptor_id'],None)
    result = []
    for entry in found.values():
        # Registry state is a snapshot: a cleared result must not be offered.
        entry['files'] = [value for value in dict.fromkeys(entry['files']) if Path(value).is_file()]
        if not entry['files']:
            continue
        entry['files'].sort(key=lambda value:(_preference(value),-entry['created'].get(value,0),len(value)))
        detail = records.get(entry['receptor_id'],{}) or records.get(entry['pdb_code'],{})
        entry.pop('created',None)
        entry.update(file=entry['files'][0],ligand=detail.get('LIGAND',''),residue=detail.get('RESNUM',''),
                     resolution=detail.get('RESOLUTION',''),rmsd=detail.get('RMSD',''),
                     label=entry['pdb_code']+(' · cadeia '+entry['chain'] if entry['chain'] else ''))
        result.append(entry)
    return sorted(result,key=lambda item:(item['pdb_code'],item['chain']))


def score_tables(store, project_id, stage_ids):
    """Consolidated score tables of the connected docking or consensus blocks."""
    from .artifact_choices import consolidated_docking_tables
    selected = []
    for stage_id in dict.fromkeys(stage_ids):
        candidates = [Path(row['storage_path']) for row in
                      _rows(store,project_id,stage_id,types=POSE_TYPES,patterns=('%.csv',),newest_only=True)
                      if row['storage_path'].lower().endswith('.csv')]
        # A job publishes its complete table at the artifacts root; the deeper
        # copies belong to single receptors or input snapshots.
        roots = [path for path in candidates if path.parent.name == 'artifacts']
        tables = consolidated_docking_tables(roots or candidates)
        tables.sort(key=lambda path:(len(path.parts),path.as_posix()))
        selected.extend(tables[:MAX_TABLES])
    return list(dict.fromkeys(selected))


def pose_records(tables):
    """Pose rows of those tables, with pose and receptor paths resolved."""
    from .pose_sources import table_records
    records = []
    for table in tables:
        records.extend(table_records(table))
    unique = {}
    for record in records:
        key = (record.get('receptor_id'),record.get('compound_id') or record.get('ligand_id'),
               record.get('engine'),record.get('pose'))
        unique.setdefault(key,record)
    return list(unique.values())


def _keys(value):
    """``4M0E_A`` and ``4M0E`` identify the same target for a selection."""
    value = str(value or '').strip()
    return {value,value.split('_')[0]} - {''}


def matches_target(selection, receptor_id, pdb_code=None):
    if not selection:
        return True
    accepted = set().union(*(_keys(value) for value in selection)) if selection else set()
    return bool(accepted & (_keys(receptor_id) | _keys(pdb_code)))


def compounds(records, selected_targets=()):
    """Molecules that produced the connected results, by compound identifier."""
    selected = [str(value) for value in selected_targets if str(value).strip()]
    result = {}
    for record in records:
        code = str(record.get('compound_id') or record.get('ligand_id') or '').strip()
        if not code:
            continue
        receptor = str(record.get('receptor_id') or '')
        if not matches_target(selected,receptor):
            continue
        entry = result.setdefault(code,{'compound_id':code,'smiles':record.get('canonical_smiles') or '',
                                        'receptors':set(),'engines':set(),'poses':0,'best_score':None,
                                        'file':record.get('pose'),'prepared_file':None})
        entry['poses'] += 1
        entry['smiles'] = entry['smiles'] or (record.get('canonical_smiles') or '')
        entry['file'] = entry['file'] or record.get('pose')
        entry['prepared_file'] = entry['prepared_file'] or next(
            (record[key] for key in ('prepared_mol2','prepared_pdbqt')
             if record.get(key) and Path(record[key]).is_file()),None)
        if receptor:
            entry['receptors'].add(receptor)
        if record.get('engine'):
            entry['engines'].add(str(record['engine']))
        try:
            score = float(record.get('score'))
        except (TypeError,ValueError):
            score = None
        if score is not None and (entry['best_score'] is None or score < entry['best_score']):
            entry['best_score'] = score
    for entry in result.values():
        entry['receptors'] = sorted(entry['receptors']); entry['engines'] = sorted(entry['engines'])
    return sorted(result.values(),key=lambda item:item['compound_id'])


def pose_identifier(record):
    """Readable pose identity; consensus publishes poses named by row index."""
    stem = Path(record.get('pose') or '').name.split('.',1)[0]
    code = str(record.get('compound_id') or '')
    if code and code.casefold() in stem.casefold():
        return stem
    parts = [str(record.get('receptor_id') or ''),code,str(record.get('engine') or '')]
    return '_'.join(part for part in parts if part) or stem


def complex_key(compound, receptor):
    """Identity of one ligand-receptor complex in the block parameters."""
    return str(compound or '')+':'+str(receptor or '')


def selected_poses(records, targets_by_id, selected_targets=(), selected_compounds=(), excluded_pairs=()):
    """Pose payloads for the chosen targets and compounds, reusing published files.

    ``excluded_pairs`` removes single complexes the editor dropped by their row,
    without removing the compound from every other target.
    """
    targets_filter = [str(v) for v in selected_targets if str(v).strip()]
    compounds_filter = {str(v) for v in selected_compounds if str(v).strip()}
    excluded = {str(v) for v in excluded_pairs if str(v).strip()}
    by_code = {}
    for identity,target in targets_by_id.items():
        by_code.setdefault(str(target.get('pdb_code') or identity.split('_')[0]),target)
    result = []
    for record in records:
        receptor = str(record.get('receptor_id') or '')
        code = str(record.get('compound_id') or record.get('ligand_id') or '')
        if not matches_target(targets_filter,receptor):
            continue
        if compounds_filter and code not in compounds_filter:
            continue
        target = targets_by_id.get(receptor) or by_code.get(receptor.split('_')[0],{})
        if excluded & {complex_key(code,receptor),
                       complex_key(code,target.get('pdb_code') or receptor.split('_')[0])}:
            continue
        payload = dict(record)
        pose = payload.get('pose')
        if not pose or not Path(pose).is_file():
            continue
        receptor_file = payload.get('receptor')
        prepared = target.get('file')
        # The connected target preparation is authoritative for this chain.
        if prepared and Path(prepared).is_file() and _preference(prepared) < _preference(receptor_file or ''):
            receptor_file = prepared
        if not receptor_file or not Path(receptor_file).is_file():
            receptor_file = prepared
        if receptor_file and Path(receptor_file).is_file():
            payload['receptor'] = str(receptor_file)
        else:
            payload.pop('receptor',None)
        # Only a raw complex still holds the crystal ligand and its solvent;
        # a prepared receptor is already clean.
        if (target.get('ligand') and payload.get('receptor')
                and _preference(payload['receptor']) >= len(RECEPTOR_PREFERENCE)-1):
            payload['exclude_residues'] = [target['ligand']]
        payload['receptor_id'] = receptor
        payload['compound_id'] = code
        payload['pdb_code'] = target.get('pdb_code') or receptor.split('_')[0]
        payload['pose_id'] = pose_identifier(payload)
        result.append(payload)
    return result
