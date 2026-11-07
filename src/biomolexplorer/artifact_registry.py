"""Transactional persistent index for execution manifests and artifacts."""
import json
import time
from pathlib import Path


def execution_manifest(stage):
    outputs={}
    for record in stage.get('artifact_records',[]):
        port=record.get('output_port') or record.get('output') or 'result'
        outputs.setdefault(port,[]).append(record['occurrence_id'])
    return {'schema_version':1,'stage_id':stage['id'],'outputs':outputs}


def register_run(db,project_id,run_id,stages):
    """Replace the index for one run in the caller's open transaction."""
    db.execute('DELETE FROM artifact_instances WHERE project_id=? AND indexed_run_id=?',(project_id,run_id))
    db.execute('DELETE FROM execution_manifests WHERE project_id=? AND run_id=?',(project_id,run_id))
    now=time.time()
    for stage in stages:
        manifest=execution_manifest(stage);stage['execution_manifest']=manifest
        db.execute('INSERT INTO execution_manifests(project_id,run_id,stage_id,schema_version,manifest,updated) VALUES (?,?,?,?,?,?)',
            (project_id,run_id,stage['id'],1,json.dumps(manifest,sort_keys=True),now))
        for record in stage.get('artifact_records',[]):
            # A cache hit in another block can reference the same immutable
            # occurrence.  The per-stage execution manifest holds that second
            # logical reference; the physical instance is stored only once.
            db.execute('''INSERT OR IGNORE INTO artifact_instances(
                occurrence_id,artifact_id,project_id,indexed_run_id,producer_run_id,stage_id,block_id,output_port,
                artifact_type,schema_version,storage_path,content_sha256,state,metadata,dependencies,provenance,created)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',(
                record['occurrence_id'],record['artifact_id'],project_id,run_id,record.get('run_id'),stage['id'],record.get('block_id'),
                record.get('output_port') or record.get('output'),record.get('artifact_type') or 'other',1,record['path'],
                record.get('content_sha256'),'available' if Path(record['path']).is_file() else 'unavailable',json.dumps(record.get('context',{}),sort_keys=True),
                json.dumps(record.get('parent_artifact_ids',[]),sort_keys=True),json.dumps(record.get('provenance',{}),sort_keys=True),now))
            ensure_entities(db,project_id,record)


def list_artifacts(db,project_id,*,stage_id=None,run_id=None):
    where=['project_id=?'];values=[project_id]
    if stage_id:where.append('stage_id=?');values.append(stage_id)
    if run_id:where.append('indexed_run_id=?');values.append(run_id)
    return [dict(row) for row in db.execute('SELECT * FROM artifact_instances WHERE '+' AND '.join(where)+' ORDER BY created,occurrence_id',values)]


def index_entities(db,project_id,artifact_id,path,artifact_type='other',occurrence_id=None):
    """Build an entity index on demand; large artifacts are never parsed at startup."""
    import csv
    path=Path(path);entities=[]
    if path.suffix.lower()=='.csv':
        with path.open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                found=[]
                for kind,columns in (('compound',('molecule_chembl_id','PubChem_CID','zinc_id','compound_id')),
                                     ('structure',('PDB_CODE','pdb_id')),('receptor',('receptor_id',)),
                                     ('pose',('pose_id','conformer_file')),
                                     ('community',('community_id','graph_component'))):
                    identifier=next((str(row.get(column)).strip() for column in columns if row.get(column) not in (None,'')),None)
                    if identifier:
                        qualifiers={key:row[key] for key in ('CHAIN','chain_id','LIGAND','ligand_id','assembly_id') if row.get(key) not in (None,'')}
                        found.append((kind,identifier,qualifiers))
                if row.get('source') not in (None,'') and row.get('target') not in (None,''):
                    found.append(('molecular_pair',str(row['source'])+'→'+str(row['target']),{}))
                entities.extend(found)
    elif path.suffix.lower() in ('.pdb','.cif','.mmcif','.pdbqt'):
        stem=path.name.split('.',1)[0];parts=stem.split('_')
        identifier=parts[0].upper() if len(parts[0])==4 else stem
        chains=set()
        header_id=None
        if path.suffix.lower() in ('.pdb','.pdbqt'):
            with path.open(encoding='utf-8',errors='replace') as stream:
                for line in stream:
                    if line.startswith('HEADER') and len(line)>=66 and line[62:66].strip():header_id=line[62:66].strip().upper()
                    if line.startswith(('ATOM  ','HETATM')) and len(line)>21 and line[21].strip():chains.add(line[21].strip())
        if header_id:identifier=header_id
        if not chains and len(parts)>1 and len(parts[1])<=4:chains.add(parts[1])
        kind=('pose' if artifact_type in ('molecular_poses','vina','dock6') or any(token in path.name.lower() for token in ('.lig.','pose','scored'))
              else 'receptor' if artifact_type in ('prepared_structures','receptor_pdb') or '.dockprep.' in path.name else 'structure')
        entities.extend((kind,identifier,{'chain_id':chain}) for chain in sorted(chains))
        if not chains:entities.append((kind,identifier,{}))
    occurrence_id=occurrence_id or artifact_id
    for kind,identifier,qualifiers in entities:
        db.execute('INSERT OR REPLACE INTO artifact_entities VALUES (?,?,?,?,?,?,?)',
            (project_id,occurrence_id,artifact_id,kind,identifier,json.dumps(qualifiers,sort_keys=True),'{}'))
    return len(entities)


def ensure_entities(db,project_id,record):
    """Index each immutable content revision once without blocking publication."""
    path=Path(record['path']);sha=record.get('content_sha256')
    if path.suffix.lower() not in ('.csv','.pdb','.cif','.mmcif','.pdbqt'):return 0
    occurrence=record['occurrence_id']
    prior=db.execute('SELECT content_sha256,state FROM artifact_entity_indexes WHERE project_id=? AND occurrence_id=?',
        (project_id,occurrence)).fetchone()
    if prior and prior['content_sha256']==sha and prior['state']=='indexed':
        return db.execute('SELECT COUNT(*) FROM artifact_entities WHERE project_id=? AND occurrence_id=?',
            (project_id,occurrence)).fetchone()[0]
    db.execute('DELETE FROM artifact_entities WHERE project_id=? AND occurrence_id=?',(project_id,occurrence))
    try:
        count=index_entities(db,project_id,record['artifact_id'],path,record.get('artifact_type','other'),occurrence)
        db.execute('INSERT OR REPLACE INTO artifact_entity_indexes VALUES (?,?,?,?,?,?,?)',
            (project_id,occurrence,record['artifact_id'],sha,'indexed',None,time.time()))
        return count
    except (OSError,UnicodeError,ValueError):
        import traceback
        db.execute('INSERT OR REPLACE INTO artifact_entity_indexes VALUES (?,?,?,?,?,?,?)',
            (project_id,occurrence,record['artifact_id'],sha,'error',traceback.format_exc(limit=1)[-1000:],time.time()))
        return 0


def list_entities(db,project_id,artifact_id):
    return [dict(row) for row in db.execute('SELECT * FROM artifact_entities WHERE project_id=? AND artifact_id=? ORDER BY entity_type,entity_id',
        (project_id,artifact_id))]


def migrate_schema(db):
    """Transactionally upgrade registry tables created by early development builds."""
    entity_columns={row['name'] for row in db.execute('PRAGMA table_info(artifact_entities)')}
    if entity_columns and 'occurrence_id' not in entity_columns:
        db.executescript('''
            ALTER TABLE artifact_entities RENAME TO artifact_entities_legacy;
            CREATE TABLE artifact_entities(
                project_id TEXT NOT NULL,occurrence_id TEXT NOT NULL,artifact_id TEXT NOT NULL,
                entity_type TEXT NOT NULL,entity_id TEXT NOT NULL,qualifiers TEXT NOT NULL,metadata TEXT NOT NULL,
                PRIMARY KEY(project_id,occurrence_id,entity_type,entity_id,qualifiers));
            INSERT INTO artifact_entities SELECT project_id,artifact_id,artifact_id,entity_type,entity_id,qualifiers,metadata
                FROM artifact_entities_legacy;
            DROP TABLE artifact_entities_legacy;
            CREATE INDEX artifact_entities_lookup ON artifact_entities(project_id,entity_type,entity_id);
        ''')
    index_columns={row['name'] for row in db.execute('PRAGMA table_info(artifact_entity_indexes)')}
    if index_columns and 'occurrence_id' not in index_columns:
        db.executescript('''
            ALTER TABLE artifact_entity_indexes RENAME TO artifact_entity_indexes_legacy;
            CREATE TABLE artifact_entity_indexes(
                project_id TEXT NOT NULL,occurrence_id TEXT NOT NULL,artifact_id TEXT NOT NULL,content_sha256 TEXT,
                state TEXT NOT NULL,error TEXT,updated REAL NOT NULL,PRIMARY KEY(project_id,occurrence_id));
            INSERT INTO artifact_entity_indexes SELECT project_id,artifact_id,artifact_id,content_sha256,state,error,updated
                FROM artifact_entity_indexes_legacy;
            DROP TABLE artifact_entity_indexes_legacy;
        ''')


class SemanticArtifactResolver:
    """Resolve logical references without inspecting filenames or extensions."""
    def __init__(self,store,project_id,pipeline,results,current_run_id=None):
        self.store,self.project_id=store,project_id
        self.pipeline={stage['id']:stage for stage in pipeline or []}
        self.results=results
        self.current_run_id=current_run_id

    @staticmethod
    def _entity_matches(row,selection):
        qualifiers=json.loads(row['qualifiers'] or '{}')
        identifiers={str(value).casefold() for value in selection.get('entity_ids',[])}
        named={'pdb_id':('structure','receptor'),'receptor_id':('receptor',),'pose_id':('pose',),
               'community_id':('community',),'ligand_id':('compound','pose')}
        if identifiers and str(row['entity_id']).casefold() not in identifiers:return False
        for key,kinds in named.items():
            value=selection.get(key)
            if value is not None and (row['entity_type'] not in kinds or str(row['entity_id']).casefold()!=str(value).casefold()):return False
        for key in ('chain_id','assembly_id','ligand_id'):
            value=selection.get(key)
            if value is not None and key not in named and str(qualifiers.get(key,qualifiers.get(key.upper(),''))).casefold()!=str(value).casefold():return False
            if key=='ligand_id' and value is not None and row['entity_type'] not in named[key] and str(qualifiers.get(key,qualifiers.get('LIGAND',''))).casefold()!=str(value).casefold():return False
        return True

    def _selected(self,reference,records):
        selection=reference.get('selection')
        if not selection:return records
        occurrences=[record.get('occurrence_id') for record in records if record.get('occurrence_id')]
        if occurrences:
            placeholders=','.join('?' for _ in occurrences)
            with self.store.connect() as db:
                rows=db.execute(f'''SELECT occurrence_id,entity_type,entity_id,qualifiers FROM artifact_entities
                    WHERE project_id=? AND occurrence_id IN ({placeholders})''',(self.project_id,*occurrences)).fetchall()
            matched={row['occurrence_id'] for row in rows if self._entity_matches(row,selection)}
            records=[record for record in records if record.get('occurrence_id') in matched]
        else:
            # Records without occurrences predate the registry.  Keep their
            # physical matching solely as an explicit compatibility path.
            from .selection_views import path_matches
            records=[record for record in records if path_matches(record['storage_path'] if 'storage_path' in record else record['path'],selection)]
        if not records:raise ValueError('A seleção não corresponde a nenhuma entidade publicada pela porta semântica.')
        return records

    def _pinned(self,reference):
        with self.store.connect() as db:
            row=db.execute('SELECT manifest FROM execution_manifests WHERE project_id=? AND run_id=? AND stage_id=?',
                (self.project_id,reference['run_id'],reference['stage'])).fetchone()
            if row is None:raise ValueError('A execução histórica não possui manifesto semântico disponível.')
            manifest=json.loads(row[0]);port=reference.get('output_port')
            occurrences=manifest.get('outputs',{}).get(port,[]) if port else [value for values in manifest.get('outputs',{}).values() for value in values]
            if not occurrences:raise ValueError('A porta não foi publicada pela execução histórica selecionada.')
            placeholders=','.join('?' for _ in occurrences)
            rows=db.execute(f'''SELECT occurrence_id,artifact_id,storage_path,state FROM artifact_instances
                WHERE project_id=? AND indexed_run_id=? AND occurrence_id IN ({placeholders})''',
                (self.project_id,reference['run_id'],*occurrences)).fetchall()
        rows=[dict(row) for row in rows if not reference.get('artifact_id') or row['artifact_id']==reference['artifact_id']]
        if reference.get('artifact_id') and len(rows)!=1:raise ValueError('O artefato histórico fixado não pertence à porta selecionada.')
        rows=self._selected(reference,rows)
        paths=[row['storage_path'] for row in rows if row['state']=='available' and Path(row['storage_path']).is_file()]
        if len(paths)!=len(rows) or not paths:raise ValueError('O artefato histórico fixado não está mais disponível.')
        return paths

    def _current(self,reference):
        stage=self.pipeline.get(reference['stage'],{})
        records=stage.get('_execution_artifact_records') or stage.get('artifact_records') or []
        available={str(path) for path in self.results.get(reference['stage'],[])}
        port=reference.get('output_port')
        if not port:
            # Explicit compatibility path for legacy bindings only.
            return list(available)
        matching=[record for record in records if record.get('output_port')==port and record.get('path') in available]
        if not matching and records and not any(record.get('output_port') for record in records):
            # Runs created before semantic manifests have descriptors but no
            # ports.  Classify their saved publication role in this isolated
            # compatibility branch; current runs never fall back to an older
            # execution's results.
            from .artifact_contracts import infer_output_port
            matching=[record for record in records if
                infer_output_port(stage.get('operation'),record.get('output'),record.get('path'),legacy=True)==port
                and record.get('path') in available]
        if not matching and not records and available:
            # Earliest project versions have only the current run's result
            # list.  Their backfill remains the sole filename-based adapter.
            matching=[{'path':path,'storage_path':path} for path in available]
        if reference.get('artifact_id'):
            matching=[record for record in matching if record.get('artifact_id')==reference['artifact_id']]
        if not matching:
            raise ValueError('A execução corrente não publicou a porta semântica selecionada.')
        matching=self._selected(reference,matching)
        return list(dict.fromkeys(record['path'] for record in matching))

    def resolve(self,reference):
        if 'asset' in reference:return None
        return self._pinned(reference) if reference.get('mode')=='pinned' else self._current(reference)


def backfill(db):
    """Index legacy runs lazily and preserve their files and configurations."""
    from .artifacts import index_run
    rows=db.execute('''SELECT r.id,r.project_id,r.stages FROM runs r
        WHERE NOT EXISTS (SELECT 1 FROM execution_manifests m WHERE m.project_id=r.project_id AND m.run_id=r.id)''').fetchall()
    for row in rows:
        stages=json.loads(row['stages'])
        try:
            index_run(stages,row['id'])
            register_run(db,row['project_id'],row['id'],stages)
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),row['id']))
        except (OSError,ValueError,KeyError,TypeError):
            # Ambiguous legacy results stay untouched and visible through their
            # original run instead of receiving a guessed semantic identity.
            continue
