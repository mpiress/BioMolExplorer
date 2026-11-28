"""Register, resolve, materialize and explain artifacts over one content store.

The catalog of PHASE-01 says which artifact an occurrence is and what produced
it; `content_store` holds the bytes.  This module is the boundary between them,
and the only place that answers the four questions a consumer actually asks:

* **register** — put the bytes of a published artifact under management, and
  record the content object the catalog already pointed at by digest;
* **resolve** — hand back content whose digest has been verified, refusing a
  divergence instead of treating it as a hit;
* **materialize** — give a job its inputs, shared and read-only by default, or
  as an isolated writable copy for an engine that rewrites what it reads;
* **provenance** — the producing operation, its inputs, its software and its
  run, read from the catalog rather than inferred from filenames.

Ordering is the invariant that makes an interruption recoverable: the bytes are
in the store before any row references them.  A publication cut short leaves
unreferenced bytes, which are inert and deduplicated by the next publication,
never a row pointing at content that is absent.  ADR-008 records this.
"""
import time
from contextlib import contextmanager
from pathlib import Path

from .stage_cache import file_digest

STORE_FOLDER=('.biomolexplorer','store')


class ArtifactManager:
    """Artifact services for one project."""

    def __init__(self,workspace,project_id):
        self.workspace=workspace
        self.project_id=project_id

    @property
    def project_root(self):
        return Path(self.workspace.project_dir(self.project_id))

    @property
    def store(self):
        from .content_store import ContentStore
        return ContentStore(self.project_root.joinpath(*STORE_FOLDER))

    @contextmanager
    def _catalog(self,db=None):
        """Join the caller's transaction, or open one of our own."""
        if db is not None:
            yield db
        else:
            with self.workspace.connect() as owned:
                yield owned

    def _instance(self,db,occurrence_id):
        rows=[dict(row) for row in db.execute('''SELECT * FROM artifact_instances
            WHERE project_id=? AND occurrence_id=? ORDER BY created''',(self.project_id,occurrence_id))]
        if not rows:
            raise LookupError('A ocorrência não está no catálogo do projeto: '+str(occurrence_id))
        # The same occurrence is indexed by every run that references it, and a
        # reuse can have materialized it elsewhere.  Prefer a row whose bytes
        # are present.
        return next((row for row in rows if Path(row['storage_path']).is_file()),rows[0])

    def _reference(self,digest):
        """Where the object lives, relative to the project folder it moves with."""
        return '/'.join((*STORE_FOLDER,'objects',digest[:2],digest))

    def register(self,occurrence_id,*,db=None,validate=True):
        """Put a published artifact under management and record its content."""
        with self._catalog(db) as catalog:
            row=self._instance(catalog,occurrence_id)
            path=Path(row['storage_path'])
            if not path.is_file():
                raise LookupError('O artefato publicado não está disponível: '+str(path))
            digest=row['content_sha256'] or file_digest(path)
            store=self.store
            check=None
            if validate:
                from .artifact_validation import validator
                check=validator(row['artifact_type'] or 'other')
            store.put(path,digest=digest,validate=check)
            self._record(catalog,digest,path,store)
            return digest

    def register_many(self,occurrence_ids,*,db=None,validate=True):
        """Register several artifacts, reporting what each one refused."""
        registered={};refused={}
        with self._catalog(db) as catalog:
            for occurrence in dict.fromkeys(occurrence_ids):
                try:
                    registered[occurrence]=self.register(occurrence,db=catalog,validate=validate)
                except (LookupError,ValueError,OSError) as error:
                    refused[occurrence]=str(error)
        return {'registered':registered,'refused':refused}

    def publish(self,occurrence_ids,*,validate=True):
        """Register artifacts under a journal that survives an interruption.

        The journal names the content before it is written, the bytes are
        written before the rows, and the rows commit with the run state.  An
        interruption on either side of that commit is recovered by
        ``result_transactions.recover`` (ADR-008).
        """
        from .result_transactions import ResultTransaction
        occurrences=list(dict.fromkeys(occurrence_ids))
        with self.workspace.connect() as db:
            transaction=ResultTransaction(self.workspace,db,self.project_id)
            planned=[]
            for occurrence in occurrences:
                try:
                    planned.append(self._digest(self._instance(db,occurrence)))
                except LookupError:
                    continue
            transaction.record_objects([digest for digest in planned if digest])
            result=self.register_many(occurrences,db=db,validate=validate)
        transaction.finish()
        return result

    def _record(self,db,digest,path,store):
        from .artifact_catalog import media_type
        size=store.size(digest)
        reference=self._reference(digest)
        db.execute('INSERT OR IGNORE INTO content_objects VALUES (?,?,?,?,?)',
            (digest,size,media_type(path),reference,time.time()))
        # An object registered by PHASE-01 has neither size nor location: the
        # phase recorded the identity, this one puts the bytes under management.
        db.execute('UPDATE content_objects SET size=?,storage_ref=? WHERE sha256=?',(size,reference,digest))

    def _digest(self,row):
        path=Path(row['storage_path'])
        if row['content_sha256']:return row['content_sha256']
        # Rows indexed before publication recorded a manifest carry no digest.
        # It is computed rather than written back: the catalog of an earlier
        # release is read, never rewritten.
        return file_digest(path) if path.is_file() else None

    def resolve(self,occurrence_id,*,db=None,validate=False):
        """Content whose digest has been verified, from the store or its origin."""
        from .artifact_validation import validate as apply_contract
        with self._catalog(db) as catalog:
            row=self._instance(catalog,occurrence_id)
        digest=self._digest(row)
        store=self.store
        if digest and store.has(digest):
            if not store.verify(digest):
                raise ValueError('O objeto armazenado divergiu do seu identificador de conteúdo: '+digest)
            return store.object_path(digest)
        path=Path(row['storage_path'])
        if not path.is_file():
            raise LookupError('O conteúdo do artefato não está disponível: '+str(path))
        if digest and file_digest(path)!=digest:
            raise ValueError('O artefato publicado divergiu da identidade registrada: '+str(path))
        if validate:apply_contract(row['artifact_type'] or 'other',path)
        return path

    def resolve_entity(self,entity_type,entity_id,*,qualifiers=None,db=None):
        """Verified content for every artifact that published one entity."""
        from .artifact_catalog import artifacts_by_entity
        with self._catalog(db) as catalog:
            rows=artifacts_by_entity(catalog,self.project_id,entity_type,entity_id,qualifiers)
            return [self.resolve(row['occurrence_id'],db=catalog) for row in rows]

    def materialize(self,occurrence_id,target,*,writable=False,db=None):
        """Give a job an input: shared and read-only, or an isolated copy.

        A shared materialization is a hard link to the read-only object, so it
        costs no bytes on a filesystem without copy-on-write and the kernel
        refuses an engine that writes in place.  An engine that rewrites its own
        inputs needs ``writable``, and then nothing it does can reach the store.
        """
        from .content_store import clone_or_copy
        import os
        with self._catalog(db) as catalog:
            row=self._instance(catalog,occurrence_id)
        source=self.resolve(occurrence_id,db=db)
        digest=self._digest(row)
        store=self.store
        target=Path(target)
        if digest and store.has(digest) and source==store.object_path(digest):
            return store.materialize(digest,target,'isolated' if writable else 'shared',verify=False)
        # Not under management yet: the isolation rules still apply to the copy
        # handed to the job.
        target.parent.mkdir(parents=True,exist_ok=True)
        clone_or_copy(source,target)
        os.chmod(target,0o600 if writable else 0o444)
        if digest and file_digest(target)!=digest:
            raise ValueError('A materialização não reproduziu o conteúdo do artefato.')
        return target

    def provenance(self,occurrence_id,*,db=None):
        """Producer operation, inputs, software and run, from the catalog."""
        from .artifact_catalog import content_object,operation_edges,producing_operation
        with self._catalog(db) as catalog:
            row=self._instance(catalog,occurrence_id)
            operation=producing_operation(catalog,self.project_id,occurrence_id)
            edges=operation_edges(catalog,operation['operation_id']) if operation else {'inputs':[],'outputs':[]}
            entities=[dict(entity) for entity in catalog.execute('''SELECT entity_type,entity_id,qualifiers
                FROM artifact_entities WHERE project_id=? AND occurrence_id=?
                ORDER BY entity_type,entity_id''',(self.project_id,occurrence_id))]
            digest=row['content_sha256']
            content=content_object(catalog,digest) if digest else None
        return {'occurrence_id':occurrence_id,'artifact_id':row['artifact_id'],
                'artifact_type':row['artifact_type'],'content_sha256':digest,
                'storage_path':row['storage_path'],'indexed_run_id':row['indexed_run_id'],
                'producer_run_id':row['producer_run_id'],'stage_id':row['stage_id'],
                'output_port':row['output_port'],'operation':operation,
                'inputs':edges['inputs'],'entities':entities,'content':content,
                'managed':bool(digest and self.store.has(digest))}

    def usage(self,*,db=None):
        """What the store holds, and what the catalog says about it."""
        store=self.store
        measured=store.usage()
        with self._catalog(db) as catalog:
            rows=catalog.execute('''SELECT COUNT(*) AS objects,COALESCE(SUM(size),0) AS bytes
                FROM content_objects WHERE storage_ref IS NOT NULL''').fetchone()
            references=catalog.execute('''SELECT COUNT(*) FROM artifact_instances
                WHERE project_id=? AND content_sha256 IS NOT NULL''',(self.project_id,)).fetchone()[0]
        return {'store':measured,'registered':{'objects':rows['objects'],'bytes':rows['bytes']},
                'logical_references':references}
