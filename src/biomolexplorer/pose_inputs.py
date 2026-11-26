"""Project-scoped targets, compounds and poses of a pose validation block.

The editor and the executor share this resolution.  Targets are the prepared
receptors of the connected redocking or preparation block, listed by PDB code;
compounds are the molecules of the connected docking or consensus results.  The
poses handed to the worker are the files those blocks already published — no
input is copied or regenerated.
"""
from pathlib import Path

from . import pose_catalog
from .bindings import sources
from .workspace import AccessDenied

TARGET_FIELD = 'base_input_path'
POSE_FIELD = 'base_selected_mols'
MAX_POSES = 20000
# Blocks whose poses complete the selection of a connected aggregate.
DOCKING_OPERATIONS = ('docking_vina','docking_dock6')


class PoseInputs:
    def __init__(self, store, project_id, results=None, pipeline=None):
        self.store, self.pid = store, project_id
        self.results = results or {}
        self.pipeline = list(pipeline or [])
        self._cache = {}

    # ------------------------------------------------------------ references
    def _asset_path(self, identifier):
        with self.store.connect() as db:
            asset = db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',
                               (self.pid,identifier)).fetchone()
        if asset is None:
            raise AccessDenied('Arquivo de pose não autorizado.')
        return self.store.scoped_path(self.pid,asset['path'])

    def _references(self, stage, field):
        stages = []; assets = []
        for reference in sources(stage.get('bindings',{}).get(field,{})):
            if reference.get('stage'):
                stages.append(reference['stage'])
            elif reference.get('asset'):
                assets.append(self._asset_path(reference['asset']))
        return stages,assets

    def signature(self, stage):
        return (tuple(self._references(stage,TARGET_FIELD)[0]),
                tuple(self._references(stage,POSE_FIELD)[0]),
                stage['parameters'].get(TARGET_FIELD),stage['parameters'].get(POSE_FIELD))

    # --------------------------------------------------------------- catalog
    def targets(self, stage):
        key = ('targets',self.signature(stage))
        if key not in self._cache:
            stages,assets = self._references(stage,TARGET_FIELD)
            found = pose_catalog.targets(self.store,self.pid,stages)
            for path in assets:
                identity,code,chain = pose_catalog._receptor_identity(path)
                if identity and path.suffix.lower() in ('.pdb','.pdbqt'):
                    found.append({'receptor_id':identity,'pdb_code':code,'chain':chain,'stage':None,
                                  'files':[str(path)],'file':str(path),'ligand':'','residue':'',
                                  'resolution':'','rmsd':'','label':code+' · arquivo enviado'})
            self._cache[key] = found
        return self._cache[key]

    def records(self, stage):
        key = ('records',self.signature(stage))
        if key not in self._cache:
            stages,assets = self._references(stage,POSE_FIELD)
            tables = pose_catalog.score_tables(self.store,self.pid,stages)
            tables += [path for path in assets if path.suffix.lower() == '.csv']
            records = pose_catalog.pose_records(tables)
            if not records and assets:
                # Uploaded pose files without a score table are still valid.
                from . import pose_sources
                records = pose_sources.collect(assets,[t['file'] for t in self.targets(stage)])
            if len(records) > MAX_POSES:
                raise ValueError(f'Selecione até {MAX_POSES} poses para validar em um bloco.')
            self._cache[key] = records
        return self._cache[key]

    def docking_stages(self, stage):
        """Docking blocks of the pipeline that this block does not read directly."""
        connected = set(self._references(stage,POSE_FIELD)[0])
        return [item['id'] for item in self.pipeline
                if item.get('operation') in DOCKING_OPERATIONS and item['id'] not in connected]

    @staticmethod
    def _complex(record):
        return (str(record.get('receptor_id') or ''),
                str(record.get('compound_id') or record.get('ligand_id') or ''),
                str(record.get('engine') or ''))

    def combined(self, stage):
        """Connected results completed by the blocks that produced each pose.

        A consensus origin republishes copies of the docking poses and keeps a
        single target per compound, so on its own it offers no pose for the
        other targets of a compound.  The docking blocks of the pipeline are
        authoritative for every complex they published; the connected results
        answer for the complexes they do not cover, such as imported poses.
        """
        extra = self.docking_stages(stage)
        key = ('combined',self.signature(stage),tuple(extra))
        if key not in self._cache:
            connected = self.records(stage)
            groups = {}
            if extra:
                tables = pose_catalog.score_tables(self.store,self.pid,extra)
                for record in pose_catalog.pose_records(tables):
                    groups.setdefault(self._complex(record),[]).append(record)
            published = dict(groups)
            for record in connected:
                identity = self._complex(record)
                if identity not in published:
                    groups.setdefault(identity,[]).append(record)
            # The catalog is not bounded here: the limit applies to the poses
            # the selection actually validates, checked when they are resolved.
            self._cache[key] = [record for items in groups.values() for record in items]
        return self._cache[key]

    def compounds(self, stage, selected_targets=None):
        if selected_targets is None:
            selected_targets = stage['parameters'].get('selected_receptors') or ()
        return pose_catalog.compounds(self.records(stage),selected_targets)

    def catalog(self, stage):
        """Targets and compounds as the editor lists them."""
        found = self.targets(stage)
        return found,self.compounds(stage,stage['parameters'].get('selected_receptors') or ())

    # ------------------------------------------------------------- execution
    def poses(self, stage):
        parameters = stage['parameters']
        by_id = {target['receptor_id']:target for target in self.targets(stage)}
        return pose_catalog.selected_poses(self.combined(stage),by_id,
            parameters.get('selected_receptors') or (),parameters.get('selected_ligands') or (),
            parameters.get('excluded_pairs') or ())

    def resolve(self, stage, item=None):
        params = dict(stage['parameters'])
        selected = self.poses(stage)
        if not selected:
            raise ValueError('Nenhuma pose disponível para validar. Conecte um bloco de docking ou consenso '
                             'e confira os targets e compostos selecionados no bloco.')
        keep = ('pose','receptor','reference','prepared_mol2','prepared_pdbqt','pose_id','receptor_id','pdb_code','ligand_id','compound_id',
                'canonical_smiles','engine','score','model','table','row','resolved_by','exclude_residues')
        params['pose_inputs'] = [{key:record[key] for key in keep if record.get(key) not in (None,'')}
                                 for record in selected]
        for record in params['pose_inputs']:
            for key in ('pose','receptor','prepared_mol2','prepared_pdbqt'):
                if record.get(key) and (not Path(record[key]).is_file() or not Path(record[key]).stat().st_size):
                    raise ValueError(f'Artefato {key} indisponível para {record.get("compound_id") or record.get("pose_id")}: {record[key]}')
            if not record.get('receptor'):
                raise ValueError(f'Receptor preparado não encontrado para {record.get("receptor_id")}.')
        poses = [Path(record['pose']) for record in selected]
        receptors = [Path(record['receptor']) for record in selected if record.get('receptor')]
        # The folders only locate the published results for the worker and the
        # cache; the exact files travel in ``pose_inputs``.
        params[POSE_FIELD] = str(poses[0].parent)
        params[TARGET_FIELD] = str(receptors[0].parent) if receptors else str(poses[0].parent)
        if item is not None:
            item.setdefault('input_files',{})[POSE_FIELD] = [str(p) for p in dict.fromkeys(poses)]
            item['input_files'][TARGET_FIELD] = [str(p) for p in dict.fromkeys(receptors or poses)]
            item['input_files']['prepared_compounds'] = list(dict.fromkeys(
                record[key] for record in params['pose_inputs']
                for key in ('prepared_mol2','prepared_pdbqt') if record.get(key)))
        return params
