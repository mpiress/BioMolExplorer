"""Independent receptor and reference-ligand preparation checkpoints."""
from pathlib import Path

from .calculation_cache import CalculationCache, identity
from .ligand_cache import preparation_identity
from .redocking_config import DEFAULTS, pair_key
from .stage_cache import file_digest
from .storage import write_dataframe


class StructureCache:
    def __init__(self, source, prepared, settings, charge_type, ph):
        self.source, self.prepared = Path(source), Path(prepared)
        self.settings, self.charge_type, self.ph = settings, charge_type, ph
        self.cache = CalculationCache('structure-preparation', self.source)
        self.software = None
        self.native_software={}
        from .paths import SOURCE_ROOT
        self.structure_adapter=file_digest(SOURCE_ROOT/'caad/docking.py')
        self.hashes = {}
        # Read effective overlays, not only the packaged templates.
        from .paths import resolve_path
        self.templates = {name: file_digest(resolve_path('src/scripts/chimera/' + name))
                          for name in ('prepare_complex.template', 'prepare_receptor.template', 'prepare_ligand.template')}

    def keys(self, record):
        code, ligand, resnum, chain = record[:4]
        from .structure_selection import source_structure
        raw = source_structure(self.source,code)
        if raw not in self.hashes:self.hashes[raw] = file_digest(raw)
        from .redocking_config import chimera_preparation
        config = chimera_preparation(self.settings.get(pair_key(record), {}))
        software=self.software
        chemistry={}
        if self.software is None:self.software=dict(preparation_identity(),structure_adapter=self.structure_adapter)
        software=self.software
        receptor = dict(source=self.hashes[raw], code=code, chain=chain, ph=self.ph, charge_type=self.charge_type,
                        options=dict(DEFAULTS, **config.get('receptor', {})),
                        cofactors=sorted(set(config.get('cofactors', []))), software=software,
                        structural_selection=config.get('structure',{}),backend=config.get('preparation_backend','auto'),
                        cofactor_templates={},
                        templates={k:v for k,v in self.templates.items() if k != 'prepare_ligand.template'})
        rkey = identity(receptor)
        if config.get('preparation_backend')=='legacy_chimera':
            from .chimera_charge_repair import local_definition
            chemistry=local_definition(self.source,record[1],config)
        ligand = dict(receptor=rkey, pair=pair_key(record),
                      options=dict(DEFAULTS, **{'charge_type':self.charge_type, **config.get('ligand', {})}),
                      ligand_protocol={name:file_digest(Path(__file__).with_name(name)) for name in
                          ('chimera_preparation.py','ligand_mol2.py','dock6_validation.py')},
                      chemical_templates=chemistry,ligand_chain=config.get('ligand_chain',chain),
                      insertion_code=config.get('ligand_insertion_code',''),template=self.templates['prepare_ligand.template'])
        return rkey, identity(ligand)

    def locked(self, record):
        return self.cache.locked(self.keys(record)[0])

    def restore(self, record):
        import pandas as pd
        rkey, lkey = self.keys(record)
        ready=[]
        for key in (rkey, lkey):
            metadata=self.cache.metadata(key)
            if not metadata:
                ready.append(False);continue
            if any(Path(name).name!=name for name in metadata.get('files',[])):
                ready.append(False);continue
            targets={name:self.prepared/name for name in metadata['files']}
            restored=self.cache.restore(key,targets)
            ready.append(restored is not None)
            if restored and 'center' in restored:
                path=self.prepared/'centers.csv'
                centers=pd.read_csv(path).to_dict('list') if path.is_file() else {}
                centers[self.prefix(record)]=restored['center']
                write_dataframe(pd.DataFrame(centers),path)
        return tuple(ready)

    @staticmethod
    def prefix(record):
        return f'{record[0]}_{record[1]}_{record[2]}{record[3]}'

    def remember_receptor(self, record):
        from .ligand_cache import validate_prepared_ligand
        validate_prepared_ligand(self.prepared/f'{record[0]}_{record[3]}.dockprep.pdbqt','pdbqt')
        self.remember(record,receptor_only=True)

    def remember(self, record, receptor_only=False):
        import pandas as pd
        rkey,lkey=self.keys(record)
        roles=[(rkey,f'{record[0]}_{record[3]}.')]
        if not receptor_only:roles.append((lkey,self.prefix(record)+'.'))
        for key,prefix in roles:
            files={p.name:p for p in self.prepared.iterdir()
                   if p.is_file() and p.name.startswith(prefix) and p.suffix in ('.pdb','.pdbqt','.mol2','.json')}
            if not files:continue
            metadata=dict(files=list(files),record=list(map(str,record[:4])))
            if key==lkey:
                metadata['center']=pd.read_csv(self.prepared/'centers.csv')[self.prefix(record)].tolist()
            self.cache.remember(key,files,metadata)
