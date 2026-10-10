"""Resolve separate graph experiments and their project-scoped compound provenance."""
import copy
from .artifact_choices import matches_selector
from pathlib import Path

from .bindings import sources
from .fingerprint_selection import filename_kind
from .input_validation import columns, contract, validate_file
from .workspace import AccessDenied

REQUIRED={'fingerprints':{'molecule_chembl_id','fingerprint'},'similarity':{'source','target','value'},
          'compounds':{'canonical_smiles','molecule_chembl_id'}}


class GraphInputs:
    def __init__(self,store,project_id,pipeline,results):
        self.store,self.pid,self.results=store,project_id,results
        self.stages={s['id']:s for s in pipeline}

    def path(self,filename):
        path=Path(filename)
        return self.store.scoped_path(self.pid,path if path.is_absolute() else self.store.project_dir(self.pid)/path)

    def files(self,reference,kind):
        if 'asset' in reference:
            with self.store.connect() as db:
                row=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(self.pid,reference['asset'])).fetchone()
            if row is None:raise AccessDenied('Arquivo não autorizado.')
            files=[self.store.scoped_path(self.pid,row[0])]
        else:
            files=[self.store.scoped_path(self.pid,p) for p in self.results.get(reference['stage'],[])]
        selector=reference.get('selector','auto')
        candidates=[p for p in files if p.is_file() and p.suffix.lower()=='.csv' and REQUIRED[kind]<=columns(p)]
        if selector!='auto':
            candidates=[p for p in candidates if matches_selector(p,selector)]
            if len(candidates)!=1:raise ValueError('Escolha um arquivo de '+kind+' disponível na origem: '+selector)
        elif kind=='compounds':
            candidates.sort(key=lambda p:(p.name!='compounds.csv',p.name!='molecules.csv',len(p.parts),p.name))
            candidates=candidates[:1]
        return sorted(candidates)

    def lineage(self,reference,filename,seen=None):
        if 'asset' in reference:
            files=self.files(reference,'fingerprints')
            return [p for p in files if REQUIRED['compounds']<=columns(p)]
        producer=self.stages.get(reference.get('stage'))
        if not producer:return []
        batch=next((b for b in producer.get('_execution_batches',[]) if str(filename) in b.get('artifacts',[])),None)
        batch_paths=batch.get('artifacts') if batch else None
        if batch:producer=batch['configuration']
        seen=set(seen or ())
        if producer['id'] in seen:return []
        seen.add(producer['id'])
        paths=[self.store.scoped_path(self.pid,p) for p in batch_paths or self.results.get(producer['id'],[])]
        structure_files=[p for p in paths if REQUIRED['compounds']<=columns(p)]
        name=Path(filename).stem
        hints=[name,name.split('_',1)[-1]]
        for hint in list(hints):
            kind=filename_kind(hint+'.csv')
            if kind:hints.append(hint.removeprefix(kind+'_'))
        exact=[p for p in structure_files if p.stem in hints]
        if exact:return exact
        if producer['operation']=='fingerprints' and structure_files:
            # Merged inputs rename the similarity table. The fingerprint
            # outputs retain the actual IDs and SMILES used in that experiment.
            selector=reference.get('selector','auto')
            algorithm=filename_kind(filename) or filename_kind(filename.split('_',1)[-1]) or filename_kind(selector)
            if algorithm:structure_files=[p for p in structure_files if filename_kind(p.name) in (None,algorithm)]
            selected=[p for p in structure_files if matches_selector(p,selector)] if selector!='auto' else structure_files
            if selected:return selected
        if producer['operation'] not in ('fingerprints','similarity') or producer.get('provided_results'):
            return structure_files
        group=producer.get('bindings',{}).get('base_input_path')
        found=[]
        for parent in sources(group) if group else []:
            if producer['operation']=='fingerprints':
                # Look at all producer tables before applying the default union.
                candidates=self.files(dict(parent,selector='auto'),'compounds')
                if 'stage' in parent and parent.get('selector','auto')=='auto':
                    all_paths=[self.path(p) for p in self.results.get(parent['stage'],[])]
                    specific=[p for p in all_paths if p.stem in hints and REQUIRED['compounds']<=columns(p)]
                    if specific:candidates=specific
                matching=[p for p in candidates if p.stem in hints]
                found+=matching or candidates
            else:
                matching=self.files(parent,'fingerprints')
                for path in matching:found+=self.lineage(parent,str(path),seen)
        if not found and producer['parameters'].get('base_input_path'):
            folder=self.path(producer['parameters']['base_input_path'])
            found=[p for p in folder.glob('*.csv') if REQUIRED['compounds']<=columns(p)]
        return list(dict.fromkeys(found))

    def resolve(self,stage,item=None):
        from .operations import validate_operation
        params=copy.deepcopy(stage['parameters']);bindings=stage.get('bindings',{})
        validate_operation('graphs',params)
        if set(bindings)-{'similarity_path','base_input_path'}:
            raise ValueError('Grafos aceitam apenas entradas de Calcular similaridade ou arquivos externos de similaridade.')
        explicit=[]
        group=bindings.get('base_input_path')
        for ref in sources(group) if group else []:
            if 'stage' in ref:
                raise ValueError('A tabela complementar de compostos deve ser um arquivo externo. Os SMILES do pipeline são identificados automaticamente.')
            candidates=self.files(ref,'compounds')
            if not candidates:raise ValueError('Inclua molecule_chembl_id e canonical_smiles na tabela externa de compostos.')
            for path in candidates:validate_file(path,'compounds',validate_rows=False)
            explicit+=candidates
        if params.get('base_input_path'):
            folder=self.path(params['base_input_path'])
            candidates=[folder] if folder.is_file() else sorted(folder.glob('*.csv'))
            explicit += [p for p in candidates if REQUIRED['compounds']<=columns(p)]
            if not explicit:raise ValueError('Inclua molecule_chembl_id e canonical_smiles na tabela externa de compostos.')
            for path in explicit:validate_file(path,'compounds',validate_rows=False)
        entries=[]
        group=bindings.get('similarity_path')
        for ref in sources(group) if group else []:
            producer=self.stages.get(ref.get('stage'),{})
            if 'stage' in ref and producer.get('operation')!='similarity':
                raise ValueError('Conecte apenas blocos Calcular similaridade à entrada dos grafos.')
            candidates=self.files(ref,'similarity')
            if not candidates:raise ValueError('A origem não possui similaridade pronta. Padrão esperado: '+contract('similarity'))
            for path in candidates:
                validate_file(path,'similarity',validate_rows=False)
                compounds=self.lineage(ref,str(path)) if 'stage' in ref else []
                entries.append({'kind':'similarity','file':str(path),'compound_files':[str(p) for p in compounds or explicit],
                    'label':producer.get('name','Arquivo próprio')+' · '+path.name,
                    'metric':producer.get('parameters',{}).get('metric','pronta'),
                    'fingerprint':filename_kind(path.name) or filename_kind(path.name.split('_',1)[-1]) or producer.get('parameters',{}).get('fingerprint','morgan')})
        if params.get('similarity_path'):
            folder=self.path(params['similarity_path'])
            paths=[folder] if folder.is_file() else sorted(folder.glob('*.csv'))
            candidates=[p for p in paths if REQUIRED['similarity']<=columns(p)]
            if not candidates:raise ValueError('Padrão esperado: '+contract('similarity'))
            for path in candidates:
                validate_file(path,'similarity',validate_rows=False)
                entries.append({'kind':'similarity','file':str(path),'compound_files':[str(p) for p in explicit],'label':path.name})
        for entry in params.get('graph_inputs',[]):
            entry['file']=str(self.path(entry['file']))
            entry['compound_files']=[str(self.path(p)) for p in entry['compound_files']]
            validate_file(entry['file'],'similarity',validate_rows=False)
            for path in entry['compound_files']:validate_file(path,'compounds',validate_rows=False)
            entries.append(entry)
        if not entries:raise ValueError('Conecte Calcular similaridade ou envie um CSV externo. Padrão esperado: '+contract('similarity'))
        # Keep different experiments separate, including their compound metadata.
        for entry in entries:
            candidates=entry['compound_files'];name=Path(entry['file']).stem
            hints=[name,name.split('_',1)[-1]]
            for hint in list(hints):
                kind=filename_kind(hint+'.csv')
                if kind:hints.append(hint.removeprefix(kind+'_'))
            exact=[p for p in candidates if Path(p).stem in hints]
            if exact:entry['compound_files']=exact
        unique=[]
        for entry in entries:
            if not any(e['file']==entry['file'] for e in unique):unique.append(entry)
        params['graph_inputs']=unique
        for key in ('base_input_path','similarity_path'):params.pop(key,None)
        if item is not None:
            item['input_files']={'similarity_path':[e['file'] for e in unique],
                'base_input_path':list(dict.fromkeys(p for e in unique for p in e['compound_files']))}
        return params
