import copy,json,time,shutil,csv,hashlib,sys
from pathlib import Path
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.pipeline import PipelineService,execution_context
from biomolexplorer.artifacts import validate_reference
import biomolexplorer.docking_data as dd
root=Path('/tmp/biomol-vina-bbb-atomic-replay');root.mkdir(exist_ok=True)
store=WorkspaceStore(root/'workspace');token=store.register('Validation','vina-bbb@local.test','temporary-validation-password');project=store.create_project(token,'BBB+ replay');pid=project['id'];dest=store.project_dir(pid)
original_root=Path('/home/michel/Downloads/Alzheimer');stages=json.loads(Path('/tmp/biomol-vina-bbb-snapshot.json').read_text())
original_docking=next(s for s in stages if s['operation']=='docking_vina');needed={r['stage'] for group in original_docking['configuration']['bindings'].values() for r in group.get('sources',[group]) if 'stage' in r}
manifest={}
for s in stages:
 if s['id'] not in needed:continue
 for value in s['artifacts']:
  source=Path(value);relative=source.relative_to(original_root);target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
  manifest[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
def rebase(value):
 if isinstance(value,str):return str(dest)+value[len(str(original_root)):] if value.startswith(str(original_root)+'/') else value
 if isinstance(value,list):return [rebase(x) for x in value]
 if isinstance(value,dict):return {rebase(k):rebase(v) for k,v in value.items()}
 return value
stages=rebase(stages);stage=next(s for s in stages if s['operation']=='docking_vina')['configuration'];results={s['id']:s['artifacts'] for s in stages if s['id'] in needed};user=store.user(token)['id']
Path(root/'configuration.json').write_text(json.dumps(dict(project=pid,token=token,stage=stage,stages=stages,results=results)))
service=PipelineService(store);counts=dict(reads=0,rows=0);reader=dd.read_compounds
def read(path):
 rows=reader(path);counts['reads']+=1;counts['rows']+=len(rows);return rows
dd.read_compounds=read
proof=dict(original_files=len(manifest),configuration_unchanged=True,cases=[])
try:
 variants=list(service._variants(pid,stage,results));proof['variants']=len(variants)
 for iteration in range(2):
  start=time.perf_counter();before=counts.copy();cases=[]
  for i,(variant,label) in enumerate(variants):
   params=service._resolve(pid,user,variant,results,cache_only=True,pipeline=execution_context(stages))
   path=Path(params['base_selected_mols'])/(params['mol_filename']+'.csv')
   with path.open() as f:rows=sum(1 for _ in csv.DictReader(f))
   cases.append(dict(label=label,rows=rows,dataset=str(path),receptor=params.get('pdb_code'),receptor_folder=params.get('base_input_path')))
   print(json.dumps(dict(iteration=iteration,batch=i+1,rows=rows,elapsed=time.perf_counter()-start)),flush=True)
  proof['cases'].append(dict(iteration=iteration,seconds=time.perf_counter()-start,reads=counts['reads']-before['reads'],rows_read=counts['rows']-before['rows'],batches=cases))
 proof['original_files_unchanged']=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in manifest.items())
 Path(root/'report.json').write_text(json.dumps(proof,indent=2));print(json.dumps({k:v for k,v in proof.items() if k!='cases'}),flush=True)
finally:service.close()
