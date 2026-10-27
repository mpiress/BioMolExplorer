"""Controlled preparation-only benchmark; no docking timings or inferred data.

Run with the scientific Python and PATH, e.g. --legacy-max 10. Repeated
chemistry (seven representative graphs) deliberately exposes charge reuse.
The legacy comparison explicitly disables minimization in BOTH protocols.
"""
import argparse
from collections import Counter,defaultdict
from contextlib import ExitStack
import os
from pathlib import Path
import subprocess
import shlex
import tempfile
import time
from unittest.mock import patch

from biomolexplorer.docking_data import write_csv,read_compounds
from biomolexplorer.docking_preparation import prepare_candidates
from biomolexplorer.ligand_cache import LigandCache
from biomolexplorer.molecular_core import MolecularPreparationCore,GasteigerEngine
from biomolexplorer.molecular_adapters import AutoDockAdapter,Dock6Adapter
from biomolexplorer.storage import write_json

SMILES=('CCO','C[NH3+]','CC(=O)[O-]','c1ccccc1','c1ccncc1','CCCl','C[C@H](O)F')


def measure(root,count,backend,warm=False,engines='vina'):
    label=backend+'-'+str(count);base=root/label;base.mkdir(exist_ok=True)
    dataset=base/'input.csv'
    write_csv(dataset,[dict(molecule_chembl_id='M'+str(i),canonical_smiles=SMILES[i%len(SMILES)]) for i in range(count)])
    commands=Counter();timings=defaultdict(float);hits=[0];observed=set()
    original_run=subprocess.run
    def launch(command,*args,**kwargs):
        name=Path(shlex.split(command)[0] if isinstance(command,str) else command[0]).name
        commands[name]+=1
        started=time.perf_counter()
        try:return original_run(command,*args,**kwargs)
        finally:
            timings['external_'+name]+=time.perf_counter()-started
            if kwargs.get('cwd'):observed.update(str(p) for p in Path(kwargs['cwd']).rglob('*') if p.is_file())
    restore=LigandCache.restore
    def reuse(*args,**kwargs):
        result=restore(*args,**kwargs);hits[0]+=int(result);return result
    def instrument(cls,name,label):
        method=getattr(cls,name)
        def timed(*args,**kwargs):
            started=time.perf_counter()
            try:return method(*args,**kwargs)
            finally:timings[label]+=time.perf_counter()-started
        return patch.object(cls,name,timed)
    with ExitStack() as stack:
        stack.enter_context(patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(base/'cache/preparations'),BIOMOL_CPU_WORKERS='1'))
        stack.enter_context(patch('subprocess.run',side_effect=launch))
        stack.enter_context(patch.object(LigandCache,'restore',reuse))
        for cls,name,label in [(MolecularPreparationCore,'load','loading'),(MolecularPreparationCore,'prepare','native_prepare_including_charge_and_cache'),
                (GasteigerEngine,'calculate','gasteiger'),(AutoDockAdapter,'export','pdbqt_export'),(Dock6Adapter,'export','mol2_export')]:
            stack.enter_context(instrument(cls,name,label))
        started=time.perf_counter()
        output=prepare_candidates(dataset,base/('warm' if warm else 'cold'),
            dict(preparation_backend=backend,ligand={'minimize':False}),engines=engines)
        elapsed=time.perf_counter()-started
    rows=read_compounds(output)
    return dict(backend=backend,engines=engines,count=count,warm=warm,successful=len(rows),wall_seconds=elapsed,
        molecules_per_second=len(rows)/elapsed,export_cache_hits=hits[0],export_cache_hit_ratio=hits[0]/count,
        external_processes=dict(commands),external_process_count=sum(commands.values()),timings_seconds=dict(timings),
        observed_tool_work_files=len(observed),output_files=sum(p.is_file() for p in output.parent.rglob('*')))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--sizes',type=int,nargs='+',default=[1,10,100,1000])
    parser.add_argument('--legacy-max',type=int,default=10);parser.add_argument('--output',required=True)
    parser.add_argument('--engines',choices=['vina','dock6','both'],default='vina')
    args=parser.parse_args();records=[]
    with tempfile.TemporaryDirectory(prefix='bme-preparation-bench-') as folder:
        for count in args.sizes:
            for backend in ('legacy_chimera','python'):
                if backend=='legacy_chimera' and count>args.legacy_max:continue
                for warm in (False,True):
                    record=measure(Path(folder),count,backend,warm,args.engines);records.append(record)
                    print(backend,count,'warm' if warm else 'cold',round(record['wall_seconds'],4),record['external_process_count'],flush=True)
    write_json(dict(description='Seven representative graphs, repeated IDs; both protocols explicitly disable minimization. Native Gasteiger is not asserted identical to Antechamber Gasteiger. Tool-work file count is a lower bound, not all OS writes.',records=records),Path(args.output))


if __name__=='__main__':main()
