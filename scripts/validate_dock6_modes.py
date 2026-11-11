"""Validate DOCK6 modes on one supplied Vina pair, without changing source files."""
import argparse
from contextlib import nullcontext
import json
import os
from pathlib import Path
import sys
import time
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biomolexplorer.docking_data import read_compounds,read_results,write_csv
from caad.docking import Dock6
from wrappers.docking import perform_docking_dock6


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-input',type=Path,required=True)
    parser.add_argument('--target',required=True)
    parser.add_argument('--poses-table',type=Path,required=True)
    parser.add_argument('--receptor-id',required=True)
    parser.add_argument('--compound-id')
    parser.add_argument('--dock6-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cache-dir',type=Path)
    parser.add_argument('--require-existing-grid',action='store_true')
    args=parser.parse_args();root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+str(args.dock6_root.resolve()/'bin')+os.pathsep+os.environ['PATH']
    os.environ['BIOMOL_PREPARATION_CACHE']=str((args.cache_dir or root/'cache').resolve()/'ligand-preparation')
    os.environ['BIOMOL_PREPARATION_WORKERS']='1'
    metadata=args.prepared_input/args.target.replace(' ','')/'pdb_codes.csv'
    record=next(row for row in read_compounds(metadata) if f"{row['PDB_CODE']}_{row['CHAIN']}"==args.receptor_id)
    record=[record[key] for key in ('PDB_CODE','LIGAND','RESNUM','CHAIN')]
    row=next(row for row in read_compounds(args.poses_table) if row.get('receptor_id')==args.receptor_id
             and (not args.compound_id or row['molecule_chembl_id']==args.compound_id))
    pose=Path(row['conformer_file'])
    row['conformer_file']=str((pose if pose.is_absolute() else args.poses_table.parent/pose).resolve())
    bound=root/'bound';bound.mkdir();write_csv(bound/'compounds.csv',[row])
    raw=root/'smiles';raw.mkdir()
    write_csv(raw/'compounds.csv',[{key:row[key] for key in ('molecule_chembl_id','canonical_smiles')}])
    original=Dock6.perform_subprocess;report={}
    for label in ('guided','score_only','restrained','denovo','guided-warm'):
        mode='guided' if label=='guided-warm' else label;calls=[]
        def track(engine,command,cwd=None):
            calls.append(str(command));return original(engine,command,cwd)
        surface=patch.object(Dock6,'prepare_surface',side_effect=AssertionError('Receptor grid was rebuilt')) if args.require_existing_grid else nullcontext()
        preparation=patch('biomolexplorer.processes.run_command',side_effect=AssertionError('Warm preparation invoked a native tool')) if label=='guided-warm' else nullcontext()
        started=time.monotonic()
        with patch.object(Dock6,'perform_subprocess',track),surface,preparation:
            details=perform_docking_dock6(str(args.prepared_input),args.target,str(root/label),
                str(raw if mode=='denovo' else bound),str(args.dock6_root),'gas','compounds',pdb_code=record,dock6_mode=mode)
        rows=read_results(root/label,'dock6')
        if len(rows)!=1:raise AssertionError((label,details))
        if mode=='score_only' and float(rows[0]['input_pose_rmsd'])>.002:raise AssertionError('Reference coordinates changed')
        if label=='guided-warm':
            if calls:raise AssertionError('Warm execution invoked DOCK6')
            first=read_results(root/'guided','dock6')[0]
            if Path(first['conformer_file']).read_bytes()!=Path(rows[0]['conformer_file']).read_bytes():raise AssertionError('Cached pose differs')
        report[label]=dict(seconds=time.monotonic()-started,native_calls=len(calls),row=rows[0],details=details)
    (root/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
