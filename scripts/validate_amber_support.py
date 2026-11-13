"""Validate local PTM libraries and optionally an exported MD pair selection.

Run with the BioMolExplorer Python environment. All products are written under
--output; the DOCK6 installation and source structures remain untouched.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from unittest.mock import patch

from biomolexplorer.md_analysis import topology
from biomolexplorer.md_residue_support import bundle
from biomolexplorer.storage import write_json


def validate_library(root,output):
    manifest=bundle();results={}
    for code,model in manifest['residues'].items():
        directory=output/'library'/code;directory.mkdir(parents=True,exist_ok=True)
        for suffix in ('.prepin','.frcmod'):
            name=code+suffix;shutil.copy2(manifest['file_sources'][name],directory/name)
        (directory/'build.leap').write_text('source leaprc.ff03\nloadAmberPrep '+code+'.prepin\nloadAmberParams '+code+'.frcmod\nx = sequence { ALA '+code+' ALA }\nsaveAmberParm x x.prmtop x.inpcrd\nsavePDB x x.pdb\nquit\n')
        command=[str(root/'bin/teLeap')]+['-I'+str(root/'parameters/leap'/name) for name in ('prep','lib','parm','cmd')]+['-f','build.leap']
        native=subprocess.run(command,cwd=directory,capture_output=True,text=True,timeout=60,check=True)
        (directory/'build.log').write_text(native.stdout+native.stderr)
        info=topology(directory/'x.prmtop');index=info['labels'].index(code)
        selected=[i for i,r in enumerate(info['residue']) if r==index]
        charge=float(info['charges'][selected].sum());names={info['names'][i] for i in selected}
        if not math.isfinite(charge) or abs(charge-model['net_charge'])>.001 or names!=set(model['atoms']):
            raise ValueError('PTM library inventory/charge mismatch: '+code)
        results[code]=dict(ok=True,charge=charge,atoms=len(names))
    return results


def validate_pairs(root,output,pairs_path):
    from wrappers.dynamics import perform_dynamics_dock6
    pairs=json.loads(pairs_path.read_text())
    if not isinstance(pairs,list) or not pairs:raise ValueError('Expected a nonempty exported md_inputs array')
    settings=dict(amber_score_before_md_minimization_cycles=100,amber_score_md_steps=1000,amber_score_after_md_minimization_cycles=10000)
    arguments=dict(dock6_app_path=str(root),md_inputs=pairs,**settings)
    started=time.monotonic();result=perform_dynamics_dock6('','',str(output/'execution'),**arguments);seconds=time.monotonic()-started
    if result['md_failures']:raise ValueError('Native validation failed; inspect dynamics_results.json and helper logs')
    started=time.monotonic()
    with patch('biomolexplorer.md_engine.run_pair',side_effect=AssertionError('Expected completed-pair cache reuse')):
        cached=perform_dynamics_dock6('','',str(output/'repeat'),**arguments)
    if cached['md_failures'] or cached['md_reused']!=len(pairs):raise ValueError('Completed-pair cache validation failed')
    return dict(result=result,seconds=seconds,cached=cached,cache_seconds=time.monotonic()-started,settings=settings)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dock6-root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--pairs-json',type=Path)
    args=parser.parse_args();root=args.dock6_root.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    os.environ['BIOMOL_PREPARATION_CACHE']=str(output/'cache/preparations')
    report=dict(library=validate_library(root,output))
    if args.pairs_json:report['pairs']=validate_pairs(root,output,args.pairs_json)
    write_json(report,output/'validation.json')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
