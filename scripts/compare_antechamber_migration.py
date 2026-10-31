"""Paired charge-method comparison: same explicit-H geometry, no preparation drift."""
import argparse
from pathlib import Path
import re
import shutil
import tempfile
import time

import numpy as np
from rdkit import Chem

from biomolexplorer.molecular_core import MolecularPreparationCore
from biomolexplorer.molecular_adapters import Dock6Adapter,amber_tool
from biomolexplorer.processes import run_command
from biomolexplorer.storage import write_text,write_json

CASES={'neutral':'CCO','cation':'C[NH3+]','anion':'CC(=O)[O-]','heterocycle':'c1ccncc1',
       'aromatic':'c1ccccc1','halogen':'CCCl','multi_hetero':'OCCO','stereocenter':'C[C@H](O)F'}
TOLERANCES=dict(max_charge=.001,mean_charge=.0002,total_charge=.001,
                total_score=.05,electrostatic=.02,vdw=.02,coordinate_rmsd=.002)


def table(path):
    block=path.read_text().split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
    return {r[1]:dict(charge=float(r[8]),type=r[5],xyz=list(map(float,r[2:5]))) for r in [s.split() for s in block.splitlines() if s.strip()]}


def run(root,dock,prefix,charge_only=False):
    core=MolecularPreparationCore(root,prefix);results=[]
    reference,_=core.load(dock/'install/test/grid_score/lig.mol2')
    center=reference.GetConformer().GetPositions().mean(axis=0)
    template=(dock/'install/test/continuous_score/continuous1.dockin').read_text()
    template=template.replace('lig.mol2','ligand.mol2').replace('rec.mol2',str(dock/'install/test/grid_generation/rec.mol2'))
    template=template.replace('../../../parameters/',str(dock/'parameters')+'/')
    for case,smiles in CASES.items():
        initial=core.prepare(Chem.MolFromSmiles(smiles))
        shift=center-initial.coordinates.mean(axis=0)
        for i,xyz in enumerate(initial.coordinates+shift):initial.molecule.GetConformer().SetAtomPosition(i,xyz)
        for method,alias in [('gasteiger','gas'),('am1bcc','am1')]:
            work=root/(case+'-'+method);work.mkdir()
            Dock6Adapter().export(initial,work/'input.mol2')
            write_text('open input.mol2\naddcharge all spec #0 chargeModel 14sb method '+alias+'\nwrite format mol2 #0 legacy.mol2\nclose all\n',work/'legacy.com')
            started=time.perf_counter();run_command(['chimera','--nogui','--silent','legacy.com'],cwd=work)
            legacy_seconds=time.perf_counter()-started
            started=time.perf_counter()
            if charge_only:
                from biomolexplorer.molecular_adapters import charge_existing_mol2
                shutil.copy2(work/'input.mol2',work/'direct.mol2')
                charge_existing_mol2(work/'direct.mol2',method,prefix)
            else:
                direct=core.prepare(initial.molecule,method=method,charge_engine='antechamber',add_hydrogens=False,remove_hydrogens=False)
                Dock6Adapter().export(direct,work/'direct.mol2')
            direct_seconds=time.perf_counter()-started
            old=table(work/'legacy.mol2');new=table(work/'direct.mol2')
            if set(old)!=set(new):raise ValueError('Atom names changed in '+case+' '+method)
            differences=np.array([new[k]['charge']-old[k]['charge'] for k in old])
            rmsd=float(np.sqrt(np.mean([np.sum((np.array(new[k]['xyz'])-old[k]['xyz'])**2) for k in old])))
            scores={}
            for backend in ('legacy','direct'):
                shutil.copy2(work/(backend+'.mol2'),work/'ligand.mol2');write_text(template,work/'score.in')
                run_command([str(dock/'bin/dock6'),'-i','score.in','-o','score.out'],cwd=work)
                text=(work/'continuous1_scored.mol2').read_text()
                scores[backend]={name:float(re.search(name+r':\s*([-+0-9.eE]+)',text)[1]) for name in ('Continuous_Score','Continuous_vdw_energy','Continuous_es_energy')}
            metrics=dict(max_charge=float(np.max(np.abs(differences))),mean_charge=float(np.mean(np.abs(differences))),
                total_charge=float(abs(differences.sum())),coordinate_rmsd=rmsd,
                total_score=abs(scores['direct']['Continuous_Score']-scores['legacy']['Continuous_Score']),
                electrostatic=abs(scores['direct']['Continuous_es_energy']-scores['legacy']['Continuous_es_energy']),
                vdw=abs(scores['direct']['Continuous_vdw_energy']-scores['legacy']['Continuous_vdw_energy']))
            record=dict(case=case,smiles=smiles,method=method,formal_charge=Chem.GetFormalCharge(initial.molecule),
                timings=dict(legacy_loading_charge_export_seconds=legacy_seconds,direct_charge_preparation_seconds=direct_seconds),
                charge_sums={label:sum(a['charge'] for a in rows.values()) for label,rows in [('legacy',old),('direct',new)]},
                atom_types_equal=all(old[k]['type']==new[k]['type'] for k in old),metrics=metrics,scores=scores,
                atoms={k:dict(legacy=old[k],direct=new[k]) for k in old},
                within_tolerances=all(metrics[k]<=limit+1e-10 for k,limit in TOLERANCES.items()) and all(old[k]['type']==new[k]['type'] for k in old))
            results.append(record);print(case,method,record['within_tolerances'],metrics,flush=True)
    return dict(tolerances=TOLERANCES,antechamber=str(amber_tool('antechamber',prefix)),charge_only=charge_only,
        protocol='Same explicit-H geometry, protonation, formal charge and atom names; no minimization. Legacy addcharge vs direct -j 5 -dr n; AM1-BCC explicitly uses AM1.',
        limitation='Fixed-pose DOCK6 continuous scoring, not a docking quality benchmark. Receptor and scoring settings identical.',results=results)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dock6',required=True);parser.add_argument('--output',required=True);parser.add_argument('--charge-only',action='store_true');args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='bme-method-equivalence-') as folder:result=run(Path(folder),Path(args.dock6),'',args.charge_only)
    write_json(result,Path(args.output))
    if not all(row['within_tolerances'] for row in result['results']):raise SystemExit(1)


if __name__=='__main__':main()
