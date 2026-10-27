"""Native, single-system comparison; not evidence of systematic equivalence."""
import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile
import time

import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem

from biomolexplorer.docking_data import write_csv,read_compounds,convert_structure
from biomolexplorer.docking_preparation import prepare_candidates
from biomolexplorer.molecular_core import MolecularPreparationCore
from biomolexplorer.processes import run_command
from biomolexplorer.storage import write_json,write_text


def atoms(text):
    return {int(line[6:11]):(line[77:].strip(),[float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))])
        for line in text.splitlines() if line.startswith(('ATOM  ','HETATM')) and line[77:].strip() not in ('H','HD','HS')}


def validate(root,dock):
    from vina import Vina
    example=dock/'install/test/amber_score_1lgu'
    source=example/'lig.mol2'
    shutil.copy2(example/'1lgu.pdb',root/'rec.pdb')
    write_text('open rec.pdb\naddh\naddcharge all spec #0 chargeModel 14sb method gas\nwrite format mol2 #0 rec.mol2\nwrite format pdb #0 rec.prepared.pdb\nclose all\n',root/'receptor.com')
    started=time.perf_counter();run_command(['chimera','--nogui','--silent','receptor.com'],cwd=root)
    receptor_seconds=time.perf_counter()-started
    run_command(['obabel','-imol2','rec.mol2','-opdbqt','-O','rec.pdbqt','-xr'],cwd=root)
    core=MolecularPreparationCore(root);crystal,_=core.load(source)
    crystal_heavy=np.array([crystal.GetConformer().GetAtomPosition(a.GetIdx()) for a in crystal.GetAtoms() if a.GetAtomicNum()!=1])
    center=crystal_heavy.mean(axis=0).tolist()
    write_csv(root/'input.csv',[dict(molecule_chembl_id='AAH',canonical_smiles=Chem.MolToSmiles(Chem.RemoveHs(crystal)),conformer_file=str(source))])
    template=(dock/'install/test/continuous_score/continuous1.dockin').read_text()
    template=template.replace('lig.mol2','ligand.mol2').replace('rec.mol2',str(root/'rec.mol2'))
    template=template.replace('../../../parameters/',str(dock/'parameters')+'/')
    results=[]
    for backend in ('legacy_chimera','python'):
        started=time.perf_counter()
        table=prepare_candidates(root/'input.csv',root/backend,dict(preparation_backend=backend,ligand={'minimize':False}),engines='both')
        rows=read_compounds(table)
        if len(rows)!=1:raise ValueError('A comparison ligand failed preparation: '+backend)
        row=rows[0];duration=time.perf_counter()-started
        initial=atoms(Path(row['prepared_pdbqt']).read_text())
        xyz=np.array([v[1] for v in initial.values()])
        if len(xyz)!=len(crystal_heavy) or np.max(cKDTree(crystal_heavy).query(xyz)[0])>.005:
            raise ValueError('Preparation moved crystallographic heavy atoms: '+backend)
        session=Vina(sf_name='vina',cpu=1,seed=2026,verbosity=0)
        session.set_receptor(str(root/'rec.pdbqt'));session.compute_vina_maps(center=center,box_size=[24,24,24])
        session.set_ligand_from_string(Path(row['prepared_pdbqt']).read_text())
        session.dock(exhaustiveness=20,n_poses=10)
        poses=session.poses(n_poses=1);final=atoms(poses)
        if set(final)!=set(initial):raise ValueError('Vina changed the heavy-atom mapping.')
        rmsd=float(np.sqrt(np.mean([np.sum((np.array(initial[k][1])-final[k][1])**2) for k in initial])))
        energies=session.energies(n_poses=1).tolist()
        shutil.copy2(row['prepared_mol2'],root/backend/'ligand.mol2')
        write_text(template,root/backend/'score.in')
        run_command([str(dock/'bin/dock6'),'-i','score.in','-o','score.out'],cwd=root/backend)
        score_text=(root/backend/'continuous1_scored.mol2').read_text()
        components={name:float(re.search(name+r':\s*([-+0-9.eE]+)',score_text)[1]) for name in ('Continuous_Score','Continuous_vdw_energy','Continuous_es_energy')}
        score=components['Continuous_Score']
        molecule,_=core.load(row['prepared_mol2'])
        block=Path(row['prepared_mol2']).read_text().split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
        charge_rows=[line.split() for line in block.splitlines() if line.strip()]
        results.append(dict(backend=backend,preparation_seconds=duration,atoms=molecule.GetNumAtoms(),
            heavy_atoms=molecule.GetNumHeavyAtoms(),formal_charge=Chem.GetFormalCharge(molecule),
            partial_charge_sum=sum(float(r[8]) for r in charge_rows),
            per_atom=[dict(name=r[1],type=r[5],charge=float(r[8]),xyz=list(map(float,r[2:5]))) for r in charge_rows],
            vina_score=energies[0][0],redocking_heavy_rmsd_angstrom=rmsd,dock6_continuous_score=score,dock6_components=components,
            preparation_metadata=json.loads(row.get('preparation_metadata') or '{}')))
    def pdb_signature(path):
        lines=[r for r in Path(path).read_text().splitlines() if r.startswith(('ATOM  ','HETATM'))]
        return dict(atoms=len(lines),residues=len({(r[21:27],r[17:20]) for r in lines}),chains=sorted({r[21] for r in lines}))
    return dict(system='DOCK6 amber_score_1lgu/AAH',minimization=False,ph=7.4,
        receptor=dict(backend='unchanged legacy_chimera',seconds=receptor_seconds,before=pdb_signature(root/'rec.pdb'),after=pdb_signature(root/'rec.prepared.pdb')),
        vina=dict(seed=2026,cpu=1,exhaustiveness=20,num_modes=10,box=[24,24,24],rmsd='First scored pose, heavy atoms, no alignment or symmetry correction; prepared coordinates checked against crystallographic input.'),
        dock6='Continuous scoring of the supplied crystallographic pose, not a DOCK6 redocking search or a grid benchmark.',
        results=results,limitations='One small system, no claim of systematic non-degradation, no receptor/cofactor/metal migration validation.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dock6',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='bme-preparation-science-') as folder:
        result=validate(Path(folder),Path(args.dock6))
    write_json(result,Path(args.output));print(json.dumps({r['backend']:{k:r[k] for k in ('vina_score','redocking_heavy_rmsd_angstrom','dock6_continuous_score')} for r in result['results']},indent=2))


if __name__=='__main__':main()
