"""Run real engines independently, exchange poses, and validate their consensus.

Uses the charged 1ABE receptor/ligand distributed with DOCK6's
ligand_sampling_demo. No provider requests or modifications to installed samples.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biomolexplorer.docking_data import (
    structure_smiles, write_csv, read_results,
    input_records, copy_input_records,
)
from wrappers.docking import perform_docking_vina, perform_docking_dock6, generate_consensus


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dock6-root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--prepare-inputs',action='store_true',help='Reuse the prepared receptor and prepare both candidate formats in one block')
    args=parser.parse_args();root=args.output.resolve();dock=args.dock6_root.resolve()
    if root.exists():raise ValueError('Use um diretório de saída novo para a validação.')
    os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+str(dock/'bin')+os.pathsep+os.environ['PATH']
    sample=dock/'tutorials/ligand_sampling_demo/1_struct'
    prepared=root/'input'/'Validation'/'Prepared';prepared.mkdir(parents=True)
    receptor=prepared/'1ABE_A.dockprep.mol2';shutil.copy2(sample/'rec_charged.mol2',receptor)
    # Receptors must be rigid in Vina.
    import subprocess
    subprocess.run(['obabel','-imol2',str(receptor),'-opdbqt','-O',str(prepared/'1ABE_A.dockprep.pdbqt'),'-xr'],check=True,capture_output=True)
    subprocess.run(['obabel','-imol2',str(receptor),'-opdb','-O',str(prepared/'1ABE_A.noH.pdb'),'-d'],check=True,capture_output=True)
    ligand=root/'reference.mol2';shutil.copy2(sample/'lig_charged.mol2',ligand)
    # Site center from the actual tutorial reference ligand.
    atoms=[];active=False
    for line in ligand.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):active=line=='@<TRIPOS>ATOM';continue
        if active and line.strip():atoms.append(list(map(float,line.split()[2:5])))
    center=[sum(a[k] for a in atoms)/len(atoms) for k in range(3)]
    (prepared/'centers.csv').write_text('1ABE_ARA_1A\n'+'\n'.join(map(str,center))+'\n')
    (prepared.parent/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABE,ARA,1,A\n')
    dataset=root/'compounds';dataset.mkdir()
    write_csv(dataset/'compounds.csv',[
        dict(molecule_chembl_id='REF',canonical_smiles=structure_smiles(ligand),conformer_file=str(ligand)),
        dict(molecule_chembl_id='ETHANOL',canonical_smiles='CCO',conformer_file=''),
    ])
    record=['1ABE','ARA',1,'A'];report={}
    receptor_input=root/'input'
    if args.prepare_inputs:
        from wrappers.redocking import prepare_structures
        prepare_structures(str(receptor_input),'Validation',str(root/'preparation'),
            base_selected_mols=str(dataset),mol_filename='compounds',receptor_prepared=True,docking_engines='both')
        receptor_input=root/'preparation';dataset=receptor_input/'Validation/Compounds'
        from biomolexplorer.docking_data import read_compounds
        candidates=read_compounds(dataset/'compounds.csv')
        assert {r['molecule_chembl_id'] for r in candidates}=={'REF','ETHANOL'}
        assert all(Path(r['prepared_pdbqt']).is_file() and Path(r['prepared_mol2']).is_file() for r in candidates)
        for name in ('1ABE_A.dockprep.mol2','1ABE_A.dockprep.pdbqt','1ABE_A.noH.pdb','centers.csv'):
            assert (receptor_input/'Validation/Prepared'/name).read_bytes()==(prepared/name).read_bytes()
        report['preparation']={'engines':'both','compounds':['REF','ETHANOL'],'receptor_reused':True}

    def run(engine,source,label):
        output=root/label
        print('Running '+label,flush=True)
        common=dict(base_input_path=str(receptor_input),target='Validation',base_output_path=str(output),
                    base_selected_mols=str(source),mol_filename='compounds',pdb_code=record)
        if engine=='vina':perform_docking_vina(**common,exhaustiveness=1,num_modes=2,sizeof_box=[18]*3)
        else:perform_docking_dock6(**common,dock6_app_path=str(dock),charge_type='gas',conformer_search_type='rigid',distance=6.)
        rows=read_results(output,engine)
        assert {r['molecule_chembl_id'] for r in rows}=={'REF','ETHANOL'},rows
        assert all(Path(r['conformer_file']).is_file() for r in rows)
        report[label]=[{k:r[k] for k in ('molecule_chembl_id','receptor_id','score')} for r in rows]
        return output

    vina=run('vina',dataset,'vina_independent')
    dock6=run('dock6',dataset,'dock6_independent')
    for engine,source,label in [('dock6',vina,'vina_to_dock6'),('vina',dock6,'dock6_to_vina')]:
        available=list(source.rglob('*'));available=[p for p in available if p.is_file()]
        rows=input_records([source/'docking_results.csv'],available)
        copied=copy_input_records(rows,root/(label+'_input'))
        run(engine,copied.parent,label)
    for label,v,d in [('independent',vina,dock6),('vina_to_dock6',vina,root/'vina_to_dock6'),
                       ('dock6_to_vina',root/'dock6_to_vina',dock6)]:
        frame=generate_consensus(str(root),str(root/('consensus_'+label)),label,
            base_vina_path=str(v),base_dock6_path=str(d))
        assert set(frame['molecule_chembl_id'])=={'REF','ETHANOL'}
        report['consensus_'+label]=frame[['molecule_chembl_id','vina','dock6','z-score','min-max']].to_dict('records')
    result=generate_consensus(str(root),str(root/'no_intersection'),'empty',base_vina_path=str(vina),base_dock6_path=str(root/'missing'))
    assert result['rows']==0 and result['skipped_reason']
    report['no_intersection']=result
    (root/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
