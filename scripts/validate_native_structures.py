"""Real 1LGU paired preparation, short minimization and docking validation."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from biomolexplorer.structure_selection import gemmi_module,read_selected,write_structure
from biomolexplorer.molecular_core import MolecularPreparationCore
from biomolexplorer.processes import run_command
from biomolexplorer.storage import write_text,write_json
from biomolexplorer.redocking_preparation import prepare_pairs
from caad.docking import Docking


def signature(path):
    structure,_=read_selected(path)
    return {(c.name,r.seqid.num,r.seqid.icode.strip(),r.name,a.name):[a.pos.x,a.pos.y,a.pos.z]
            for c in structure[0] for r in c for a in r if not a.element.is_hydrogen}


def clash_count(path):
    from scipy.spatial import cKDTree
    xyz=np.array(list(signature(path).values()))
    return len(cKDTree(xyz).query_pairs(.8))


def run(root,dock6,amber):
    gemmi=gemmi_module();example=dock6/'install/test/amber_score_1lgu'
    source=root/'source';source.mkdir();(source/'chemistry').mkdir()
    shutil.copy2(example/'lig.mol2',source/'chemistry/AAH.mol2')
    ligand,_=MolecularPreparationCore(root,amber).load(example/'lig.mol2')
    for atom in ligand.GetAtoms():
        info=atom.GetPDBResidueInfo();info.SetResidueNumber(300);info.SetChainId('A')
    structure=gemmi.read_structure(str(example/'1lgu.pdb'))
    ligand_structure=gemmi.read_pdb_string(Chem.MolToPDBBlock(ligand))
    structure[0][0].add_residue(ligand_structure[0][0][0])
    write_structure(structure,source/'1LGU.pdb')
    record=['1LGU','AAH',300,'A'];config={'ambertools_path':amber,'receptor':{'minimize':False},'ligand':{'minimize':False}}
    native=source/'Prepared';native.mkdir();legacy=root/'legacy';legacy.mkdir()
    # Selection reference only; no protonation/minimization differences here.
    write_text('open '+str(source/'1LGU.pdb')+'\nselect #0:.A\nselect invert\ndelete selected\nwrite format pdb #0 complex.pdb\nclose all\n',legacy/'select.com')
    run_command(['chimera','--nogui','--silent','select.com'],cwd=legacy)
    old=signature(legacy/'complex.pdb');new=signature(source/'1LGU.pdb')
    if old!=new:raise ValueError('Native/legacy selected heavy atoms, identities or coordinates differ.')
    write_text('open '+str(example/'1lgu.pdb')+'\naddh\naddcharge all spec #0 chargeModel 14sb method gas\nwrite format mol2 #0 rec.mol2\nwrite format pdb #0 rec.pdb\nclose all\n',legacy/'rec.com')
    started=time.perf_counter();run_command(['chimera','--nogui','--silent','rec.com'],cwd=legacy)
    legacy_receptor_seconds=time.perf_counter()-started
    from biomolexplorer.native_preparation import receptor_pdbqt
    receptor_pdbqt(legacy/'rec.mol2',legacy/'rec.pdbqt')
    calls=[]
    import biomolexplorer.processes as processes
    actual=processes.subprocess.run
    def traced(argv,*args,**kwargs):
        calls.append(Path(str(argv[0])).name);return actual(argv,*args,**kwargs)
    timings={}
    warm_source=root/'warm-source';warm_source.mkdir();shutil.copy2(source/'1LGU.pdb',warm_source/'1LGU.pdb')
    shutil.copytree(source/'chemistry',warm_source/'chemistry')
    for label,destination in [('cold',native),('warm',warm_source/'Prepared')]:
        destination.mkdir(exist_ok=True)
        instance=Docking(complex_input_path=str(destination.parent),output_path=str(destination))
        before=len(calls);started=time.perf_counter()
        with patch('biomolexplorer.processes.subprocess.run',side_effect=traced):
            prepared,failures=prepare_pairs(instance,[record],{'1LGU|AAH|300|A':config},'gas',7.4)
        if failures or not prepared:raise ValueError(str(failures))
        timings[label]=dict(seconds=time.perf_counter()-started,external_calls=calls[before:])
    if any('chimera' in name for name in calls) or timings['warm']['external_calls']:
        raise ValueError('Native preparation called Chimera or warm reuse repeated a calculation.')
    legacy_full=root/'legacy-full';legacy_full.mkdir()
    legacy_config=dict(config,preparation_backend='legacy_chimera')
    instance=Docking(complex_input_path=str(source),output_path=str(legacy_full))
    started=time.perf_counter();instance.prepare_for_docking([record],'gas',7.4,True,
        {'1LGU|AAH|300|A':legacy_config},_prepared_receptors=set())
    legacy_full_seconds=time.perf_counter()-started
    legacy_min=root/'legacy-min';legacy_min.mkdir()
    instance=Docking(complex_input_path=str(source),output_path=str(legacy_min))
    instance.prepare_for_docking([record],'gas',7.4,True,
        {'1LGU|AAH|300|A':{'preparation_backend':'legacy_chimera'}},_prepared_receptors=set())
    receptor_report=json.loads((native/'1LGU_A.preparation.json').read_text())
    ligand_report=json.loads((native/'1LGU_AAH_300A.preparation.json').read_text())
    old_atoms=signature(legacy/'rec.pdb');new_atoms=signature(native/'1LGU_A.prepared.pdb')
    # Original HIS labels become explicit HIE/HID/HIP; compare immutable IDs.
    old_map={(c,n,i,a):v for (c,n,i,r,a),v in old_atoms.items()}
    new_map={(c,n,i,a):v for (c,n,i,r,a),v in new_atoms.items()}
    shared=old_map.keys() & new_map.keys()
    if max(np.linalg.norm(np.array(old_map[k])-new_map[k]) for k in shared)>.002:
        raise ValueError('Receptor preparation moved original heavy atoms before minimization.')
    # Fixed pose scores plus actual Vina redocking. Receptor changes are visible,
    # not hidden behind a common legacy receptor in both arms.
    from vina import Vina
    from biomolexplorer.amber_preparation import minimize_topology,minimize_ligand,export_receptor
    from biomolexplorer.structure_selection import subset,protein_residue
    from biomolexplorer.molecular_adapters import Dock6Adapter,AutoDockAdapter
    rec_work=native/'amber/1LGU_A'
    parameters={name:rec_work/name for name in ('receptor.prmtop','receptor.inpcrd','receptor.amber.pdb','receptor.parameters.json')}
    xyz,min_report=minimize_topology(parameters['receptor.prmtop'],parameters['receptor.inpcrd'],rec_work,amber)
    prepared=MolecularPreparationCore(root,amber).prepare(ligand,charge_engine='antechamber')
    minimized=minimize_ligand(prepared,root/'ligand-min',amber)
    minimization_cases={}
    for label,smiles in [('neutral','CCO'),('cation','C[NH3+]'),('anion','CC(=O)[O-]'),('heterocycle','c1ccncc1')]:
        small=MolecularPreparationCore(root,amber).prepare(Chem.MolFromSmiles(smiles),charge_engine='antechamber')
        small=minimize_ligand(small,root/('min-'+label),amber)
        minimization_cases[label]=small.metadata['minimization']
    minimized_dir=root/'minimized';minimized_dir.mkdir()
    receptor_structure=subset(structure,lambda c,r:protein_residue(r))
    export_receptor(parameters,receptor_structure,minimized_dir/'rec.mol2',minimized_dir/'rec.pdb',
                    minimize=True,prefix=amber,directory=rec_work)
    receptor_pdbqt(minimized_dir/'rec.mol2',minimized_dir/'rec.pdbqt')
    Dock6Adapter().export(minimized,minimized_dir/'lig.mol2');AutoDockAdapter().export(minimized,minimized_dir/'lig.pdbqt')
    poses={};center=ligand.GetConformer().GetPositions()[:9].mean(axis=0).tolist()
    continuous=(dock6/'install/test/continuous_score/continuous1.dockin').read_text()
    for label,rec_mol2,rec_pdbqt,lig_mol2,lig_pdbqt in [('legacy',legacy/'rec.mol2',legacy/'rec.pdbqt',native/'1LGU_AAH_300A.lig.mol2',native/'1LGU_AAH_300A.lig.pdbqt'),
        ('native',native/'1LGU_A.dockprep.mol2',native/'1LGU_A.dockprep.pdbqt',native/'1LGU_AAH_300A.lig.mol2',native/'1LGU_AAH_300A.lig.pdbqt'),
        ('legacy_minimized',legacy_min/'1LGU_A.dockprep.mol2',legacy_min/'1LGU_A.dockprep.pdbqt',legacy_min/'1LGU_AAH_300A.lig.mol2',legacy_min/'1LGU_AAH_300A.lig.pdbqt'),
        ('native_minimized',minimized_dir/'rec.mol2',minimized_dir/'rec.pdbqt',minimized_dir/'lig.mol2',minimized_dir/'lig.pdbqt')]:
        work=root/('scores-'+label);work.mkdir()
        shutil.copy2(lig_mol2,work/'ligand.mol2')
        template=continuous.replace('lig.mol2','ligand.mol2').replace('rec.mol2',str(rec_mol2))
        template=template.replace('../../../parameters/',str(dock6/'parameters')+'/')
        write_text(template,work/'score.in');run_command([str(dock6/'bin/dock6'),'-i','score.in','-o','score.out'],cwd=work)
        text=(work/'continuous1_scored.mol2').read_text()
        components={name:float(re.search(name+r':\s*([-+0-9.eE]+)',text)[1]) for name in ('Continuous_Score','Continuous_vdw_energy','Continuous_es_energy')}
        session=Vina(sf_name='vina',cpu=1,seed=2026,verbosity=0);session.set_receptor(str(rec_pdbqt))
        session.compute_vina_maps(center=center,box_size=[24,24,24]);session.set_ligand_from_file(str(lig_pdbqt))
        session.dock(exhaustiveness=20,n_poses=10)
        input_atoms={r[12:16].strip():np.array([float(r[k:k+8]) for k in (30,38,46)]) for r in (native/'1LGU_AAH_300A.lig.pdbqt').read_text().splitlines() if r.startswith(('ATOM  ','HETATM')) and r.split()[-1]!='HD'}
        docked={r[12:16].strip():np.array([float(r[k:k+8]) for k in (30,38,46)]) for r in session.poses(n_poses=1).splitlines() if r.startswith(('ATOM  ','HETATM')) and r.split()[-1]!='HD'}
        rmsd=float(np.sqrt(np.mean([np.sum((input_atoms[k]-docked[k])**2) for k in input_atoms])))
        poses[label]=dict(vina_score=float(session.energies(n_poses=1)[0][0]),heavy_pose_rmsd=rmsd,dock6=components)
    # Native SES -> spheres -> box -> grid -> DOCK6 flexible redocking.
    from caad.docking import Dock6
    ligand_dir=root/'dock6-ligands';ligand_dir.mkdir();shutil.copy2(native/'1LGU_AAH_300A.lig.mol2',ligand_dir/'AAH.lig.mol2')
    dock6_redocking={}
    min_prepared=root/'minimized-prepared';min_prepared.mkdir()
    shutil.copy2(minimized_dir/'rec.mol2',min_prepared/'1LGU_A.dockprep.mol2')
    min_structure,_=read_selected(minimized_dir/'rec.pdb');write_structure(min_structure,min_prepared/'1LGU_A.noH.pdb',remove_hydrogens=True)
    min_ligands=root/'minimized-ligands';min_ligands.mkdir();shutil.copy2(minimized_dir/'lig.mol2',min_ligands/'AAH.lig.mol2')
    for label,receptor_dir,ligands in [('legacy',legacy_full,ligand_dir),('native',native,ligand_dir),('native_minimized',min_prepared,min_ligands)]:
        destination=root/('dock6-redocking-'+label)
        engine=Dock6(dock6_path=str(dock6),ligand_input_path=str(ligands),receptor_input_path=str(receptor_dir),
                     base_output_path=str(destination),pdb_code='1LGU_A',binding_site_center=center)
        try:
            engine.validate_setup();engine.prepare_site();engine.prepare_minimization();engine.perform_dock6_evaluation()
        finally:engine.export_results()
        scored=list((destination/'flex').glob('*_scored.mol2'))
        if not scored:raise ValueError('DOCK6 flexible redocking produced no scored poses.')
        text=scored[0].read_text()
        dock6_redocking[label]=dict(scored_files=[str(path.relative_to(destination)) for path in scored],
            grid_score=float(re.search(r'Grid_Score:\s*([-+0-9.eE]+)',text)[1]),
            dock6_reported_rmsd=float(re.search(r'HA_RMSDs:\s*([-+0-9.eE]+)',text)[1]))
    return dict(system='DOCK6 1LGU/AAH; authoritative local MOL2 chemistry; receptor + ligand recombined without altering coordinates',
        selection=dict(heavy_atoms=len(new),same_atoms_identifiers_coordinates=True),
        receptor=dict(legacy_heavy=len(old_atoms),native_heavy=len(new_atoms),shared_heavy=len(shared),
                      native_net_charge=receptor_report['parameterization']['net_charge'],
                      added_atoms=receptor_report['parameterization']['added_atoms'],
                      protonation_states=receptor_report['parameterization']['states']),
        ligand=dict(extraction_heavy_rmsd=ligand_report['extraction_heavy_rmsd'],atom_mapping=ligand_report['original_atom_mapping']),
        docking=poses,dock6_flexible_redocking=dock6_redocking,
        minimization=dict(receptor=json.loads(min_report.read_text()),ligand=minimized.metadata['minimization'],representative_ligands=minimization_cases),
        benchmark=dict(legacy_receptor_only_seconds=legacy_receptor_seconds,legacy_full_preparation_seconds=legacy_full_seconds,native_full_preparation=timings),
        limitations='One small system. Paired continuous DOCK6 scores isolate receptor preparation; the separate flexible run checks native SES/spheres/grid/docking. Vina RMSD is measured against the original unminimized crystal pose, first scored pose, without alignment/symmetry. Both full-preparation timing arms disable minimization and emit both formats/centers. Short vacuum minimization is not equilibration. No assertion of systematic non-regression from one target.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dock6',required=True);parser.add_argument('--ambertools',required=True);parser.add_argument('--output',required=True);parser.add_argument('--work-dir')
    args=parser.parse_args()
    if args.work_dir:
        root=Path(args.work_dir);root.mkdir(parents=True,exist_ok=True)
        if any(root.iterdir()):raise ValueError('Validation work directory must be empty.')
        os.environ['BIOMOL_PREPARATION_CACHE']=str(root/'shared/preparations')
        data=run(root,Path(args.dock6).resolve(),str(Path(args.ambertools).resolve()))
        write_json(data,Path(args.output));print(json.dumps({k:v for k,v in data.items() if k in ('selection','docking','benchmark')},indent=2));return
    with tempfile.TemporaryDirectory(prefix='bme-native-validation-') as folder:
        root=Path(folder);os.environ['BIOMOL_PREPARATION_CACHE']=str(root/'shared/preparations')
        data=run(root,Path(args.dock6).resolve(),str(Path(args.ambertools).resolve()))
        write_json(data,Path(args.output));print(json.dumps({k:v for k,v in data.items() if k in ('selection','docking','benchmark')},indent=2))


if __name__=='__main__':main()
