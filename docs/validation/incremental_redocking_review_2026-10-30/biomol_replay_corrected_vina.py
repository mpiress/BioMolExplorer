from pathlib import Path
import json,shutil,time
from caad.docking import Docking
from biomolexplorer.vina_batch import run_resilient
from biomolexplorer.redocking_dock6 import crystal_reference,heavy_atom_rmsd
from biomolexplorer.docking_data import score
import pandas as pd
root=Path('/home/michel/Downloads/Alzheimer/runs/384ffe70b0484994ab021a823aa78efb')
for code,ligand,resnum,folder in [('6ZWI','QRH',604,'/tmp/biomol-graph-regression-6ZWI-exact'),('6QAA','HUN',613,'/tmp/biomol-aromatic-current-regression')]:
 work=Path(folder);prefix=f'{code}_{ligand}_{resnum}A';source=next(root.rglob(code+'_A.dockprep.pdbqt')).parent
 prepared=work/'prepared.mol2' if code=='6ZWI' else work/(prefix+'.lig.mol2')
 dest=work/(prefix+'.lig.mol2')
 if prepared!=dest:
  shutil.copyfile(prepared,dest)
  shutil.copyfile(work/'prepared.pdb',work/(prefix+'.lig.pdb'))
 shutil.copyfile(source/(code+'_A.dockprep.pdbqt'),work/(code+'_A.dockprep.pdbqt'))
 dock=Docking(complex_input_path=str(source.parent),output_path=str(work));dock.prepare_on_obabel(dest.name,prefix+'.lig.pdbqt',[('p','7.4')],input_format='mol2',output_format='pdbqt')
 if (source/'centers.csv').is_file() and prefix in pd.read_csv(source/'centers.csv').columns:center=pd.read_csv(source/'centers.csv')[prefix].tolist()
 else:center=dock.calculate_ligand_centerofmass(prefix+'.lig.pdb',ligand)
 config=work/'test.vina';config.write_text(f'receptor = {code}_A.dockprep.pdbqt\nligand = {prefix}.lig.pdbqt\nout = docked.pdbqt\n'+''.join(f'center_{axis} = {value}\n' for axis,value in zip('xyz',center))+''.join(f'size_{axis} = 24\n' for axis in 'xyz')+'exhaustiveness = 20\nnum_modes = 10\nseed = 73\nverbosity = 0\nscoring = vina\ncpu = 1\n')
 started=time.monotonic();outcome=run_resilient([config],work)
 crystal_reference(dest,source.parent/(code+'.pdb'),[code,ligand,resnum,'A'],work/'crystal.mol2')
 print(json.dumps(dict(pair=prefix,completed=len(outcome['completed']),failed=outcome['failed'],rmsd=heavy_atom_rmsd(work/'crystal.mol2',work/'docked.pdbqt'),score=score(work/'docked.pdbqt','vina'),seconds=time.monotonic()-started)),flush=True)
