from pathlib import Path
import shutil,json,os,time
from biomolexplorer.vina_batch import run_resilient
from biomolexplorer.redocking_dock6 import crystal_reference,heavy_atom_rmsd
from biomolexplorer.docking_data import score
root=Path('/home/michel/Downloads/Alzheimer/runs/384ffe70b0484994ab021a823aa78efb')
work=Path('/tmp/biomol-vina-current-regression');work.mkdir(exist_ok=True)
reference=next(p for p in root.rglob('4M0E_1YL_604A.lig.pdbqt') if p.parent.name=='Prepared');source=reference.parent
for name in ('4M0E_1YL_604A.lig.pdbqt','4M0E_A.dockprep.pdbqt','4M0E_1YL_604A.lig.mol2','centers.csv'):shutil.copyfile(source/name,work/name)
import pandas as pd
center=pd.read_csv(work/'centers.csv')['4M0E_1YL_604A'].tolist()
config=work/'4M0E.vina';config.write_text('receptor = 4M0E_A.dockprep.pdbqt\nligand = 4M0E_1YL_604A.lig.pdbqt\nout = docked.pdbqt\n'+''.join(f'center_{axis} = {value}\n' for axis,value in zip('xyz',center))+''.join(f'size_{axis} = 24\n' for axis in 'xyz')+'exhaustiveness = 20\nnum_modes = 10\nseed = 73\nverbosity = 0\nscoring = vina\ncpu = 1\n')
started=time.monotonic();outcome=run_resilient([config],work)
raw=next(p for p in source.parent.glob('4M0E.pdb'))
crystal_reference(work/'4M0E_1YL_604A.lig.mol2',raw,['4M0E','1YL',604,'A'],work/'crystal.mol2')
print(json.dumps(dict(outcome=outcome,rmsd=heavy_atom_rmsd(work/'crystal.mol2',work/'docked.pdbqt'),score=score(work/'docked.pdbqt','vina'),seconds=time.monotonic()-started),default=str),flush=True)
