from pathlib import Path
import shutil,json,time,os
from wrappers.redocking import perform_redocking
root=Path('/home/michel/Downloads/Alzheimer/runs/384ffe70b0484994ab021a823aa78efb');work=Path('/tmp/biomol-dock6-current-regression');source=work/'source'/'Target';prepared=source/'Prepared';prepared.mkdir(parents=True,exist_ok=True)
original=next(p for p in root.rglob('4M0E_1YL_604A.lig.mol2') if p.parent.name=='Prepared');old=original.parent
for p in old.iterdir():
 if p.is_file() and p.suffix in ('.pdb','.mol2','.pdbqt','.csv','.json'):shutil.copyfile(p,prepared/p.name)
shutil.copyfile(old.parent/'4M0E.pdb',source/'4M0E.pdb')
os.environ['BIOMOL_CPU_WORKERS']='1'
started=time.monotonic();details=perform_redocking(str(source.parent),'Target',str(work/'results'),[['4M0E','1YL',604,'A',2.0]],prepare_complex=False,redocking_engine='dock6',dock6_app_path='/home/michel/progs/dock6',conformer_search_type='rigid')
print(json.dumps(dict(details=details,metadata=(source/'pdb_codes.csv').read_text(),seconds=time.monotonic()-started)),flush=True)
