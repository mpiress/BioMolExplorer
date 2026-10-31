from pathlib import Path
import shutil,json,os,time
from caad.docking import Docking
from biomolexplorer.dock6_validation import validate_mol2
from rdkit import Chem
root=Path('/home/michel/Downloads/Alzheimer/runs/384ffe70b0484994ab021a823aa78efb');original=next(root.rglob('prepare_ligand_6QAA_HUN_613A.com'));work=Path('/tmp/biomol-aromatic-current-regression');work.mkdir(exist_ok=True)
for name in ('6QAA_A.complex.pdb','biomol_ligand_cleanup.py'):shutil.copyfile(original.parent/name,work/name)
script=work/original.name;script.write_text(original.read_text().replace(str(original.parent),str(work)))
os.environ['BIOMOL_CACHE_DIR']=str(work/'cache');os.environ['BIOMOL_PREPARATION_CACHE']=str(work/'cache'/'ligand-preparation')
started=time.monotonic();dock=Docking(complex_input_path=str(work),output_path=str(work));dock.prepare_on_chimera(script.name)
molecule=validate_mol2(work/'6QAA_HUN_613A.lig.mol2')
print(json.dumps(dict(pair='6QAA/HUN/613/A',atoms=molecule.GetNumAtoms(),formal_charge=Chem.GetFormalCharge(molecule),preparation=json.loads((work/'6QAA_HUN_613A.lig.preparation.json').read_text()),seconds=time.monotonic()-started)),flush=True)
