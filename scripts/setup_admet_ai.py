"""Install the optional predictor in an isolated local CPU environment."""
import argparse
from pathlib import Path
import subprocess
import venv

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--directory',type=Path,default=Path(__file__).resolve().parents[1]/'.admet-env')
args=parser.parse_args()
venv.EnvBuilder(with_pip=True).create(args.directory)
python=args.directory/('Scripts/python.exe' if __import__('os').name=='nt' else 'bin/python')
subprocess.run([str(python),'-m','pip','install','torch==2.8.0','--index-url','https://download.pytorch.org/whl/cpu'],check=True)
subprocess.run([str(python),'-m','pip','install','admet-ai==2.0.1'],check=True)
subprocess.run([str(python),'-c','from admet_ai import ADMETModel; model=ADMETModel(include_physchem=False,drugbank_path=None,num_workers=0); print(model.predict(smiles="CCO"))'],check=True)
(args.directory/'.ready').write_text('admet-ai==2.0.1\n')
print('Predictor ready:',python)
