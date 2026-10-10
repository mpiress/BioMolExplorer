"""Optional local trained predictions, kept separate from descriptor heuristics."""
from importlib.metadata import version
import math
import csv
import os
import json
import subprocess
import sys
from importlib.util import find_spec
from pathlib import Path


def isolated_python():
    from .paths import SOURCE_ROOT
    configured=os.environ.get('BIOMOL_ADMET_PYTHON')
    if configured:
        path=Path(configured).expanduser()
        return path if path.is_file() else None
    root=SOURCE_ROOT.parent/'.admet-env'
    python=root/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
    return python if python.is_file() and (root/'.ready').is_file() else None


def predictor_available():
    if isolated_python() is not None:return True
    try:return find_spec('admet_ai') is not None and version('admet-ai').startswith('2.')
    except Exception:return False


def predict_admet_ai(smiles):
    python=isolated_python()
    if python is not None and python.resolve()!=Path(sys.executable).resolve():
        from .paths import SOURCE_ROOT
        env=os.environ.copy();env['PYTHONPATH']=str(SOURCE_ROOT)
        result=subprocess.run([str(python),str(Path(__file__).with_name('admet_ai_runner.py'))],
            input=json.dumps(smiles),text=True,capture_output=True,env=env)
        if result.returncode:
            raise RuntimeError('Falha no ambiente ADMET-AI isolado: '+result.stderr[-4000:])
        payload=json.loads(result.stdout)
        return payload['records'],payload['provider']
    return predict_local(smiles)


def predict_local(smiles):
    try:
        from admet_ai import ADMETModel
    except ImportError as exc:
        raise RuntimeError('ADMET-AI não está instalado no ambiente do worker. Instale o extra admet-ml nesse ambiente ou escolha RDKit local.') from exc
    # No DrugBank percentiles (these are not probabilities); no repeated descriptors.
    import admet_ai
    if not version('admet-ai').startswith('2.'):
        raise RuntimeError('A integração requer ADMET-AI 2.x; versões anteriores usam modelos diferentes.')
    model=ADMETModel(include_physchem=False,drugbank_path=None,num_workers=0)
    unique=list(dict.fromkeys(smiles));records={};columns=[]
    # Bounded batches prevent a large library from building one enormous graph dataset.
    for offset in range(0,len(unique),256):
        batch=unique[offset:offset+256]
        predicted=model.predict(smiles=batch)
        if not hasattr(predicted,'index') or list(predicted.index)!=batch or not predicted.columns.is_unique:
            raise ValueError('ADMET-AI retornou moléculas ausentes, fora de ordem ou colunas duplicadas.')
        columns=list(predicted.columns)
        for smile,row in predicted.iterrows():
            values={}
            for key,value in row.items():
                if not math.isfinite(float(value)):
                    raise ValueError('ADMET-AI retornou uma predição não finita.')
                values['ML_'+str(key)]=float(value)
            records[smile]=values
    # Preserve raw endpoint names/units: Pgp_Broccatelli is inhibition, not substrate.
    definitions={}
    metadata=Path(admet_ai.__file__).parent/'resources'/'data'/'admet.csv'
    if metadata.is_file():
        with metadata.open() as stream:
            definitions={row['id']:{k:row[k] for k in ('name','category','task_type','units','species','url')}
                for row in csv.DictReader(stream) if row['id'] in columns}
    return records,dict(endpoint_definitions=definitions,provider='ADMET-AI',version=version('admet-ai'),
        endpoints=columns,domain_assessment='not_assessed',uncertainty='not_provided',
        source='https://github.com/swansonk14/admet_ai',
        interpretation='Raw model endpoint outputs; see TDC definitions and provider version. No generic probability threshold is applied.')
