"""Prepare external candidates once and export engine-specific, reusable inputs."""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .docking_data import read_compounds, write_csv, convert_structure
from .redocking_config import DEFAULTS, validate_preparation_settings
from .visualizations import molecule_sdf


def validate_receptor_outputs(prepared, records, engines):
    """Check the reusable receptor bundle before exporting candidate inputs."""
    from .docking_data import binding_center
    prepared=Path(prepared)
    if not records:
        raise ValueError('A preparação não produziu receptores utilizáveis.')
    suffixes=['.dockprep.pdbqt']  # Stable receptor selector for either engine.
    if engines in ('dock6','both'):
        suffixes += ['.dockprep.mol2','.noH.pdb']
    for record in records:
        receptor=f'{record[0]}_{record[3]}'
        for suffix in suffixes:
            file=prepared/(receptor+suffix)
            if not file.is_file() or not file.stat().st_size:
                raise ValueError('Arquivo do receptor necessário para '+engines+': '+file.name)
        try:
            binding_center(prepared,record)
        except (OSError,KeyError,TypeError,ValueError) as error:
            raise ValueError('Centro do sítio de docking ausente ou inválido para '+receptor+'.') from error


def prepare_candidates(dataset, destination, config=None, ph=7.4, engines='both'):
    if engines not in ('vina','dock6','both'):
        raise ValueError('Selecione Vina, DOCK6 ou ambos para a saída de preparação.')
    config=config or {}
    validate_preparation_settings(config)
    options=dict(DEFAULTS,**config.get('ligand',{}))
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    rows=[]
    for original in read_compounds(dataset):
        row=dict(original);code=row['molecule_chembl_id']
        if not re.fullmatch(r'[A-Za-z0-9_.+-]{1,100}',code) or '..' in code:
            raise ValueError('Código molecular inválido para o docking.')
        pose=row.get('conformer_file') or row.get('pose_file')
        for key in ('prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin'):row.pop(key,None)
        with tempfile.TemporaryDirectory(prefix='bme-prepare-') as folder:
            work=Path(folder)
            if pose:
                source=Path(pose)
                if not source.is_absolute():source=Path(dataset).parent/source
                convert_structure(source,work/'input.mol2',ph=ph if options['add_hydrogens'] else None,index=int(row.get('conformer_index') or 1))
            else:
                (work/'input.sdf').write_text(molecule_sdf(row['canonical_smiles']),encoding='utf-8')
                convert_structure(work/'input.sdf',work/'input.mol2',ph=ph if options['add_hydrogens'] else None)
            # Apply the same role-scoped operations used for reference-ligand preparation.
            commands=['open input.mol2']
            if options['remove_solvent']:commands.append('delete solvent')
            if options['remove_hydrogens']:commands.append('delete element.H')
            if options['add_hydrogens']:commands.append('addh')
            commands.append('addcharge all spec #0 chargeModel 14sb method '+options['charge_type'])
            if options['minimize']:commands.append('minimize spec #0')
            commands += ['write format mol2 #0 prepared.mol2','close all']
            (work/'prepare.com').write_text('\n'.join(commands)+'\n',encoding='utf-8')
            result=subprocess.run(['chimera','--nogui','--silent','prepare.com'],cwd=work,capture_output=True,text=True,timeout=300)
            prepared=work/'prepared.mol2'
            if result.returncode or not prepared.is_file() or not prepared.stat().st_size:
                raise ValueError('Não foi possível preparar o composto para docking: '+code)
            # Export both selected formats from this one conformation; no second minimization.
            formats=('pdbqt','mol2') if engines=='both' else ('pdbqt',) if engines=='vina' else ('mol2',)
            for format in formats:
                output=destination/(code+'.lig.'+format)
                if format=='mol2':shutil.copy2(prepared,output)
                else:convert_structure(prepared,output,ph=ph if options['add_hydrogens'] else None)
                row['prepared_'+format]=str(output.resolve())
            row['conformer_file']=row['prepared_mol2'] if engines!='vina' else row['prepared_pdbqt']
            row.pop('pose_file',None);row.pop('conformer_index',None)
            row['prepared_origin']='library' if pose and row.get('conformer_origin')=='library' else 'pose' if pose else 'smiles'
            row['docking_engines']=engines
            rows.append(row)
    write_csv(destination/'compounds.csv',rows)
    return destination/'compounds.csv'
