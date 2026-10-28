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
    suffixes=['.dockprep.pdbqt'] if engines in ('vina','both') else []
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


def export_receptor_pdbqt(prepared,records):
    """Add only Vina's rigid format to a working copy, preserving charges/geometry."""
    from .processes import run_command
    prepared=Path(prepared)
    for receptor in sorted({f'{r[0]}_{r[3]}' for r in records}):
        destination=prepared/(receptor+'.dockprep.pdbqt')
        if destination.is_file() and destination.stat().st_size:continue
        source=prepared/(receptor+'.dockprep.mol2')
        if not source.is_file():continue  # Normal bundle validation reports the missing format.
        run_command(['obabel','-imol2',str(source),'-opdbqt','-O',str(destination),'-xr'],cwd=str(prepared))


def prepare_candidates(dataset, destination, config=None, ph=7.4, engines='both'):
    if engines not in ('vina','dock6','both'):
        raise ValueError('Selecione Vina, DOCK6 ou ambos para a saída de preparação.')
    config=config or {}
    validate_preparation_settings(config)
    options=dict(DEFAULTS,**config.get('ligand',{}))
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    from .preparation_failures import CURRENT_REPORT,PreparationFailures
    report=CURRENT_REPORT.get() or PreparationFailures()
    inputs=read_compounds(dataset);rows=[]
    from .diagnostics import event,get_logger
    from .progress import report_progress
    logger=get_logger('scientific_tools')
    for row_number,original in enumerate(inputs,2):
        owned=[]
        with report.compound(dataset,row_number,original,destination,owned):
            row=dict(original);code=row['molecule_chembl_id']
            formats=('pdbqt','mol2') if engines=='both' else ('pdbqt',) if engines=='vina' else ('mol2',)
            if all(row.get('prepared_'+format) for format in formats):
                for format in formats:
                    source=Path(row['prepared_'+format])
                    if not source.is_absolute():source=Path(dataset).parent/source
                    if not source.is_file() or not source.stat().st_size:raise ValueError('Arquivo preparado ausente ou vazio: '+source.name)
                    row['prepared_'+format]=str(source.resolve())
                rows.append(row)
                event(logger,'ligand.preparation.reused','Reutilizando composto já preparado.',compound_identifier=code)
                continue
            if not re.fullmatch(r'[A-Za-z0-9_.+-]{1,100}',code) or '..' in code:
                raise ValueError('Código molecular inválido para o docking.')
            owned.extend(destination/(code+'.lig.'+format) for format in formats)
            pose=row.get('conformer_file') or row.get('pose_file')
            for key in ('prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin'):row.pop(key,None)
            with tempfile.TemporaryDirectory(prefix='bme-prepare-') as folder:
                event(logger,'ligand.preparation.started','Preparando ligante externo (não o receptor).',compound_identifier=code)
                report_progress(None,row_number-2,len(inputs),key='preparation.ligand',arguments={'code':code,'index':row_number-1,'total':len(inputs)})
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
                result=subprocess.run(['chimera','--nogui','--silent','prepare.com'],cwd=work,capture_output=True,text=True)
                prepared=work/'prepared.mol2'
                details='\n'.join((getattr(result,'stdout','') or '',getattr(result,'stderr','') or ''))
                hidden_error=re.search(r'Error while sourcing|Traceback \(most recent call last\)|^Error(?:\s|:)',details,re.M)
                if result.returncode or hidden_error or not prepared.is_file() or not prepared.stat().st_size:
                    raise ValueError('Não foi possível preparar o composto para docking: '+code+' · '+details[-2000:])
                # Export both selected formats from this one conformation; no second minimization.
                formats=('pdbqt','mol2') if engines=='both' else ('pdbqt',) if engines=='vina' else ('mol2',)
                for format in formats:
                    output=destination/(code+'.lig.'+format);owned.append(output)
                    if format=='mol2':shutil.copy2(prepared,output)
                    else:convert_structure(prepared,output,ph=ph if options['add_hydrogens'] else None)
                    row['prepared_'+format]=str(output.resolve())
                row['conformer_file']=row['prepared_mol2'] if engines!='vina' else row['prepared_pdbqt']
                row.pop('pose_file',None);row.pop('conformer_index',None)
                row['prepared_origin']='library' if pose and row.get('conformer_origin')=='library' else 'pose' if pose else 'smiles'
                row['docking_engines']=engines
                rows.append(row)
                event(logger,'ligand.preparation.completed','Preparação do ligante concluída.',compound_identifier=code)
    fields=list(dict.fromkeys([k for row in inputs for k in row]+['molecule_chembl_id','canonical_smiles','prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin','conformer_file']))
    write_csv(destination/'compounds.csv',rows,fields)
    report.write(destination)
    return destination/'compounds.csv'
