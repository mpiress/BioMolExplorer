"""Prepare external candidates once and export engine-specific, reusable inputs."""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .docking_data import read_compounds, convert_structure
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


def _prepare_one(original,dataset,destination,options,ph,engines,formats,cache,owned):
    import time
    from .diagnostics import event,get_logger
    logger=get_logger('scientific_tools')
    started=time.monotonic()
    row=dict(original);code=row['molecule_chembl_id']
    if not re.fullmatch(r'[A-Za-z0-9_.+-]{1,100}',code) or '..' in code:
        raise ValueError('Código molecular inválido para o docking.')
    compatible=True
    if row.get('preparation_metadata'):
        import json
        try:metadata=json.loads(row['preparation_metadata'])
        except (TypeError,ValueError):metadata={}
        compatible=isinstance(metadata,dict) and metadata.get('backend') in ('chimera','legacy_chimera')
        if compatible:
            compatible=float(metadata.get('ph',ph))==float(ph) and bool(metadata.get('minimization',options['minimize']))==bool(options['minimize'])
            method={'gasteiger':'gas','am1bcc':'am1'}.get(metadata.get('charge_method'),metadata.get('charge_method'))
            requested={'gasteiger':'gas','am1bcc':'am1'}.get(options['charge_type'],options['charge_type'])
            compatible=compatible and method==requested
    if compatible and row.get('docking_engines') and any(not row.get('prepared_'+fmt) for fmt in formats):
        raise ValueError('A saída de preparação não inclui o formato solicitado pelo docking.')
    if compatible and all(row.get('prepared_'+fmt) for fmt in formats):
        for fmt in formats:
            source=Path(row['prepared_'+fmt])
            if not source.is_absolute():source=Path(dataset).parent/source
            if not source.is_file() or not source.stat().st_size:
                raise ValueError('Arquivo preparado ausente ou vazio: '+source.name)
            from .ligand_cache import validate_prepared_ligand
            validate_prepared_ligand(source,fmt)
            row['prepared_'+fmt]=str(source.resolve())
        event(logger,'ligand.preparation.reused','Reutilizando composto já preparado.',compound_identifier=code)
        return row,True
    owned.extend(destination/(code+'.lig.'+fmt) for fmt in formats)
    pose=row.get('prepared_mol2') or row.get('conformer_file') or row.get('pose_file')
    for key in ('prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin'):row.pop(key,None)
    key=cache.key(row,dataset,options,ph,engines)
    with cache.locked(key):
        reused=cache.restore(key,row,destination,formats)
        if not reused:
            with tempfile.TemporaryDirectory(prefix='bme-prepare-') as folder:
                event(logger,'ligand.preparation.started','Preparando ligante externo (não o receptor).',compound_identifier=code)
                work=Path(folder)
                step_started=time.monotonic()
                if pose:
                    source=Path(pose)
                    if not source.is_absolute():source=Path(dataset).parent/source
                    convert_structure(source,work/'input.mol2',ph=ph if options['add_hydrogens'] else None,index=int(row.get('conformer_index') or 1))
                else:
                    (work/'input.sdf').write_text(molecule_sdf(row['canonical_smiles']),encoding='utf-8')
                    convert_structure(work/'input.sdf',work/'input.mol2',ph=ph if options['add_hydrogens'] else None)
                event(logger,'ligand.preparation.conformer','Conformação inicial disponível.',compound_identifier=code,duration_ms=round((time.monotonic()-step_started)*1000))
                # Apply the same role-scoped operations used for reference-ligand preparation.
                commands=['open input.mol2']
                if options['remove_solvent']:commands.append('delete solvent')
                if options['remove_hydrogens']:commands.append('delete element.H')
                if options['add_hydrogens']:commands.append('addh')
                method={'gasteiger':'gas','am1bcc':'am1'}.get(options['charge_type'],options['charge_type'])
                commands.append('addcharge all spec #0 chargeModel 14sb method '+method)
                if options['minimize']:commands.append('minimize spec #0 prep false')
                commands += ['write format mol2 #0 prepared.mol2','close all']
                (work/'prepare.com').write_text('\n'.join(commands)+'\n',encoding='utf-8')
                step_started=time.monotonic()
                from .processes import run_command
                run_command(['chimera','--nogui','prepare.com'],cwd=work)
                event(logger,'ligand.preparation.chimera','Preparação pelo Chimera finalizada.',compound_identifier=code,duration_ms=round((time.monotonic()-step_started)*1000))
                prepared=work/'prepared.mol2'
                if not prepared.is_file() or not prepared.stat().st_size:
                    raise ValueError('Não foi possível preparar o composto para docking: '+code)
                from .chimera_mol2 import normalize_carboxylates
                normalize_carboxylates(prepared)
                # Export both selected formats from this one conformation; no second minimization.
                for format in formats:
                    output=destination/(code+'.lig.'+format)
                    if format=='mol2':shutil.copy2(prepared,output)
                    else:convert_structure(prepared,output)
                    from .ligand_cache import validate_prepared_ligand
                    validate_prepared_ligand(output,format)
                    row['prepared_'+format]=str(output.resolve())
                import json
                row['preparation_metadata']=json.dumps(dict(backend='legacy_chimera',charge_method={'gas':'gasteiger','am1':'am1bcc'}.get(method,method),ph=ph,minimization=options['minimize']),sort_keys=True)
                cache.save(key,row,formats)
        row['conformer_file']=row['prepared_mol2'] if engines!='vina' else row['prepared_pdbqt']
        row.pop('pose_file',None);row.pop('conformer_index',None)
        row['prepared_origin']='library' if pose and row.get('conformer_origin')=='library' else 'pose' if pose else 'smiles'
        row['docking_engines']=engines
        event(logger,'ligand.preparation.reused' if reused else 'ligand.preparation.completed',
              'Preparação reutilizada do cache verificado.' if reused else 'Preparação do ligante concluída.',
              compound_identifier=code,duration_ms=round((time.monotonic()-started)*1000))
        return row,reused


def prepare_candidates(dataset, destination, config=None, ph=7.4, engines='both'):
    """Prepare independent ligands concurrently; checkpoint every successful result."""
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    from contextvars import copy_context
    from .diagnostics import event,get_logger
    from .ligand_cache import LigandCache
    from .molecule_quality import CURRENT_REPORT as QUALITY, QualityReport
    from .paths import worker_count
    from .preparation_failures import CURRENT_REPORT,PreparationFailures
    from .progress import report_progress
    from .storage import write_text
    import csv
    import io
    import time

    if engines not in ('vina','dock6','both'):
        raise ValueError('Selecione Vina, DOCK6 ou ambos para a saída de preparação.')
    from .redocking_config import chimera_preparation
    config=chimera_preparation(config)
    validate_preparation_settings(config)
    options=dict(DEFAULTS,**config.get('ligand',{}))
    options['preparation_backend']='legacy_chimera'
    backend='legacy_chimera'
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    report=CURRENT_REPORT.get() or PreparationFailures()
    quality=QUALITY.get()
    inputs=read_compounds(dataset);rows={}
    codes=[r.get('molecule_chembl_id') for r in inputs]
    if len(set(codes))!=len(codes):
        raise ValueError('Códigos moleculares duplicados na entrada da preparação.')
    formats=('pdbqt','mol2') if engines=='both' else ('pdbqt',) if engines=='vina' else ('mol2',)
    fields=list(dict.fromkeys([k for row in inputs for k in row]+['molecule_chembl_id','canonical_smiles','prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin','conformer_file','preparation_metadata']))
    cache=LigandCache(destination,backend)
    workers=min(worker_count(),len(inputs) or 1)
    logger=get_logger('scientific_tools')
    event(logger,'ligand.preparation.batch','Preparando ligantes antes do docking.',total=len(inputs),workers=workers)
    report_progress('Preparando ligantes antes do docking (0/'+str(len(inputs))+').',0,len(inputs))

    last_checkpoint=[0.]
    def checkpoint(force=False):
        now=time.monotonic()
        last_checkpoint[0]=now
        stream=io.StringIO(newline='')
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader();writer.writerows(rows[i] for i in sorted(rows))
        write_text(stream.getvalue(),destination/'compounds.csv')
        report.write(destination)

    def prepare(index):
        local=PreparationFailures();local_quality=QualityReport()
        token=QUALITY.set(local_quality)
        result=None;owned=[]
        try:
            with local.compound(dataset,index+2,inputs[index],destination,owned):
                result=_prepare_one(inputs[index],dataset,destination,options,ph,engines,formats,cache,owned)
        finally:QUALITY.reset(token)
        return result,local.records,local_quality.records

    completed=reused=0
    checkpoint()
    # Bound the queue to the configured CPU budget. Each Chimera process uses
    # its own working directory; correlation contexts are copied per task.
    with ThreadPoolExecutor(max_workers=workers,thread_name_prefix='ligand-preparation') as pool:
        pending={};next_index=0
        while next_index<len(inputs) or pending:
            while next_index<len(inputs) and len(pending)<workers:
                index=next_index;next_index+=1
                pending[pool.submit(copy_context().run,prepare,index)]=index
            finished,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in finished:
                index=pending.pop(future)
                try:result,failures,exclusions=future.result()
                except BaseException:
                    checkpoint(force=True);raise
                report.records.extend(failures)
                if quality:
                    for entry in exclusions:
                        quality.exclude(entry['source'],entry['row'],entry['identifiers'],entry['reason'])
                if result:
                    row,was_reused=result;rows[index]=row;reused+=int(was_reused)
                completed+=1
                checkpoint()
                report_progress('Preparação de ligantes: '+str(completed)+'/'+str(len(inputs))+
                                ' avaliados; '+str(len(rows))+' prontos; '+str(reused)+
                                ' reutilizados; '+str(completed-len(rows))+' falhas. Docking iniciará após a preparação.',
                                completed,len(inputs))
    checkpoint(force=True)
    return destination/'compounds.csv'
