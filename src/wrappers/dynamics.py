"""Batch AMBER refinement with the pipeline's verified per-pair checkpoints."""
import hashlib
import json
from pathlib import Path
import subprocess
import time

from biomolexplorer.md_config import validate_settings,execution_settings
from biomolexplorer.md_cache import PairCache
from biomolexplorer import md_engine
from biomolexplorer.stage_cache import file_digest
from biomolexplorer.storage import write_json,write_text
from biomolexplorer.progress import report_progress
from biomolexplorer.diagnostics import event,get_logger


def perform_dynamics_dock6(base_input_path,base_selected_mols,base_output_path,
        dock6_app_path='',ambertools_path='',ambertools_env='',md_inputs=None,md_seed_path=None,
        selected_pairs=None,excluded_pairs=None,
        amber_score_type='dock6',amber_score_movable_region='ligand',
        amber_score_minimization_rmsgrad=.01,amber_score_before_md_minimization_cycles=100,
        amber_score_md_steps=3000,amber_score_after_md_minimization_cycles=100,
        amber_score_gb_model=5,amber_score_nonbonded_cutoff=18.,amber_score_temperature=300.,
        amber_score_movable_distance_cutoff=10.,amber_score_receptor_movable_atom_expr='::',
        amber_score_complex_movable_atom_expr='::',amber_score_abort_on_unprepped_ligand=True,
        use_existing_ligand_charges=False,md_time_step=.001,trajectory_stride=10,
        md_heating_steps=1000,md_equilibration_steps=1000,
        random_seed=2026,md_movable_mask=':1',analyze_trajectory=True):
    settings=validate_settings(locals())
    pairs=md_inputs or []
    if not pairs:raise ValueError('Nenhum par ativo de receptor e ligante foi selecionado.')
    root,found=md_engine.tools(dock6_app_path,ambertools_path,settings,ambertools_env)
    identity=md_engine.tool_identity(root,found)
    output=Path(base_output_path);output.mkdir(parents=True,exist_ok=True)
    cache=PairCache(output,md_seed_path);records=[]
    logger=get_logger('scientific_tools')
    def checkpoint():
        write_json(dict(version=1,settings=settings,pairs=records),output/'dynamics_results.json')
        import csv,io
        stream=io.StringIO(newline='')
        fields=['id','receptor','ligand','status','reused','amber_score','complex_energy','receptor_contribution','ligand_contribution','error','directory']
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(records)
        write_text(stream.getvalue(),output/'dynamics_results.csv')
    for index,pair in enumerate(pairs,1):
        key=pair['id'];directory=output/'pairs'/key
        record=dict(id=key,receptor=pair['receptor']['name'],ligand=pair['ligand']['name'],
            status='failed',reused=False,directory='pairs/'+key)
        report_progress('Dinâmica: '+str(index)+'/'+str(len(pairs))+' · '+record['receptor']+' / '+record['ligand'],index-1,len(pairs))
        started=time.monotonic()
        try:
            from biomolexplorer.md_inputs import molecular_identity
            ligand_identity=molecular_identity(pair['ligand'])
            request=dict(version=2,settings=execution_settings(settings),software=identity,
                receptor=file_digest(pair['receptor']['path']),ligand=ligand_identity)
            signature=hashlib.sha256(json.dumps(request,sort_keys=True,allow_nan=False).encode()).hexdigest()
            reused=cache.restore(signature,directory)
            if reused:
                # Labels and identities belong to the current selection.
                record.update({k:v for k,v in reused.items() if k not in ('id','receptor','ligand','directory')})
            else:
                record.update(md_engine.run_pair(pair,directory,root,found,settings,ambertools_path,ambertools_env),status='succeeded')
            cache.remember(signature,directory,record)
            event(logger,'dynamics.pair.reused' if reused else 'dynamics.pair.completed',
                'Dinâmica do par concluída.',compound_identifier=pair['ligand']['name'],duration_ms=round((time.monotonic()-started)*1000))
        except (ValueError,RuntimeError,OSError,subprocess.SubprocessError) as error:
            if getattr(error,'error_code',None) in ('TOOL_NOT_FOUND','TOOL_START_FAILED'):raise
            record.update(status='failed',error=str(error))
            progress=directory/'progress.json'
            if progress.is_file():record['failed_stage']=json.loads(progress.read_text()).get('stage')
            logger.warning('Falha na dinâmica do par %s / %s: %s',record['receptor'],record['ligand'],error,exc_info=True)
        records.append(record);checkpoint()
    failures=sum(r['status']=='failed' for r in records)
    report_progress('Dinâmica concluída: '+str(len(records)-failures)+' sucessos; '+str(failures)+' falhas.',len(records),len(records))
    return dict(md_pairs=len(records),md_failures=failures,md_reused=sum(r['reused'] for r in records))
