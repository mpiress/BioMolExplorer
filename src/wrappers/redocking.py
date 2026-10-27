from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Redocking analysis",

    module_description=(
    "Wrapper module for managing and integrating redocking "
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from typing import Optional, List, Tuple, Literal
from pathlib import Path
from pandas import DataFrame, isna
import os

from caad.docking import DockVina, Docking
from kernel.loggers import LoggerManager
from kernel.utilities import fileHandling
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
logger     = LoggerManager.get_logger('wrapper_docking', log_file='logs/docking.log')
ChargeType = Literal['gas', 'am1']
#----------------------------------------------------------------------------------------------


def perform_redocking(base_input_path:str, target:str, base_output_path:str, pdb_codes:Optional[List[Tuple[str, str, str, str]]]=None,
                      pH:Optional[float]=7.4, sizeof_box:Optional[List]=[24,24,24], exhaustiveness:Optional[int]=20,
                      num_modes:Optional[int]=10, prepare_complex:Optional[bool]=True, charge_type:Optional[ChargeType]='gas', preparation_pairs:Optional[dict]=None,
                      redocking_engine:str='vina', dock6_app_path:str='', conformer_search_type:str='flex') -> None:

    logger = LoggerManager.get_logger('wrapper_redocking', log_file='logs/redocking.log')

    try:

        base_prepared_complexes = f'{base_input_path}/{target.replace(' ','')}/Prepared/'
        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}/'
        path                    = directory(base_input_path)

        from biomolexplorer.redocking_config import validate_structure_pairs, validate_redocking_tools
        records = [list(r) for r in pdb_codes.to_records(index=False)] if isinstance(pdb_codes, DataFrame) else pdb_codes
        from biomolexplorer.redocking_config import chimera_pairs
        preparation_pairs=chimera_pairs(records or [],preparation_pairs)
        validate_structure_pairs(base_input_path, records, preparation_pairs or {}, prepared=not prepare_complex)
        if redocking_engine not in ('vina','dock6'):raise ValueError('Motor de redocking inválido.')
        if conformer_search_type not in ('flex','rigid'):raise ValueError('Busca conformacional inválida.')
        validate_redocking_tools(prepare_complex,preparation_pairs=preparation_pairs,engine=redocking_engine)
        if redocking_engine=='dock6':
            from biomolexplorer.docking_tools import dock6_root
            root=dock6_root(dock6_app_path or None)
            if root is None:raise ValueError('Instalação DOCK6 não encontrada. Configure a pasta bin/ e parameters/.')
            dock6_app_path=str(root)
        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)

        if pdb_codes is None:
            raise ValueError('Selecione pelo menos um par receptor / ligante e sua cadeia na aba Input Data.')
        from biomolexplorer.redocking_config import validate_pairs, pair_key
        validate_pairs([list(record) for record in pdb_codes.to_records(index=False)] if isinstance(pdb_codes, DataFrame) else pdb_codes, preparation_pairs or {})
        columns=['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION']
        if isinstance(pdb_codes,DataFrame):
            pdb_codes=pdb_codes.reindex(columns=columns).copy()
        else:
            pdb_codes=DataFrame([list(record[:4])+[record[4] if len(record)>4 else None] for record in pdb_codes],columns=columns)


        if f1.isFile('pdb_codes')[0]:
            metadata = f1.csv_to_dataframe('pdb_codes')
            if 'RESOLUTION' in metadata:
                resolutions = {(str(r.PDB_CODE), str(r.LIGAND), str(r.RESNUM), str(r.CHAIN)): r.RESOLUTION for r in metadata.itertuples()}
                for idx, row in pdb_codes.iterrows():
                    if row['RESOLUTION'] is None or isna(row['RESOLUTION']):
                        key = (str(row['PDB_CODE']), str(row['LIGAND']), str(row['RESNUM']), str(row['CHAIN']))
                        pdb_codes.at[idx, 'RESOLUTION'] = resolutions.get(key)

        selected=pdb_codes.copy();failures=[]
        from biomolexplorer.redocking_cache import PairCache
        cache=PairCache(base_input_path,base_output_path,
            [list(r) for r in selected.to_records(index=False)],preparation_pairs or {},
            dict(pH=pH,sizeof_box=sizeof_box,exhaustiveness=exhaustiveness,num_modes=num_modes,
                 prepare_complex=prepare_complex,charge_type=charge_type,
                 redocking_engine=redocking_engine,dock6_app_path=dock6_app_path,conformer_search_type=conformer_search_type))
        ready=cache.restore()
        # Seeds are transport files, not additional result datasets.
        import shutil
        for seed in (Path(base_input_path)/'cache').glob('redocking-*'):
            if seed.is_dir():shutil.rmtree(seed)
        pdb_codes=pdb_codes[[pair_key(list(r)) not in cache.restored for r in pdb_codes.to_records(index=False)]]
        if prepare_complex and not pdb_codes.empty:
            dock = Docking(complex_input_path=base_input_path, output_path=base_prepared_complexes)
            from biomolexplorer.redocking_preparation import prepare_pairs,write_failures
            pdb_codes,failures=prepare_pairs(dock,[list(r) for r in pdb_codes.to_records(index=False)],
                preparation_pairs or {},charge_type,pH,ready=ready)
            write_failures(failures,base_output_path)
            pdb_codes = DataFrame(pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION'])
            f1.dataframe_to_csv('pdb_codes', pdb_codes)

        if not pdb_codes.empty:
            keys={pair_key(list(r)) for r in pdb_codes.to_records(index=False)}
            settings={k:v for k,v in (preparation_pairs or {}).items() if k in keys}
            validate_structure_pairs(base_input_path, [list(r) for r in pdb_codes.to_records(index=False)], settings, prepared=True)
            vina = DockVina(ligand_input_path=base_prepared_complexes, receptor_input_path=base_prepared_complexes,
                        output_path=base_output_path, complex_input_path=base_input_path,
                        pdb_codes=pdb_codes, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box,
                        exhaustiveness=exhaustiveness, num_modes=num_modes)
            checkpoint_rows=[]
            def checkpoint(record,rmsd):
                match=next((r for r in selected.to_records(index=False)
                            if pair_key(list(r))==pair_key(record)),record+[None])
                row=dict(zip(columns,match));row.update(RMSD=float(rmsd),STATUS='succeeded',FAILED_STEP='',ERROR='')
                checkpoint_rows.append(row)
                cache.publish(checkpoint_rows)
            if redocking_engine=='vina':
                vina.pair_checkpoint=checkpoint
                vina.redocking(pH=pH)
            else:
                from biomolexplorer.dock6_redocking import redock
                redock([list(r) for r in pdb_codes.to_records(index=False)],
                    base_prepared_complexes,base_output_path,dock6_app_path,conformer_search_type,checkpoint,settings)

        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)

        completed=f1.csv_to_dataframe('pdb_codes') if not pdb_codes.empty else DataFrame(columns=columns+['RMSD'])
        outcomes={pair_key(list(r)[:4]):r for r in completed.to_records(index=False)}
        failed={pair_key([r[k] for k in columns[:4]]):r for r in failures}
        rows=[]
        for record in selected.to_records(index=False):
            key=pair_key(list(record));row=dict(zip(columns,record))
            failure=failed.get(key)
            row.update(RMSD=None,STATUS=failure['STATUS'] if failure else 'succeeded',
                       FAILED_STEP=failure['FAILED_STEP'] if failure else '',ERROR=failure['ERROR'] if failure else '')
            if not failure and key in cache.restored:
                row['RMSD']=cache.restored[key]['row']['RMSD']
            elif not failure and key in outcomes:
                row['RMSD']=outcomes[key]['RMSD']
            rows.append(row)
        f1.dataframe_to_csv('pdb_codes',DataFrame(rows))
        from biomolexplorer.redocking_preparation import write_failures
        write_failures(failures,base_output_path)
        cache.publish(rows)
        from biomolexplorer.progress import report_progress
        report_progress(f'Redocking concluído: {len(rows)-len(failures)} pares preparados; {len(failures)} falhas de preparação.')
        details=dict(redocking_failures=len(failures),redocking_succeeded=len(rows)-len(failures))
        if cache.restored:details['redocking_reused']=len(cache.restored)
        return details


    except Exception as e:
        logger.error(f'Error during to perform the {target} in redocking wrapper function', exc_info=True)
        raise

    finally:
        if os.path.exists(directory(base_output_path)):
            [os.remove(directory(base_output_path) + file) for file in os.listdir(directory(base_output_path)) if file.endswith('.vina')]


def prepare_structures(base_input_path:str, target:str, base_output_path:str,
                       pdb_codes=None, pH:float=7.4, charge_type:str='gas', preparation_options=None,
                       base_selected_mols=None, mol_filename:str='compounds',
                       receptor_prepared:bool=False, docking_engines:str='both'):
    """Prepare or reuse PDB receptors and prepare independent compound sources."""
    import shutil
    if docking_engines not in ("vina","dock6","both"):
        raise ValueError("Selecione Vina, DOCK6 ou ambos para a saída de preparação.")
    from biomolexplorer.storage import write_dataframe
    source = resolve_path(base_input_path) / target.replace(' ', '')
    if pdb_codes is None:
        import pandas as pd
        metadata = pd.read_csv(source / 'pdb_codes.csv')
        if 'STATUS' in metadata:metadata=metadata[metadata['STATUS'].isin(('succeeded','SUCCESS'))]
        pdb_codes = metadata[['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN']].values.tolist()
    output = resolve_path(base_output_path) / target.replace(' ', '') / 'Prepared'
    if receptor_prepared:
        if base_selected_mols is not None:
            from biomolexplorer.docking_preparation import validate_receptor_outputs
            validate_receptor_outputs(source/'Prepared',pdb_codes,docking_engines)
        output.mkdir(parents=True,exist_ok=True)
        receptors={f'{record[0]}_{record[3]}' for record in pdb_codes}
        for file in (source/'Prepared').iterdir():
            if file.name=='centers.csv' or file.name.split('.',1)[0] in receptors:
                shutil.copy2(file,output/file.name)
        prepared=pdb_codes
    else:
        docking = Docking(complex_input_path=str(source), output_path=str(output))
    settings=None
    if preparation_options is not None:
        from biomolexplorer.redocking_config import pair_key,validate_pairs
        settings={pair_key(record):preparation_options for record in pdb_codes}
        if not receptor_prepared:validate_pairs([list(r) for r in pdb_codes],settings)
    if not receptor_prepared:
        from biomolexplorer.redocking_preparation import prepare_pairs,write_failures
        prepared,failures=prepare_pairs(docking,pdb_codes,settings or {},charge_type,pH)
        write_failures(failures,output.parent)
    records = [list(record[:4]) + [record[4] if len(record) > 4 else None] for record in prepared]
    write_dataframe(DataFrame(records, columns=['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION']),
                    output.parent / 'pdb_codes.csv')
    if base_selected_mols is not None:
        from biomolexplorer.docking_preparation import prepare_candidates,validate_receptor_outputs
        validate_receptor_outputs(output,prepared,docking_engines)
        prepare_candidates(resolve_path(base_selected_mols)/(mol_filename+'.csv'),
                           output.parent/'Compounds',preparation_options,pH,docking_engines)
        receptors={f'{record[0]}_{record[3]}' for record in prepared}
        for file in output.iterdir():
            if file.suffix in ('.pdb','.pdbqt','.mol2') and file.name.split('.',1)[0] not in receptors:
                file.unlink()
