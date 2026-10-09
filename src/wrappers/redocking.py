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
                      num_modes:Optional[int]=10, prepare_complex:Optional[bool]=True, charge_type:Optional[ChargeType]='gas', preparation_pairs:Optional[dict]=None) -> None:

    logger = LoggerManager.get_logger('wrapper_redocking', log_file='logs/redocking.log')

    try:

        base_prepared_complexes = f'{base_input_path}/{target.replace(' ','')}/Prepared/'
        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}/'
        path                    = directory(base_input_path)

        from biomolexplorer.redocking_config import validate_structure_pairs, validate_redocking_tools
        records = [list(r) for r in pdb_codes.to_records(index=False)] if isinstance(pdb_codes, DataFrame) else pdb_codes
        validate_structure_pairs(base_input_path, records, preparation_pairs or {}, prepared=not prepare_complex)
        validate_redocking_tools(prepare_complex)
        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)

        if pdb_codes is None:
            raise ValueError('Selecione pelo menos um par receptor / ligante e sua cadeia na aba Input Data.')
        from biomolexplorer.redocking_config import validate_pairs
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

        if prepare_complex:
            dock = Docking(complex_input_path=base_input_path, output_path=base_prepared_complexes)
            pdb_codes = dock.prepare_for_docking(pdb_codes=pdb_codes.to_records(index=False), charge_type=charge_type, pH=pH, redefine_centerofmass=True, preparation_pairs=preparation_pairs or {})
            if not pdb_codes:raise ValueError('A preparação não produziu complexos utilizáveis. Consulte o log de preparação.')
            pdb_codes = DataFrame(pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION'])
            f1.dataframe_to_csv('pdb_codes', pdb_codes)


        validate_structure_pairs(base_input_path, [list(r) for r in pdb_codes.to_records(index=False)], preparation_pairs or {}, prepared=True)
        vina = DockVina(ligand_input_path=base_prepared_complexes, receptor_input_path=base_prepared_complexes,
                        output_path=base_output_path, complex_input_path=base_input_path,
                        pdb_codes=pdb_codes, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box,
                        exhaustiveness=exhaustiveness, num_modes=num_modes)

        vina.redocking(pH=pH)

        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)

        pdb  = f1.csv_to_dataframe('pdb_codes')

        # Selecting a subset must not delete other structures from the input collection.
        pdb = pdb[pdb['RMSD'].notnull()]
        f1.dataframe_to_csv('pdb_codes', pdb)


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
        prepared = docking.prepare_for_docking(pdb_codes, charge_type, pH, True, preparation_pairs=settings)
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
