from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Wrapper module for managing and integrating redocking " 
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from typing import Optional, List, Tuple, Literal
from pathlib import Path
from pandas import DataFrame
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
                      num_modes:Optional[int]=10, prepare_complex:Optional[bool]=True, charge_type:Optional[ChargeType]='gas') -> None:
    
    logger = LoggerManager.get_logger('wrapper_redocking', log_file='logs/redocking.log')

    try:

        base_prepared_complexes = f'{base_input_path}/{target.replace(' ','')}/Prepared/'
        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}/'
        path                    = str(Path.cwd())+'/' + base_input_path

        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)
            
        if pdb_codes is None:
            pdb_codes = f1.csv_to_dataframe('pdb_codes')
            pdb_codes = pdb_codes[['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION']].to_records(index=False)
        
        
        if prepare_complex:
            dock = Docking(complex_input_path=base_input_path, output_path=base_prepared_complexes)
            pdb_codes = dock.prepare_for_docking(pdb_codes=pdb_codes, charge_type=charge_type, pH=pH, redefine_centerofmass=True)
            pdb_codes = DataFrame(pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION'])
            f1.dataframe_to_csv('pdb_codes', pdb_codes)

        
        vina = DockVina(ligand_input_path=base_prepared_complexes, receptor_input_path=base_prepared_complexes,
                        output_path=base_output_path, complex_input_path=base_input_path,
                        pdb_codes=pdb_codes, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box, 
                        exhaustiveness=exhaustiveness, num_modes=num_modes)
        
        vina.redocking(pH=pH)
        
        f1 = fileHandling(input_path=base_input_path, output_path=base_input_path)
        
        pdb  = f1.csv_to_dataframe('pdb_codes')
        
        complexes = set([file + '.pdb' for file in pdb[pdb['RMSD'].notnull()]['PDB_CODE'].to_list()])
        files     = set([f for f in os.listdir(path) if f.endswith('.pdb')])
        files     = files - complexes
        [os.remove(path + f) for f in files]

        prefix   = [f.split('.')[0] for f in files]
        files    = os.listdir(path+'/Prepared/')
        [os.remove(path +'/Prepared/' + f) for f in files if any(f.startswith(p) for p in prefix)]
       
        pdb = pdb[pdb['RMSD'].notnull()]
        f1.dataframe_to_csv('pdb_codes', pdb)
        

    except Exception as e:
        logger.error(f'Error during to perform the {target} in redocking wrapper function', exc_info=True)

    finally:
        if os.path.exists(base_output_path[1:]):
            [os.remove(base_output_path[1:] + file) for file in os.listdir(base_output_path[1:]) if file.endswith('.vina')]