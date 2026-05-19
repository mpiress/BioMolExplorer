from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Wrapper module for managing and integrating crawler " 
    "implementations available in the src/crawlers directory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
import json
from pathlib import Path
from typing import Optional, List
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from crawlers.targets import Targets
from crawlers.bioactivities import Bioactivity
from crawlers.molecules import Molecule
from crawlers.molecules import SimilarMols
from crawlers.molecules import ZincMols
from crawlers.complex import PDBComplex, PolymerEntityType, ExperimentalMethod

from kernel.loggers import LoggerManager
from kernel.utilities import fileHandling
#----------------------------------------------------------------------------------------------


logger = LoggerManager.get_logger('crawlers', log_file='logs/loaders.log')


def read_filters(path:str):

    try:
        path = str(Path.cwd()) + path
        with open(path, 'r') as fp:
            filters = json.load(fp)
        return filters
    
    except Exception as e:
        logger.error(f'Error during to perform {path} in read_filters wrapper function', exc_info=True)



def load_chembl(target_name:str, base_output_path:str):

    try:
        target_output_path = f'{base_output_path}/ChEMBL/targets/'
        bioactivity_output_path = f'{base_output_path}/ChEMBL/bioactivity/{target_name.replace(' ','')}/'
        molecule_output_path = f'{base_output_path}/ChEMBL/molecules/{target_name.replace(' ','')}/'
        similar_output_path = f'{base_output_path}/ChEMBL/similars/{target_name.replace(' ','')}/'
        

        target = Targets()
        bioact = Bioactivity()
        mols = Molecule()
        sims = SimilarMols()

        script_path = '/src/scripts/crawlers/target.json'
        filters = read_filters(script_path)
        target.set_outputpath(target_output_path)
        target.search(target_name, filters)

        script_path = '/src/scripts/crawlers/bioactivity.json'
        filters = read_filters(script_path)
        bioact.set_outputpath(bioactivity_output_path)
        bioact.set_targetpath(target_output_path)
        bioact.search(target_name, filters)

        script_path = '/src/scripts/crawlers/molecules.json'
        filters = read_filters(script_path)
        mols.set_outputpath(molecule_output_path)
        mols.set_bioactivitypath(bioactivity_output_path)
        mols.search(filters)

        script_path = '/src/scripts/crawlers/similarmols.json'
        filters = read_filters(script_path)
        sims.set_outputpath(similar_output_path)
        sims.set_bioactivitypath(bioactivity_output_path)
        sims.search(filters)


        drugbank_output_path = f'{base_output_path}/ChEMBL/DrugBank/'
        molecules = fileHandling(output_path=drugbank_output_path)

        molecules.prepare_datamols(target=target_name,
                                   inputpath_mols=molecule_output_path,
                                   inputpath_similars=similar_output_path)
    
    except Exception as e:
        logger.error(f'Error during to perform {target_name} in load_chembl wrapper function', exc_info=True)



def is_valid(value):
    return value is not None and (not isinstance(value, str) or value.strip() != '')


def load_pdb(target:str, base_output_path:str, pdb_ec:Optional[str]=None, organism:Optional[List[str]]=None,
             PolymerEntityTypeID:Optional[List[PolymerEntityType]]=None,
             ExperimentalMethodID:Optional[List[ExperimentalMethod]]=None,
             max_resolution:Optional[float]=None, must_have_ligand:Optional[bool]=True):

    try:
        
        pdb_output_path = f'{base_output_path}/PDB/{target.replace(' ','')}/'
        pdb = PDBComplex(output_path=pdb_output_path)

        filters = {
            key: value
            for key, value in {
                'PolymerEntityTypeID': PolymerEntityTypeID,
                'ExperimentalMethodID': ExperimentalMethodID,
                'ec_target': pdb_ec,
                'organism': organism,
                'max_resolution': max_resolution,
                'must_have_ligand': must_have_ligand
            }.items()
            if is_valid(value)
}
        pdb.get_pdb_files(filters=filters)
    
    except Exception as e:
        logger.error(f'Error during to perform {target} in load_pdb wrapper function', exc_info=True)




def load_zinc(base_output_path:str, filename:str, verbose=False):

    try:
        
        zinc = ZincMols()
        output = filename.split('.')[0]

        zinc_output_path = f'{base_output_path}/'
        zinc.set_uri_inputpath(f'{base_output_path}/{filename}')
        zinc.set_outputpath(zinc_output_path)
        zinc.search(output_filename=output, verbose=verbose)
            
        
    except Exception as e:
        logger.error(f'Error during to perform the wrapper load_zinc function', exc_info=True)