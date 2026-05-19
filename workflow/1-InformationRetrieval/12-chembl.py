#----------------------------------------------------------------------------------------------
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Core functions for managing and extracting chemical " 
    "information from the ChEMBL dataset"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.crawlers import load_chembl
#----------------------------------------------------------------------------------------------


if __name__ == "__main__":
    

    #----------------------------------------------------------------------------------------------
    # Example 1: Retrival information from ChEMBL database for monoamine oxidase
    # @param target_name: str = 'monoamine oxidase' - specific target name defined by ChEMBL
    # @param base_output_path: str = '/datasets' - base path to save the output files
    # @obs: Filters to compose retrieval information from ChEMBL database are defined by the
    # scripts in the scripts folder located in the src > scripts > crawlers folder.
    #----------------------------------------------------------------------------------------------
    load_chembl(target_name='monoamine oxidase',
                base_output_path='/datasets') 
    

    #----------------------------------------------------------------------------------------------
    # Example 2: Retrival information from ChEMBL database for Amine oxidase [flavin-containing] B
    # @param target_name: str = 'Amine oxidase [flavin-containing] B' - specific target name 
    # defined by ChEMBL @param base_output_path: str = '/datasets' - base path to save the output 
    # files @obs: Filters to compose retrieval information from ChEMBL database are defined by the
    # scripts in the scripts folder located in the src > scripts > crawlers folder.
    #----------------------------------------------------------------------------------------------
    load_chembl(target_name='Amine oxidase [flavin-containing] B',
                base_output_path='/datasets') 
    

    