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
    # @param search_term: str = specific target name defined by ChEMBL or ChEMBL_ID reference
    # @param base_output_path: str = '/datasets' - base path to save the output files
    # @obs: Filters to compose retrieval information from ChEMBL database are defined by the
    # scripts in the scripts folder located in the src > scripts > crawlers folder.
    #----------------------------------------------------------------------------------------------
    load_chembl(search_term='CHEMBL1914',
                base_output_path='/datasets') 
    


    