#----------------------------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Compound retrieval",

    module_description=(
    "Core functions for managing and extracting chemical "
    "information from ChEMBL and expanding the compound set with PubChem"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.crawlers import retrieve_compounds
#----------------------------------------------------------------------------------------------


if __name__ == "__main__":


    #----------------------------------------------------------------------------------------------
    # Example 1: Retrieve target compounds and expand with structural similars
    # @param search_term: str = specific target name defined by ChEMBL or ChEMBL_ID reference
    # @param base_output_path: str = '/datasets' - base path to save the output files
    # @obs: Filters to compose retrieval information from ChEMBL database are defined by the
    # scripts in the scripts folder located in the src > scripts > crawlers folder.
    #----------------------------------------------------------------------------------------------
    retrieve_compounds(search_term='CHEMBL220',
                base_output_path='/datasets',
                include_pubchem=True,
                pubchem_threshold=75,
                pubchem_max_records=1000)




