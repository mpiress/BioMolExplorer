#----------------------------------------------------------------------------------------------
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Core functions for managing and extracting chemical " 
    "information from the ZINC dataset"
),

    module_version="1.0.0"
)


#----------------------------------------------------------------------------------------------
from wrappers.crawlers import load_chembl, load_pdb, load_zinc
#----------------------------------------------------------------------------------------------



if __name__ == "__main__":
    

    #----------------------------------------------------------------------------------------------
    # Example 1: Load ZINC database
    # @param base_output_path: str = '/datasets' - base output path to save the data and find the 
    #                                              URIs files
    # @param zinc2d: bool  = False - if True, download 2D structures if URI is available
    # @param zinc3d: bool  = True  - if True, download 3D structures if URI is available
    # @param verbose: bool = True  - if True, print the progress in the console
    #----------------------------------------------------------------------------------------------
    load_zinc(base_output_path='/datasets/ZINC',
              filename='ZINC2D.uri',
              verbose=True) 
    
    
    
   