#----------------------------------------------------------------------------------------------
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Main functions for managing and integrating redocking " 
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.redocking import perform_redocking
#----------------------------------------------------------------------------------------------



if __name__ == "__main__":

    #----------------------------------------------------------------------------------------------
    # REDOCKING EXPERIMENTS USING AUTODOCK VINA WITH COMPLEXES FROM PDB - RMSD CALCULATION
    #----------------------------------------------------------------------------------------------
    perform_redocking(base_input_path='/datasets/PDB',
                      target='MonoamineOxidaseB',
                      base_output_path='/resultados/redocking',
                      prepare_complex=True, charge_type='am1')
    