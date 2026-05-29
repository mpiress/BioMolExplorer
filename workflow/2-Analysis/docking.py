#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Main functions for managing and integrating docking " 
    "analysis available in the wrappers folder (docking.py)"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
from wrappers.docking import perform_consensus
#----------------------------------------------------------------------------------------------


if __name__ == "__main__":


    perform_consensus(base_input_path='/datasets/PDB',
                      target='Butyrylcholinesterase',
                      pdb_code=('7AWG', 'S6Q', 605, 'A'),
                      base_output_path='/resultados/docking',
                      base_selected_mols='/datasets/ChEMBL/DrugBank/ADMET/Molecules',
                      mol_filename='molecules',
                      dock6_app_path='/home/michel/progs/dock6/',
                      prepare_complex=False, charge_type='am1')
    
   
    