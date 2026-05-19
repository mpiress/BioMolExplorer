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
    "information from the PDB dataset"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.crawlers import load_pdb
from crawlers.complex import PolymerEntityType, ExperimentalMethod
#----------------------------------------------------------------------------------------------



if __name__ == "__main__":
    

    #----------------------------------------------------------------------------------------------
    #EXAMPLE 1: Load PDBs for the target Acetylcholinesterase
    # @param target: A reference name for the target, used to create the output folder
    # @param base_output_path: The base path to save the output files
    # @param pdb_ec: The Enzyme Commission number for the target
    # @param organism: The organism name for the target
    # @param PolymerEntityTypeID: The type of polymer entity (DNA, NA_HYBRID, PROTEIN, RNA, OTHER)
    # @param ExperimentalMethodID: The experimental method used to obtain the structure (X_RAY_DIFFRACTION,
    # ELECTRON_MICROSCOPY, SOLID_STATE_NMR, SOLUTION_NMR, NEUTRON_DIFFRACTION, FIBER_DIFFRACTION,
    # POWDER_DIFFRACTION, ELECTRON_CRYSTALLOGRAPHY, OTHER)
    # @param max_resolution: The maximum resolution for the structure
    # @param must_have_ligand: If the structure must have a ligand as part of the complex
    # @obs: The parameters are flexible and can be adjusted according to the research needs.
    #----------------------------------------------------------------------------------------------
    load_pdb(target='Monoamine Oxidase B', 
             base_output_path='/datasets', 
             pdb_ec='1.4.3.1',
             PolymerEntityTypeID=[PolymerEntityType.PROTEIN],
             ExperimentalMethodID=[ExperimentalMethod.X_RAY_DIFFRACTION],
             max_resolution=2.0, must_have_ligand=True)
    
    #----------------------------------------------------------------------------------------------
    #EXAMPLE 2: Load PDBs for the target Butyrylcholinesterase
    # @param target: A reference name for the target, used to create the output folder
    # @param base_output_path: The base path to save the output files
    # @param pdb_ec: The Enzyme Commission number for the target
    # @param organism: The organism name for the target
    # @param PolymerEntityTypeID: The type of polymer entity (DNA, NA_HYBRID, PROTEIN, RNA, OTHER)
    # @param ExperimentalMethodID: The experimental method used to obtain the structure (X_RAY_DIFFRACTION,
    # ELECTRON_MICROSCOPY, SOLID_STATE_NMR, SOLUTION_NMR, NEUTRON_DIFFRACTION, FIBER_DIFFRACTION,
    # POWDER_DIFFRACTION, ELECTRON_CRYSTALLOGRAPHY, OTHER)
    # @param max_resolution: The maximum resolution for the structure
    # @param must_have_ligand: If the structure must have a ligand as part of the complex
    # @obs: The parameters are flexible and can be adjusted according to the research needs.
    #----------------------------------------------------------------------------------------------
    load_pdb(target='Monoamine Oxidase B', 
             base_output_path='/datasets', 
             pdb_ec='1.4.3.21',
             PolymerEntityTypeID=[PolymerEntityType.PROTEIN],
             ExperimentalMethodID=[ExperimentalMethod.X_RAY_DIFFRACTION],
             max_resolution=2.0, must_have_ligand=True)
    