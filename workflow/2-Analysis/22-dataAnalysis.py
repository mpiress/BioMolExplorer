#----------------------------------------------------------------------------------------------
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Main functions for managing and integrating data " 
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.molecular_analyzer import compute_similarity, analyze_graphs
from kernel.descriptors import similarityFunctions, fingerprints
from wrappers.molecular_analyzer import generate_fingerprints
#----------------------------------------------------------------------------------------------


if __name__ == "__main__":
    
    
    #----------------------------------------------------------------------------------------------
    # PREPARE FINGERPRINTS FOR EACH MOL AND SIM FILES AVAILABLE IN THE DRUGBANK FOLDER.
    # @param base_input_path: The path to the folder containing the MOL and SIM files.
    # @param amputate: If True, the amputate function is used to remove columns with no variation.
    # @param morgan: If True, the Morgan fingerprint is generated.
    # @param maccs: If True, the MACCS fingerprint is generated.
    # @param pharmacophore: If True, the pharmacophore fingerprint is generated.
    #----------------------------------------------------------------------------------------------
    generate_fingerprints(base_input_path='/datasets/ChEMBL/DrugBank',
                          morgan=True, maccs=True, pharmacophore=True
    )
    
    #----------------------------------------------------------------------------------------------
    # COMPUTE THE SIMILARITY BETWEEN THE MOLECULES BASED ON THE FINGERPRINTS GENERATED.
    # @param base_input_path: The path to the folder containing the fingerprints.
    # @param base_output_path: The path to the folder where the similarity computations will be saved.
    # @param metric: The similarity metric to be used (Tanimoto, Dice, Cosine, etc).
    # @param fingerprint: The fingerprint to be used (Morgan, MACCS, and Pharmacophore).
    #----------------------------------------------------------------------------------------------
    compute_similarity(base_input_path='/datasets/ChEMBL/DrugBank/Fingerprints',
                       base_output_path='/datasets/ChEMBL/DrugBank',
                       metric=similarityFunctions.TanimotoSimilarity,
                       fingerprint=fingerprints.Morgan                  
    )

    #----------------------------------------------------------------------------------------------
    # COMPUTE THE RELATIONSHIP BETWEEN THE MOLECULES BASED ON THE SIMILARITY VALUES COMPUTED. THE GRAPH
    # REPRESENTATION IS USED TO IDENTIFY RELEVANT CHARACTERISTICS IN MOLECULES AND ANALOGS STRUCTURES.
    # @param base_input_path: The path to the folder containing the Similarity folder computed.
    # @param base_output_path: The path to the folder where the graph representations will be saved.
    # @param metric: The similarity metric to be considered in this step (Tanimoto, Dice, Cosine, etc).
    # @param fingerprint: The fingerprint to be considered in this step (Morgan, MACCS, and Pharmacophore).
    #----------------------------------------------------------------------------------------------
    analyze_graphs(base_input_path='/datasets/ChEMBL/DrugBank',
                    base_output_path='/resultados/grafos',
                    metric=similarityFunctions.TanimotoSimilarity,
                    fingerprint=fingerprints.Morgan
    )
    