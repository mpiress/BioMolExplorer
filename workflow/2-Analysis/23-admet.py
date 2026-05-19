#----------------------------------------------------------------------------------------------
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Main functions for managing and integrating ADMET " 
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from wrappers.admet import ADMETWrapper
#----------------------------------------------------------------------------------------------
    
if __name__ == "__main__":

    adme_pipeline = ADMETWrapper(base_path='/datasets/ChEMBL/DrugBank/Molecules',
                                 input_file='molecules.csv')
    
    adme_pipeline.run_pipeline()
    
    adme_pipeline.export_results(
        output_based_path = '/resultados/admet',
        csv_path='all_properties.csv',
        bbb_pos_path='bbb_positive.smi',
        bbb_neg_path='bbb_negative.smi',
        hia_pos_path='hia_positive.smi'
    )
    
    adme_pipeline.generate_plot(output_path='/resultados/admet/plots',
                                output_image_file='adme_boiled_egg.png')
    
    adme_pipeline.print_summary()