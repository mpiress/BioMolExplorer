from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Wrapper module for managing and integrating molecular " 
    "analysis available in the src/caad directory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
import os
from pathlib import Path
from typing import Optional, List
import json
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.utilities import MolExplorer
from kernel.descriptors import MolSimilarity
from caad.complex_network import GraphAnalysis
from kernel.descriptors import similarityFunctions, fingerprints
from kernel.loggers import LoggerManager
from kernel.descriptors import Descriptors
from kernel.utilities import fileHandling
#----------------------------------------------------------------------------------------------


logger = LoggerManager.get_logger('molecular_analyzer', log_file='logs/analyzer.log')

def read_filters(path:str):

    try:
        path = str(Path.cwd()) + path
        with open(path, 'r') as fp:
            filters = json.load(fp)
        return filters
    
    except Exception as e:
        logger.error(f'Error during to perform {path} in read_filters wrapper function', exc_info=True)


def compute_similarity(base_input_path:str, base_output_path:str, 
                       metric:Optional[similarityFunctions]=similarityFunctions.TanimotoSimilarity,
                       fingerprint:Optional[fingerprints]=fingerprints.Morgan,
                       filename:Optional[str]=None, threshold:Optional[int]=None):

    try:
        
        base_input_path  = f'{base_input_path}/'
        base_output_path = f'{base_output_path}/Similarity/'
        
        script_path = '/src/scripts/crawlers/similarmols.json'
        filters = read_filters(script_path)
        threshold = int(filters['similarity']) / 100 if threshold is None else threshold / 100
        similarity = MolSimilarity(threshold=threshold)
        similarity.set_inputpath(base_input_path)
        similarity.set_outputpath(base_output_path)
        similarity.perform_similarity(filename=filename, metric=metric, fp=fingerprint.value)
    
    except Exception as e:
        logger.error(f'Error during to perform the wrapper compute_similarity function', exc_info=True)
    
    

def analyze_graphs(base_input_path:str, base_output_path:str, metric, fingerprint):
    
    try:
        
        base_input_path  = f'{base_input_path}/'
        base_output_path = f'{base_output_path}/'
        
        ga = GraphAnalysis()
        ga.set_inputpath(base_input_path)
        ga.set_outputpath(base_output_path)
        ga.prepare_graph_analysis(metric=metric, fp=fingerprint)
    
    except Exception as e:
        logger.error(f'Error during to perform the wrapper analyze_graphs function', exc_info=True)




def filter_mutagenic_tumorigenic(base_input_path:str, base_output_path:str, datawarrior_filename:str, delimiter:str,
                                 mutagenic:Optional[List[str]]=['high', 'low'], tumorigenic:Optional[List[str]]=['high', 'low'],
                                 druglikeness:Optional[List[float]]=[-1.0, 2.0]):
    
    try:
        
        base_input_path  = f'{base_input_path}/'
        base_output_path = f'{base_output_path}/'

        mutagenic = [m.lower() for m in mutagenic]
        tumorigenic = [t.lower() for t in tumorigenic]
        
        file = str(Path.cwd()) + base_input_path + datawarrior_filename
        if not os.path.isfile(file):
            logger.error(f'Error during to perform the filter_mutagenic_tumorigenic function', exc_info=True)
            logger.error(f'File {file} not found!!', exc_info=True)
            return
        
        mols = MolExplorer(input_path=base_input_path, output_path=base_output_path)
        mols.extract_mutagenic_tumorigenic(datawarrior_filename=datawarrior_filename,
                                           delimiter=delimiter, mutagenic=mutagenic,
                                           tumorigenic=tumorigenic, druglikeness=druglikeness)

    except Exception as e:
        logger.error(f'Error during to perform the wrapper filter_mutagenic_tumorigenic function', exc_info=True)



def generate_fingerprints(base_input_path:str, morgan_n_bits:Optional[int]=2048, radius:Optional[int]=2, files:Optional[List]=None,
                          morgan:Optional[bool]=True, maccs:Optional[bool]=True, pharmacophore:Optional[bool]=True):

    try:

        base_output_path  = f'{base_input_path}/Fingerprints/'
        base_input_path  = f'{base_input_path}/'

        f1 = fileHandling(input_path=base_input_path, output_path=base_output_path)
        ds = Descriptors(inputpath=base_input_path, outputpath=base_output_path)
                
        if files is None:
            files = [f for f in os.listdir(base_input_path[1:]) if f.endswith('_MOLS.csv') or f.endswith('_SIMS.csv')]
            
        for filename in files:
            filename = filename.split('.')[0]
            data = f1.csv_to_dataframe(filename)
            data.rename(columns={'canonical_smiles': 'smiles'}, inplace=True)
            data = data[['molecule_chembl_id','smiles']]
            
            if morgan:
                fingerprints = ds.get_fingerprints(smiles_df=data, morgan=True, morgan_n_bits=morgan_n_bits, radius=radius,
                                                   maccs=False, pharmacophore=False)
                
                f1.dataframe_to_csv('morgan_'+filename, fingerprints)
            
            if maccs:
                fingerprints = ds.get_fingerprints(smiles_df=data, morgan=False, maccs=True, pharmacophore=False)
                
                f1.dataframe_to_csv('maccs_'+filename, fingerprints)
            
            if pharmacophore:
                fingerprints = ds.get_fingerprints(smiles_df=data, morgan=False, maccs=False, pharmacophore=True)
                
                f1.dataframe_to_csv('pharmacophore_'+filename, fingerprints)

    except Exception as e:
        logger.error(f'Error during to perform the wrapper generate_fingerprints function', exc_info=True)
        logger.error(f'[ERROR]: {e}', exc_info=True)