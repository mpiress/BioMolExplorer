from biomolexplorer.paths import directory, resolve_path, worker_count
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
        path = resolve_path(path)
        with open(path, 'r') as fp:
            filters = json.load(fp)
        return filters

    except Exception as e:
        logger.error(f'Error during to perform {path} in read_filters wrapper function', exc_info=True)
        raise


def compute_similarity(base_input_path:str, base_output_path:str,
                       metric:Optional[similarityFunctions]=similarityFunctions.TanimotoSimilarity,
                       fingerprint:Optional[fingerprints]=fingerprints.Morgan,
                       filename:Optional[str]=None, threshold:Optional[int]=None, approximate:bool=True):

    try:

        base_input_path  = f'{base_input_path}/'
        base_output_path = f'{base_output_path}/Similarity/'

        script_path = '/src/scripts/crawlers/similarmols.json'
        filters = read_filters(script_path)
        threshold = int(filters['similarity']) / 100 if threshold is None else threshold / 100
        similarity = MolSimilarity(threshold=threshold, approximate=approximate)
        similarity.set_inputpath(base_input_path)
        similarity.set_outputpath(base_output_path)
        similarity.perform_similarity(filename=filename, metric=metric, fp=fingerprint.value)

    except Exception as e:
        logger.error(f'Error during to perform the wrapper compute_similarity function', exc_info=True)
        raise



def analyze_graphs(base_input_path:Optional[str]=None, base_output_path:Optional[str]=None,
                   metric:Optional[similarityFunctions]=similarityFunctions.TanimotoSimilarity,
                   fingerprint:Optional[fingerprints]=fingerprints.Morgan,
                   similarity_path:Optional[str]=None, fingerprints_path:Optional[str]=None,
                   threshold:int=70, mcs_timeout:int=30, mcs_ring_matches_ring_only:bool=True,
                   mcs_complete_rings_only:bool=False, graph_inputs:Optional[List]=None):

    try:

        from caad.graph_results import analyze_inputs, folder_inputs
        if base_output_path is None:raise ValueError('Informe a pasta de saída dos grafos.')
        base_output_path=resolve_path(base_output_path)
        base_input_path=resolve_path(base_input_path) if base_input_path is not None else None
        fingerprints_path=resolve_path(fingerprints_path) if fingerprints_path is not None else None
        similarity_path=resolve_path(similarity_path) if similarity_path is not None else None
        entries=graph_inputs or folder_inputs(base_input_path,fingerprints_path,similarity_path,metric,fingerprint)
        analyze_inputs(entries,base_output_path,metric=metric,fingerprint=fingerprint,threshold=threshold,
            mcs_timeout=mcs_timeout,ring_matches_ring_only=mcs_ring_matches_ring_only,complete_rings_only=mcs_complete_rings_only)

    except Exception as e:
        logger.error(f'Error during to perform the wrapper analyze_graphs function', exc_info=True)
        raise




def filter_mutagenic_tumorigenic(base_input_path:str, base_output_path:str, datawarrior_filename:str, delimiter:str,
                                 mutagenic:Optional[List[str]]=['high', 'low'], tumorigenic:Optional[List[str]]=['high', 'low'],
                                 druglikeness:Optional[List[float]]=[-1.0, 2.0]):

    try:

        base_input_path  = f'{base_input_path}/'
        base_output_path = f'{base_output_path}/'

        mutagenic = [m.lower() for m in mutagenic]
        tumorigenic = [t.lower() for t in tumorigenic]

        file = directory(base_input_path) + datawarrior_filename
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
        raise



def generate_fingerprints(base_input_path:str, morgan_n_bits:Optional[int]=2048, radius:Optional[int]=2,
                          files:Optional[List]=None, morgan:Optional[bool]=True, maccs:Optional[bool]=True,
                          pharmacophore:Optional[bool]=True, base_output_path:Optional[str]=None,
                          chunk_size:int=100):
    """Generate fingerprints in bounded CSV chunks, with atomic final outputs."""
    import pandas as pd
    from biomolexplorer.storage import write_dataframe_chunks
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    output = resolve_path(base_output_path or f'{base_input_path}/Fingerprints/')
    input_dir = resolve_path(base_input_path)
    if files is None:
        files = [p.name for p in sorted(input_dir.glob('*.csv'))
                 if {'molecule_chembl_id', 'canonical_smiles'}.issubset(pd.read_csv(p, nrows=0).columns)]
    if not files:
        raise ValueError('No compound CSV files found for fingerprint generation')
    enabled = {'morgan': morgan, 'maccs': maccs, 'pharmacophore': pharmacophore}
    if not any(enabled.values()):
        raise ValueError('Select at least one fingerprint algorithm')
    ds = Descriptors()
    for filename in files:
        stem = Path(filename).stem
        source = input_dir / (stem + '.csv')
        for kind, selected in enabled.items():
            if not selected:
                continue
            def frames():
                metadata_columns=[c for c in pd.read_csv(source,nrows=0).columns if c!='fingerprint']
                yield pd.DataFrame(columns=[*metadata_columns,'fingerprint'])
                for frame in pd.read_csv(source, chunksize=chunk_size,dtype={'molecule_chembl_id':str}):
                    from biomolexplorer.molecule_quality import clean_dataframe
                    frame = clean_dataframe(frame,source=str(source)).rename(columns={'canonical_smiles': 'smiles'})
                    if frame.empty:
                        continue
                    values=ds.get_fingerprints(smiles_df=frame, morgan=kind == 'morgan', maccs=kind == 'maccs',
                        pharmacophore=kind == 'pharmacophore', radius=radius, morgan_n_bits=morgan_n_bits)
                    metadata=frame.rename(columns={'smiles':'canonical_smiles'}).drop(columns=['fingerprint'],errors='ignore')
                    metadata['fingerprint']=values['fingerprint'].to_numpy()
                    yield metadata
            write_dataframe_chunks(frames(), output / f'{kind}_{stem}.csv')
