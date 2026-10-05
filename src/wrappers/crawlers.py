from biomolexplorer.paths import directory, resolve_path, worker_count
from biomolexplorer.storage import write_dataframe
from biomolexplorer.progress import report_progress
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Information retrieval",

    module_description=(
    "Wrapper module for managing and integrating crawler "
    "implementations available in the src/crawlers directory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
import json
from pathlib import Path
from typing import Optional, List
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from crawlers.targets import Targets
from crawlers.bioactivities import Bioactivity
from crawlers.molecules import Molecule
from crawlers.molecules import SimilarMols
from crawlers.molecules import ZincMols
from crawlers.pubchem import PubChemSimilarMols
from crawlers.complex import PDBComplex, PolymerEntityType, ExperimentalMethod

from kernel.loggers import LoggerManager
from kernel.utilities import fileHandling
#----------------------------------------------------------------------------------------------


logger = LoggerManager.get_logger('crawlers', log_file='logs/loaders.log')


def read_filters(path:str):

    try:
        path = resolve_path(path)
        with open(path, 'r') as fp:
            filters = json.load(fp)
        return filters

    except Exception as e:
        logger.error(f'Error during to perform {path} in read_filters wrapper function', exc_info=True)
        raise



def retrieve_compounds(search_term:str, base_output_path:str, include_pubchem:bool=False,
                pubchem_threshold:int=75, pubchem_max_records:int=1000,
                chembl_filters:Optional[dict]=None):
    """Retrieve target compounds, bioactivities and optional structural similars.

    ChEMBL supplies the target data; include_pubchem expands the compound set
    through PubChem. Source-specific names identify providers and their filters.
    """

    try:
        target_output_path = f'{base_output_path}/ChEMBL/targets/'
        bioactivity_output_path = f'{base_output_path}/ChEMBL/bioactivity/{search_term.replace(' ','')}/'
        molecule_output_path = f'{base_output_path}/ChEMBL/molecules/{search_term.replace(' ','')}/'
        similar_output_path = f'{base_output_path}/ChEMBL/similars/{search_term.replace(' ','')}/'


        target = Targets()
        bioact = Bioactivity()
        mols = Molecule()
        sims = SimilarMols()

        script_path = '/src/scripts/crawlers/target.json'
        filters = dict((chembl_filters or {}).get('target', read_filters(script_path)))
        target.set_outputpath(target_output_path)
        report_progress(f'Consultando o alvo {search_term} no ChEMBL…')
        target.search(search_term, filters)

        script_path = '/src/scripts/crawlers/bioactivity.json'
        filters = dict((chembl_filters or {}).get('bioactivity', read_filters(script_path)))
        bioact.set_outputpath(bioactivity_output_path)
        bioact.set_targetpath(target_output_path)
        report_progress(f'Baixando as bioatividades de {search_term} e aplicando os filtros…')
        bioact.search(search_term, filters)

        script_path = '/src/scripts/crawlers/molecules.json'
        filters = dict((chembl_filters or {}).get('molecules', read_filters(script_path)))
        mols.set_outputpath(molecule_output_path)
        mols.set_bioactivitypath(bioactivity_output_path)
        mols.search(filters)

        script_path = '/src/scripts/crawlers/similarmols.json'
        filters = dict((chembl_filters or {}).get('similars', read_filters(script_path)))
        sims.set_outputpath(similar_output_path)
        sims.set_bioactivitypath(bioactivity_output_path)
        report_progress('Buscando compostos similares no ChEMBL…')
        sims.search(filters)


        drugbank_output_path = f'{base_output_path}/ChEMBL/DrugBank/'
        molecules = fileHandling(output_path=drugbank_output_path)

        report_progress('Consolidando os compostos recuperados…')
        molecules.prepare_datamols(target=search_term,
                                   inputpath_mols=molecule_output_path,
                                   inputpath_similars=similar_output_path)

        if include_pubchem:
            report_progress('Buscando compostos similares no PubChem…')
            return expand_similar_compounds(search_term, base_output_path, pubchem_threshold, pubchem_max_records)
        chembl = PubChemSimilarMols.read_chembl([
            resolve_path(molecule_output_path), resolve_path(similar_output_path)])
        if chembl.empty:
            raise ValueError('No ChEMBL molecules matched the configured filters')
        return _save_compound_dataset(chembl, None, resolve_path(base_output_path), search_term)

    except Exception as e:
        logger.error(f'Error during to perform {search_term} in retrieve_compounds wrapper function', exc_info=True)
        raise


def expand_similar_compounds(search_term:str, base_output_path:str, threshold:int=75, max_records:int=1000,
                             base_input_path:Optional[str]=None, input_file:Optional[str]=None):
    """Expand existing ChEMBL downloads without repeating ChEMBL retrieval.

    Paths follow the existing wrappers: /datasets means <cwd>/datasets.
    """
    root = resolve_path(base_output_path)
    input_root = resolve_path(base_input_path) if base_input_path else root
    target = search_term.replace(' ', '')
    if input_file:
        from biomolexplorer.input_validation import validate_file
        if Path(input_file).name!=input_file:
            raise ValueError('Selecione um arquivo da pasta de entrada.')
        validate_file(input_root/input_file,'compounds')
        from pandas import read_csv
        chembl=read_csv(input_root/input_file)
    else:
        chembl = PubChemSimilarMols.read_chembl([
            input_root / 'ChEMBL' / 'molecules' / target,
            input_root / 'ChEMBL' / 'similars' / target,
        ])
    crawler = PubChemSimilarMols(root / 'PubChem' / 'similars' / target, threshold, max_records)
    try:
        compounds = crawler.search(chembl)
    finally:
        crawler.session.close()
    return _save_compound_dataset(chembl, compounds, root, search_term)


def _save_compound_dataset(chembl, compounds, root, search_term):
    from biomolexplorer.molecule_quality import clean_dataframe
    chembl=clean_dataframe(chembl,source='ChEMBL')
    if compounds is not None:compounds=clean_dataframe(compounds,source='PubChem')
    original = chembl.reindex(columns=['molecule_chembl_id', 'canonical_smiles', 'molecule_properties']).copy()
    original['source'] = 'ChEMBL'
    from pandas import concat
    combined = concat([original, compounds], ignore_index=True)
    combined['_structure'] = combined['canonical_smiles'].map(lambda s: PubChemSimilarMols._identity(s)[0])
    combined = combined.dropna(subset=['_structure']).drop_duplicates('_structure').drop(columns='_structure')
    output = root / 'compounds' / search_term.replace(' ', '')
    output.mkdir(parents=True, exist_ok=True)
    write_dataframe(combined, output / 'compounds.csv')
    return combined



def is_valid(value):
    return value is not None and (not isinstance(value, str) or value.strip() != '')


def load_pdb(target:str, base_output_path:str, pdb_ec:Optional[str]=None, organism:Optional[List[str]]=None,
             PolymerEntityTypeID:Optional[List[PolymerEntityType]]=None,
             ExperimentalMethodID:Optional[List[ExperimentalMethod]]=None,
             max_resolution:Optional[float]=None, must_have_ligand:Optional[bool]=True):

    try:

        pdb_output_path = f'{base_output_path}/PDB/{target.replace(' ','')}/'
        pdb = PDBComplex(output_path=pdb_output_path)

        filters = {
            key: value
            for key, value in {
                'PolymerEntityTypeID': PolymerEntityTypeID,
                'ExperimentalMethodID': ExperimentalMethodID,
                'ec_target': pdb_ec,
                'organism': organism,
                'max_resolution': max_resolution,
                'must_have_ligand': must_have_ligand
            }.items()
            if is_valid(value)
}
        pdb.get_pdb_files(filters=filters)

    except Exception as e:
        logger.error(f'Error during to perform {target} in load_pdb wrapper function', exc_info=True)
        raise




def load_zinc(base_output_path:str, filename:str, verbose=False, base_input_path:Optional[str]=None):

    try:

        zinc = ZincMols()
        output = filename.split('.')[0]

        zinc_output_path = f'{base_output_path}/'
        zinc.set_uri_inputpath(str(resolve_path(base_input_path or base_output_path) / filename))
        zinc.set_outputpath(zinc_output_path)
        zinc.search(output_filename=output, verbose=verbose)


    except Exception as e:
        logger.error(f'Error during to perform the wrapper load_zinc function', exc_info=True)
        raise
