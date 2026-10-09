from biomolexplorer.paths import directory, resolve_path, worker_count
from biomolexplorer.storage import write_dataframe, write_json
from biomolexplorer.retrieval import collection_name, target_criteria, identifiers, MODES
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
        return _clean_filters(filters)

    except Exception as e:
        logger.error(f'Error during to perform {path} in read_filters wrapper function', exc_info=True)
        raise



def _clean_filters(filters):
    return {key: value for key, value in filters.items() if value is not None and value != '' and value != []}


def retrieve_compounds(search_term:str, base_output_path:str, include_pubchem:bool=False,
                pubchem_threshold:int=75, pubchem_max_records:int=1000,
                chembl_filters:Optional[dict]=None, search_mode:str="target",
                max_targets:int=25, max_records:int=1000, expand_chembl:bool=True,
                similarity_threshold:int=70):
    """Retrieve target compounds, bioactivities and optional structural similars.

    ChEMBL supplies the target data; include_pubchem expands the compound set
    through PubChem. Source-specific names identify providers and their filters.
    """

    if not isinstance(search_term,str) or not search_term.strip():
        raise ValueError('Informe uma consulta ChEMBL.')
    for value in (max_targets,max_records):
        if type(value) is not int or value<1:raise ValueError('Os limites de recuperação devem ser inteiros positivos.')
    if type(similarity_threshold) is not int or not 1<=similarity_threshold<=100:
        raise ValueError('A similaridade ChEMBL deve estar entre 1 e 100.')
    if search_mode not in MODES:
        raise ValueError('Modo de busca ChEMBL inválido.')
    if search_mode in ('molecule_id', 'molecule_name', 'similarity', 'substructure'):
        return _retrieve_direct_compounds(search_term, base_output_path, search_mode, max_records,
                                         similarity_threshold, chembl_filters, include_pubchem,
                                         pubchem_threshold, pubchem_max_records)
    label = collection_name(search_term)
    try:
        target_output_path = f'{base_output_path}/ChEMBL/targets/'
        bioactivity_output_path = f'{base_output_path}/ChEMBL/bioactivity/{label}/'
        molecule_output_path = f'{base_output_path}/ChEMBL/molecules/{label}/'
        similar_output_path = f'{base_output_path}/ChEMBL/similars/{label}/'


        effective_filters = {}
        target = Targets()
        bioact = Bioactivity()
        mols = Molecule()
        sims = SimilarMols()

        script_path = '/src/scripts/crawlers/target.json'
        filters = _clean_filters((chembl_filters or {}).get('target', read_filters(script_path)))
        effective_filters['target'] = dict(filters)
        target.set_outputpath(target_output_path)
        report_progress(f'Consultando o alvo {search_term} no ChEMBL…')
        target.search(search_term, filters, search_mode=search_mode, max_targets=max_targets)

        script_path = '/src/scripts/crawlers/bioactivity.json'
        filters = _clean_filters((chembl_filters or {}).get('bioactivity', read_filters(script_path)))
        bioact.set_outputpath(bioactivity_output_path)
        bioact.set_targetpath(target_output_path)
        report_progress(f'Baixando as bioatividades de {search_term} e aplicando os filtros…')
        filters.setdefault("max_records", max_records)
        effective_filters['bioactivity'] = dict(filters)
        bioact.search(search_term, filters)

        script_path = '/src/scripts/crawlers/molecules.json'
        filters = _clean_filters((chembl_filters or {}).get('molecules', read_filters(script_path)))
        mols.set_outputpath(molecule_output_path)
        mols.set_bioactivitypath(bioactivity_output_path)
        effective_filters['molecules'] = dict(filters)
        mols.search(filters)

        script_path = '/src/scripts/crawlers/similarmols.json'
        filters = _clean_filters((chembl_filters or {}).get('similars', read_filters(script_path)))
        sims.set_outputpath(similar_output_path)
        sims.set_bioactivitypath(bioactivity_output_path)
        if expand_chembl:
            filters.setdefault('max_records', max_records)
            report_progress('Buscando compostos similares no ChEMBL…')
            effective_filters['similars'] = dict(filters)
            sims.search(filters)


        drugbank_output_path = f'{base_output_path}/ChEMBL/DrugBank/'
        molecules = fileHandling(output_path=drugbank_output_path)

        report_progress('Consolidando os compostos recuperados…')
        molecules.prepare_datamols(target=label,
                                   inputpath_mols=molecule_output_path,
                                   inputpath_similars=similar_output_path)

        write_json({'provider': 'ChEMBL', 'search_term': search_term, 'search_mode': search_mode,
                    'max_targets': max_targets, 'max_records_per_target': max_records,
                    'filters': effective_filters, 'expand_chembl': expand_chembl},
                   resolve_path(base_output_path)/'retrieval_report.json')
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


def _retrieve_direct_compounds(term, output_path, mode, maximum, threshold, filters,
                               include_pubchem, pubchem_threshold, pubchem_max_records):
    """Search molecules independently of targets or activity evidence."""
    from crawlers.chembl_client import ChEMBLClient
    from crawlers.molecules import _select_molecules
    from rdkit import Chem
    root = resolve_path(output_path)
    label = collection_name(term)
    client = ChEMBLClient()
    parameters = _clean_filters((filters or {}).get('molecules', {}))
    natural = parameters.pop('natural_product', None)
    if natural is not None:
        natural=int(natural)
        if natural not in (0,1):raise ValueError('Produto natural deve ser 0, 1 ou não especificado.')
    moltype = parameters.pop('molecule_type', None)
    weight = parameters.pop('molecule_weight', None)
    resource = client.molecule
    if mode == 'molecule_id':
        parameters['molecule_chembl_id__in'] = identifiers(term, 'chembl')
    elif mode == 'molecule_name':
        parameters['pref_name__icontains'] = term.strip()
    else:
        reference = term.strip()
        if mode == 'similarity' and reference.upper().startswith('CHEMBL'):
            reference = identifiers(reference, 'chembl')[0]
        elif Chem.MolFromSmiles(reference) is None:
            raise ValueError('Informe um SMILES válido para a consulta estrutural.')
        if mode == 'similarity':
            resource = client.similarity
            parameters.update(chembl_id=reference, similarity=threshold)
        else:
            resource = client.substructure
            parameters['smiles'] = reference
    report_progress('Consultando compostos diretamente no ChEMBL…')
    records = list(resource.filter(**parameters).take(maximum))
    molecules = _select_molecules(records, Molecule().str_to_dict, natural,
                                 moltype.lower() if moltype else None, float(weight) if weight else None)
    if molecules.empty:
        raise ValueError('Nenhum composto corresponde à busca e aos filtros selecionados.')
    folder = root/'ChEMBL/molecules'/label
    write_dataframe(molecules, folder/'results.csv')
    chembl = PubChemSimilarMols.read_chembl([folder])
    write_json({'provider': 'ChEMBL', 'search_term': term, 'search_mode': mode,
                'max_records': maximum, 'limit_reached': len(records)>=maximum, 'fetched': len(records), 'selected': len(chembl),
                'filters': filters or {}, 'similarity_threshold': threshold if mode=='similarity' else None,
                'activity_evidence': False}, root/'retrieval_report.json')
    if include_pubchem:
        return expand_similar_compounds(term, output_path, pubchem_threshold, pubchem_max_records)
    return _save_compound_dataset(chembl, None, root, term)


def expand_similar_compounds(search_term:str, base_output_path:str, threshold:int=75, max_records:int=1000,
                             base_input_path:Optional[str]=None, input_file:Optional[str]=None):
    """Expand existing ChEMBL downloads without repeating ChEMBL retrieval.

    Paths follow the existing wrappers: /datasets means <cwd>/datasets.
    """
    root = resolve_path(base_output_path)
    input_root = resolve_path(base_input_path) if base_input_path else root
    target = collection_name(search_term)
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
    output = root / 'compounds' / collection_name(search_term)
    output.mkdir(parents=True, exist_ok=True)
    write_dataframe(combined, output / 'compounds.csv')
    return combined



def is_valid(value):
    return value is not None and (not isinstance(value, str) or value.strip() != '')


def load_pdb(target:str="Estruturas", base_output_path:str="datasets", pdb_ec:Optional[str]=None, organism:Optional[List[str]]=None,
             PolymerEntityTypeID:Optional[List[PolymerEntityType]]=None,
             ExperimentalMethodID:Optional[List[ExperimentalMethod]]=None,
             max_resolution:Optional[float]=None, must_have_ligand:Optional[bool]=True,
             pdb_query:Optional[str]=None, pdb_ids:Optional[List[str]]=None,
             uniprot_ids:Optional[List[str]]=None, ligand_ids:Optional[List[str]]=None,
             max_records:int=100):

    try:

        if type(max_records) is not int or max_records<1:
            raise ValueError('O limite PDB deve ser um inteiro positivo.')
        pdb_output_path = f'{base_output_path}/PDB/{collection_name(target)}/'
        pdb = PDBComplex(output_path=pdb_output_path)

        filters = {
            key: value
            for key, value in {
                'pdb_query': pdb_query, 'pdb_ids': pdb_ids, 'uniprot_ids': uniprot_ids,
                'ligand_ids': ligand_ids, 'max_records': max_records,
                'PolymerEntityTypeID': PolymerEntityTypeID,
                'ExperimentalMethodID': ExperimentalMethodID,
                'ec_target': pdb_ec,
                'organism': organism,
                'max_resolution': max_resolution,
                'must_have_ligand': must_have_ligand
            }.items()
            if is_valid(value)
}
        if not any(filters.get(k) for k in ('pdb_query', 'pdb_ids', 'uniprot_ids', 'ligand_ids', 'ec_target',
                                                    'organism', 'PolymerEntityTypeID', 'ExperimentalMethodID', 'max_resolution')):
            if target not in ('Estruturas', 'MeuAlvo'):
                filters['pdb_query'] = target
            else:
                raise ValueError('Escolha uma busca PDB por texto, identificadores ou filtros; o nome da coleção é opcional.')
        return pdb.get_pdb_files(filters=filters)

    except Exception as e:
        logger.error(f'Error during to perform {target} in load_pdb wrapper function', exc_info=True)
        raise




def load_zinc(base_output_path:str, filename:str='zinc_urls.txt', verbose=0, base_input_path:Optional[str]=None, download_workers:int=4):

    try:

        from biomolexplorer.zinc_retrieval import retrieve_tranches
        return retrieve_tranches(resolve_path(base_input_path or base_output_path)/filename,base_output_path,download_workers=download_workers)


    except Exception as e:
        logger.error(f'Error during to perform the wrapper load_zinc function', exc_info=True)
        raise
