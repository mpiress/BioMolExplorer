from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Bioactivities analysis by ChEMBL",

    module_description=(
        "Module responsible for extracting PDBs by target, "
        "using EC reference and filtering information."
    ),

    module_version="1.4.0"
)

#----------------------------------------------------------------------------------------------
import os

from enum import Enum
from pathlib import Path
from Bio.PDB import PDBParser, Polypeptide
from rcsbsearchapi.const import STRUCTURE_ATTRIBUTE_SEARCH_SERVICE
from rcsbsearchapi.search import AttributeQuery, TextQuery
from biomolexplorer.retrieval import identifiers
from biomolexplorer.storage import write_dataframe, write_json
from biomolexplorer.progress import report_progress
from pandas import DataFrame
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.loggers import LoggerManager
#----------------------------------------------------------------------------------------------


class PolymerEntityType(Enum):
    DNA = "DNA"
    NA_HYBRID = "NA-hybrid"
    OTHER = "Other"
    RNA = "RNA"
    PROTEIN = "Protein"


class ExperimentalMethod(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    POWDER_DIFFRACTION = "POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    SOLUTION_SCATTERING = "SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"



class PDBComplex():

    def __init__(self, output_path=None):
        self.__path = str(Path.cwd())
        self.set_outputpath(output_path) if output_path != None else None
        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/complex.log')



    def set_outputpath(self, output_path:str):
        self.__outputpath = output_path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)



    def get_pdb_ids_with_filters(self, filter_params:dict) -> list:

        try:

            filter_params = dict(filter_params)
            for key in ('organism','PolymerEntityTypeID','ExperimentalMethodID'):
                if isinstance(filter_params.get(key),str):filter_params[key]=[filter_params[key]]
            queries = []
            if filter_params.get('pdb_query'):
                queries.append(TextQuery(filter_params['pdb_query']))
            if filter_params.get('pdb_ids'):
                queries.append(AttributeQuery('rcsb_id', 'in', identifiers(filter_params['pdb_ids'], 'pdb')))
            if filter_params.get('uniprot_ids'):
                queries.append(AttributeQuery('rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
                                              'in', identifiers(filter_params['uniprot_ids'], 'uniprot')))
                queries.append(AttributeQuery('rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name',
                                              'exact_match', 'UniProt'))
            if filter_params.get('ligand_ids'):
                queries.append(AttributeQuery('rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id',
                                              'in', identifiers(filter_params['ligand_ids'], 'ligand')))

            if 'ec_target' in filter_params and filter_params['ec_target']:
                queries.append(AttributeQuery("rcsb_polymer_entity.rcsb_ec_lineage.id", "exact_match", filter_params['ec_target'], STRUCTURE_ATTRIBUTE_SEARCH_SERVICE))
            if 'PolymerEntityTypeID' in filter_params and filter_params['PolymerEntityTypeID']:
                polymer = AttributeQuery("entity_poly.rcsb_entity_polymer_type", "exact_match", getattr(filter_params['PolymerEntityTypeID'][0], 'value', filter_params['PolymerEntityTypeID'][0]), STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                for polymer_type in filter_params['PolymerEntityTypeID'][1:]:
                    polymer |= AttributeQuery("entity_poly.rcsb_entity_polymer_type", "exact_match", getattr(polymer_type, "value", polymer_type), STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                queries.append(polymer)
            if 'organism' in filter_params and filter_params['organism']:
                organism = AttributeQuery("rcsb_entity_source_organism.taxonomy_lineage.name", "exact_match", filter_params['organism'][0], STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                for org in filter_params['organism'][1:]:
                    organism |= AttributeQuery("rcsb_entity_source_organism.taxonomy_lineage.name", "exact_match", org, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                queries.append(organism)
            if 'ExperimentalMethodID' in filter_params and filter_params['ExperimentalMethodID']:
                ExpMethod = AttributeQuery("exptl.method", "exact_match", getattr(filter_params['ExperimentalMethodID'][0], 'value', filter_params['ExperimentalMethodID'][0]), STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                for method in filter_params['ExperimentalMethodID'][1:]:
                    ExpMethod |= AttributeQuery("exptl.method", "exact_match", getattr(method, "value", method), STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
                queries.append(ExpMethod)
            if 'max_resolution' in filter_params and filter_params['max_resolution'] != None:
                queries.append(AttributeQuery("rcsb_entry_info.resolution_combined", "less_or_equal", filter_params['max_resolution'], STRUCTURE_ATTRIBUTE_SEARCH_SERVICE))
            if 'must_have_ligand' in filter_params and filter_params['must_have_ligand']:
                queries.append(AttributeQuery("rcsb_entry_info.nonpolymer_entity_count", "greater", 0, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE))

            if queries:
                combined_query = queries[0]
                for query in queries[1:]:
                    combined_query &= query

                maximum = filter_params.get('max_records', 100)
                results = self._search(combined_query, maximum)

                return results
            else:
                raise ValueError('Informe texto, IDs PDB, UniProt, EC ou outro critério de busca.')

        except Exception as e:
            self.logger.exception('PDB search failed: %s', filter_params)
            raise


    def _search(self, query, maximum):
        """Bounded RCSB pagination with finite timeouts and transient retries."""
        import json
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        if type(maximum) is not int or maximum < 1:
            raise ValueError('O limite PDB deve ser um inteiro positivo.')
        node = json.loads(query.to_json())
        session = requests.Session()
        session.mount('https://', HTTPAdapter(max_retries=Retry(total=2, backoff_factor=.5,
                      status_forcelist=(429, 500, 502, 503, 504), allowed_methods=('POST',))))
        results, seen, start, total = [], set(), 0, None
        try:
            while len(results) < maximum:
                payload = {'query': node, 'return_type': 'entry', 'request_options': {
                    'paginate': {'start': start, 'rows': min(1000, maximum-len(results))},
                    'results_content_type': ['experimental']}}
                response = session.post('https://search.rcsb.org/rcsbsearch/v2/query', json=payload, timeout=(5, 30))
                if response.status_code == 204:
                    break
                response.raise_for_status()
                data = response.json()
                batch = data.get('result_set', [])
                total = data.get('total_count')
                if not isinstance(batch, list) or not isinstance(total, int):
                    raise ValueError('Resposta inválida da busca RCSB.')
                if not batch:
                    break
                added = 0
                for item in batch:
                    code = identifiers([item.get('identifier', '')], 'pdb')[0]
                    if code not in seen:
                        seen.add(code); results.append(code); added += 1
                        if len(results) >= maximum: break
                if not added:
                    raise RuntimeError('A paginação RCSB repetiu resultados sem avançar.')
                start += len(batch)
                if start >= total:
                    break
        except requests.RequestException as exc:
            raise RuntimeError('Falha na consulta RCSB após as tentativas automáticas. Confira a conexão e os critérios.') from exc
        finally:
            session.close()
        self._last_search = {'query': node, 'total_count': total, 'limit_reached': total is not None and total>len(results)}
        return results

    def __identify_ligands(self, pdb_file):

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('structure', pdb_file)

            resolution = structure.header.get('resolution', None)

            ligands = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if Polypeptide.is_aa(residue, standard=True):
                            continue
                        if residue.id[0] != ' ' and residue.resname != 'HOH':
                            ligands.append((residue.resname, residue.id[1], chain.id))

            return set(ligands), resolution

        except Exception:
            self.logger.exception('Invalid downloaded PDB: %s', pdb_file)
            raise


    def get_pdb_files(self, filters:dict):
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from biomolexplorer.storage import write_text
        root = Path(directory(self.__outputpath))
        pdb_codes = self.get_pdb_ids_with_filters(filters)
        if not pdb_codes:
            write_json({'filters': self._serializable(filters), 'matched': [], 'downloads': []}, root/'retrieval_report.json')
            raise ValueError('Nenhuma estrutura PDB corresponde aos critérios. Remova filtros ou altere a busca.')

        def download(code):
            session = requests.Session()
            session.mount('https://', HTTPAdapter(max_retries=Retry(total=2, backoff_factor=.5,
                          status_forcelist=(429, 500, 502, 503, 504), allowed_methods=('GET',))))
            try:
                response = session.get(f'https://files.rcsb.org/download/{code}.pdb', timeout=(5, 30))
                response.raise_for_status()
                if not any(line.startswith(('ATOM  ', 'HETATM')) for line in response.text.splitlines()):
                    raise ValueError('Resposta sem coordenadas PDB; a entrada pode exigir mmCIF.')
                text = ''.join(line+'\n' for line in response.text.splitlines() if not line.startswith(('LINK', 'SSBOND')))
                path = root/(code+'.pdb')
                from io import StringIO
                ligands, resolution = self.__identify_ligands(StringIO(text))
                write_text(text, path)
                return {'pdb_id': code, 'status': 'downloaded', 'ligands': len(ligands),
                        'resolution': resolution}, sorted(ligands)
            except Exception as exc:
                self.logger.exception('PDB download failed: %s', code)
                return {'pdb_id': code, 'status': 'failed', 'error': str(exc)}, []
            finally:
                session.close()

        outcomes = {}
        with ThreadPoolExecutor(max_workers=min(4, worker_count())) as executor:
            futures = {executor.submit(download, code): code for code in pdb_codes}
            for completed, future in enumerate(as_completed(futures), 1):
                outcomes[futures[future]] = future.result()
                report_progress('Baixando estruturas PDB…', completed, len(pdb_codes))
        records = []
        for code in pdb_codes:
            info, ligands = outcomes[code]
            for name, number, chain in ligands:
                records.append([code, name, number, chain, info['resolution']])
        write_dataframe(DataFrame(records, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION']), root/'pdb_codes.csv')
        report = {'filters': self._serializable(filters), 'matched': pdb_codes,
                  'downloads': [outcomes[code][0] for code in pdb_codes],
                  'limit': filters.get('max_records', 100), 'search': getattr(self, '_last_search', {}),
                  'downloaded': sum(info['status']=='downloaded' for info,_ in outcomes.values()),
                  'failed': sum(info['status']=='failed' for info,_ in outcomes.values())}
        write_json(report, root/'retrieval_report.json')
        report_progress(f"PDB: {report['downloaded']} estruturas baixadas; {report['failed']} falhas de download.")
        if not any(info['status']=='downloaded' for info, _ in outcomes.values()):
            raise RuntimeError('Nenhum PDB foi baixado. Consulte retrieval_report.json para os erros individuais.')
        return report

    @staticmethod
    def _serializable(filters):
        return {key: [item.value if isinstance(item, Enum) else item for item in value]
                if isinstance(value, list) else value for key, value in filters.items()}
