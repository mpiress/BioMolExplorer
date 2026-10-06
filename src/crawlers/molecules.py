from biomolexplorer.paths import directory, resolve_path, worker_count
from biomolexplorer.progress import report_progress
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Bioactivities analysis by ChEMBL",

    module_description=(
        "Module responsible for extracting molecules by target "
        "on ChEMBL by a target bioactivity reference."
    ),

    module_version="1.4.0"
)

#----------------------------------------------------------------------------------------------
import os
import requests

from pandas import DataFrame, concat
from concurrent import futures
from threading import Lock
from pathlib import Path
from typing import Optional
import ast
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.utilities import fileHandling, fileReading
from kernel.loggers import LoggerManager
from crawlers.settings import CrawlerSettings
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
#----------------------------------------------------------------------------------------------


def _monitor_queries(pool, message, sizes=None):
    completed = 0
    total = sum(sizes.values()) if sizes else len(pool)
    for future in futures.as_completed(pool):
        try:
            future.result()
        except Exception:
            # Avoid starting every queued request after a provider has failed.
            for pending in pool:
                pending.cancel()
            raise
        completed += sizes[future] if sizes else 1
        report_progress(message, completed, total)


def _select_molecules(records, convert_properties, np_filter, mol_filter, mwt_filter):
    molecules = DataFrame.from_records(records)
    if molecules.empty:
        return molecules
    molecules.drop_duplicates(subset='molecule_chembl_id', inplace=True, ignore_index=True)
    molecules['molecule_type'] = molecules['molecule_type'].astype(str).str.lower()
    molecules['natural_product'] = molecules['natural_product'].fillna(-1).astype(int)
    molecules['molecule_properties'] = molecules['molecule_properties'].apply(convert_properties)
    if np_filter is not None:
        molecules = molecules[molecules['natural_product'] == np_filter]
    if mol_filter:
        molecules = molecules[molecules['molecule_type'] == mol_filter]
    if mwt_filter is not None:
        molecules = molecules[molecules['molecule_properties'].apply(
            lambda value: value.get('full_mwt', float('inf')) <= mwt_filter)]
    return molecules


class MyMolecules():

    def __init__(self) -> None:
        self.path   = str(Path.cwd())
        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/molecules.log')

    def get_path(self):
        return self.path



class Molecule(CrawlerSettings, MyMolecules):

    def __init__(self, path=None, bioactivity_path=None, extension='csv') -> None:
        CrawlerSettings.__init__(self)
        MyMolecules.__init__(self)
        self.__molecule    = self.get_client_connection().molecule
        self.__extension   = extension
        self.set_outputpath(path) if path != None else None
        self.set_bioactivitypath(bioactivity_path) if bioactivity_path != None else None



    def set_bioactivitypath(self, path:str):
        self.__bioactivitypath = path
        if not os.path.exists(directory(self.__bioactivitypath)):
            print('[ERROR]: The bioactivity path needs to be informed before!')
            raise FileNotFoundError('Required input directory does not exist')


    def set_outputpath(self, path:str):
        self.__outputpath = path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)


    def str_to_dict(self, value):
        if isinstance(value, dict):
            converted = value
        elif isinstance(value, str):
            try:
                converted = ast.literal_eval(value)
            except (SyntaxError, ValueError):
                return {}
        else:
            return {}

        if 'full_mwt' in converted:
            try:
                converted['full_mwt'] = float(converted['full_mwt'])
            except (ValueError, TypeError):
                converted['full_mwt'] = float('inf')

        return converted



    def __search_mol(self, molecule_id:str, filter_params:dict, np_filter:int, mol_filter:str, mwt_filter:float) -> None:
        filter_params = dict(filter_params or {})

        try:
            files   = fileHandling(input_path=self.__outputpath, ext=self.__extension)
            infile  =  files.isFile(molecule_id)[0]

            filter_params['molecule_chembl_id'] = molecule_id

            molecule = files.csv_to_dataframe(molecule_id) if infile else self.__molecule.filter(**filter_params)
            if len(molecule) == 0:
                return

            try:
                molecule = _select_molecules(molecule, self.str_to_dict, np_filter, mol_filter, mwt_filter)

            except Exception as e:
                self.logger.error(f'Error during to perform {molecule_id} molecule in __search_mol function', exc_info=True)
                raise

            self.save_molecule(molecule, molecule_id) if molecule.shape[0] > 0 else None

        except Exception as e:
            self.logger.error(f'Error during to perform {molecule_id} molecule in __search_mol function', exc_info=True)
            raise




    def save_molecule(self, molecule:DataFrame, file_name:str) -> None:
        files   = fileHandling(output_path=self.__outputpath, ext=self.__extension)
        infile  =  files.isFile(file_name)[1]
        if molecule.shape[0] > 0 and not infile:
            files.dataframe_to_csv(file_name, molecule)




    def search(self, filter_params:dict):
        filter_params = dict(filter_params or {})

        f1    = fileHandling(input_path=self.__bioactivitypath, ext=self.__extension)
        files = [f.rsplit('.')[0] for f in os.listdir(directory(self.__bioactivitypath)) if f.endswith('.csv')]

        mols = []
        for file in files:
            tmp  = f1.csv_to_dataframe(file)
            mols  = mols + tmp['molecule_chembl_id'].tolist()

        np_filter  = filter_params.pop('natural_product', None)
        np_filter  = int(np_filter) if np_filter != None else None

        mol_filter = filter_params.pop('molecule_type', None)
        mol_filter = mol_filter.lower() if mol_filter != None else None

        mwt_filter = filter_params.pop('molecule_weight', None)
        mwt_filter = float(mwt_filter) if mwt_filter != None else None


        identifiers = list(dict.fromkeys(mol for mol in mols if isinstance(mol, str) and mol))
        report_progress('Baixando os registros moleculares do ChEMBL…', 0, len(identifiers))
        with futures.ThreadPoolExecutor(max_workers=min(2, worker_count())) as executor:
            batches = [identifiers[offset:offset + 25] for offset in range(0, len(identifiers), 25)]
            pool = {executor.submit(self.__search_batch, batch, filter_params, np_filter, mol_filter, mwt_filter): batch
                    for batch in batches}
            _monitor_queries(pool, 'Baixando os registros moleculares do ChEMBL…',
                             {future: len(batch) for future, batch in pool.items()})

    def __search_batch(self, identifiers, filter_params, np_filter, mol_filter, mwt_filter):
        """Fetch a small set at once, retaining the existing per-compound exports."""
        files = fileHandling(input_path=self.__outputpath, ext=self.__extension)
        pending = []
        for identifier in identifiers:
            if files.isFile(identifier)[0]:
                self.__search_mol(identifier, filter_params, np_filter, mol_filter, mwt_filter)
            else:
                pending.append(identifier)
        if not pending:
            return
        parameters = dict(filter_params)
        parameters['molecule_chembl_id__in'] = pending
        records = list(self.__molecule.filter(**parameters))
        if any(row.get('molecule_chembl_id') not in pending for row in records):
            raise RuntimeError('Resposta ChEMBL contém um composto não solicitado.')
        molecules = _select_molecules(records, self.str_to_dict, np_filter, mol_filter, mwt_filter)
        if molecules.empty:
            return
        for identifier, rows in molecules.groupby('molecule_chembl_id', sort=False):
            self.save_molecule(rows, identifier)







class SimilarMols(CrawlerSettings, MyMolecules):

    def __init__(self, path:Optional[str]=None,  bioactivity_path:Optional[str]=None, extension:Optional[str]='csv') -> None:
        CrawlerSettings.__init__(self)
        MyMolecules.__init__(self)
        self.__similarity  = self.get_client_connection().similarity
        self.__extension   = extension
        self.set_outputpath(path) if path != None else None
        self.set_bioactivitypath(bioactivity_path) if bioactivity_path != None else None



    def set_bioactivitypath(self, path:str):
        self.__bioactivitypath = path
        if not os.path.exists(directory(self.__bioactivitypath)):
            print('[ERROR]: The bioactivity path needs to be informed before!')
            raise FileNotFoundError('Required input directory does not exist')


    def set_outputpath(self, path:str):
        self.__outputpath = path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)


    def str_to_dict(self, value):
        if isinstance(value, dict):
            converted = value
        elif isinstance(value, str):
            try:
                converted = ast.literal_eval(value)
            except (SyntaxError, ValueError):
                return {}
        else:
            return {}

        if 'full_mwt' in converted:
            try:
                converted['full_mwt'] = float(converted['full_mwt'])
            except (ValueError, TypeError):
                converted['full_mwt'] = float('inf')

        return converted


    def __search_similar_mols(self, molecule_id:str, filter_params:dict, np_filter:int, mol_filter:str, mwt_filter:float) -> None:
        filter_params = dict(filter_params or {})

        try:
            files   = fileHandling(input_path=self.__outputpath, ext=self.__extension)
            infile  =  files.isFile(molecule_id)[0]

            filter_params['chembl_id'] = molecule_id

            maximum = filter_params.pop('max_records', 1000)
            molecules = self.__similarity.filter(**filter_params).take(maximum)
            if len(molecules) == 0:
                return

            try:

                molecules = _select_molecules(molecules, self.str_to_dict, np_filter, mol_filter, mwt_filter)

            except Exception:
                self.logger.exception('Invalid ChEMBL similar molecule response for %s', molecule_id)
                raise

            self.save_molecule(molecules, molecule_id) if molecules.shape[0] > 0 else None

        except Exception as e:
            self.logger.error(f'Error during to perform {molecule_id} molecule in __search_similar_mols function', exc_info=True)
            raise




    def save_molecule(self, molecule:DataFrame, file_name:str) -> None:
        files   = fileHandling(output_path=self.__outputpath, ext=self.__extension)
        infile  =  files.isFile(file_name)[1]
        if molecule.shape[0] > 0 and not infile:
            files.dataframe_to_csv(file_name, molecule)



    def search(self, filter_params:dict) -> None:
        filter_params = dict(filter_params or {})

        f1    = fileHandling(input_path=self.__bioactivitypath, ext=self.__extension)
        files = [f.rsplit('.')[0] for f in os.listdir(directory(self.__bioactivitypath)) if f.endswith('.csv')]

        mols = []
        for file in files:
            tmp  = f1.csv_to_dataframe(file)
            mols  = mols + tmp['molecule_chembl_id'].tolist()

        np_filter  = filter_params.pop('natural_product', None)
        np_filter  = int(np_filter) if np_filter != None else None

        mol_filter = filter_params.pop('molecule_type', None)
        mol_filter = mol_filter.lower() if mol_filter != None else None

        mwt_filter = filter_params.pop('molecule_weight', None)
        mwt_filter = float(mwt_filter) if mwt_filter != None else None

        identifiers = list(dict.fromkeys(mol for mol in mols if isinstance(mol, str) and mol))
        report_progress('Buscando compostos similares no ChEMBL…', 0, len(identifiers))
        with futures.ThreadPoolExecutor(max_workers=min(2, worker_count())) as executor:
            pool = {executor.submit(self.__search_similar_mols, mol, filter_params, np_filter, mol_filter, mwt_filter): mol for mol in identifiers}
            _monitor_queries(pool, 'Buscando compostos similares no ChEMBL…')





class ZincMols(MyMolecules):

    def __init__(self, uri_filename:Optional[str]=None, outputpath:Optional[str]=None) -> None:
        super().__init__()
        self.set_uri_inputpath(uri_filename) if uri_filename != None else None
        self.set_outputpath(outputpath) if outputpath != None else None
        self.lock = Lock()



    def set_uri_inputpath(self, path:str):
        self.__uri_inputpath = path
        if not os.path.exists(resolve_path(self.__uri_inputpath)):
            print('[ERROR]: The uri path needs to be informed before!')
            raise FileNotFoundError('Required input directory does not exist')


    def set_outputpath(self, path:str):
        self.__outputpath = path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)



    def __search_in_zinc(self, idx, url, verbose) -> DataFrame:

        url = url.strip()
        mol = DataFrame(columns=['smile', 'zinc_id'])

        # Configuração de Retentativas
        retry_strategy = Retry(
            total=5, # Tenta 5 vezes antes de desistir
            backoff_factor=2, # Espera: 2s, 4s, 8s, 16s, 32s...
            status_forcelist=[429, 500, 502, 503, 504], # Erros que disparam a retentativa
            allowed_methods=["GET"]
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("https://", adapter)
        session.mount("http://", adapter)

        # Adicionando Headers para evitar o 403
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        }

        try:
            response = session.get(url, headers=headers, timeout=120)
            response.raise_for_status()

            if response.status_code == 200:
                conteudo = response.text.splitlines()[1:]
                conteudo = [token.split() for token in conteudo if token.strip()]
                if any(len(row)!=2 for row in conteudo):
                    raise ValueError('Tabela ZINC inválida: esperado SMILES e identificador por linha.')
                mol = DataFrame(conteudo, columns=['smile', 'zinc_id'])
                print('File number:', idx, ' URL:', url) if verbose else None

            return mol

        except Exception as e:
            self.logger.error(f'Error during to perform {idx} molecule in {url} url in __search_in_zinc function', exc_info=True)
            raise
        finally:
            session.close()




    def search(self, output_filename='zinc', verbose=False) -> None:

        try:

            files = fileHandling(output_path=self.__outputpath, ext='csv')
            uri   = self.__uri_inputpath[:self.__uri_inputpath.rfind('/')+1]
            urls  = fileReading(inputpath=uri, file=self.__uri_inputpath.rsplit('/')[-1])
            mols  = DataFrame(columns=['smile', 'zinc_id'])
            files.dataframe_to_csv(output_filename, mols)

            chunk_size   = 100
            chunk_number = 0
            max_itens    = urls.get_size()
            processed    = 0

            with futures.ThreadPoolExecutor(max_workers=1) as executor:
                while(processed < max_itens):
                    data = urls.get_chunk(chunk_number, chunk_size)
                    chunk_number += 1
                    idx = processed
                    processed += len(data)


                    pool = {executor.submit(self.__search_in_zinc, item[0]+idx, item[1], verbose) : item for item in enumerate(data)}
                    for future in pool:
                        tmp = future.result()
                        if isinstance(tmp, DataFrame):
                            with self.lock:
                                files.dataframe_to_csv(output_filename, tmp, mode='a')


        except Exception as e:
            self.logger.error(f'Error during to perform the search function', exc_info=True)
            raise
