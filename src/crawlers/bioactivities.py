from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Bioactivities analysis by ChEMBL",

    module_description=(
        "Module responsible for extracting bioactivities by target "
        "on ChEMBL by a target reference."
    ),

    module_version="1.4.0"
)

#----------------------------------------------------------------------------------------------
import os

from pandas import DataFrame, to_numeric
from pathlib import Path
from concurrent import futures
from typing import Optional
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from crawlers.settings import CrawlerSettings
from kernel.utilities import fileHandling
from kernel.loggers import LoggerManager
#----------------------------------------------------------------------------------------------



class Bioactivity(CrawlerSettings):


    def __init__(self, target_path=None, path=None, extension='csv') -> None:
        super().__init__()
        self.__bioactivity = super().get_client_connection().activity
        self.__path = str(Path.cwd())
        self.__extension   = extension
        self.set_outputpath(path) if path != None else None
        self.set_targetpath(target_path) if target_path != None else None
        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/bioactivities.log')



    def set_targetpath(self, path:str):
        self.__targetpath = path
        if not os.path.exists(directory(self.__targetpath)):
            print('[ERROR]: The target path needs to be informed before!')
            raise FileNotFoundError('Required input directory does not exist')



    def set_outputpath(self, path:str):
        self.__outputpath = path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)




    def __search_bioactivity(self, target_id:str, filter_params:dict) -> None:
        filter_params = dict(filter_params or {})

        try:
            max_value_ref = filter_params.pop('max_value_ref', None)
            max_records = filter_params.pop('max_records', 1000)
            filter_params.pop('molecule_type', None)
            if max_value_ref is not None and (not isinstance(filter_params.get('standard_units'),str) or not filter_params['standard_units'].strip()):
                raise ValueError('Informe a unidade padrão para aplicar um limite numérico de atividade.')
            if max_value_ref is not None:
                import math
                max_value_ref=float(max_value_ref)
                if not math.isfinite(max_value_ref) or max_value_ref<0:
                    raise ValueError('O limite de atividade deve ser finito e não negativo.')
            filter_params['target_chembl_id'] = target_id

            columns = ['activity_id', 'activity_properties', 'canonical_smiles', 'molecule_chembl_id', 'molecule_pref_name',
                    'parent_molecule_chembl_id', 'pchembl_value', 'qudt_units', 'target_organism', 'target_pref_name', 'type', 'units', 'value',
                    'standard_type','standard_units','standard_value','standard_relation']

            bioact = self.__bioactivity.filter(**filter_params).only(columns).take(max_records)

            if len(bioact) > 0:
                bioact = DataFrame.from_records(bioact,columns=columns)
                bioact.dropna(subset=['canonical_smiles'], inplace=True)
                value_column='standard_value' if 'standard_value' in bioact.columns else 'value'
                bioact[value_column]=to_numeric(bioact[value_column],errors='coerce')
                if max_value_ref is not None:
                    bioact=bioact.dropna(subset=[value_column])
            else:
                bioact = DataFrame(columns=columns)

            value_column='standard_value' if 'standard_value' in bioact.columns else 'value'
            if max_value_ref is not None:
                bioact = bioact[bioact[value_column] <= float(max_value_ref)]
            self.save_bioactivity(bioact, target_id)

        except Exception as e:
            self.logger.error(f'Error during to perform {target_id} in __search_bioactivity function', exc_info=True)
            raise




    def search(self, target:str, filter_params:dict) -> None:
        filter_params = dict(filter_params or {})

        files       = fileHandling(input_path=self.__targetpath, ext=self.__extension)
        from biomolexplorer.retrieval import collection_name
        target_ids  = files.csv_to_dataframe(collection_name(target).upper())
        target_ids  = target_ids['target_chembl_id'].tolist()


        with futures.ThreadPoolExecutor(max_workers=min(2, worker_count())) as executor:
            pool = {executor.submit(self.__search_bioactivity, target, filter_params) : target for target in target_ids}

        [future.result() for future in pool]



    def save_bioactivity(self, bioactivity:DataFrame, file_name:str) -> None:
        from biomolexplorer.storage import write_dataframe
        write_dataframe(bioactivity, Path(directory(self.__outputpath))/(file_name+'.csv'))
