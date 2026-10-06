from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Bioactivities analysis by ChEMBL",

    module_description=(
        "Module responsible for extracting targets "
        "on ChEMBL by a target name reference."
    ),

    module_version="1.4.0"
)

#----------------------------------------------------------------------------------------------
import os

from pandas import DataFrame
from pathlib import Path
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from crawlers.settings import CrawlerSettings
from kernel.utilities import fileHandling
from kernel.loggers import LoggerManager
#----------------------------------------------------------------------------------------------


class Targets(CrawlerSettings):

    def __init__(self, path=None, extension='csv') -> None:
        super().__init__()
        self.__target    = super().get_client_connection().target
        self.__path      = str(Path.cwd())
        self.__extension = extension
        self.logger      = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/targets.log')
        self.set_outputpath(path) if path != None else None


    def set_outputpath(self, path:str):
        self.__outputpath = path
        if not os.path.exists(directory(self.__outputpath)):
            os.makedirs(directory(self.__outputpath), exist_ok=True)



    def search(self, search_term:str, filter_params:dict, search_mode='target', max_targets=25) -> None:
        from biomolexplorer.retrieval import target_criteria, collection_name
        from biomolexplorer.storage import write_dataframe
        filters = dict(filter_params or {})
        if 'type__in' in filters:
            filters['target_type__in'] = filters.pop('type__in')
        filters.pop('relationship_type', None)
        filters.update(target_criteria(search_term, search_mode))
        columns = ['pref_name', 'target_chembl_id', 'target_components', 'target_type', 'organism']
        try:
            query = self.__target.filter(**filters).only(columns).take(max_targets)
            target = DataFrame.from_records(list(query), columns=columns)
            if target.empty:
                raise ValueError(f'Nenhum alvo ChEMBL corresponde a "{search_term}". Revise o modo e os filtros.')
            target = target.drop_duplicates('target_chembl_id', ignore_index=True)
            write_dataframe(target, Path(directory(self.__outputpath))/(collection_name(search_term).upper()+'.csv'))
        except Exception:
            self.logger.exception('Error searching ChEMBL target %s', search_term)
            raise

    def save_target(self, targets:DataFrame, file_name:str) -> None:
        files   = fileHandling(output_path=self.__outputpath, ext=self.__extension)
        infile  =  files.isFile(file_name)[1]
        if targets.shape[0] > 0 and not infile:
            files.dataframe_to_csv(file_name.upper(), targets)





