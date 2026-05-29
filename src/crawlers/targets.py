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
        if not os.path.exists(self.__path + self.__outputpath):
            os.makedirs(self.__path + self.__outputpath, exist_ok=True)
      
     
       
    def search(self, search_term:str, filter_params:dict) -> None:
        try:
            if filter_params is None:
                filter_params = {}
                
            files = fileHandling(input_path=self.__outputpath, ext=self.__extension)
            search_term_upper = search_term.upper()
            
            # O nome do arquivo salvo ainda pode ser baseado no termo de busca
            infile = files.isFile(search_term_upper)[0]
            columns = ['pref_name', 'target_chembl_id', 'target_components', 'target_type']

            # Se o arquivo local já existe, lê do CSV, senão busca na API do ChEMBL
            if infile:
                target = files.csv_to_dataframe(search_term_upper)
            else:
                # 1. Se parecer um ChEMBL ID (ex: CHEMBL240)
                if search_term_upper.startswith("CHEMBL"):
                    filter_params["target_chembl_id"] = search_term_upper
                
                # 2. Se parecer um ID de Acesso do UniProt (ex: P00533 - geralmente 6 ou 10 caracteres alfanuméricos)
                elif len(search_term) in [6, 10] and any(char.isdigit() for char in search_term):
                    filter_params["target_components__accession"] = search_term_upper
                
                # 3. Se for um nome de texto, usamos "__icontains" para busca parcial (ou mantém iexact se preferir)
                else:
                    filter_params["pref_name__icontains"] = search_term

                # Executa o filtro na API do ChEMBL
                target = self.__target.filter(**filter_params).only(columns)

            # Processamento dos dados retornados
            if len(target) > 0:
                # Se veio da API, convertemos o query result para DataFrame
                if not infile:
                    target = DataFrame.from_records(target)
                
                target.drop_duplicates(subset='target_chembl_id', inplace=True, ignore_index=True)
                
                # Garante que o DataFrame final tenha as colunas desejadas (se existirem)
                available_cols = [col for col in columns if col in target.columns]
                target = target[available_cols]
            else:
                target = DataFrame()
            
            if self.__outputpath is not None:
                self.save_target(target, search_term_upper)
                
        except Exception as e:
            self.logger.error(f'Error during search for target "{search_term}"', exc_info=True)
            
    
     

    def save_target(self, targets:DataFrame, file_name:str) -> None:
        files   = fileHandling(output_path=self.__outputpath, ext=self.__extension)
        infile  =  files.isFile(file_name)[1]
        if targets.shape[0] > 0 and not infile:
            files.dataframe_to_csv(file_name.upper(), targets)
            
            
            
        
    
