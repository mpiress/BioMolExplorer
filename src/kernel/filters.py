from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="Molecular Relevant Filters",

    module_description=(
        "Module responsible for filtering relevant molecules."
    ),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
from pandas import DataFrame
from pathlib import Path

from rdkit import Chem
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.loggers import LoggerManager
#----------------------------------------------------------------------------------------------


class MyFilters():
    
    def __init__(self) -> None:
        self.path   = str(Path.cwd())
        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/filters.log')
        
    def get_path(self):
        return self.path
        


class LookupWizard(MyFilters):
    
    def __init__(self) -> None:
        super(MyFilters, self).__init__()
    
        
    def filter(self, input:DataFrame, query:str) -> DataFrame:
        
        try:

            result = input.query(query)
            return result
        
        except Exception as e:
            self.logger.error(f'Error during to perform {query} in filter function', exc_info=True)
        



class Molecule(MyFilters):
    
    def __init__(self) -> None:
        super(MyFilters, self).__init__()
    
    
    def clean_fragments(self, smiles:list) -> list:
        
        try:

            molecules  = [Chem.MolFromSmiles(s) for s in smiles]
            molecules = [smi for smi, mol in zip(smiles, molecules) if Chem.GetMolFrags(mol) == 1]
            return molecules
        
        except Exception as e:
            self.logger.error(f'Error during to perform the clean_fragments function', exc_info=True)
    
    
    
    def remove_duplicates(self, smiles:list) -> list:
        
        try:
            
            molecules = set([Chem.CanonSmiles(smi) for smi in smiles])
            return list(molecules)
        
        except Exception as e:
            self.logger.error(f'Error during to perform the remove_duplicates function', exc_info=True)
        
            
        