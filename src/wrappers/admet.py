from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

   module_description=(
    "Wrapper module for ADMET-based molecular evaluation "
    "functionalities available in the src/caad diretory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
from pathlib import Path
import pandas as pd
import os
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from caad.admet import MoleculeEvaluator, BoiledEggPlotter
from kernel.utilities import fileHandling
#----------------------------------------------------------------------------------------------

class ADMETWrapper:
    
    def __init__(self, base_output_path, base_input_path=None, input_file=None, verbose=False):
        self.base_output_path = base_output_path 
        self.set_output_path()

        self.base_input_path = base_input_path
        self.csv_files = self.set_input_path(input_file)

        self.evaluator = MoleculeEvaluator()
        self.df_results = None
        self.excluded_count = 0

        self.verbose = verbose
       

    def set_input_path(self, input_file) -> list:
        
        base_input_path = Path(str(Path.cwd()) + self.base_input_path)
        
        if not base_input_path.is_dir():
            raise FileNotFoundError(f"O caminho especificado não existe ou não é um diretório: {self.base_input_path}")
        
        csv_files = [str(arquivo.name) for arquivo in base_input_path.glob('*.csv')] if input_file is None else [input_file]
        return csv_files
    

    def set_output_path(self) -> None:
        base_output_path = Path(str(Path.cwd()) + self.base_output_path)
        if not base_output_path.exists() or not self.base_output_path is None:
           os.makedirs(base_output_path, exist_ok=True)


    def _load_data(self, filename):
        f1 = fileHandling(input_path=self.base_input_path, output_path=self.base_output_path)
        df = f1.csv_to_dataframe(filename)
        
        df = df[['molecule_chembl_id','canonical_smiles']]

        return df


    def run_pipeline(self):
        
        for filename in self.csv_files:

            filename = filename.split('.')[0]
            df_input = self._load_data(filename)
            
            rows = []
            self.excluded_count = 0

            for _, row in df_input.iterrows():
                smiles = row['canonical_smiles']
                name = row['molecule_chembl_id']
                
                props = self.evaluator.calculate_properties(smiles)
                
                if props is None:
                    continue
                
                # Filtro Toxicológico usando o Kernel
                if self.evaluator.is_toxic(props['Mol']):
                    self.excluded_count += 1
                    
                # Extração de descritores e predições
                tpsa, logp, mw = props['TPSA'], props['WLOGP'], props['MW']
                hbd, hba, rb = props['HBD'], props['HBA'], props['RB']

                bbb = self.evaluator.classify_bbb(tpsa, logp, mw, hbd, rb)
                hia = self.evaluator.classify_hia(tpsa, logp)
                pgp = self.evaluator.predict_pgp(mw, tpsa, logp, hbd)

                rows.append({
                    'canonical_smiles': smiles, 'molecule_chembl_id': name,
                    'TPSA': round(tpsa, 2), 'WLOGP': round(logp, 2), 'MW': round(mw, 2),
                    'HBD': hbd, 'HBA': hba, 'RB': rb,
                    'BBB': bbb, 'HIA': hia, 'PGP': pgp
                })

            results = pd.DataFrame(rows)
            self.export_results(results, filename)
            self.generate_plot(results, filename)
            self.print_summary(results) if self.verbose else None
            

    def export_results(self, results, filename):
        f1 = fileHandling(input_path=self.base_input_path, output_path=self.base_output_path)

        if results is None:
            raise ValueError("O pipeline precisa ser executado (.run_pipeline()) antes da exportação.")

        f1.dataframe_to_csv(filename, results)

        # Filtros e exportações parciais (SMILES e Name sem header/índice)
        for criteria, data in [( 'BBB==\"BBB+\"', filename+'_BBB+'), 
                               ( 'BBB==\"BBB-\"', filename+'_BBB-'), 
                               ( 'HIA==\"HIA+\"', filename+'_HIA+')]:
            
            tmp = results.query(criteria)[['canonical_smiles', 'molecule_chembl_id']]
            f1.dataframe_to_csv(data, tmp)


    def generate_plot(self, results, filename):
        if results is None:
            raise ValueError("Nenhum dado processado para gerar o gráfico.")
        BoiledEggPlotter.plot(results, self.base_output_path, filename+'_egg.png')


    def print_summary(self, results):
        if results is None:
            return
        total_valid = len(results)
        print('\n========== SUMMARY ==========')
        print(f'Total molecules input: {total_valid + self.excluded_count}')
        print(f'Excluded toxic compounds: {self.excluded_count}')
        print(f'Molecules processed: {total_valid}')
        print(f'BBB+: {len(results[results["BBB"]=="BBB+"])}')
        print(f'BBB-: {len(results[results["BBB"]=="BBB-"])}')
        print(f'HIA+: {len(results[results["HIA"]=="HIA+"])}')
        print('=============================\n')