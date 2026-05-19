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
from src.caad.admet import MoleculeEvaluator, BoiledEggPlotter
#----------------------------------------------------------------------------------------------

class ADMETWrapper:
    
    def __init__(self, base_output_path, base_input_path=None, input_file=None):
        self.set_path(str(Path.cwd()) + base_output_path)
        self.evaluator = MoleculeEvaluator()
        self.df_results = None
        self.excluded_count = 0
       

    def set_path(self, path:str) -> None:
        self.path = path if path != None else None
        if path != None and not os.path.exists(path):
           os.makedirs(path, exist_ok=True)


    def _load_data(self):
        input = self.path + '/' + self.input_file
        df = pd.read_csv(input, sep='\t', header=None)
        if df.shape[1] == 1:
            df.columns = ['SMILES']
            df['Name'] = [f'Cmpd_{i+1}' for i in range(len(df))]
        else:
            df = df.iloc[:, :2]
            df.columns = ['SMILES', 'Name']
        return df

    def run_pipeline(self):
        df_input = self._load_data()
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
                return True if False else None
                continue

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

        self.df_results = pd.DataFrame(rows)

    def export_results(self, output_path, csv_path, bbb_pos_path, bbb_neg_path, hia_pos_path):
        if self.df_results is None:
            raise ValueError("O pipeline precisa ser executado (.run_pipeline()) antes da exportação.")

        output = str(Path.cwd()) + output_path + '/' + csv_path
        
        self.df_results.to_csv(output, index=False)

        # Filtros e exportações parciais (SMILES e Name sem header/índice)
        for criteria, path in [( 'BBB==\"BBB+\"', bbb_pos_path), 
                               ( 'BBB==\"BBB-\"', bbb_neg_path), 
                               ( 'HIA==\"HIA+\"', hia_pos_path)]:
            
            output = str(Path.cwd()) + output_path + '/' + path    
            self.df_results.query(criteria)[['canonical_smiles', 'molecule_chembl_id']].to_csv(output, sep='\t', index=False, header=False)

    def generate_plot(self, output_path, output_image_file):
        if self.df_results is None:
            raise ValueError("Nenhum dado processado para gerar o gráfico.")
        BoiledEggPlotter.plot(self.df_results, output_path, output_image_file)

    def print_summary(self):
        if self.df_results is None:
            return
        total_valid = len(self.df_results)
        print('\n========== SUMMARY ==========')
        print(f'Total molecules input: {total_valid + self.excluded_count}')
        print(f'Excluded toxic compounds: {self.excluded_count}')
        print(f'Molecules processed: {total_valid}')
        print(f'BBB+: {len(self.df_results[self.df_results["BBB"]=="BBB+"])}')
        print(f'BBB-: {len(self.df_results[self.df_results["BBB"]=="BBB-"])}')
        print(f'HIA+: {len(self.df_results[self.df_results["HIA"]=="HIA+"])}')
        print('=============================\n')