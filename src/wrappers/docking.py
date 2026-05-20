#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Wrapper module for managing and integrating molecular " 
    "analysis available by consensus docking strategy"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
from typing import Optional, List, Tuple, Literal
from pathlib import Path
from pandas import DataFrame
import os
import re
import time
import shutil

from Bio.PDB import PDBParser, Polypeptide

from caad.docking import DockVina, Dock6, Docking
from kernel.loggers import LoggerManager
from kernel.utilities import fileHandling

import matplotlib.pyplot as plt
import seaborn as sns
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
logger     = LoggerManager.get_logger('wrapper_docking', log_file='logs/docking.log')
ChargeType = Literal['gas', 'am1']
#----------------------------------------------------------------------------------------------



def identify_ligand(pdb_file, ligand, chain_id):
    parser = PDBParser(QUIET=True)

    structure = parser.get_structure('structure', pdb_file)
    lig = None

    chain_id = [item for item in chain_id if item]
    
    for model in structure:
        for chain in model:
            if chain.id not in chain_id: 
                continue
            for residue in chain:
                if Polypeptide.is_aa(residue, standard=True):
                    continue
                if residue.id[0] != ' ' and residue.resname != 'HOH' and residue.resname == ligand:
                    lig = residue.id[1] 

    return lig




def get_better_complex(base_input_path:str) -> list:
    
    try:

        f2 = fileHandling(input_path=base_input_path)
        pdb_codes  = f2.csv_to_dataframe('pdb_codes')

        pdb_codes['Score'] = pdb_codes['RESOLUTION'] + pdb_codes['RMSD']
        min_index = pdb_codes['Score'].idxmin()
            
        best_pdb = pdb_codes.loc[min_index]
        pdb_codes = [(best_pdb['PDB_CODE'], best_pdb['LIGAND'], best_pdb['RESNUM'], best_pdb['CHAIN'])]

        return pdb_codes
    
    except Exception as e:
        logger.error(f'Error during to perform the get_better_complex wrapper function', exc_info=True)




def plot_scatter_comparison(df:DataFrame, output_path:str):
    df_minmax = df.copy()
    df_zscore = df.copy()

    # Normalização Min-Max
    df_minmax['vina'] = df_minmax['vina'].abs()
    df_minmax['dock6'] = df_minmax['dock6'].abs()
    df_minmax['vina'] = (df_minmax['vina'] - df_minmax['vina'].min()) / (df_minmax['vina'].max() - df_minmax['vina'].min())
    df_minmax['dock6'] = (df_minmax['dock6'] - df_minmax['dock6'].min()) / (df_minmax['dock6'].max() - df_minmax['dock6'].min())

    # Normalização Z-score
    df_zscore['vina'] = df_zscore['vina'].abs()
    df_zscore['dock6'] = df_zscore['dock6'].abs()
    df_zscore['vina'] = (df_zscore['vina'] - df_zscore['vina'].mean()) / df_zscore['vina'].std()
    df_zscore['dock6'] = (df_zscore['dock6'] - df_zscore['dock6'].mean()) / df_zscore['dock6'].std()

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    sns.scatterplot(x=df_minmax['vina'], y=df_minmax['dock6'], ax=axes[0], color='blue')
    axes[0].set_title("Min-Max Correlaction")

    sns.scatterplot(x=df_zscore['vina'], y=df_zscore['dock6'], ax=axes[1], color='red')
    axes[1].set_title("Z-score Correlaction")

    plt.tight_layout()
    plt.savefig(output_path[1:] + '/correlation.png')
    



def generate_consensus(base_input_path:str, base_output_path:str, target:str, repulsion_weight: float = 1.0) -> None:
    
    vinapath  = f'{base_input_path}/Vina/' 
    dock6path = f'{base_input_path}/Dock6/flex/' 
    dock6path = dock6path if os.path.isdir(dock6path[1:]) else f'{base_input_path}/Dock6/rigid/'
    path      = str(Path.cwd())

    if not os.path.isdir(vinapath[1:]):
        print(f'[ERROR] The path {vinapath} is not valid!')
        exit(1)

    if not os.path.isdir(dock6path[1:]):
        print(f'[ERROR] The path {dock6path} is not valid!')
        exit(1)

    dock6  = [f.split('_scored')[0] for f in os.listdir(dock6path[1:]) if f.endswith('.mol2')]
    
    molecules = {}
    for file in dock6:
    
        with open(path + vinapath + file + '.lig.pdbqt', 'r') as fp:
            content = fp.read()
        
        # Expressão regular para capturar o resultado do Vina
        vina_score = re.search(r'REMARK VINA RESULT:\s+(-?\d+\.\d+)', content)
        vina_score = float(vina_score.group(1)) if vina_score else 0.0
        
        with open(path + dock6path + file + '_scored.mol2', 'r') as fp:
            content = fp.read() 
        
        # Expressão regular para capturar o Grid_Score
        dock6_score = re.search(r'Grid_Score:\s+(-?\d+\.\d+)', content)
        dock6_score = float(dock6_score.group(1)) if dock6_score else 0.0

        # Captura da energia repulsiva interna
        repulsion_energy = re.search(r'Internal_energy_repulsive:\s+(-?\d+\.\d+)', content)
        repulsion_energy = float(repulsion_energy.group(1)) if repulsion_energy else 0.0

        # Penalização do Grid_Score pela energia repulsiva com peso (λ)
        dock6_score = (dock6_score + (repulsion_weight * repulsion_energy)) if (dock6_score + (repulsion_weight * repulsion_energy)) < 0 else 0.0

        molecules[file] = (vina_score, dock6_score)
    
    df = DataFrame.from_dict(molecules, orient='index', columns=['vina', 'dock6'])
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'molecule'}, inplace=True)
    

    # Z-score
    df['z-score'] = ((df['vina'].abs() - df['vina'].abs().mean()) / df['vina'].abs().std() + 
                    (df['dock6'].abs() - df['dock6'].abs().mean()) / df['dock6'].abs().std()) / 2


    # Min-Max
    df['min-max'] = ((df['vina'].abs() - df['vina'].abs().min()) / (df['vina'].abs().max() - df['vina'].abs().min()) + 
                 (df['dock6'].abs() - df['dock6'].abs().min()) / (df['dock6'].abs().max() - df['dock6'].abs().min())) / 2

    # Arredondamento para duas casas decimais
    df = df.round(2)

    plot_scatter_comparison(df, base_output_path)
    
    # Ordenação por ambas as colunas (primeiro por Z-score, depois por Min-Max)
    df = df.sort_values(by=['z-score', 'min-max'], ascending=[False, False])

    f1 = fileHandling(input_path=base_output_path+'/', output_path=base_output_path+'/')
    f1.dataframe_to_csv(target, df)    

    


def perform_docking_vina(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str, 
                         mol_filename:str, pdb_code:Optional[Tuple[str, str, str, str]]=None, 
                         pH:Optional[float]=7.4, sizeof_box:Optional[List]=[24,24,24], 
                         exhaustiveness:Optional[int]=20, num_modes:Optional[int]=10) -> None:
    
    try:

        base_prepared_complexes = f'{base_input_path}/{target.replace(' ','')}/Prepared/'
        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_input_mols         = f'{base_output_path}/Molecules/'
        base_selected_mols      = f'{base_selected_mols}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}'

        
        if os.path.exists(f'{str(Path.cwd())+base_output_path}/Vina/'):
            print(f'[INFO] The path {base_output_path}/Vina/ already exists!')
            return       
     
        vina = DockVina(ligand_input_path=base_selected_mols, receptor_input_path=base_prepared_complexes,
                        output_path=base_input_mols, complex_input_path=base_prepared_complexes,
                        pdb_codes=pdb_code, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box, 
                        exhaustiveness=exhaustiveness, num_modes=num_modes, mol_filename=mol_filename)
        
        
        if len(os.listdir(str(Path.cwd())+'/' + base_input_mols)) == 0:
            vina.prepare_compounds_for_vina(pH=pH)
            
        vina.set_ligandpath(base_input_mols)
        vina.set_outputpath(f'{base_output_path}/Vina/')
        vina.docking(base_selected_mols)
        del vina

    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_docking wrapper function', exc_info=True)
    
    
            


def perform_docking_dock6(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,  
                          dock6_app_path:str, charge_type:str, mol_filename:str,
                          pdb_code:Optional[Tuple[str, str, str, str]]=None, density:Optional[float]=0.5,
                          radius:Optional[float]=1.4, distance:Optional[float]=10.0,
                          conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex', 
                          plot_max_residues:Optional[int]=50) -> None:
    
    try:

        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_selected_mols      = f'{base_selected_mols}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}'
        
        pdb_code=f'{pdb_code[0]}_{pdb_code[3]}'
        
        if os.path.exists(f'{str(Path.cwd())+base_output_path}/Dock6/'):
            print(f'[INFO] The path {str(Path.cwd())+base_output_path}/Dock6/ already exists!')
            return
        
        gbvp = Dock6(ligand_input_path=f'{base_output_path}/Vina/', base_output_path=f'{base_output_path}/Molecules')
        gbvp.recover_better_conforms_of_vina(charge_type=charge_type)
        del gbvp
        
        dock6 = Dock6(dock6_path=dock6_app_path, ligand_input_path=f'{base_output_path}/Molecules/',
                    receptor_input_path=f'{base_input_path}Prepared/', base_output_path=f'{base_output_path}/Dock6',
                    pdb_code=pdb_code, density=density, radius=radius, distance=distance,
                    max_residues=plot_max_residues, conformer_search_type=conformer_search_type, mol_filename=mol_filename)
        
        dock6.prepare_surface()
        dock6.prepare_showbox()
        dock6.prepare_gridbox()
        dock6.prepare_minimization()
        dock6.prepare_footprint()
        dock6.plot_footprint_results()
        dock6.perform_dock6_evaluation()
        del dock6
        
    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_docking wrapper function', exc_info=True)


    


def perform_consensus(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,  dock6_app_path:str,
                    pdb_code:Optional[Tuple[str, str, str, str]]=None, density:Optional[float]=0.5, radius:Optional[float]=1.4,
                    distance:Optional[float]=10.0, conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex', 
                    plot_max_residues:Optional[int]=50, pH:Optional[float]=7.4, sizeof_box:Optional[List]=[24,24,24], 
                    exhaustiveness:Optional[int]=20, num_modes:Optional[int]=10, prepare_complex:Optional[bool]=True,
                    charge_type:Optional[ChargeType]='gas', mol_filename:Optional[str]='molecules') -> None:
     
                         
    try:

        input_path     = f'{base_input_path}/{target.replace(' ','')}/'
        output_path    = f'{base_output_path}/{target.replace(' ','')}' 
        dock6_app_path = dock6_app_path + '/' if not dock6_app_path.endswith('/') else dock6_app_path

        if pdb_code is None:    
            pdb_code = get_better_complex(input_path)
            print(f'[INFO] To perform docking with the better prepared complex for {target} - {pdb_code[0][0]}_{pdb_code[0][3]}.complex.pdb is necessary to refine loops fist')
            print(f'[INFO] Please, using chimera in Tools -> Structure Editing -> Model/Refine Loops')
            print(f'[INFO] After refine loops:')
            print(f'    1. Save the refined complex in the folder {input_path} as {pdb_code[0][0]}.pdb')
            print(f'    2. Execute the perform_docking again with pdb_code={pdb_code[0]} as input of function')
            exit(1)

        resnum = identify_ligand(f'{input_path[1:]}{pdb_code[0]}.pdb', pdb_code[1], pdb_code[3])
        pdb_codes=[(pdb_code[0], pdb_code[1], resnum, pdb_code[3])]

        if prepare_complex:
            dock = Docking(complex_input_path=input_path, output_path=f'{input_path}/Prepared/')
            dock.prepare_for_docking(pdb_codes=pdb_codes, charge_type=charge_type, pH=pH, redefine_centerofmass=True)
            del dock


        perform_docking_vina(base_input_path=base_input_path, target=target, base_output_path=base_output_path,
                             base_selected_mols=base_selected_mols, pdb_code=pdb_codes, pH=pH,
                             sizeof_box=sizeof_box, exhaustiveness=exhaustiveness, num_modes=num_modes,
                             mol_filename=mol_filename)
        

        time.sleep(10)
        dir_fd = os.open(output_path[1:] + '/Vina/', os.O_DIRECTORY)
        os.fsync(dir_fd)
        time.sleep(10)
        
        
        perform_docking_dock6(base_input_path=base_input_path, target=target, base_output_path=base_output_path,
                              base_selected_mols=base_selected_mols,  dock6_app_path=dock6_app_path, charge_type=charge_type,
                              pdb_code=pdb_code, density=density, radius=radius, distance=distance, 
                              conformer_search_type=conformer_search_type, plot_max_residues=plot_max_residues,
                              mol_filename=mol_filename)
        
        
        
        time.sleep(1)
        dir_fd = os.open(output_path[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)
        
    
        if os.path.exists(f'{str(Path.cwd())+output_path}/Vina/') and os.path.exists(f'{str(Path.cwd())+output_path}/Dock6/{conformer_search_type}/'):
            generate_consensus(base_input_path=output_path, base_output_path=output_path, target=target)
        
    
    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_consensus wrapper function', exc_info=True)

    finally:
        if os.path.exists(f'{str(Path.cwd())+output_path}/Molecules/'):
            shutil.rmtree(output_path[1:] + '/Molecules')
        
            