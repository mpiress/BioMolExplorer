from biomolexplorer.storage import sync_directory
from biomolexplorer.paths import directory, resolve_path, worker_count
#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
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
import math
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
        raise




def normalized_score(values, method):
    values=values.abs()
    denominator=values.std() if method=='z-score' else values.max()-values.min()
    if not math.isfinite(denominator) or denominator==0:
        return values*0
    return (values-(values.mean() if method=='z-score' else values.min()))/denominator


def plot_scatter_comparison(df:DataFrame, output_path:str):
    df_minmax = df.copy()
    df_zscore = df.copy()

    # Normalização Min-Max
    df_minmax['vina'] = df_minmax['vina'].abs()
    df_minmax['dock6'] = df_minmax['dock6'].abs()
    df_minmax['vina'] = normalized_score(df_minmax['vina'],'min-max')
    df_minmax['dock6'] = normalized_score(df_minmax['dock6'],'min-max')

    # Normalização Z-score
    df_zscore['vina'] = df_zscore['vina'].abs()
    df_zscore['dock6'] = df_zscore['dock6'].abs()
    df_zscore['vina'] = normalized_score(df_zscore['vina'],'z-score')
    df_zscore['dock6'] = normalized_score(df_zscore['dock6'],'z-score')

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    sns.scatterplot(x=df_minmax['vina'], y=df_minmax['dock6'], ax=axes[0], color='blue')
    axes[0].set_title("Min-Max Correlaction")

    sns.scatterplot(x=df_zscore['vina'], y=df_zscore['dock6'], ax=axes[1], color='red')
    axes[1].set_title("Z-score Correlaction")

    plt.tight_layout()
    plt.savefig(directory(output_path) + '/correlation.png')
    plt.close(fig)




def generate_consensus(base_input_path:str, base_output_path:str, target:str, repulsion_weight: float = 1.0,
                       base_vina_path:str=None, base_dock6_path:str=None):
    from biomolexplorer.docking_data import read_results, consensus_rows, write_csv, NO_INTERSECTION
    output=Path(base_output_path);output.mkdir(parents=True,exist_ok=True)
    rows=consensus_rows(read_results(base_vina_path or Path(base_input_path)/'Vina','vina'),
                        read_results(base_dock6_path or Path(base_input_path)/'Dock6','dock6'),repulsion_weight)
    if not rows:
        import json
        (output/(target+'.csv')).unlink(missing_ok=True)
        (output/'consensus_report.json').write_text(json.dumps({'skipped_reason':NO_INTERSECTION},ensure_ascii=False),encoding='utf-8')
        return {'skipped_reason':NO_INTERSECTION,'rows':0}
    (output/'consensus_report.json').unlink(missing_ok=True)
    for index,row in enumerate(rows):
        for key in ('vina_pose','dock6_pose'):
            source=Path(row[key]);destination=output/'poses'/f'{index}_{key}{source.suffix}'
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
            row[key]=destination.relative_to(output).as_posix()
        row['conformer_file']=row['dock6_pose']
    df=DataFrame(rows)
    df['z-score']=(normalized_score(df['vina'],'z-score')+normalized_score(df['dock6'],'z-score'))/2
    df['min-max']=(normalized_score(df['vina'],'min-max')+normalized_score(df['dock6'],'min-max'))/2
    df[['vina','dock6','z-score','min-max']]=df[['vina','dock6','z-score','min-max']].round(3)
    df=df.sort_values(['z-score','min-max'],ascending=False)
    plot_scatter_comparison(df,str(output))
    df.to_csv(output/(target+'.csv'),index=False)
    return df


def perform_docking_vina(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,
                         mol_filename:str, pdb_code:Optional[Tuple[str, str, str, str]]=None,
                         pH:Optional[float]=7.4, sizeof_box:Optional[List]=[24,24,24],
                         exhaustiveness:Optional[int]=20, num_modes:Optional[int]=10) -> None:

    try:

        output_root=Path(base_output_path)
        dataset=Path(base_selected_mols)/(mol_filename+'.csv')
        base_prepared_complexes = f'{base_input_path}/{target.replace(' ','')}/Prepared/'
        base_input_path         = f'{base_input_path}/{target.replace(' ','')}/'
        base_input_mols         = f'{base_output_path}/Molecules/'
        base_selected_mols      = f'{base_selected_mols}/'
        base_output_path        = f'{base_output_path}/{target.replace(' ','')}'

        if pdb_code is None:
            metadata = fileHandling(input_path=base_input_path, output_path=base_input_path).csv_to_dataframe('pdb_codes')
            pdb_code = metadata[['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN']].to_records(index=False)
        elif pdb_code and not isinstance(pdb_code[0], (list, tuple)):
            pdb_code = [pdb_code]


        if os.path.exists(f'{directory(base_output_path)}/Vina/'):
            print(f'[INFO] The path {base_output_path}/Vina/ already exists!')
            return

        vina = DockVina(ligand_input_path=base_selected_mols, receptor_input_path=base_prepared_complexes,
                        output_path=base_input_mols, complex_input_path=base_prepared_complexes,
                        pdb_codes=pdb_code, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box,
                        exhaustiveness=exhaustiveness, num_modes=num_modes, mol_filename=mol_filename)


        from biomolexplorer.docking_data import prepare_ligands, write_results
        prepare_ligands(dataset,base_input_mols,'pdbqt',ph=pH)

        vina.set_ligandpath(base_input_mols)
        vina.set_outputpath(f'{base_output_path}/Vina/')
        vina.docking(base_selected_mols)
        records=list(pdb_code)
        write_results('vina',Path(base_output_path,'Vina').glob('*.pdbqt'),dataset,records,output_root)
        del vina

    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_docking wrapper function', exc_info=True)
        raise





def perform_docking_dock6(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,
                          dock6_app_path:str, charge_type:str, mol_filename:str,
                          pdb_code:Optional[Tuple[str, str, str, str]]=None, density:Optional[float]=0.5,
                          radius:Optional[float]=1.4, distance:Optional[float]=10.0,
                          conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex',
                          plot_max_residues:Optional[int]=50, base_vina_path:Optional[str]=None) -> None:
    from biomolexplorer.docking_data import prepare_ligands,read_compounds,read_results,write_csv,write_results,binding_center
    output=Path(base_output_path);data=Path(base_input_path)/target.replace(' ','')
    dataset=Path(base_selected_mols)/(mol_filename+'.csv')
    if pdb_code is None:
        import pandas as pd
        records=pd.read_csv(data/'pdb_codes.csv')[['PDB_CODE','LIGAND','RESNUM','CHAIN']].values.tolist()
    elif isinstance(pdb_code[0],(list,tuple)):records=list(pdb_code)
    else:records=[pdb_code]
    for record in records:
        receptor=f'{record[0]}_{record[3]}'
        root=output/target.replace(' ','')/receptor
        prepared_dataset=dataset
        if base_vina_path:
            poses=read_results(base_vina_path,'vina');rows=read_compounds(dataset)
            for row in rows:
                matches=[p for p in poses if p['molecule_chembl_id']==row['molecule_chembl_id'] and p['receptor_id'] in ('',receptor)]
                # Legacy pose names include receptor and reference-ligand prefixes.
                matches+= [p for p in poses if p['molecule_chembl_id']==f'{record[0]}_{record[1]}_{record[2]}{record[3]}_'+row['molecule_chembl_id']]
                if matches:
                    row['conformer_file']=min(matches,key=lambda p:p['score'])['conformer_file']
                    for key in ('prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin'):row.pop(key,None)
                elif not row.get('conformer_file'):raise ValueError('As poses Vina não correspondem ao receptor e aos compostos selecionados.')
            root.mkdir(parents=True,exist_ok=True);prepared_dataset=root/'input_compounds.csv';write_csv(prepared_dataset,rows)
        center=binding_center(data/'Prepared',record)
        prepare_ligands(prepared_dataset,root/'Molecules','mol2',charge_type=charge_type,center=center)
        dock6=Dock6(dock6_path=dock6_app_path,ligand_input_path=str(root/'Molecules'),
                    receptor_input_path=str(data/'Prepared'),base_output_path=str(root/'Dock6'),
                    pdb_code=receptor,density=density,radius=radius,distance=distance,
                    max_residues=plot_max_residues,conformer_search_type=conformer_search_type,mol_filename=mol_filename,
                    binding_site_center=center)
        dock6.prepare_surface();dock6.prepare_showbox();dock6.prepare_gridbox()
        dock6.prepare_minimization();dock6.prepare_footprint();dock6.plot_footprint_results();dock6.perform_dock6_evaluation()
        dock6.export_results()
        write_results('dock6',(root/'Dock6'/conformer_search_type).glob('*_scored.mol2'),dataset,[record],root)
    # One table represents all receptors so automatic consumers keep every result.
    rows=[]
    for table in sorted(output.rglob('docking_results.csv')):
        for row in read_compounds(table):
            row['conformer_file']=(table.parent/row['conformer_file']).relative_to(output).as_posix();rows.append(row)
    write_csv(output/'docking_results.csv',rows)


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
            raise ValueError('Required inputs or directories are missing; see logs')

        resnum = identify_ligand(f'{directory(input_path)}{pdb_code[0]}.pdb', pdb_code[1], pdb_code[3])
        pdb_codes=[(pdb_code[0], pdb_code[1], resnum, pdb_code[3])]

        if prepare_complex:
            dock = Docking(complex_input_path=input_path, output_path=f'{input_path}/Prepared/')
            dock.prepare_for_docking(pdb_codes=pdb_codes, charge_type=charge_type, pH=pH, redefine_centerofmass=True)
            del dock


        #perform_docking_vina(base_input_path=base_input_path, target=target, base_output_path=base_output_path,
        #                     base_selected_mols=base_selected_mols, pdb_code=pdb_codes, pH=pH,
        #                     sizeof_box=sizeof_box, exhaustiveness=exhaustiveness, num_modes=num_modes,
        #                     mol_filename=mol_filename)

        sync_directory(directory(output_path) + '/Vina/')

        perform_docking_dock6(base_input_path=base_input_path, target=target, base_output_path=base_output_path,
                              base_selected_mols=base_selected_mols,  dock6_app_path=dock6_app_path, charge_type=charge_type,
                              pdb_code=pdb_code, density=density, radius=radius, distance=distance,
                              conformer_search_type=conformer_search_type, plot_max_residues=plot_max_residues,
                              mol_filename=mol_filename)

        sync_directory(directory(output_path))


        if os.path.exists(f'{directory(output_path)}/Vina/') and os.path.exists(f'{directory(output_path)}/Dock6/{conformer_search_type}/'):
            generate_consensus(base_input_path=output_path, base_output_path=output_path, target=target)


    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_consensus wrapper function', exc_info=True)
        raise

    finally:
        if os.path.exists(f'{directory(output_path)}/Molecules/'):
            shutil.rmtree(directory(output_path) + '/Molecules')
