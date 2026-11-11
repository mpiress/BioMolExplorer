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
    if 'consensus_model' in df:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = Figure(figsize=(7, 6))
        canvas = FigureCanvasAgg(fig)
        ax = fig.subplots()
        ax.scatter(df['vina'], df['dock6'], c=df['normalized_score'], cmap='viridis', s=14)
        ax.set(xlabel='Relevância Vina', ylabel='Relevância DOCK6',
               title='Consenso normalizado por composto', xlim=(0, 1.02), ylim=(0, 1.02))
        fig.tight_layout()
        canvas.print_png(str(Path(output_path)/'correlation.png'))
        fig.clear()
        return
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
    from biomolexplorer.docking_data import NO_INTERSECTION
    from biomolexplorer.consensus_io import read_scores, write_table
    from biomolexplorer.calculation_cache import copy_artifact
    from biomolexplorer.progress import report_progress
    from biomolexplorer.consensus_model import consensus_rows
    output=Path(base_output_path);output.mkdir(parents=True,exist_ok=True)
    report_progress(key='consensus.reading')
    rows=consensus_rows(read_scores(base_vina_path or Path(base_input_path)/'Vina','vina'),
                        read_scores(base_dock6_path or Path(base_input_path)/'Dock6','dock6'),repulsion_weight)
    if not rows:
        import json
        (output/(target+'.csv')).unlink(missing_ok=True)
        (output/'consensus_report.json').write_text(json.dumps({'skipped_reason':NO_INTERSECTION},ensure_ascii=False),encoding='utf-8')
        return {'skipped_reason':NO_INTERSECTION,'rows':0}
    (output/'consensus_report.json').unlink(missing_ok=True)
    copied_context = {}
    for index,row in enumerate(rows):
        if index % 100 == 0:
            report_progress(key='consensus.exporting',completed=index,total=len(rows))
        for key in ('vina_pose','dock6_pose'):
            if not row.get(key):continue
            source=Path(row[key])
            if not source.is_file():
                row.pop(key,None);continue
            destination=output/'poses'/f'{index}_{key}{source.suffix}'
            destination.parent.mkdir(parents=True,exist_ok=True);copy_artifact(source,destination)
            row[key]=destination.relative_to(output).as_posix()
        for field in ('vina_receptor_file','dock6_receptor_file','vina_reference_file','dock6_reference_file','footprint_file'):
            if row.get(field):
                import hashlib
                source=Path(row[field])
                if not source.is_file():
                    row.pop(field,None);continue
                key=hashlib.sha256(str(source).encode()).hexdigest()[:16]
                destination=output/'context'/(key+'_'+source.name)
                destination.parent.mkdir(exist_ok=True)
                if source not in copied_context:
                    copy_artifact(source,destination)
                    copied_context[source] = destination
                row[field]=destination.relative_to(output).as_posix()
        row['conformer_file']=row.get('dock6_pose') or row.get('vina_pose') or ''
    df=DataFrame(rows)
    df=df.sort_values('normalized_score',ascending=False,kind='stable')
    write_table(df,output/(target+'.csv'))
    try:
        plot_scatter_comparison(df,str(output))
    except (OSError,ValueError,RuntimeError) as error:
        import json
        (output/'correlation.png').unlink(missing_ok=True)
        (output/'consensus_report.json').write_text(json.dumps(
            {'rows':len(df),'visualization_error':str(error)},ensure_ascii=False),encoding='utf-8')
        logger.warning('Consensus visualization failed; score table is complete.',exc_info=True)
    report_progress(key='consensus.done',completed=len(rows),total=len(rows))
    return df


def prepare_docking_receptor(base_input_path, target, output, pdb_code, ph, options, prepared, engine='vina'):
    """Reuse prepared receptors; raw PDBs use the shared preparation workflow."""
    from biomolexplorer.redocking_config import validate_preparation_settings,chimera_preparation
    options=chimera_preparation(options)
    validate_preparation_settings(options)
    from biomolexplorer.diagnostics import event
    if prepared:
        from biomolexplorer.docking_preparation import validate_receptor_outputs
        source=resolve_path(base_input_path)/target.replace(' ','')
        records=pdb_code
        if records is None:
            import pandas as pd
            records=pd.read_csv(source/'pdb_codes.csv')[['PDB_CODE','LIGAND','RESNUM','CHAIN']].values.tolist()
        else:
            records=list(records)
            if records and isinstance(records[0],str):records=[records]
        if engine=='vina' and any(not (source/'Prepared'/f'{r[0]}_{r[3]}.dockprep.pdbqt').is_file() and
                (source/'Prepared'/f'{r[0]}_{r[3]}.dockprep.mol2').is_file() for r in records):
            from biomolexplorer.docking_preparation import export_receptor_pdbqt
            destination=Path(output)/'Receptor'/target.replace(' ','')
            shutil.copytree(source,destination,dirs_exist_ok=True)
            export_receptor_pdbqt(destination/'Prepared',records)
            source=destination;base_input_path=str(destination.parent)
        validate_receptor_outputs(source/'Prepared',records,engine)
        event(logger,'receptor.reused','Receptor preparado reutilizado; sem nova preparação.',engine=engine)
        return base_input_path
    from wrappers.redocking import prepare_structures
    records=pdb_code
    if records and not isinstance(records[0],(list,tuple)):records=[records]
    destination=Path(output)/'Receptor'
    event(logger,'receptor.preparation.started','Preparando receptor bruto uma vez para este lote.',engine=engine)
    prepare_structures(base_input_path,target,str(destination),pdb_codes=records,
                       pH=ph,preparation_options=options)
    return str(destination)


def perform_docking_vina(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,
                         mol_filename:str, pdb_code:Optional[Tuple[str, str, str, str]]=None,
                         pH:Optional[float]=7.4, sizeof_box:Optional[List]=[24,24,24],
                         exhaustiveness:Optional[int]=20, num_modes:Optional[int]=10,
                         preparation_options=None, receptor_prepared:bool=True) -> None:

    try:

        output_root=Path(base_output_path)
        output_root.mkdir(parents=True,exist_ok=True)
        base_input_path=prepare_docking_receptor(base_input_path,target,output_root,pdb_code,pH,preparation_options,receptor_prepared)
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


        vina = DockVina(ligand_input_path=base_selected_mols, receptor_input_path=base_prepared_complexes,
                        output_path=base_input_mols, complex_input_path=base_prepared_complexes,
                        pdb_codes=pdb_code, centerofmasspath=base_prepared_complexes, sizeof_box=sizeof_box,
                        exhaustiveness=exhaustiveness, num_modes=num_modes, mol_filename=mol_filename)


        from biomolexplorer.docking_data import prepare_ligands, write_results
        from biomolexplorer.docking_accounting import reconcile
        from biomolexplorer.preparation_failures import CURRENT_REPORT,PreparationFailures
        failures=CURRENT_REPORT.get() or PreparationFailures()
        token=CURRENT_REPORT.set(failures)
        try:
            ligands=prepare_ligands(dataset,base_input_mols,'pdbqt',ph=pH,
                            preparation_options=preparation_options or {})
        finally:CURRENT_REPORT.reset(token)
        records=list(pdb_code)
        if not ligands:
            from biomolexplorer.docking_data import write_csv
            write_csv(output_root/'docking_results.csv',[],['molecule_chembl_id','canonical_smiles','receptor_id','engine','score','conformer_file'])
            reconcile('vina',dataset,records,[],output_root,failures.records if failures else ())
            return dict(skipped_reason='Nenhum composto pôde ser preparado para docking. Consulte a lista de falhas na preparação.')

        vina.set_ligandpath(base_input_mols)
        vina.selected_ligand_files=ligands
        vina.set_outputpath(f'{base_output_path}/Vina/')
        vina_summary=vina.docking(base_selected_mols)
        records=list(pdb_code)
        # A reused output directory can contain historical poses from removed
        # references. Only this request's completed calculations enter its table.
        poses=[Path(base_output_path,'Vina',Path(entry['configuration']).stem+'.lig.pdbqt')
               for entry in vina_summary.get('completed',[])]
        if poses:
            rows=write_results('vina',poses,dataset,records,output_root,prepared=base_prepared_complexes)
        else:
            from biomolexplorer.docking_data import write_csv
            rows=[]
            write_csv(output_root/'docking_results.csv',rows,['molecule_chembl_id','canonical_smiles','receptor_id','engine','score','conformer_file'])
        reconcile('vina',dataset,records,rows,output_root,failures.records if failures else (),vina_summary.get('failed',[]))
        del vina
        return dict(vina_completed=len(vina_summary.get('completed',[])),
                    docking_failures=len(vina_summary.get('failed',[])),
                    vina_reused=sum(bool(row.get('reused')) for row in vina_summary.get('completed',[])),
                    vina_failures=len(vina_summary.get('failed',[])),
                    vina_attempts=len(vina_summary.get('attempts',[])))

    except Exception as e:
        logger.error(f'Error during to perform the {target} in perform_docking wrapper function', exc_info=True)
        raise





def perform_docking_dock6(base_input_path:str, target:str, base_output_path:str, base_selected_mols:str,
                          dock6_app_path:str, charge_type:str, mol_filename:str,
                          pdb_code:Optional[Tuple[str, str, str, str]]=None, density:Optional[float]=0.5,
                          radius:Optional[float]=1.4, distance:Optional[float]=10.0,
                          conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex',
                          plot_max_residues:Optional[int]=50, base_vina_path:Optional[str]=None,
                          pH:float=7.4, preparation_options=None, receptor_prepared:bool=True, dock6_mode:str='auto') -> None:
    from biomolexplorer.docking_data import prepare_ligands,read_compounds,read_results,write_csv,write_results,binding_center
    from biomolexplorer.preparation_failures import CURRENT_REPORT,PreparationFailures
    from biomolexplorer.dock6_protocol import resolve_mode,mode_for_row
    failures=CURRENT_REPORT.get() or PreparationFailures()
    dataset=Path(base_selected_mols)/(mol_filename+'.csv')
    input_rows=read_compounds(dataset)
    resolve_mode(dock6_mode,True)
    input_modes=[mode_for_row(dock6_mode,row) for row in input_rows] if not base_vina_path else []
    output=Path(base_output_path)
    output.mkdir(parents=True,exist_ok=True)
    base_input_path=prepare_docking_receptor(base_input_path,target,output,pdb_code,pH,preparation_options,receptor_prepared,engine='dock6')
    data=Path(base_input_path)/target.replace(' ','')
    if pdb_code is None:
        import pandas as pd
        records=pd.read_csv(data/'pdb_codes.csv')[['PDB_CODE','LIGAND','RESNUM','CHAIN']].values.tolist()
    elif isinstance(pdb_code[0],(list,tuple)):records=list(pdb_code)
    else:records=[pdb_code]
    pose_input=bool(base_vina_path) or any(mode in ('guided','score_only','restrained') for mode in input_modes)
    pose_index={}
    def source_energy(row):
        value=row.get('score')
        return float(value) if value not in (None,'') else float('inf')
    if pose_input:
        candidates=read_results(base_vina_path,'vina') if base_vina_path else input_rows
        for original in candidates:
            row=dict(original);key=(row.get('receptor_id',''),row['molecule_chembl_id'])
            previous=pose_index.get(key)
            if previous is None or source_energy(row)<source_energy(previous):pose_index[key]=row
    shared_ligands=None;shared_rows={}
    if not pose_input:
        token=CURRENT_REPORT.set(failures)
        try:
            shared_ligands=prepare_ligands(dataset,output/'Molecules','mol2',charge_type=charge_type,ph=pH,
                preparation_options=preparation_options or {})
        finally:CURRENT_REPORT.reset(token)
        prepared_index=output/'Molecules'/'compounds.csv'
        shared_rows={row['molecule_chembl_id']:row for row in read_compounds(
            prepared_index if prepared_index.is_file() else dataset)}
    current_tables=[]
    for record in records:
        receptor=f'{record[0]}_{record[3]}'
        root=output/target.replace(' ','')/receptor
        prepared_dataset=dataset
        if pose_input:
            from biomolexplorer.docking_data import best_vina_model
            selected_rows=[]
            compounds={row['molecule_chembl_id']:row for row in input_rows
                       if not row.get('receptor_id') or row['receptor_id']==receptor or base_vina_path}
            for row_number,original in enumerate(compounds.values(),2):
                row=dict(original)
                if not base_vina_path and not (row.get('conformer_file') or row.get('pose_file')):
                    resolve_mode(dock6_mode,False)
                    selected_rows.append(row);continue
                with failures.compound(dataset,row_number,row,root,[],step='Seleção da pose Vina'):
                    code=row['molecule_chembl_id']
                    matches=[pose_index[key] for key in ((receptor,code),('',code),
                        ('',f'{record[0]}_{record[1]}_{record[2]}{record[3]}_'+code)) if key in pose_index]
                    if matches:
                        selected=min(matches,key=source_energy)
                        pose=Path(selected.get('conformer_file') or selected.get('pose_file',''))
                        if not pose.is_absolute():pose=dataset.parent/pose
                        if not pose.is_file():raise ValueError('Arquivo da pose selecionada ausente: '+str(pose))
                        if selected.get('engine')=='vina' or (pose.suffix.lower()=='.pdbqt' and 'REMARK VINA RESULT:' in pose.read_text()):
                            index,value=best_vina_model(pose)
                        else:
                            index,value=int(selected.get('conformer_index') or 1),selected.get('score')
                        row.update(conformer_file=str(pose.resolve()),conformer_index=index,receptor_id=receptor)
                        if value not in (None,''):row['score']=float(value)
                        for key in ('prepared_pdbqt','prepared_mol2','docking_engines','prepared_origin'):row.pop(key,None)
                    else:raise ValueError('As poses Vina não correspondem ao receptor e aos compostos selecionados.')
                    selected_rows.append(row)
            root.mkdir(parents=True,exist_ok=True);prepared_dataset=root/'input_compounds.csv'
            fields=list(dict.fromkeys(k for row in selected_rows for k in row)) or ['molecule_chembl_id','canonical_smiles']
            write_csv(prepared_dataset,selected_rows,fields);failures.write(root)
        center=binding_center(data/'Prepared',record)
        token=CURRENT_REPORT.set(failures)
        try:
            if shared_ligands is None:
                rows=read_compounds(prepared_dataset)
                for row in rows:row['dock6_operation']=mode_for_row(dock6_mode,row)
                write_csv(prepared_dataset,rows,list(dict.fromkeys(k for row in rows for k in row)) or ['molecule_chembl_id','canonical_smiles'])
                ligands=prepare_ligands(prepared_dataset,root/'Molecules','mol2',charge_type=charge_type,center=center,ph=pH,
                                preparation_options=preparation_options or {})
            else:
                from biomolexplorer.calculation_cache import copy_artifact
                from biomolexplorer.docking_data import center_mol2
                ligands=[]
                for source in shared_ligands:
                    source=Path(source);destination=root/'Molecules'/source.name
                    code=source.name.removesuffix('.lig.mol2');row=shared_rows.get(code,{})
                    copy_artifact(source,destination)
                    if row.get('prepared_origin') in ('smiles','library'):
                        center_mol2(destination,center)
                    ligands.append(destination)
        finally:CURRENT_REPORT.reset(token)
        if not ligands:continue
        dock6=Dock6(dock6_path=dock6_app_path,ligand_input_path=str(root/'Molecules'),
                    receptor_input_path=str(data/'Prepared'),base_output_path=str(root/'Dock6'),
                    pdb_code=receptor,density=density,radius=radius,distance=distance,
                    max_residues=plot_max_residues,conformer_search_type=conformer_search_type,mol_filename=mol_filename,
                    binding_site_center=center)
        dock6.set_compound_report(prepared_dataset,failures)
        from biomolexplorer.dock6_validation import compound_code
        dock6.selected_compound_ids={compound_code(path) for path in ligands}
        candidate_rows={row['molecule_chembl_id']:row for row in read_compounds(prepared_dataset)}
        dock6.compound_modes={code:mode_for_row(dock6_mode,candidate_rows[code])
                             for code in dock6.selected_compound_ids}
        try:
            dock6.validate_setup()
            if not dock6._ligand_files(dock6.ligandpath):continue
            dock6.prepare_site()
            dock6.prepare_minimization();dock6.perform_dock6_evaluation()
            completed_compounds=dock6.selected_compound_ids-dock6.excluded_compounds
            if dock6._ligand_files(dock6.ligandpath):
                dock6.prepare_footprint(docked=True);dock6.plot_footprint_results()
        finally:
            failures.write(output)
            dock6.export_results()
        selected=completed_compounds
        poses=[path for path in (root/'Dock6'/conformer_search_type).glob('*_scored.mol2')
               if path.name.removesuffix('_scored.mol2') in selected]
        if not poses:continue
        scored_rows=write_results('dock6',poses,prepared_dataset,[record],root,prepared=data/'Prepared')
        from biomolexplorer.redocking_dock6 import heavy_atom_rmsd
        for row in scored_rows:
            code=row['molecule_chembl_id'];mode=dock6.compound_modes[code]
            row['dock6_mode']=mode
            if mode in ('guided','score_only','restrained'):
                row['input_pose_rmsd']=heavy_atom_rmsd(root/'Molecules'/(code+'.lig.mol2'),root/row['conformer_file'])
                if mode=='score_only' and row['input_pose_rmsd']>0.002:
                    raise ValueError('O modo score-only alterou as coordenadas da referência: '+code)
        write_csv(root/'docking_results.csv',scored_rows)
        current_tables.append(root/'docking_results.csv')
    # One table represents all receptors so automatic consumers keep every result.
    rows=[]
    for table in sorted(current_tables):
        for row in read_compounds(table):
            for field in ('conformer_file','receptor_file','reference_file','footprint_file'):
                if row.get(field):row[field]=(table.parent/row[field]).relative_to(output).as_posix()
            rows.append(row)
    fields=list(dict.fromkeys(k for row in rows for k in row)) or ['molecule_chembl_id','canonical_smiles','receptor_id','engine','score','conformer_file']
    write_csv(output/'docking_results.csv',rows,fields)
    from biomolexplorer.docking_accounting import reconcile
    reconcile('dock6',dataset,records,rows,output,failures.records)
    failures.write(output)
    if not rows:return dict(skipped_reason='Nenhum composto produziu resultado utilizável no DOCK6. Consulte a lista de compostos não utilizados.')
    return dict(docking_failures=len(failures.records))


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
