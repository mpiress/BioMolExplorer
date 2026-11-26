"""Pose validation wrapper: physical and geometric screening of docked poses.

Criteria follow PoseBusters (Buttenschoen et al., Chem. Sci. 15, 3130, 2024):
chemical validity, internal geometry against a relaxed conformer, planarity,
relative internal energy and intermolecular steric clash.

Scientific dependencies are imported inside the functions; importing this
module costs nothing beyond the standard library.
"""
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from biomolexplorer.paths import worker_count

REPORT = 'pose_validation.csv'
SUMMARY = 'pose_validation.json'
FIELDS = ('pose_id','receptor_id','ligand_id','compound_id','engine','score','pose_valid',
          'valid_chemistry','valid_connectivity','valid_geometry','valid_planarity','valid_energy','valid_clash',
          'min_vdw_ratio','min_vdw_overlap','min_heavy_distance','clashing_contacts','delta_energy','energy_method',
          'max_bond_deviation','max_angle_deviation','max_planarity_deviation','planar_groups',
          'atoms','bond_orders','resolved_by','rejection_reasons','rejection_details','error',
          'pose_file','receptor_file','source_table')


def _logger():
    from biomolexplorer.diagnostics import get_logger
    return get_logger('pose_validation')


def validate_poses(base_input_path: str, base_output_path: str, target: str = 'MeuAlvo',
                   base_selected_mols: str = None, pose_inputs: list = None,
                   selected_receptors: list = None, selected_ligands: list = None,
                   excluded_pairs: list = None,
                   strict_filter: bool = True, vdw_clash_scale: float = 0.75,
                   delta_energy_limit: float = 100.0, bond_tolerance: float = 0.25,
                   angle_tolerance: float = 25.0, planarity_tolerance: float = 0.25,
                   minimization_steps: int = 400):
    from biomolexplorer import pose_sources
    from biomolexplorer.diagnostics import event
    from biomolexplorer.pose_geometry import describe, evaluate_pose
    from biomolexplorer.progress import report_progress
    logger = _logger()
    output = Path(base_output_path); output.mkdir(parents=True,exist_ok=True)
    records = [dict(record) for record in (pose_inputs or [])]
    if not records:
        # Standalone and CLI execution resolve the poses from the input folder.
        records = pose_sources.from_directory(base_selected_mols or base_input_path,base_input_path)
    receptors = {str(value) for value in (selected_receptors or []) if str(value).strip()}
    ligands = {str(value) for value in (selected_ligands or []) if str(value).strip()}
    if receptors:
        records = [r for r in records if str(r.get('receptor_id') or '') in receptors]
    if ligands:
        records = [r for r in records if {str(r.get('compound_id') or ''),str(r.get('ligand_id') or '')} & ligands]
    excluded = {str(value) for value in (excluded_pairs or []) if str(value).strip()}
    if excluded:
        # One complex the editor dropped by its row, as compound:receptor.  The
        # block receives the selection already filtered; a standalone or CLI
        # execution applies the same exclusions over the resolved poses.
        from biomolexplorer.pose_catalog import complex_key
        records = [r for r in records if complex_key(r.get('compound_id') or r.get('ligand_id'),
                                                     r.get('receptor_id')) not in excluded]
    records = [r for r in records if r.get('pose')]
    if not records:
        raise ValueError('Nenhuma pose foi encontrada para validar. Conecte resultados de Vina, DOCK6 ou Consenso, '
                         'ou selecione os arquivos de pose nas configurações do bloco.')
    event(logger,'pose_validation.started',f'Pose validation started for {len(records)} poses',
          poses=len(records),configuration=('strict' if strict_filter else 'annotate'))
    report_progress(key='pose.validating',completed=0,total=len(records))
    options = dict(vdw_clash_scale=vdw_clash_scale,delta_energy_limit=delta_energy_limit,
                   bond_tolerance=bond_tolerance,angle_tolerance=angle_tolerance,
                   planarity_tolerance=planarity_tolerance,minimization_steps=minimization_steps)

    def evaluate(record):
        try:
            return record,evaluate_pose(record,**options)
        except Exception as error:                      # noqa: BLE001 - one pose must not stop the batch
            logger.warning('Pose evaluation failed for %s',record.get('pose'),exc_info=True)
            return record,{'pose_valid':False,'rejection_reasons':['unreadable'],'error':str(error)[:400]}

    rows = []; completed = 0
    from rdkit import rdBase
    # Reading a PDBQT tries the explicit-hydrogen molecule before the heavy-atom
    # fallback; RDKit prints those expected attempts as errors on stderr.  The
    # outcome of every pose is recorded in the report (bond_orders, error,
    # rejection reasons), so the parser diagnostics are muted for the batch.
    muted = rdBase.BlockLogs()
    with ThreadPoolExecutor(max_workers=max(1,min(worker_count(),len(records)))) as pool:
        for record,result in pool.map(evaluate,records):
            completed += 1
            if completed % 10 == 0 or completed == len(records):
                report_progress(key='pose.validating',completed=completed,total=len(records))
            reasons = result.get('rejection_reasons') or []
            rows.append({**{key:'' for key in FIELDS},
                         'pose_id':record.get('pose_id') or Path(record['pose']).name.split('.',1)[0],
                         'receptor_id':record.get('receptor_id') or '','ligand_id':record.get('ligand_id') or '',
                         'compound_id':record.get('compound_id') or '','engine':record.get('engine') or '',
                         'score':record.get('score') if record.get('score') not in (None,'') else '',
                         'resolved_by':record.get('resolved_by') or 'binding',
                         'pose_file':str(record['pose']),'receptor_file':str(record.get('receptor') or ''),
                         'source_table':str(record.get('table') or ''),
                         **{key:value for key,value in result.items() if key in FIELDS and key != 'rejection_reasons'},
                         'rejection_reasons':';'.join(reasons),'rejection_details':describe(reasons)})
    del muted
    rows.sort(key=lambda row:(not row['pose_valid'],str(row['receptor_id']),str(row['compound_id']),str(row['pose_id'])))
    with (output/REPORT).open('w',encoding='utf-8',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(FIELDS),extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    published = _publish(rows,output,strict_filter)
    counts = {}
    for row in rows:
        for reason in (row['rejection_reasons'].split(';') if row['rejection_reasons'] else []):
            counts[reason] = counts.get(reason,0)+1
    valid = sum(1 for row in rows if row['pose_valid'])
    summary = {'target':target,'poses':len(rows),'valid_poses':valid,'invalid_poses':len(rows)-valid,
               'published_poses':len(published),'strict_filter':bool(strict_filter),
               'rejection_counts':counts,'criteria':'posebusters_inspired_v1','tolerances':options,
               'unresolved_inputs':sum(1 for row in rows if row['resolved_by'] == 'unresolved'),
               'identifier_resolved_inputs':sum(1 for row in rows if row['resolved_by'] == 'identifier')}
    (output/SUMMARY).write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True),encoding='utf-8')
    event(logger,'pose_validation.completed',
          f'Pose validation completed: {valid}/{len(rows)} valid, {len(published)} published',
          **{key:summary[key] for key in ('poses','valid_poses','invalid_poses','published_poses')})
    report_progress(key='pose.done',completed=len(rows),total=len(rows))
    return summary


def _publish(rows, output, strict_filter):
    """Copy the poses handed to the next blocks, keeping their native names."""
    from biomolexplorer.calculation_cache import copy_artifact
    poses = output/'poses'; context = output/'receptors'; published = []; owners = {}
    for row in rows:
        if strict_filter and not row['pose_valid']:
            continue
        source = Path(row['pose_file'])
        if not source.is_file():
            continue
        poses.mkdir(parents=True,exist_ok=True)
        # Keep the producer's filename when it already identifies the compound;
        # a consensus pose named by row index gets the readable pose identity.
        name = (source.name if str(row['compound_id'] or '').casefold() in source.name.casefold()
                else row['pose_id']+''.join(source.suffixes[-2:] if '.lig' in source.suffixes else source.suffixes[-1:]))
        # DOCK6 names a pose after its compound, so one compound docked to
        # three receptors arrives three times under the same name.  Every pose
        # keeps its own file: equal sizes are not the same pose, and
        # overwriting one with another would hand the next block, and the
        # validation table, a complex that was never evaluated.
        owner = (str(row['receptor_id'] or ''),str(row['pose_id'] or ''),str(source))
        qualified = [name]+[prefix+'_'+name for prefix in
                            (str(row['receptor_id'] or ''),str(row['pose_id'] or '')) if prefix]
        qualified.append(str(len(published))+'_'+name)
        destination = next((poses/candidate for candidate in qualified
                            if owners.setdefault(candidate,owner) == owner),poses/name)
        copy_artifact(source,destination)
        published.append(destination)
        receptor = Path(row['receptor_file']) if row['receptor_file'] else None
        if receptor and receptor.is_file() and receptor.suffix.lower() == '.pdb':
            context.mkdir(parents=True,exist_ok=True)
            copy = context/receptor.name
            if not copy.exists():
                copy_artifact(receptor,copy)
    return published
