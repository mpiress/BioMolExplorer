"""Reconcile requested receptor/compound pairs without inventing successful scores."""
from pathlib import Path

from .docking_data import read_compounds, write_csv
from .storage import write_json


def reconcile(engine, dataset, records, rows, output, failures=(), execution_failures=()):
    output = Path(output)
    compounds = read_compounds(dataset)
    identifiers = {row['molecule_chembl_id'] for row in compounds}
    input_identities={(row.get('receptor_id','') if row.get('conformer_file') or row.get('pose_file') else '',
                       row['molecule_chembl_id']) for row in compounds}
    receptors = {f'{record[0]}_{record[3]}' for record in records}
    successful = set()
    for row in rows:
        pair = (row['receptor_id'], row['molecule_chembl_id'])
        if pair[0] not in receptors or pair[1] not in identifiers:
            raise ValueError('Resultado de docking fora das entradas selecionadas: '+str(pair))
        if pair in successful:
            raise ValueError('Resultado de docking duplicado: '+str(pair))
        successful.add(pair)
    preparation = {}
    for failure in failures:
        code=failure['compound_identifier']
        context=Path(failure.get('preparation_context','')).parts
        affected=[receptor for receptor in receptors if receptor in context] or receptors
        status='docking_failed' if 'DOCK6' in failure.get('preparation_step','') else 'preparation_failed'
        for receptor in affected:
            preparation[(receptor,code)]=(status,failure['reason'])
    engine_errors = {}
    for failure in execution_failures:
        code = failure.get('compound_identifier', '')
        configuration = Path(failure.get('configuration', '')).stem
        for receptor in receptors:
            pdb, chain = receptor.rsplit('_', 1)
            if configuration.startswith(pdb+'_') and ('_'+chain+'_' in configuration or chain+'_'+code in configuration):
                engine_errors[(receptor, code)] = failure.get('error', 'Falha do motor de docking.')
    missing = []
    for receptor in sorted(receptors):
        for code in sorted(identifiers):
            if (receptor, code) in successful:
                continue
            failure=preparation.get((receptor,code))
            reason = failure[1] if failure else engine_errors.get((receptor, code))
            missing.append(dict(receptor_id=receptor, molecule_chembl_id=code,
                status=failure[0] if failure else 'docking_failed' if reason else 'unaccounted',
                reason=reason or 'A entrada selecionada não produziu uma pose; consulte o log do motor.'))
    report = dict(schema_version=1, engine=engine, input_rows=len(compounds),
        unique_compounds=len(identifiers), duplicate_input_rows=len(compounds)-len(input_identities),
        receptors=len(receptors), expected_pairs=len(identifiers)*len(receptors),
        successful_pairs=len(successful), unsuccessful_pairs=len(missing),
        unaccounted_pairs=sum(row['status']=='unaccounted' for row in missing),
        per_receptor={receptor: dict(expected=len(identifiers),
            successful=sum(pair[0]==receptor for pair in successful),
            unsuccessful=sum(row['receptor_id']==receptor for row in missing)) for receptor in sorted(receptors)})
    write_json(report, output/'docking_accounting.json')
    write_csv(output/'docking_unprocessed.csv', missing,
        ['receptor_id', 'molecule_chembl_id', 'status', 'reason'])
    return report
