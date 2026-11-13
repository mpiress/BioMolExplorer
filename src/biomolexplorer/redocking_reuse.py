"""Verify completed redocking outcomes, including recorded preparation exclusions."""
import csv
import math
from pathlib import Path


def accounted_preparation_exclusions(entry,root):
    from .result_state import request,input_digest
    from .redocking_config import pair_key
    if entry.get('status')!='succeeded' or any(entry.get(k) for k in
            ('outputs_removed','dependency_obsolete','partial_results')):
        return False
    root=Path(root).resolve()
    total=0
    try:
        for unit in entry.get('batches') or [entry]:
            if unit.get('status')!='succeeded' or any(unit.get(k) for k in
                    ('outputs_removed','dependency_obsolete','partial_results')):
                return False
            config=unit.get('configuration',entry.get('configuration',{}))
            if config.get('operation')!='redocking':return False
            expected=set(request(config)['axes'].get('pairs',[]))
            if not expected:return False
            engine=config.get('parameters',{}).get('redocking_engine','vina')
            manifest=unit.get('artifact_manifest') or entry.get('artifact_manifest',{})
            outcomes={};failures={}
            for value in unit.get('artifacts',[]):
                path=Path(value)
                if path.name not in ('pdb_codes.csv','redocking_failures.csv'):continue
                path=path.resolve()
                if not path.is_relative_to(root) or input_digest(path)!=manifest.get(str(path)):
                    return False
                with path.open(encoding='utf-8-sig',newline='') as stream:
                    reader=csv.DictReader(stream)
                    if not {'PDB_CODE','LIGAND','RESNUM','CHAIN','STATUS'}<=set(reader.fieldnames or []):continue
                    target=failures if path.name=='redocking_failures.csv' else outcomes
                    for row in reader:
                        key=pair_key([row[k] for k in ('PDB_CODE','LIGAND','RESNUM','CHAIN')])
                        if key in target:return False
                        target[key]=row
            if set(outcomes)!=expected:return False
            failed=set()
            for key,row in outcomes.items():
                if row['STATUS']=='succeeded':
                    rmsd=float(row['RMSD'])
                    if row.get('ENGINE',engine)!=engine or not math.isfinite(rmsd) or rmsd<0:return False
                elif row['STATUS'] in ('preparation_failed','UNSUPPORTED'):
                    if (row.get('RMSD','').strip() or not row.get('ERROR','').strip() or
                            row.get('FAILED_STEP') not in ('Preparação do complexo','Classificação do complexo')):
                        return False
                    failed.add(key)
                else:return False
            if set(failures)!=failed or len(failed)!=unit.get('redocking_failures',0):return False
            for key,row in failures.items():
                if (row['STATUS'] not in ('preparation_failed','UNSUPPORTED') or
                        any(row.get(k)!=outcomes[key].get(k) for k in ('ERROR','FAILED_STEP'))):return False
            total+=len(failed)
        return total>0 and total==entry.get('redocking_failures',total)
    except (OSError,ValueError,TypeError,KeyError,AttributeError,RuntimeError):
        return False
