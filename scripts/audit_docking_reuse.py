"""Audit saved docking output and optionally consolidate verified legacy caches.

No engine is invoked and no original scientific output is rewritten or deleted.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sqlite3
import time

from biomolexplorer.artifacts import annotate_files
from biomolexplorer.calculation_cache import CalculationCache,copy_artifact
from biomolexplorer.ligand_cache import validate_prepared_ligand
from biomolexplorer.stage_cache import file_digest
from biomolexplorer.storage import write_json


def audit(database,run_id,destination,consolidate=False):
    started=time.monotonic()
    with sqlite3.connect('file:'+str(Path(database).resolve())+'?mode=ro',uri=True) as db:
        run=db.execute('SELECT id,project_id,status FROM runs WHERE id=?',(run_id,)).fetchone()
        if run is None:raise ValueError('Run not found')
        query="""SELECT json_remove(stage.value,'$.artifact_records','$.execution_manifest',
            '$.processing_runtime','$.artifact_manifest') FROM runs r,json_each(r.stages) stage
            WHERE r.id=? AND json_extract(stage.value,'$.operation') IN ('docking_vina','docking_dock6')"""
        stages=[json.loads(row[0]) for row in db.execute(query,(run_id,))]
    report=dict(run_id=run_id,run_status=run[2],stages=[],cache_imports=[])
    cache_roots=set();project_root=None
    for stage in stages:
        artifacts=stage.get('artifacts',[])
        if artifacts:
            parts=Path(artifacts[0]).parts
            project_root=Path(*parts[:parts.index('runs')])
        summaries=[Path(p) for p in artifacts if Path(p).name=='docking_results.csv']
        from biomolexplorer.docking_data import result_tables
        summaries=result_tables(summaries)
        molecules=set();pairs=Counter();engines=Counter();rows=0
        for table in summaries:
            with table.open(newline='',encoding='utf-8-sig') as stream:
                for row in csv.DictReader(stream):
                    rows+=1;molecules.add(row['molecule_chembl_id'])
                    pairs[(row.get('receptor_id'),row['molecule_chembl_id'])]+=1
                    engines[row.get('engine')]+=1
        public=[p for p in artifacts if 'Molecules' not in Path(p).parts and not Path(p).name.startswith('.')]
        for p in artifacts:
            parts=Path(p).parts
            if '.biomolexplorer' in parts:
                cache_roots.add(Path(*parts[:parts.index('.biomolexplorer')+1])/'cache')
        before=time.monotonic()
        annotated=annotate_files(stage,[dict(path=str(p),name=p.name) for p in summaries],run_id)
        report['stages'].append(dict(stage_id=stage['id'],operation=stage['operation'],
            batches=len(stage.get('batches',[])),artifact_count=len(artifacts),
            public_artifact_count=len(public),intermediate_artifacts=len(artifacts)-len(public),
            summary_tables=len(summaries),rows=rows,unique_molecules=len(molecules),
            unique_receptor_compound_pairs=len(pairs),duplicate_pair_rows=sum(n-1 for n in pairs.values()),
            engines=dict(engines),annotation_seconds=round(time.monotonic()-before,4),
            annotation_count=len(annotated)))
    if consolidate and project_root:
        target=project_root/'.biomolexplorer'/'cache'
        preparation=target/'ligand-preparation';preparation.mkdir(parents=True,exist_ok=True)
        seen=set()
        for root in sorted(cache_roots):
            count=0
            for manifest_path in (root/'ligand-preparation').glob('*/manifest.json'):
                manifest=json.loads(manifest_path.read_text());key=manifest['key']
                if len(key)!=64 or any(c not in '0123456789abcdef' for c in key):continue
                if key in seen:continue
                files=manifest.get('files',{})
                if not files or set(files)-{'mol2','pdbqt'}:continue
                for fmt,sha in files.items():
                    source=manifest_path.parent/('ligand.'+fmt)
                    if file_digest(source)!=sha:raise ValueError('Corrupt preparation: '+str(source))
                    validate_prepared_ligand(source,fmt)
                entry=preparation/key;entry.mkdir(exist_ok=True)
                for fmt,sha in files.items():
                    output=entry/('ligand.'+fmt)
                    if not output.exists() or file_digest(output)!=sha:
                        copy_artifact(manifest_path.parent/('ligand.'+fmt),output)
                write_json(manifest,entry/'manifest.json');seen.add(key);count+=1
            calculations=root/'calculations'/'index.sqlite3';imported=0
            if calculations.is_file():
                with sqlite3.connect('file:'+str(calculations)+'?mode=ro',uri=True) as db:
                    experiments=db.execute('SELECT namespace,key,metadata FROM experiments').fetchall()
                    for namespace,key,metadata in experiments:
                        rows=db.execute('SELECT role,sha256 FROM artifacts WHERE namespace=? AND experiment_key=?',(namespace,key)).fetchall()
                        files={role:calculations.parent/'objects'/sha[:2]/sha for role,sha in rows}
                        if not files:continue
                        for role,sha in rows:
                            if file_digest(files[role])!=sha:raise ValueError('Corrupt calculation object')
                        cache=CalculationCache(namespace,target,root=target/'calculations')
                        if cache.metadata(key) is None:
                            cache.remember(key,files,json.loads(metadata));imported+=1
            report['cache_imports'].append(dict(source=str(root),preparations=count,calculations=imported))
    report['elapsed_seconds']=round(time.monotonic()-started,3)
    write_json(report,destination)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database',type=Path,required=True)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--consolidate-caches',action='store_true')
    args=parser.parse_args()
    print(json.dumps(audit(args.database,args.run_id,args.output,args.consolidate_caches),indent=2))


if __name__=='__main__':main()
