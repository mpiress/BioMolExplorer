"""Compare CSV rescans with indexed pages; optionally audit a real run read-only."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import time
import zlib

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biomolexplorer.result_index import ResultIndex
from biomolexplorer.result_catalog import _project_blob,proof_stage
from biomolexplorer.table_sort import ordered


def digest(stream):
    value=hashlib.sha256()
    for chunk in iter(lambda:stream.read(1048576),b''):value.update(chunk)
    stream.seek(0)
    return value.hexdigest()


def benchmark(rows):
    with tempfile.TemporaryDirectory(prefix='biomol-index-benchmark-') as directory:
        root=Path(directory);source=root/'results.csv'
        with source.open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['molecule_chembl_id','canonical_smiles','score','receptor_id','description'])
            for index in range(rows):writer.writerow([f'M{index:08}', 'CCO',-float(index%997),'TARGET_A','measurement '+('x'*100)])
        baseline=[];before=time.monotonic()
        for offset in (0,25,50):
            with source.open('rb') as stream:version=digest(stream)
            with source.open(newline='') as stream:values=list(enumerate(csv.DictReader(stream)))
            values=ordered(values,lambda item:item[1]['score'],True,True)
            baseline.append([index for index,_ in values[offset:offset+25]])
        baseline_seconds=time.monotonic()-before
        index=ResultIndex(root,source,digest);before=time.monotonic()
        first=index.page(0,25,'','score',True,True)
        cold=time.monotonic()-before
        timings=[];actual=[]
        for offset in (0,25,50):
            before=time.monotonic();page=index.page(offset,25,'','score',True,True)
            timings.append(time.monotonic()-before);actual.append([ordinal for ordinal,_ in page['visible']])
        if actual!=baseline or first['version']!=version:raise AssertionError('Indexed result differs from CSV')
        with source.open('rb') as stream:
            if digest(stream)!=version:raise AssertionError('Source changed')
        return dict(rows=rows,csv_bytes=source.stat().st_size,legacy_equivalent_three_pages_seconds=baseline_seconds,
            indexed_cold_build_and_sort_seconds=cold,indexed_warm_page_seconds=timings,
            indexed_database_bytes=index.database.stat().st_size,identical_rows_and_version=True)


def audit_run(database,run_id):
    with sqlite3.connect('file:'+str(database.resolve())+'?mode=ro',uri=True) as db:
        row=db.execute('SELECT rowid,length(CAST(stages AS BLOB)) FROM runs WHERE id=?',(run_id,)).fetchone()
        if row is None:raise ValueError('Run not found')
        before=time.monotonic()
        with db.blobopen('runs','stages',row[0],readonly=True) as blob:stages=_project_blob(blob)
        duration=time.monotonic()-before
        encoded=json.dumps(stages,separators=(',',':')).encode()
        compressed=zlib.compress(encoded,1)
        before=time.monotonic();decoded=json.loads(zlib.decompress(compressed));warm=time.monotonic()-before
        if len(decoded)!=len(stages):raise AssertionError('Catalog differs')
        result=dict(run_id=run_id,original_bytes=row[1],review_bytes=len(encoded),projection_seconds=duration,
            compressed_review_bytes=len(compressed),cached_catalog_decode_seconds=warm,
            proof_bytes=len(json.dumps([proof_stage(stage) for stage in stages],separators=(',',':')).encode()),
            artifacts=sum(len(stage.get('artifacts',[])) for stage in stages),
            stages=len(stages),original_database_unchanged=True)
        from biomolexplorer.artifact_naming import native_name
        tables=[Path(value) for stage in stages if stage.get('operation')=='docking_vina'
                for value in stage.get('artifacts',[]) if str(value).endswith('.csv') and native_name(value)=='docking_results.csv']
        if tables:
            with tempfile.TemporaryDirectory(prefix='biomol-real-result-index-') as directory:
                index=ResultIndex(Path(directory),tables[0],digest)
                before=time.monotonic();page=index.page(0,25,'','score',True,True);cold=time.monotonic()-before
                before=time.monotonic();warm_page=index.page(25,25,'','score',True,True);warm=time.monotonic()-before
                with tables[0].open(encoding='utf-8-sig',newline='') as stream:original=list(enumerate(csv.DictReader(stream)))
                expected=ordered(original,lambda item:item[1]['score'],True,True)
                if page['visible']!=expected[:25] or warm_page['visible']!=expected[25:50]:raise AssertionError('Real scores differ')
                result['vina_table']=dict(rows=page['total'],cold_seconds=cold,warm_page_seconds=warm,identical_sorted_rows=True)
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,default=100000)
    parser.add_argument('--database',type=Path)
    parser.add_argument('--run-id')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report={'tables':benchmark(args.rows)}
    if args.database and args.run_id:report['real_run']=audit_run(args.database,args.run_id)
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
