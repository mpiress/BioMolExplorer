"""Exercise real, verified cache restores without starting scientific engines."""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
import tempfile
import time

from biomolexplorer.calculation_cache import CalculationCache
from biomolexplorer.ligand_cache import LigandCache
from biomolexplorer.stage_cache import file_digest
from biomolexplorer.storage import write_json
from biomolexplorer.vina_batch import validate_pose


class ReadOnlyCalculations(CalculationCache):
    def __init__(self,root,namespace):
        self.root=Path(root);self.database=self.root/'index.sqlite3'
        self.namespace=namespace;self.reuse=True

    @contextmanager
    def connect(self):
        db=sqlite3.connect('file:'+str(self.database)+'?mode=ro',uri=True)
        try:yield db
        finally:db.close()


class ReadOnlyPreparations(LigandCache):
    def __init__(self,root):
        self.root=Path(root);self.reuse=True


def validate(root,output,limit=36):
    root=Path(root);started=time.monotonic();prepared=poses=0
    calculation=ReadOnlyCalculations(root/'calculations','vina')
    preparation=ReadOnlyPreparations(root/'ligand-preparation')
    with tempfile.TemporaryDirectory(prefix='biomol-warm-cache-') as directory:
        destination=Path(directory)
        for manifest in sorted(preparation.root.glob('*/manifest.json'))[:limit]:
            data=json.loads(manifest.read_text());formats=tuple(data['files'])
            row={'molecule_chembl_id':'validation'}
            if not preparation.restore(data['key'],row,destination,formats):
                raise ValueError('Prepared ligand was not restored')
            for fmt,sha in data['files'].items():
                if file_digest(row['prepared_'+fmt])!=sha:raise ValueError('Prepared coordinates changed')
            prepared+=1
        with calculation.connect() as db:
            count=db.execute('SELECT count(*) FROM experiments WHERE namespace=?',('vina',)).fetchone()[0]
            keys=[row[0] for row in db.execute('SELECT key FROM experiments WHERE namespace=? ORDER BY key LIMIT ?',('vina',limit))]
        for key in keys:
            pose=destination/'pose.pdbqt'
            if calculation.restore(key,{'pose':pose},validate_pose) is None:
                raise ValueError('Completed Vina pose was not restored')
            poses+=1
    report=dict(preparations_restored=prepared,poses_restored=poses,
        indexed_vina_calculations=count,scientific_engine_invocations=0,
        unchanged_coordinates_and_charges=True,elapsed_seconds=round(time.monotonic()-started,4))
    write_json(report,output);return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(validate(args.cache_root,args.output),indent=2))
