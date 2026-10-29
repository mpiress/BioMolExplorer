"""Durable, compressed per-molecule fingerprint memoization without pickle."""
import hashlib
import inspect
import json
import os
from pathlib import Path
import sqlite3
import zlib
import ast
import textwrap


def fingerprints(frame,compute,kind,radius,n_bits):
    base=os.environ.get('BIOMOL_CACHE_DIR')
    if not base or os.environ.get('BIOMOL_REUSE_RESULTS')=='0':return compute(frame)
    from rdkit import Chem,rdBase
    from kernel.descriptors import Descriptors
    from .molecular_parsing import mol_from_smiles
    from .stage_cache import file_digest
    from . import molecular_parsing
    directory=Path(base)/'calculations/molecular';directory.mkdir(parents=True,exist_ok=True)
    worker=getattr(Descriptors,{'morgan':'morgan_worker','maccs':'maccs_worker','pharmacophore':'pharmacophore_worker'}[kind])
    protocol=json.dumps(dict(schema=1,kind=kind,radius=radius if kind=='morgan' else None,
        bits=n_bits if kind=='morgan' else None,rdkit=rdBase.rdkitVersion,
        implementation=hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(worker))),include_attributes=False).encode()).hexdigest(),
        parsing=file_digest(molecular_parsing.__file__)),sort_keys=True)
    keys=[];labels={}
    for smile in frame['smiles']:
        molecule=mol_from_smiles(smile)
        canonical=Chem.MolToSmiles(molecule,isomericSmiles=True) if molecule is not None else str(smile)
        key=hashlib.sha256((protocol+'\n'+canonical).encode()).hexdigest();keys.append(key);labels[key]=canonical
    with sqlite3.connect(directory/'fingerprints.sqlite3',timeout=60) as db:
        db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=FULL')
        db.execute('CREATE TABLE IF NOT EXISTS entries (key TEXT PRIMARY KEY,label TEXT,protocol TEXT,payload BLOB,sha TEXT)')
        values={};unique=list(dict.fromkeys(keys))
        for start in range(0,len(unique),400):
            chunk=unique[start:start+400]
            for key,stored_protocol,payload,sha in db.execute('SELECT key,protocol,payload,sha FROM entries WHERE key IN ('+','.join('?' for _ in chunk)+')',chunk):
                try:
                    if stored_protocol!=protocol:continue
                    if hashlib.sha256(payload).hexdigest()!=sha:continue
                    vector=json.loads(zlib.decompress(payload))
                    if not isinstance(vector,list) or not vector or any(type(v) is not int or v not in (0,1) for v in vector):continue
                    values[key]=vector
                except (ValueError,zlib.error):continue
        missing={key:index for index,key in reversed(list(enumerate(keys))) if key not in values}
        if missing:
            selected=list(missing);computed=compute(frame.iloc[[missing[key] for key in selected]])
            for key,vector in zip(selected,computed['fingerprint']):
                if vector is None:continue
                vector=vector if isinstance(vector,list) else list(vector)
                payload=zlib.compress(json.dumps(vector,separators=(',',':'),default=int).encode())
                db.execute('INSERT OR REPLACE INTO entries VALUES (?,?,?,?,?)',(key,labels[key],protocol,payload,hashlib.sha256(payload).hexdigest()))
                values[key]=vector
    import pandas as pd
    return pd.DataFrame({'fingerprint':[values.get(key) for key in keys]})
