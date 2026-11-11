"""Bounded score transport for consensus; producer artifacts remain read-only."""
import csv
import json
import os
import tempfile
from pathlib import Path


def read_scores(directory, engine):
    from .docking_data import RESULT_FILE, read_results
    from .input_validation import ALIASES
    directory = Path(directory)
    summary = directory if directory.is_file() else directory/RESULT_FILE
    if not summary.is_file():
        yield from read_results(directory, engine)
        return
    resolved = {}
    with summary.open(encoding='utf-8-sig', newline='') as stream:
        for original in csv.DictReader(stream):
            row = {ALIASES.get(k, k): value for k, value in original.items()}
            if row.get('engine') and row['engine'] != engine:
                continue
            row['engine'] = engine
            for field in ('conformer_file', 'receptor_file', 'reference_file', 'footprint_file'):
                value = row.get(field)
                if value:
                    if value not in resolved:
                        path = Path(value)
                        resolved[value] = str(path if path.is_absolute() else summary.parent/path)
                    row[field] = resolved[value]
            yield row


def materialize_scores(rows, selected, destination):
    """Keep authorized absolute references, without cloning every input pose."""
    from .docking_data import write_csv
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    attachments = {}
    for row in rows:
        for field in ('conformer_file', 'pose_file', 'prepared_pdbqt', 'prepared_mol2',
                      'receptor_file', 'reference_file', 'footprint_file'):
            value = row.get(field)
            if value and value not in attachments:
                try:
                    stat = Path(value).stat()
                    attachments[value] = [stat.st_dev, stat.st_ino, stat.st_size,
                                          stat.st_mtime_ns, stat.st_ctime_ns]
                except FileNotFoundError:
                    attachments[value] = None
    fields = list(dict.fromkeys(key for row in rows for key in row)) or [
        'molecule_chembl_id', 'canonical_smiles', 'receptor_id', 'engine', 'score', 'conformer_file']
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination, delete=False) as stream:
            temporary = Path(stream.name)
        write_csv(temporary, rows, fields)
        temporary.replace(destination/'docking_results.csv')
    finally:
        if temporary:
            temporary.unlink(missing_ok=True)
    (destination/'consensus_input_files.json').write_text(json.dumps({
        'sources': [str(path) for path in selected], 'attachments': attachments},
        sort_keys=True, ensure_ascii=False), encoding='utf-8')
    return destination


def write_table(frame, path):
    """Publish complete CSV bytes atomically before optional visualization."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='',
                                         dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            frame.to_csv(stream, index=False)
        os.replace(temporary, path)
    finally:
        if temporary:
            temporary.unlink(missing_ok=True)
