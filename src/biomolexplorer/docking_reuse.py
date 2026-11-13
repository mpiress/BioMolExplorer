"""Accept complete docking experiments with explicitly accounted exclusions."""
import csv
import json
import math
from collections import Counter
from pathlib import Path


def accounted_preparation_exclusions(entry, operation, root):
    """Never adopt interrupted jobs, unknown pairs or failed engine executions.

    Chemical preparation exclusions are recorded outcomes of unchanged inputs.
    Retrying them requires an explicit forced execution or changed inputs.
    All accompanying files must still match the producer's original hashes.
    """
    if operation not in ('docking_vina', 'docking_dock6'):
        return False
    if entry.get('status') != 'succeeded' or any(entry.get(key) for key in (
            'outputs_removed', 'dependency_obsolete', 'partial_results')):
        return False
    from .result_state import input_digest
    engine = 'vina' if operation == 'docking_vina' else 'dock6'
    root = Path(root).resolve()
    excluded = 0
    units = entry.get('batches') or [entry]
    try:
        for unit in units:
            if unit.get('status') != 'succeeded' or any(unit.get(key) for key in (
                    'outputs_removed', 'dependency_obsolete', 'partial_results')):
                return False
            manifest = unit.get('artifact_manifest') or entry.get('artifact_manifest', {})
            files = {Path(path).resolve() for path in unit.get('artifacts', [])}
            reports = sorted(path for path in files if path.name == 'docking_accounting.json')
            if not reports:
                return False
            for path in reports:
                missing = path.with_name('docking_unprocessed.csv')
                scores = path.with_name('docking_results.csv')
                for source in (path, missing, scores):
                    if source not in files or not source.is_relative_to(root):
                        return False
                    if input_digest(source) != manifest.get(str(source)):
                        return False
                report = json.loads(path.read_text())
                fields = ('expected_pairs', 'successful_pairs', 'unsuccessful_pairs',
                          'unaccounted_pairs', 'unique_compounds', 'receptors')
                if any(type(report.get(key)) is not int or report[key] < 0 for key in fields):
                    return False
                if (report.get('schema_version') != 1 or report.get('engine') != engine or
                        report['unaccounted_pairs'] or report['successful_pairs'] == 0 or
                        report['expected_pairs'] != report['unique_compounds']*report['receptors'] or
                        report['expected_pairs'] != report['successful_pairs']+report['unsuccessful_pairs']):
                    return False
                unsuccessful = set()
                with missing.open(encoding='utf-8-sig', newline='') as stream:
                    for row in csv.DictReader(stream):
                        pair = row['receptor_id'], row['molecule_chembl_id']
                        if (not all(pair) or pair in unsuccessful or
                                row['status'] != 'preparation_failed' or not row['reason'].strip()):
                            return False
                        unsuccessful.add(pair)
                successful = set()
                with scores.open(encoding='utf-8-sig', newline='') as stream:
                    for row in csv.DictReader(stream):
                        pair = row['receptor_id'], row['molecule_chembl_id']
                        if (not all(pair) or pair in successful or pair in unsuccessful or
                                row.get('engine') != engine or not math.isfinite(float(row['score']))):
                            return False
                        successful.add(pair)
                if (len(unsuccessful) != report['unsuccessful_pairs'] or
                        len(successful) != report['successful_pairs'] or
                        len({code for receptor, code in successful | unsuccessful}) != report['unique_compounds']):
                    return False
                good = Counter(receptor for receptor, code in successful)
                bad = Counter(receptor for receptor, code in unsuccessful)
                per_receptor = report['per_receptor']
                if set(per_receptor) != set(good) | set(bad) or len(per_receptor) != report['receptors']:
                    return False
                for receptor, counts in per_receptor.items():
                    if counts != dict(expected=report['unique_compounds'],
                                      successful=good[receptor], unsuccessful=bad[receptor]):
                        return False
                excluded += len(unsuccessful)
        return excluded > 0
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return False
