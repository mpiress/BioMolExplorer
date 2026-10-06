"""Differential check against an independently built official DMS executable."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import tempfile

from biomolexplorer.molecular_surface import generate_surface


def canonical(path):
    return Counter(row.replace('-0.000', ' 0.000') for row in path.read_text().splitlines())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dms-executable', required=True)
    args = parser.parse_args()
    executable = str(Path(args.dms_executable).resolve())
    fixtures = Path(__file__).resolve().parents[1]/'tests/fixtures/dms'
    failures = []
    with tempfile.TemporaryDirectory(prefix='biomol-dms-validation-') as folder:
        reference, native = Path(folder)/'reference.dms', Path(folder)/'native.dms'
        for case in json.loads((fixtures/'manifest.json').read_text()):
            source = fixtures/f"{case['name']}.pdb"
            subprocess.run([executable, str(source), '-d', str(case['density']), '-w',
                            str(case['probe_radius']), '-n', '-o', str(reference)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            summary = generate_surface(source, native, density=case['density'], probe_radius=case['probe_radius'])
            old, new = canonical(reference), canonical(native)
            missing, extra = sum((old-new).values()), sum((new-old).values())
            print(f"{case['name']}: missing={missing}, extra={extra}; {summary}")
            if missing or extra:
                failures.append(case['name'])
    if failures:
        raise SystemExit('Differential check failed: '+', '.join(failures))


if __name__ == '__main__':
    main()
