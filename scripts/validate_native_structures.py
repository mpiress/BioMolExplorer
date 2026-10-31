"""Smoke-test native DOCK6 6.11 Amber Score on its bundled 1LGU example.

All work is isolated in a temporary directory; no installation/project edits.
This replaces the retired experimental external ff14SB/sander validator.
"""
import argparse
import json
from pathlib import Path
import tempfile
import time

from biomolexplorer.md_config import validate_settings
from biomolexplorer.md_engine import tools,run_pair
from biomolexplorer.storage import write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dock6',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--steps',type=int,default=2)
    args=parser.parse_args()
    settings=validate_settings(dict(amber_score_md_steps=args.steps,
        amber_score_before_md_minimization_cycles=2,amber_score_after_md_minimization_cycles=2))
    root,found=tools(args.dock6,'',settings)
    example=root/'install/test/amber_score_1lgu'
    pair=dict(receptor=dict(path=str(example/'1lgu.pdb'),name='1LGU'),
        ligand=dict(path=str(example/'lig.mol2'),name='1LGU ligand'))
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='biomol-dock6-amber-') as folder:
        result=run_pair(pair,Path(folder),root,found,settings)
        data=dict(energies=result,analysis=json.loads((Path(folder)/'analysis.json').read_text()),
            protocol=json.loads((Path(folder)/'protocol.json').read_text()),
            elapsed_seconds=time.monotonic()-started,
            limitation='One short native example; not full scientific acceptance or conventional MD.')
        write_json(data,Path(args.output))
    print('Native DOCK6 6.11 Amber Score completed; report: '+args.output)


if __name__=='__main__':main()
