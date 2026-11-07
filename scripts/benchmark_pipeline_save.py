"""Reproducible configuration-save benchmark; no scientific work is executed."""
import argparse
import statistics
import tempfile
import time
from pathlib import Path

from biomolexplorer.catalog import new_stage
from biomolexplorer.workspace import WorkspaceStore


def measure(blocks,files,repeats):
    with tempfile.TemporaryDirectory(prefix='biomol-save-benchmark-') as folder:
        store=WorkspaceStore(folder);token=store.register('Benchmark','benchmark@example.org','benchmark-password')
        project=store.create_project(token,'Save benchmark');root=store.project_dir(project['id'])/'results';root.mkdir()
        for index in range(files):(root/f'{index:04d}.dat').write_bytes(b'x'*256)
        pipeline=[]
        for index in range(blocks):
            stage=new_stage('fingerprints');stage['name']=f'Fingerprint {index+1}';pipeline.append(stage)
        project=store.save_pipeline(token,project['id'],pipeline,project['revision'])
        unchanged=[];presentation=[];changed=[]
        for index in range(repeats):
            started=time.perf_counter();store.save_pipeline(token,project['id'],project['pipeline'],project['revision']);unchanged.append(time.perf_counter()-started)
            project=store.project(token,project['id'])
            project['pipeline'][index%blocks]['position']={'x':index+1,'y':index+2}
            started=time.perf_counter();project=store.save_pipeline(token,project['id'],project['pipeline'],project['revision'],presentation_only=True);presentation.append(time.perf_counter()-started)
            project['pipeline'][index%blocks]['parameters']['radius']=2+(index%3)
            started=time.perf_counter();project=store.save_pipeline(token,project['id'],project['pipeline'],project['revision']);changed.append(time.perf_counter()-started)
        return {'blocks':blocks,'files':files,'unchanged_median_ms':statistics.median(unchanged)*1000,
            'presentation_median_ms':statistics.median(presentation)*1000,
            'changed_median_ms':statistics.median(changed)*1000}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=3);parser.add_argument('--files',type=int,default=300)
    args=parser.parse_args()
    for blocks in (1,20,80):print(measure(blocks,args.files,args.repeats))
