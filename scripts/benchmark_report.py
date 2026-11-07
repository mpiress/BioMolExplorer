"""Synthetic report benchmark; never opens or changes user projects."""
import json
import tempfile
import time
from statistics import median
from uuid import uuid4

from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.catalog import new_stage
from biomolexplorer.experimental_report import build_report,build_report_preview


def main():
    with tempfile.TemporaryDirectory() as directory:
        store=WorkspaceStore(directory)
        token=store.register('Audit','report@example.org','audit-password')
        project=store.create_project(token,'Report benchmark');pid=project['id']
        stages=[new_stage('admet') for _ in range(4)]
        store.save_pipeline(token,pid,stages,project['revision'])
        for index,stage in enumerate(stages):
            paths=[]
            for table in range(8):
                path=store.project_dir(pid)/f'metrics-{index}-{table}.csv'
                path.write_text('compound,ADMET_selected,MW\n'+''.join(f'C{i},true,{100+i%300}\n' for i in range(1000)))
                paths.append(str(path))
            stage.update(status='succeeded',artifacts=paths)
        with store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(uuid4().hex,pid,store.user(token)['id'],'succeeded',json.dumps(stages),time.time(),time.time(),None))
        start=time.perf_counter();report=build_report(store,token,pid);full=time.perf_counter()-start
        previews=[]
        for _ in range(10):
            start=time.perf_counter();build_report_preview(store,token,pid);previews.append(time.perf_counter()-start)
        print(json.dumps({'scenario':{'stages':4,'tables':32,'rows_per_table':1000},'full_seconds':full,
            'preview_median_seconds':median(previews),'preview_max_seconds':max(previews),
            'charts':sum(len(analysis['charts']) for analysis in report['analyses']),
            'datasets_with_metrics':len({metric['dataset'] for analysis in report['analyses'] for metric in analysis['metrics']})},indent=2))


if __name__=='__main__':main()
