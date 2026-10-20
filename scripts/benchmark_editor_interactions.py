"""Controlled Python editor benchmark (no Flutter rendering or real docking).

PYTHONPATH=src python scripts/benchmark_editor_interactions.py --baseline HEAD
"""
import argparse
import asyncio
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import time
import types
from types import SimpleNamespace
from unittest.mock import patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.workspace import WorkspaceStore

ROOT = Path(__file__).resolve().parents[1]


class Page:
    width=1440; height=1000; web=False
    def __init__(self):self.controls=[];self.dialogs=[]
    def update(self,*controls):pass
    def show_dialog(self,dialog):dialog.open=True;self.dialogs.append(dialog)
    def pop_dialog(self):return self.dialogs.pop() if self.dialogs else None


def previous(revision, name):
    module=types.ModuleType('biomolexplorer.ui.benchmark_'+name)
    module.__package__='biomolexplorer.ui'
    source=subprocess.check_output(['git','show',f'{revision}:src/biomolexplorer/ui/{name}.py'],cwd=ROOT,text=True)
    exec(compile(source,f'{revision}/{name}.py','exec'),module.__dict__)
    return module


def summary(values):
    ordered=sorted(values)
    return dict(median_ms=round(statistics.median(values),3),
                p95_ms=round(ordered[math.ceil(.95*len(values))-1],3))


async def measure(ui_type,editor_type,repeats):
    results=[]
    with tempfile.TemporaryDirectory() as tmp:
        store=WorkspaceStore(tmp)
        token=store.register('Bench','bench@example.org','benchmark-password')
        project=store.create_project(token,'Bench')
        ui=ui_type(Page(),store,SimpleNamespace(dock6_path=None));ui.token=token;ui.current=project
        ui.schedule_save=lambda:None
        # Match the existing UI tests: real SQLite calls, synchronous dispatch.
        # Executor scheduling and socket/Flutter latency are excluded.
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=call
        try:
            for count in (1,10,25,50,100):
                samples={key:[] for key in ('add','open','drag')};created=[]
                for repeat in range(repeats+1):
                    project['pipeline']=[new_stage('retrieve_structures') for _ in range(count-1)];ui.selected=None
                    editor=editor_type(ui,True);ui.flow_editor=editor
                    ui.page.controls=[editor.stack,editor.selection]
                    with patch.object(editor,'node',wraps=editor.node) as build:
                        start=time.perf_counter();await editor.add('docking_vina')
                        elapsed=(time.perf_counter()-start)*1000
                        if repeat:samples['add'].append(elapsed);created.append(build.call_count)
                    stage=project['pipeline'][-1]
                    start=time.perf_counter();await ui.open_stage_dialog(stage['id'])
                    elapsed=(time.perf_counter()-start)*1000
                    if repeat:samples['open'].append(elapsed)
                    ui.page.pop_dialog();ui.editing_stage=False
                    header=editor.nodes[stage['id']].content.controls[0]
                    start=time.perf_counter();header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1)))
                    if repeat:samples['drag'].append((time.perf_counter()-start)*1000)
                results.append(dict(blocks=count,nodes_created=statistics.median(created),
                                    **{key:summary(values) for key,values in samples.items()}))
            # Unconnected Vina in a project with unrelated historical compounds.
            folder=store.project_dir(project['id']);files=[]
            for i in range(5):
                path=folder/f'compounds_{i}.csv'
                path.write_text('molecule_chembl_id,canonical_smiles\n'+''.join(f'CHEMBL{j},CCO\n' for j in range(1000)))
                files.append(str(path))
            item=dict(new_stage('retrieve_compounds'),status='succeeded',artifacts=files)
            durations=[];reads=[]
            from biomolexplorer.pubchem_retrieval import reference_choices
            with patch.object(store,'list_runs',return_value=[dict(id='history',stages=[item])]):
                for _ in range(repeats):
                    with patch('biomolexplorer.pubchem_retrieval.reference_choices',wraps=reference_choices) as reader:
                        start=time.perf_counter();await ui.open_stage_dialog(stage['id'])
                        durations.append((time.perf_counter()-start)*1000);reads.append(reader.call_count)
                    ui.page.pop_dialog();ui.editing_stage=False
            return dict(interactions=results,unconnected_vina_history=dict(files=5,rows_per_file=1000,
                        reads=statistics.median(reads),**summary(durations)))
        finally:ui.structure_viewer.close()


async def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',default='HEAD');parser.add_argument('--repeats',type=int,default=21)
    args=parser.parse_args()
    if args.repeats<2:parser.error('use at least two repeats')
    before_app=previous(args.baseline,'app');before_canvas=previous(args.baseline,'flow_canvas')
    baseline=subprocess.check_output(['git','rev-parse',args.baseline],cwd=ROOT,text=True).strip()
    before=await measure(before_app.WorkspaceUI,before_canvas.FlowCanvas,args.repeats)
    after=await measure(WorkspaceUI,FlowCanvas,args.repeats)
    print(json.dumps(dict(baseline=baseline,repeats=args.repeats,before=before,after=after),indent=2))


if __name__=='__main__':asyncio.run(main())
