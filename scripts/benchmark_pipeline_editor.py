"""Compare Python drag handling against a Git revision; no browser or network.

PYTHONPATH=src python scripts/benchmark_pipeline_editor.py --baseline <revision>
"""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import time
import types
from types import SimpleNamespace

from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage
import biomolexplorer.ui.localization as localization

ROOT=Path(__file__).resolve().parents[1]


class Page:
    def __init__(self):self.controls=[]
    def update(self,*controls):pass


def baseline_module(revision,name):
    module=types.ModuleType('biomolexplorer.ui.benchmark_'+name)
    module.__package__='biomolexplorer.ui'
    module.__file__=str(Path(localization.__file__).with_name(name+'.py'))
    source=subprocess.check_output(
        ['git','show',f'{revision}:src/biomolexplorer/ui/{name}.py'],cwd=ROOT,text=True)
    exec(compile(source,module.__file__,'exec'),module.__dict__)
    return module


def measure(editor_type,page_type,count,frames):
    raw_page=Page()
    ui=SimpleNamespace(current={'pipeline':[new_stage('retrieve_structures') for _ in range(count)]},
        selected=None,dirty=False,page=page_type(raw_page),notify=lambda message:None)
    editor=editor_type(ui,True)
    raw_page.controls=[editor.stack,editor.selection]
    ui.page.update()  # Warm localization caches; initial construction is excluded.
    header=editor.stack.controls[2].content.controls[0]
    event=SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1))
    durations=[]
    for _ in range(frames):
        start=time.perf_counter()
        header.on_pan_update(event)
        durations.append(time.perf_counter()-start)
    return round(statistics.median(durations)*1000,3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--frames',type=int,default=50)
    args=parser.parse_args()
    if args.frames<1:parser.error('--frames must be positive')
    previous=baseline_module(args.baseline,'flow_canvas')
    previous_localization=baseline_module(args.baseline,'localization')
    results=[]
    for count in (1,20,80):
        results.append(dict(blocks=count,
            before_ms=measure(previous.FlowCanvas,previous_localization.LocalizedPage,count,args.frames),
            after_ms=measure(FlowCanvas,LocalizedPage,count,args.frames)))
    print(json.dumps(dict(baseline=args.baseline,frames=args.frames,results=results),indent=2))


if __name__=='__main__':main()
