"""Compare connected canvas CPU work; transport and Flutter rendering excluded."""
import argparse
import json
import math
from pathlib import Path
import statistics
import time
from types import SimpleNamespace

from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage
from benchmark_pipeline_editor import baseline_module


class Page:
    def __init__(self):self.controls=[];self.overlay=[]
    def update(self,*controls):self.updated=len(controls)


def measure(canvas_type,page_type,count,files,frames):
    stages=[new_stage('retrieve_compounds')]+[new_stage('admet') for _ in range(count-1)]
    for index,stage in enumerate(stages[1:],1):
        refs=[dict(stage=stages[index-1]['id'],selector=f'file-{j}.csv') for j in range(files)]
        stage['bindings']['base_input_path']={'sources':refs}
    raw=Page();ui=SimpleNamespace(current={'pipeline':stages},selected=None,dirty=False,
        page=page_type(raw),notify=lambda value:None)
    start=time.perf_counter();editor=canvas_type(ui,True);initial=(time.perf_counter()-start)*1000
    raw.controls=[editor.stack,editor.selection];ui.page.update()
    header=editor.nodes[stages[count//2]['id']].content.controls[0]
    times=[]
    for _ in range(frames):
        start=time.perf_counter()
        header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1)))
        times.append((time.perf_counter()-start)*1000)
    return dict(blocks=count,files_per_connection=files,visual_edges=len(editor.edges.shapes),
        initial_ms=round(initial,3),drag_median_ms=round(statistics.median(times),3),
        drag_p95_ms=round(sorted(times)[math.ceil(len(times)*.95)-1],3))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',default='HEAD');parser.add_argument('--frames',type=int,default=40)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.frames<1:parser.error('--frames must be positive')
    old_canvas=baseline_module(args.baseline,'flow_canvas')
    old_page=baseline_module(args.baseline,'localization')
    report=dict(scope='Python event handling with production localization; no Flutter FPS measurement',
        baseline=args.baseline,frames=args.frames,before=[],after=[])
    for count,files in ((10,1),(50,1),(100,1),(100,25)):
        report['before'].append(measure(old_canvas.FlowCanvas,old_page.LocalizedPage,count,files,args.frames))
        report['after'].append(measure(FlowCanvas,LocalizedPage,count,files,args.frames))
    output=json.dumps(report,indent=2)
    if args.output:args.output.write_text(output+'\n')
    print(output)
