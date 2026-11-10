"""Serve the production canvas with synthetic data, without a backend pipeline."""
import argparse
import inspect
from types import SimpleNamespace
import time

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage,verbatim

STATE={}


def main(raw_page):
    page=LocalizedPage(raw_page,'en')
    page.title='BioMolExplorer canvas performance audit';page.padding=8
    page.theme=ft.Theme(color_scheme_seed='#2552E8')
    ui=SimpleNamespace(page=page,current={'id':'canvas-audit','pipeline':[]},selected=None,
        dirty=False,tab='Pipeline',service=SimpleNamespace(dock6_path=None),sidebar_expanded={})
    ui.tr=page.translator
    info=verbatim(ft.Text('Canvas audit'),'value')
    async def guard(action):
        try:
            result=action()
            if inspect.isawaitable(result):await result
        except Exception as error:
            STATE.setdefault('errors',[]).append(str(error));raise
    ui.guard=guard
    ui.notify=lambda message:STATE.setdefault('notifications',[]).append(str(message))
    ui.selected_stage=lambda:next((stage for stage in ui.current['pipeline'] if stage['id']==ui.selected),None)
    def save():STATE['autosaves']=STATE.get('autosaves',0)+1
    ui.schedule_save=save
    async def configure(stage_id):
        stage=next(stage for stage in ui.current['pipeline'] if stage['id']==stage_id)
        dialog=ft.AlertDialog(title=ft.Text('Configure '+stage['name']),content=ft.Text(stage['operation']),
            actions=[ft.TextButton('Close',on_click=lambda e:page.pop_dialog())])
        page.show_dialog(dialog)
    ui.open_stage_dialog=configure;ui.open_stage_results=configure

    def route(e=None):
        size=raw_page.route.strip('/')
        count=int(size) if size in ('10','50','100') else 100
        stages=[new_stage('retrieve_compounds')]+[new_stage('admet') for _ in range(count-1)]
        for index,stage in enumerate(stages):
            stage['name']=f'Block {index+1:03d}'
            if index:
                stage['bindings']['base_input_path']={'sources':[
                    dict(stage=stages[index-1]['id'],selector=f'selected-{j}.csv') for j in range(25)]}
        ui.current['pipeline']=stages;ui.selected=stages[0]['id']
        STATE.clear();STATE.update(errors=[],drag_ms=[],autosaves=0)
        editor=FlowCanvas(ui,True);ui.flow_editor=editor
        for node in editor.nodes.values():
            header=node.content.controls[0];move=header.on_pan_update
            def measured(event,move=move):
                start=time.perf_counter();move(event)
                STATE['drag_ms'].append(round((time.perf_counter()-start)*1000,4))
            header.on_pan_update=measured
        info.value=f'Canvas ready: {count} blocks, {len(editor.edges.shapes)} connections'
        raw_page.on_keyboard_event=editor.keyboard
        raw_page.controls=[ft.Column([info,editor.build(False)],expand=True)]
        STATE['ui']=ui
        page.update()
    raw_page.on_route_change=route
    route()


def metrics():
    ui=STATE.get('ui')
    data={key:value for key,value in STATE.items() if key!='ui'}
    if ui:
        editor=ui.flow_editor
        data.update(blocks=len(editor.stages),connections=len(editor.edges.shapes),
            positions={stage['name']:stage.get('position') for stage in editor.stages},
            viewport_scale=editor.viewport_scale,
            viewport_offset=[editor.viewport_offset.x,editor.viewport_offset.y],
            selected=ui.selected,history=len(editor.history))
    return data


if __name__=='__main__':
    import flet_web.fastapi as flet_fastapi
    import uvicorn
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=8578);args=parser.parse_args()
    app=flet_fastapi.FastAPI()
    app.get('/canvas-audit')(metrics)
    app.mount('/',flet_fastapi.app(main,no_cdn=True,route_url_strategy=ft.RouteUrlStrategy.HASH))
    uvicorn.run(app,host='127.0.0.1',port=args.port,access_log=False,ws='websockets-sansio')
