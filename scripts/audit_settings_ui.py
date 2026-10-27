"""Render production configuration forms for browser regression inspection.

Run with PYTHONPATH=src python scripts/audit_settings_ui.py, then visit
http://127.0.0.1:8573/#/docking_vina (append /disabled for read-only controls).
No project files or scientific calculations are created.
"""
import argparse
from types import SimpleNamespace
import flet as ft
from biomolexplorer.catalog import TITLES,new_stage
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.settings_layout import dialog_size


def main(page):
    page.title='BioMolExplorer · configuration audit'
    page.theme=ft.Theme(color_scheme_seed='#1A5E9A')
    page.padding=0
    stages=[new_stage(op) for op in TITLES]
    ui=SimpleNamespace(current={'pipeline':stages},service=SimpleNamespace(dock6_path=None),
        page=page,artifact_choices={},redocking_records={})
    dialog=None
    body=None
    def route(e=None):
        nonlocal dialog,body
        path=page.route.strip('/').split('/')
        operation=path[0] if path[0] in TITLES else 'docking_vina'
        stage=next(s for s in stages if s['operation']==operation)
        if operation=='redocking':stage['parameters']['pdb_codes']=[['1LGU','AAH',163,'A',1.8]]
        form=GuidedForm(ui,stage,[],len(path)<2 or path[1]!='disabled')
        if dialog:page.pop_dialog()
        body=ft.Container(**dialog_size(page),padding=16,content=form.layout())
        dialog=ft.AlertDialog(title=ft.Text(TITLES[operation][0],size=22),content=body,
            actions=[ft.TextButton('Validar configurações',on_click=lambda e:page.show_dialog(
                ft.AlertDialog(title=ft.Text('Configuração lida'),content=ft.Text(str(form.read()['parameters'])))))])
        page.show_dialog(dialog)
    def resize(e):
        if body:
            for key,value in dialog_size(page).items():setattr(body,key,value)
            page.update()
    page.on_route_change=route
    page.on_resize=resize
    page.add(ft.Text('Configuration audit · use #/<operation> in the URL'))
    route()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8573)
    args=parser.parse_args()
    import flet_web.fastapi as flet_fastapi
    import uvicorn
    app=flet_fastapi.FastAPI()
    app.mount('/',flet_fastapi.app(main,no_cdn=True,route_url_strategy=ft.RouteUrlStrategy.HASH))
    uvicorn.run(app,host='127.0.0.1',port=args.port,access_log=False,ws='websockets-sansio')
