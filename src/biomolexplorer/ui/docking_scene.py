"""Shared docking overlay and geometric contact actions for result tables."""
import csv
import io
import flet as ft
from biomolexplorer.docking_scene import scene_payload
from .feedback import close_dialog
from .localization import verbatim


async def residue_dialog(ui,pid,spec_loader,valid,open_3d):
    if not valid():return
    cutoff=ft.TextField(label='Distância de contato (Å)',value='4',width=190)
    # A scrollable Column still grows to its content unless its parent gives it
    # the remaining height. Keep both scroll directions inside a clipped viewport.
    listing=ft.Column(expand=True,scroll=ft.ScrollMode.AUTO,spacing=8)
    state={}
    async def load():
        spec=await ui.call(spec_loader)
        data=await ui.call(scene_payload,ui.store,ui.token,pid,spec,float(cutoff.value))
        if not valid():return
        state['data']=data
        rows={}
        for role,key in (('pose','contacts'),('reference','reference_contacts')):
            for item in data[key]:
                identity=tuple(item[k] for k in ('chain','resi','icode','resn'))
                entry=rows.setdefault(identity,dict(item,pose=None,reference=None))
                entry[role]=item['distance']
                entry[role+'_atoms']=item['receptor_atom']+' / '+item['ligand_atom']
        state['rows']=list(rows.values())
        table=ft.DataTable(column_spacing=20,data_row_min_height=40,data_row_max_height=48,
            columns=[ft.DataColumn(ft.Text(label)) for label in
            ('Resíduo','Cadeia','Pose (Å)','Referência (Å)','Átomos na pose','Átomos na referência')],rows=[])
        for row in sorted(state['rows'],key=lambda r:r['pose'] if r['pose'] is not None else r['reference']):
            values=[f"{row['resn']} {row['resi']}{row['icode']}",row['chain'] or '—',
                f"{row['pose']:.3f}" if row['pose'] is not None else '—',
                f"{row['reference']:.3f}" if row['reference'] is not None else '—',row.get('pose_atoms','—'),row.get('reference_atoms','—')]
            table.rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(value))) for value in values]))
        listing.controls=[ft.Row([table],scroll=ft.ScrollMode.AUTO)] if state['rows'] else [ft.Text('Nenhum contato nesta distância.')]
        if data['reference_kind']:
            listing.controls.insert(0,ft.Text('Referência: ligante cristalográfico.' if data['reference_kind']=='crystal' else 'Referência: ligante preparado.'))
        ui.page.update()
    async def apply(e):await ui.guard(load)
    async def export():
        if not valid() or not state.get('rows'):return
        await ui.call(spec_loader)
        if not valid():return
        stream=io.StringIO();fields=['resn','resi','icode','chain','pose','reference','pose_atoms','reference_atoms']
        writer=csv.DictWriter(stream,fields,extrasaction='ignore');writer.writeheader();writer.writerows(state['rows'])
        await ui.picker.save_file(file_name='residue_contacts.csv',src_bytes=stream.getvalue().encode())
    async def download(e):await ui.guard(export)
    async def view(e):
        if valid():await ui.guard(open_3d)
    note='Contatos geométricos por distância mínima entre átomos pesados. Não classifica ligações de hidrogênio ou outras interações.'
    content=ft.Container(clip_behavior=ft.ClipBehavior.HARD_EDGE,
        content=ft.Column([ft.Text(note,max_lines=3,overflow=ft.TextOverflow.ELLIPSIS,tooltip=note),
            ft.Row([cutoff,ft.TextButton('Aplicar',on_click=apply)],wrap=True),listing],expand=True,spacing=12))
    previous_resize=getattr(ui.page,'on_resize',None)
    active=True
    def resize():
        width=getattr(ui.page,'width',None) or 1000
        height=getattr(ui.page,'height',None) or 800
        content.width=max(180,min(980,width-80))
        content.height=max(100,min(620,height-(300 if width<620 else 220)))
    async def resized(e):
        if previous_resize:
            import inspect
            result=previous_resize(e)
            if inspect.isawaitable(result):await result
        if active:
            resize();ui.page.update()
    def dismissed(e=None):
        nonlocal active
        active=False
        if getattr(ui.page,'on_resize',None) is resized:ui.page.on_resize=previous_resize
    def close(e):
        dismissed();close_dialog(ui.page,dialog)
    dialog=ft.AlertDialog(title=ft.Text('Resíduos próximos ao ligante'),content=content,
        inset_padding=16,content_padding=ft.Padding.symmetric(horizontal=24,vertical=12),
        clip_behavior=ft.ClipBehavior.HARD_EDGE,on_dismiss=dismissed,
        actions=[ft.TextButton('Abrir comparação 3D',on_click=view),ft.TextButton('Baixar contatos (CSV)',on_click=download),
            ft.TextButton('Fechar',on_click=close)])
    resize()
    await load()
    if valid():
        ui.page.on_resize=resized
        ui.page.show_dialog(dialog)
