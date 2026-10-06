"""Structure-specific actions and audited ligand curation."""
import flet as ft
from .feedback import close_dialog
from .localization import verbatim


class PDBActions:
    def __init__(self,results,editable):
        self.results=results;self.ui=results.ui;self.editable=editable

    def actions(self,item):
        if not item['name'].lower().endswith('.pdb'):return []
        async def ligands(e):
            if self.results.valid():await self.ui.guard(lambda:self.open_ligands(item))
        async def structure(e):
            if self.results.valid():await self.ui.guard(lambda:self.ui.preview_artifact(self.results.project_id,item['path']))
        code=item['name'][:-4].upper()
        url=self.ui.pdb_view_url(self.results.project_id,item['path']) if hasattr(self.ui,'pdb_view_url') else None
        return [ft.IconButton(ft.Icons.SCIENCE_OUTLINED,tooltip='Ligantes · '+code,on_click=ligands),
                ft.IconButton(ft.Icons.VISIBILITY_OUTLINED,tooltip='Visualizar estrutura 3D',url=ft.Url(url,target=ft.UrlTarget.BLANK) if url else None,on_click=None if url else structure),
                ft.IconButton(ft.Icons.OPEN_IN_NEW,tooltip='Abrir no RCSB PDB',url='https://www.rcsb.org/structure/'+code)]

    async def open_ligands(self,item):
        result=self.results
        args=(result.token,result.project_id,result.run_id,result.stage['id'],item['path'])
        context=await self.ui.call(result.service.pdb_ligands,*args)
        if not result.valid():return
        previous_editing=getattr(self.ui,'editing_stage',False)
        def dismissed(e=None):self.ui.editing_stage=previous_editing
        def finish(e=None):
            dismissed();close_dialog(self.ui.page,dialog)
        rows=[];listing=ft.Column(spacing=12)
        def append(record=None):
            record=record or {}
            name=ft.TextField(label='Ligante',value=record.get('LIGAND',''),width=120,disabled=not self.editable)
            number=ft.TextField(label='Resíduo',value=str(record.get('RESNUM','')),width=130,disabled=not self.editable)
            chain=ft.TextField(label='Cadeia',value=record.get('CHAIN',''),width=90,disabled=not self.editable)
            row=ft.Row(wrap=True,spacing=12)
            entry=(row,name,number,chain);rows.append(entry)
            def remove(e):
                rows.remove(entry);listing.controls.remove(row);self.ui.page.update()
            row.controls=[name,number,chain,ft.IconButton(ft.Icons.DELETE_OUTLINE,tooltip='Remover ligante',
                            on_click=remove,disabled=not self.editable)]
            listing.controls.append(row)
        for record in context['ligands']:append(record)
        candidates=context.get('available_ligands',[])
        choice=ft.Dropdown(label='Ligantes presentes na estrutura',hint_text='Selecione um resíduo para adicionar',
            options=[verbatim(ft.DropdownOption(key=str(i),text=f"{r['LIGAND']} · {r['RESNUM']} · {r['CHAIN']}"),'text')
                     for i,r in enumerate(candidates)],enable_filter=True,enable_search=True,disabled=not self.editable)
        def add(e):
            append(candidates[int(choice.value)] if choice.value is not None else None)
            self.ui.page.update()
        async def save(e):
            async def action():
                if not result.valid():return
                records=[{'LIGAND':name.value,'RESNUM':number.value,'CHAIN':chain.value} for _,name,number,chain in rows]
                await self.ui.call(result.service.set_pdb_ligands,*args,records,context['revision'])
                if not result.valid():return
                finish()
                self.ui.notify('Ligantes atualizados. A alteração foi registrada no histórico do projeto.')
            await self.ui.guard(action)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Ligantes · '+context['pdb_id']),on_dismiss=dismissed,
            content=ft.Container(width=640,height=380,content=ft.Column([
                ft.Text('Mantenha os ligantes de interesse. Para adicionar, informe o código, resíduo e cadeia presentes na estrutura.',size=12),
                listing,choice,ft.TextButton('Adicionar ligante',icon=ft.Icons.ADD,on_click=add,disabled=not self.editable)],
                scroll=ft.ScrollMode.AUTO,spacing=20)),
            actions=[ft.TextButton('Fechar',on_click=finish),
                     ft.Button('Salvar ligantes',on_click=save,disabled=not self.editable)])
        self.ui.editing_stage=True
        self.ui.page.show_dialog(dialog)
