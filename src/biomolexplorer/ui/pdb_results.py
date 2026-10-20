"""Structure-specific actions and audited ligand curation."""
import copy
import flet as ft
from biomolexplorer.pdb_ligands import REASONS
from .feedback import close_dialog
from .localization import verbatim


class PDBActions:
    def __init__(self,results,editable):
        self.results=results;self.ui=results.ui;self.editable=editable

    async def refresh_pipeline(self):
        """Adopt the audited plan without overwriting pending local canvas edits."""
        project=getattr(self.ui,'current',None)
        if not project or project['id']!=self.results.project_id:return
        base=copy.deepcopy(getattr(self.ui,'base_pipeline',project['pipeline']))
        fresh=await self.ui.call(self.results.service.store.project,self.results.token,project['id'])
        if not self.results.valid() or self.ui.current is not project:return
        if fresh['revision']<project['revision']:return
        from biomolexplorer.project_merge import merge_pipeline
        merged=merge_pipeline(base,project['pipeline'],fresh['pipeline'])
        project.update(dict(fresh,pipeline=merged))
        self.ui.base_pipeline=copy.deepcopy(fresh['pipeline'])
        self.ui.dirty=merged!=fresh['pipeline']
        if self.ui.dirty and hasattr(self.ui,'schedule_save'):self.ui.schedule_save()
        if getattr(self.ui,'tab',None)=='Pipeline' and hasattr(self.ui,'flow_editor'):
            self.ui.flow_editor.refresh()

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
        recommended=context.get('recommended_ligands',[])
        choice=ft.Dropdown(label='Ligantes presentes na estrutura',hint_text='Selecione um resíduo para adicionar',
            enable_filter=True,enable_search=True,disabled=not self.editable)
        show_all=ft.Switch(label='Mostrar também os demais resíduos',
            value=not (context.get('triage_available') or recommended))
        evidence=verbatim(ft.Text(size=12))
        def filter_candidates(e=None):
            choice.options=[verbatim(ft.DropdownOption(key=str(i),text=f"{r['LIGAND']} · {r['RESNUM']}{r.get('ICODE','')} · {r['CHAIN']}"),'text')
                for i,r in enumerate(candidates) if show_all.value or r.get('RECOMMENDED')]
            choice.value=None;evidence.value=''
            if e:self.ui.page.update()
        def describe(e):
            if choice.value is not None:
                row=candidates[int(choice.value)]
                reason=REASONS.get(row.get('REASON'),'')
                if hasattr(self.ui,'tr'):reason=self.ui.tr(reason)
                evidence.value=f"{reason} · {row.get('CONTACT_RESIDUES',0)} " + (
                    self.ui.tr('resíduos da proteína a até 4 Å') if hasattr(self.ui,'tr') else 'resíduos da proteína a até 4 Å')
                chain_label=self.ui.tr('Cadeias em contato:') if hasattr(self.ui,'tr') else 'Cadeias em contato:'
                evidence.value+=f" · {chain_label} {row.get('RECEPTOR_CHAINS','')}"
                if row.get('LIGAND_NAME'):evidence.value+='\n'+row['LIGAND_NAME']
                evidence.tooltip=row.get('CONTACT_RESIDUE_IDS','')
                self.ui.page.update()
        show_all.on_change=filter_candidates;choice.on_select=describe
        filter_candidates()
        def add(e):
            if choice.value is None:return
            append(candidates[int(choice.value)])
            self.ui.page.update()
        def use_recommended(e):
            rows.clear();listing.controls.clear()
            for record in recommended:append(record)
            self.ui.page.update()
        async def save(e):
            async def action():
                if not result.valid():return
                records=[{'LIGAND':name.value,'RESNUM':number.value,'CHAIN':chain.value} for _,name,number,chain in rows]
                if hasattr(self.ui,'save_draft'):await self.ui.save_draft(silent=True)
                if not result.valid():return
                saved=await self.ui.call(result.service.set_pdb_ligands,*args,records,context['revision'])
                await self.refresh_pipeline()
                if not result.valid():return
                finish()
                if saved.get('updated_redocking_stages'):
                    self.ui.notify('Ligantes e pares de redocking atualizados no histórico do projeto. Execute o pipeline para processar as novas combinações.')
                else:self.ui.notify('Ligantes atualizados. A alteração foi registrada no histórico do projeto.')
            await self.ui.guard(action)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Ligantes · '+context['pdb_id']),on_dismiss=dismissed,
            content=ft.Container(width=640,height=380,content=ft.Column([
                ft.Text('Mantenha os ligantes de interesse. Para adicionar, informe o código, resíduo e cadeia presentes na estrutura.',size=12),
                listing,choice,ft.TextButton('Adicionar ligante',icon=ft.Icons.ADD,on_click=add,disabled=not self.editable),
                show_all,evidence,
                ft.Text('A lista recomendada usa contatos com a proteína e anotações do PDB. Contato geométrico não confirma sítio ativo.',size=12),
                ft.Text('Nenhum candidato recomendado. Revise os demais resíduos ou selecione outra estrutura.',size=12,visible=not bool(recommended)),
                ft.TextButton('Usar candidatos recomendados',on_click=use_recommended,disabled=not self.editable or not recommended)],
                scroll=ft.ScrollMode.AUTO,spacing=20)),
            actions=[ft.TextButton('Fechar',on_click=finish),
                     ft.Button('Salvar ligantes',on_click=save,disabled=not self.editable)])
        self.ui.editing_stage=True
        self.ui.page.show_dialog(dialog)
