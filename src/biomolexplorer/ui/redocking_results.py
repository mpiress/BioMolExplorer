"""RMSD table with a simulation-specific file and molecular preview dialog."""
import flet as ft
import math
from biomolexplorer.i18n import RuntimeMessages
from .runtime_messages import session_language
from biomolexplorer.redocking_results import RedockingResults, best_rmsd_by_target
from .removal import RemovalProgress
from .processing import save_file
from .feedback import close_dialog
from .file_table import FileTable
from .localization import verbatim
from .table_sort import TableSort


class RedockingResultsTable:
    def __init__(self, results, simulations):
        self.results = results
        self.ui = results.ui
        self.simulations = simulations
        self.active = True
        self.service = RedockingResults(self.ui.store)
        self.file_tables = []
        self.offset = 0
        self.limit = 25
        self.sort_key = 'target'
        self.ascending = True
        self.sorting=TableSort('target')
        self.messages = RuntimeMessages(session_language(self.ui))
        targets=sorted({t for row in simulations for t in row.get('targets',[])},key=str.casefold)
        self.target_filter=ft.Dropdown(label=self.messages('redocking.filter.target'),value='all',
            options=[ft.DropdownOption(key='all',text=self.messages('redocking.filter.all'))]+
                [verbatim(ft.DropdownOption(key='target:'+t,text=t)) for t in targets],on_select=lambda e:self.filter_changed(),
            width=280)
        self.search=ft.TextField(label=self.messages('redocking.filter.search'),on_change=lambda e:self.filter_changed(),width=280)
        self.highlight=ft.Checkbox(label=self.messages('redocking.highlight'),value=True,on_change=lambda e:self.filter_changed())
        self.table = ft.DataTable(columns=[ft.DataColumn(ft.Text(title), numeric=title == 'RMSD (Å)')
            for title in ('Target', 'PDB', 'Ligante', 'Resíduo', 'Cadeia', 'Estado', 'RMSD (Å)', '3D', 'Resíduos', 'Ações')], rows=[],
            column_spacing=28, heading_row_color='#F1F5F9')
        self.table.columns[0]=verbatim(ft.DataColumn('Target'))
        self.sorting.bind(self.table,['target','pdb','ligand','residue','chain','status','rmsd',None,None,None],self.sorted,numeric={'residue','rmsd'})
        self.count = ft.Text(size=13, color='#64748B')
        self.warning = ft.Text(size=13,color='#B91C1C')
        self.previous = ft.TextButton('Página anterior', on_click=lambda e: self.move(-1))
        self.next = ft.TextButton('Próxima página', on_click=lambda e: self.move(1))
        self.export = ft.TextButton(self.messages('redocking.export'),icon=ft.Icons.DOWNLOAD,on_click=self.download_summary)
        self.redraw()

    def valid(self):
        return self.active and self.results.valid()

    def close(self):
        self.active = False
        for table in self.file_tables:
            table.active = False

    def filtered(self):
        target=self.target_filter.value
        query=(self.search.value or '').strip().casefold()
        rows=[r for r in self.simulations if (target=='all' or target[7:] in r.get('targets',[]))
            and (not query or query in ' '.join(str(r.get(k,'')) for k in ('target','pdb','ligand','residue','chain')).casefold())]
        rows=TableSort('rmsd').order(rows,lambda r,k:None if r.get('status')=='preparation_failed' else r.get(k),numeric={'rmsd'})
        return self.sorting.order(rows,lambda r,k:None if k=='rmsd' and r.get('status')=='preparation_failed' else r.get(k),numeric={'residue','rmsd'})

    def filter_changed(self):
        if not self.valid():return
        self.offset=0;self.redraw();self.ui.page.update()

    def sort(self,key):
        if not self.valid():return
        self.sorting.ascending=not self.sorting.ascending if self.sorting.key==key else True
        self.sorting.key=key
        self.sorting.indicator(self.table,['target','pdb','ligand','residue','chain','status','rmsd',None,None,None])
        self.sorted()

    def sorted(self):
        if not self.valid():return
        self.sort_key=self.sorting.key;self.ascending=self.sorting.ascending;self.offset=0
        self.redraw();self.ui.page.update()

    def redraw(self):
        rows = []
        filtered=self.filtered()
        last=max(0,(len(filtered)-1)//self.limit*self.limit)
        self.offset=min(max(0,self.offset),last)
        minima=best_rmsd_by_target(self.simulations)
        for simulation in filtered[self.offset:self.offset+self.limit]:
            unsupported=simulation.get('status')=='UNSUPPORTED'
            failed=simulation.get('status','succeeded')!='succeeded'
            async def view(e, key=simulation['id']):
                if self.valid():
                    await self.ui.guard(lambda: self.open(key))
            async def overlay(e,key=simulation['id']):
                if self.valid():await self.ui.guard(lambda:self.overlay(key))
            async def contacts(e,key=simulation['id']):
                if self.valid():await self.ui.guard(lambda:self.contacts(key))
            async def remove(e,key=simulation['id']):
                if self.valid():self.confirm_remove(key)
            values = [simulation.get('target') or '—']+[simulation[k] for k in ('pdb', 'ligand', 'residue', 'chain')]
            action=ft.TextButton('Ver simulação', icon=ft.Icons.VISIBILITY, on_click=view)
            if unsupported:action=ft.TextButton('Ver classificação',on_click=view)
            elif failed:action=ft.Row([ft.TextButton('Ver erro',on_click=view),
                ft.IconButton(ft.Icons.DELETE_OUTLINE,tooltip='Excluir falha e complexo de origem',on_click=remove,disabled=not self.results.writable)])
            status_label='Não suportado neste backend' if unsupported else 'Falha na execução' if failed and simulation.get('failed_step','').startswith(('Redocking','Docking','RMSD')) else 'Falha na preparação' if failed else 'Concluído'
            best=not failed and any(math.isclose(simulation['rmsd'],minima[t],rel_tol=0,abs_tol=1e-9) for t in simulation.get('targets',[]) if t in minima)
            rows.append(ft.DataRow(color='#FEF2F2' if failed else '#DCFCE7' if self.highlight.value and best else None,cells=[ft.DataCell(verbatim(ft.Text(value))) for value in values] + [
                ft.DataCell(ft.Text(status_label,
                    color='#B91C1C' if failed else '#15803D',tooltip=simulation.get('error') or None)),
                ft.DataCell(ft.Text('—' if failed else f"{simulation['rmsd']:.3f}")),
                ft.DataCell(self.view_button(simulation,overlay,failed)),
                ft.DataCell(ft.TextButton('Resíduos',icon=ft.Icons.SCIENCE_OUTLINED,on_click=contacts,disabled=failed)),
                ft.DataCell(action)]))
        self.table.rows = rows
        self.count.value = self.messages('redocking.count',shown=len(rows),total=len(self.simulations),
            start=self.offset+1 if rows else 0,end=self.offset+len(rows),filtered=len(filtered),
            page=self.offset//self.limit+1,pages=max(1,(len(filtered)+self.limit-1)//self.limit))
        self.previous.disabled = self.offset == 0
        self.next.disabled = self.offset+self.limit >= len(filtered)
        failures=sum(s.get('status','succeeded')!='succeeded' for s in self.simulations)
        self.warning.visible=bool(failures)
        self.warning.value=self.messages('redocking.failures',count=failures)

    def move(self, direction):
        if not self.valid():
            return
        self.offset = max(0, self.offset+direction*self.limit)
        self.redraw()
        self.ui.page.update()

    def build(self):
        # Keep navigation above long pages and inside the available width. The
        # wide scientific table scrolls internally rather than pushing Next
        # beyond the workspace's right edge.
        from .settings_layout import grid
        return ft.Column([grid([self.target_filter,self.search,self.highlight,self.export]),
            self.count,self.warning,
            ft.Row([self.previous,self.next],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
            ft.Row([self.table],scroll=ft.ScrollMode.ALWAYS)],spacing=16,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    async def download_summary(self,e):
        async def action():
            if not self.valid():return
            r=self.results
            data=await self.ui.call(self.service.rmsd_csv,r.token,r.project_id,r.run_id,r.stage['id'],
                [row['id'] for row in self.filtered()])
            if self.valid():await self.ui.picker.save_file(file_name='redocking_rmsd_by_target.csv',src_bytes=data)
        await self.ui.guard(action)

    def confirm_remove(self,key):
        if not self.valid() or not self.results.writable:return
        simulation=next(s for s in self.simulations if s['id']==key)
        previous_editing=getattr(self.ui,'editing_stage',False)
        self.ui.editing_stage=True
        def dismissed(e=None):self.ui.editing_stage=previous_editing
        def finish(e=None):
            dismissed();close_dialog(self.ui.page,dialog)
        async def remove(e):
            async def action():
                if not self.valid():return
                r=self.results
                await self.ui.call(self.service.remove_failure,r.token,r.project_id,r.run_id,r.stage['id'],key)
                if not self.valid():return
                finish()
                self.simulations=await self.ui.call(self.service.simulations,r.token,r.project_id,r.run_id,r.stage['id'])
                self.offset=min(self.offset,max(0,(len(self.simulations)-1)//self.limit*self.limit))
                self.redraw();self.ui.page.update()
                self.ui.notify('Falha e complexo de origem removidos. Os demais resultados foram preservados.')
            await progress.run(action)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Excluir falha e complexo de origem?'),on_dismiss=dismissed,
            content=ft.Column([verbatim(ft.Text(simulation['pdb']+'.pdb')),
                ft.Text('O PDB de origem será apagado e seus pares removidos da próxima execução deste redocking. Outros pares do mesmo PDB perderão essa entrada. Resultados já concluídos serão preservados. A recuperação original deixará de baixar esse PDB enquanto ele estiver na lista de exclusões.')],tight=True),
            actions=[ft.TextButton('Cancelar',on_click=finish),
                ft.TextButton('Excluir complexo',on_click=remove)])
        progress=RemovalProgress(self.ui,dialog,'Removendo complexo…')
        self.ui.page.show_dialog(dialog)

    def view_url(self,key):
        if not self.valid() or not getattr(self.ui.page,'web',False):return None
        r=self.results
        return self.ui.docking_view_url(r.project_id,r.run_id,r.stage['id'],'redocking',key,r.token)

    def view_button(self,simulation,callback,failed=False):
        url=self.view_url(simulation['id']) if not failed else None
        return ft.TextButton('3D',icon=ft.Icons.VIEW_IN_AR,disabled=failed,
            on_click=None if url else callback,url=ft.Url(url,target=ft.UrlTarget.BLANK) if url else None)

    async def overlay(self,key):
        if self.valid():
            await self.ui.preview_docking(self.results.project_id,self.results.run_id,self.results.stage['id'],
                'redocking',key,self.results.token)

    async def contacts(self,key):
        from .docking_scene import residue_dialog
        r=self.results
        await residue_dialog(self.ui,r.project_id,
            lambda:self.service.scene(r.token,r.project_id,r.run_id,r.stage['id'],key),
            self.valid,lambda:self.overlay(key),open_3d_url=self.view_url(key))

    async def open(self, key):
        results = self.results
        args = (results.token, results.project_id, results.run_id, results.stage['id'], key)
        simulation = await self.ui.call(self.service.simulation, *args)
        if not self.valid():
            return
        def preview_actions(item):
            if not item['molecular']:
                return []
            async def preview(e):
                if self.valid() and files.valid():
                    await self.ui.guard(lambda: self.ui.preview_artifact(results.project_id, item['path']))
            return [ft.IconButton(ft.Icons.VIEW_IN_AR, tooltip='Visualizar estrutura 3D', on_click=preview)]
        files = FileTable(self.ui, results.project_id, simulation['files'], extra_actions=preview_actions,run_id=results.run_id,stage_id=results.stage['id'])
        self.file_tables.append(files)
        previous_editing = getattr(self.ui, 'editing_stage', False)
        self.ui.editing_stage = True
        def dismissed(e=None):
            files.active = False
            if files in self.file_tables:
                self.file_tables.remove(files)
            self.ui.editing_stage = previous_editing
        def finish(e):
            dismissed()
            close_dialog(self.ui.page, dialog)
        async def download_all(e):
            async def action():
                if not self.valid() or not files.valid():
                    return
                data = await self.ui.call(self.service.archive, *args)
                if self.valid() and files.valid():
                    name = f"{simulation['pdb']}_{simulation['ligand']}_{simulation['residue']}{simulation['chain']}_redocking.zip"
                    await save_file(self.ui,file_name=name, src_bytes=data)
            await self.ui.guard(action)
        failed=simulation.get('status')!='succeeded';unsupported=simulation.get('status')=='UNSUPPORTED'
        failure_label='Não suportado neste backend' if unsupported else self.ui.tr('Falha na preparação') if hasattr(self.ui,'tr') else 'Falha na preparação'
        suffix=failure_label if failed else f"RMSD {simulation['rmsd']:.3f} Å"
        title = f"{simulation['pdb']} / {simulation['ligand']} / {simulation['residue']} / {simulation['chain']} · {suffix}"
        dialog = ft.AlertDialog(title=verbatim(ft.Text(title, size=20)), on_dismiss=dismissed,
            content=ft.Container(width=max(280, min(1100, (self.ui.page.width or 1440)-100)),
                height=max(240, min(650, (getattr(self.ui.page,'height',None) or 900)-240)),
                content=ft.Column([ft.Text('Consulte os arquivos da simulação. Use a opção 3D para examinar o receptor, o ligante de referência ou as poses.'),
                    *([verbatim(ft.Text(simulation.get('error',''),color='#B91C1C'))] if failed else []),
                    files.build()], scroll=ft.ScrollMode.AUTO, spacing=16)),
            actions=[ft.TextButton('Baixar todos (ZIP)', icon=ft.Icons.DOWNLOAD, on_click=download_all),
                ft.TextButton('Fechar', on_click=finish)])
        self.ui.page.show_dialog(dialog)
