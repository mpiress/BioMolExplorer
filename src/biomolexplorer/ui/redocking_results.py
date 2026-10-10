"""RMSD table with a simulation-specific file and molecular preview dialog."""
import flet as ft
from biomolexplorer.redocking_results import RedockingResults
from .feedback import close_dialog
from .file_table import FileTable
from .localization import verbatim


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
        self.table = ft.DataTable(columns=[ft.DataColumn(ft.Text(title), numeric=title == 'RMSD (Å)')
            for title in ('PDB', 'Ligante', 'Resíduo', 'Cadeia', 'RMSD (Å)', '3D', 'Resíduos', 'Ações')], rows=[],
            column_spacing=28, heading_row_color='#F1F5F9')
        self.count = ft.Text(size=13, color='#64748B')
        self.previous = ft.TextButton('Página anterior', on_click=lambda e: self.move(-1))
        self.next = ft.TextButton('Próxima página', on_click=lambda e: self.move(1))
        self.redraw()

    def valid(self):
        return self.active and self.results.valid()

    def close(self):
        self.active = False
        for table in self.file_tables:
            table.active = False

    def redraw(self):
        rows = []
        for simulation in self.simulations[self.offset:self.offset+self.limit]:
            async def view(e, key=simulation['id']):
                if self.valid():
                    await self.ui.guard(lambda: self.open(key))
            async def overlay(e,key=simulation['id']):
                if self.valid():await self.ui.guard(lambda:self.overlay(key))
            async def contacts(e,key=simulation['id']):
                if self.valid():await self.ui.guard(lambda:self.contacts(key))
            values = [simulation[k] for k in ('pdb', 'ligand', 'residue', 'chain')]
            rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(value))) for value in values] + [
                ft.DataCell(ft.Text(f"{simulation['rmsd']:.3f}")),
                ft.DataCell(ft.TextButton('3D',icon=ft.Icons.VIEW_IN_AR,on_click=overlay)),
                ft.DataCell(ft.TextButton('Resíduos',icon=ft.Icons.SCIENCE_OUTLINED,on_click=contacts)),
                ft.DataCell(ft.TextButton('Ver simulação', icon=ft.Icons.VISIBILITY, on_click=view))]))
        self.table.rows = rows
        self.count.value = f'{len(self.simulations)} simulações · Página {self.offset//self.limit+1}'
        self.previous.disabled = self.offset == 0
        self.next.disabled = self.offset+self.limit >= len(self.simulations)

    def move(self, direction):
        if not self.valid():
            return
        self.offset = max(0, self.offset+direction*self.limit)
        self.redraw()
        self.ui.page.update()

    def build(self):
        return ft.Column([self.count, ft.Row([self.table], scroll=ft.ScrollMode.AUTO),
            ft.Row([self.previous, self.next], alignment=ft.MainAxisAlignment.SPACE_BETWEEN)], spacing=16)

    async def overlay(self,key):
        if self.valid():
            await self.ui.preview_docking(self.results.project_id,self.results.run_id,self.results.stage['id'],
                'redocking',key,self.results.token)

    async def contacts(self,key):
        from .docking_scene import residue_dialog
        r=self.results
        await residue_dialog(self.ui,r.project_id,
            lambda:self.service.scene(r.token,r.project_id,r.run_id,r.stage['id'],key),
            self.valid,lambda:self.overlay(key))

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
        files = FileTable(self.ui, results.project_id, simulation['files'], extra_actions=preview_actions)
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
                    await self.ui.picker.save_file(file_name=name, src_bytes=data)
            await self.ui.guard(action)
        title = f"{simulation['pdb']} / {simulation['ligand']} / {simulation['residue']} / {simulation['chain']} · RMSD {simulation['rmsd']:.3f} Å"
        dialog = ft.AlertDialog(title=verbatim(ft.Text(title, size=20)), on_dismiss=dismissed,
            content=ft.Container(width=max(280, min(1100, (self.ui.page.width or 1440)-100)),
                height=max(240, min(650, (getattr(self.ui.page,'height',None) or 900)-240)),
                content=ft.Column([ft.Text('Consulte os arquivos da simulação. Use a opção 3D para examinar o receptor, o ligante de referência ou as poses.'),
                    files.build()], scroll=ft.ScrollMode.AUTO, spacing=16)),
            actions=[ft.TextButton('Baixar todos (ZIP)', icon=ft.Icons.DOWNLOAD, on_click=download_all),
                ft.TextButton('Fechar', on_click=finish)])
        self.ui.page.show_dialog(dialog)
