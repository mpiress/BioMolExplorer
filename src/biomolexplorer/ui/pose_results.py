"""Validated poses as the block's own result table, with its files beside it."""
import flet as ft
from .feedback import close_dialog
from .localization import verbatim

APPROVED = '#166534'
REJECTED = '#B91C1C'
MUTED = '#64748B'
PAGE = 25
DETAILS = (('score','Score do motor',''),('min_vdw_ratio','Contato mais próximo · fração da soma dos raios vdW',''),
           ('min_vdw_overlap','Sobreposição de van der Waals','Å'),('min_heavy_distance','Distância mínima entre átomos pesados','Å'),
           ('clashing_contacts','Contatos avaliados no sítio',''),('delta_energy','Energia interna acima do confôrmero relaxado','kcal/mol'),
           ('max_bond_deviation','Desvio máximo relativo nas ligações',''),('max_angle_deviation','Desvio máximo nos ângulos','°'),
           ('max_planarity_deviation','Desvio máximo de planaridade','Å'),('planar_groups','Grupos planares avaliados',''),
           ('atoms','Átomos pesados da pose',''))
TOLERANCES = {'vdw_clash_scale':'Fração mínima da soma dos raios de van der Waals',
              'delta_energy_limit':'Energia interna máxima (kcal/mol)',
              'bond_tolerance':'Desvio relativo máximo nas ligações',
              'angle_tolerance':'Desvio máximo nos ângulos (°)',
              'planarity_tolerance':'Desvio máximo de planaridade (Å)',
              'minimization_steps':'Iterações da minimização de referência'}


def quantity(value, digits=4):
    if value is None or value == '':
        return '—'
    if isinstance(value,(int,float)):
        return str(value) if isinstance(value,int) else f'{value:.{digits}g}'
    return str(value)


def card(label, value, color=None):
    return ft.Container(ft.Column([ft.Text(label,size=12,color=MUTED),
        ft.Text(str(value),size=22,weight=ft.FontWeight.W_600,color=color)],spacing=4),
        padding=14,col={'xs':6,'sm':3})


class PoseValidationTable:
    """Paged validation table; the block's files stay collapsed beneath it."""

    def __init__(self, results, report, files=None):
        self.results = results
        self.ui = results.ui
        self.report = report
        self.files = files
        self.active = True
        self.offset = 0
        self.filter = 'all'
        self.count = ft.Text(size=13,color=MUTED)
        self.table = ft.DataTable(columns=[])
        self.pager = ft.Row(spacing=8)
        self.choice = ft.Dropdown(label='Poses exibidas',value='all',width=260,
            options=[ft.DropdownOption(key='all',text='Todas as poses'),
                     ft.DropdownOption(key='valid',text='Somente aprovadas'),
                     ft.DropdownOption(key='invalid',text='Somente reprovadas')],on_select=self.change_filter)
        self.redraw()

    def valid(self):
        return self.active and self.results.valid()

    def close(self):
        self.active = False
        if self.files is not None:
            self.files.active = False

    def rows(self):
        rows = self.report['rows']
        if self.filter == 'valid':
            return [row for row in rows if row['pose_valid']]
        if self.filter == 'invalid':
            return [row for row in rows if row['pose_valid'] is False]
        return rows

    async def change_filter(self, event):
        async def action():
            if not self.valid():
                return
            self.filter = self.choice.value or 'all'
            self.offset = 0
            self.redraw()
            self.ui.page.update()
        await self.ui.guard(action)

    def move(self, direction):
        rows = self.rows()
        self.offset = max(0,min(max(0,len(rows)-1),self.offset+direction*PAGE))
        self.redraw()
        self.ui.page.update()

    def redraw(self):
        rows = self.rows()
        self.table.columns = [ft.DataColumn(ft.Text(label),numeric=numeric) for label,numeric in
            (('Pose',False),('Alvo',False),('Composto',False),('Motor',False),('Score',True),
             ('Situação',False),('Critérios reprovados',False),('Contato mín. · fração vdW',True),
             ('ΔE interna (kcal/mol)',True),('Ações',False))]
        self.table.rows = []
        for row in rows[self.offset:self.offset+PAGE]:
            approved = bool(row['pose_valid'])
            failed = ', '.join(row['failed']) or ('—' if approved else 'Não avaliada')
            actions = [ft.IconButton(ft.Icons.INFO_OUTLINE,tooltip='Métricas e tolerâncias desta pose',
                                     on_click=self.detail_handler(row))]
            if row.get('pose_path') and row.get('receptor_path'):
                actions.append(ft.IconButton(ft.Icons.VIEW_IN_AR,tooltip='Abrir o complexo receptor–pose em 3D',
                                             on_click=self.complex_handler(row)))
            elif row.get('pose_path'):
                actions.append(ft.IconButton(ft.Icons.VISIBILITY_OUTLINED,tooltip='Abrir a pose em 3D',
                                             on_click=self.pose_handler(row)))
            self.table.rows.append(ft.DataRow(color='#FEF2F2' if not approved else None,cells=[
                ft.DataCell(verbatim(ft.Text(row['pose_id'] or '—',selectable=True))),
                ft.DataCell(verbatim(ft.Text(row['receptor_id'] or '—',selectable=True))),
                ft.DataCell(verbatim(ft.Text(row['compound_id'] or '—',selectable=True))),
                ft.DataCell(verbatim(ft.Text((row['engine'] or '—').upper()))),
                ft.DataCell(ft.Text(quantity(row['score']))),
                ft.DataCell(ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE if approved else ft.Icons.ERROR_OUTLINE,
                        color=APPROVED if approved else REJECTED,size=18),
                    ft.Text('Aprovada' if approved else 'Reprovada',color=APPROVED if approved else REJECTED)],spacing=6)),
                ft.DataCell(verbatim(ft.Text(failed,width=220,max_lines=2,
                    tooltip=row['rejection_details'] or failed))),
                ft.DataCell(ft.Text(quantity(row['min_vdw_ratio']))),
                ft.DataCell(ft.Text(quantity(row['delta_energy']))),
                ft.DataCell(ft.Row(actions,spacing=2))]))
        total = len(self.report['rows'])
        counts = self.report['counts']
        shown = f'{min(len(rows),self.offset+PAGE)-min(self.offset,len(rows))} de {len(rows)}'
        self.count.value = (f"{total} poses avaliadas · {counts['valid']} aprovadas · {counts['invalid']} reprovadas"
                            + (f' · exibindo {shown}' if len(rows) > PAGE else ''))
        self.pager.controls = [ft.TextButton('Anterior',on_click=lambda e:self.move(-1),disabled=self.offset <= 0),
            ft.TextButton('Próxima',on_click=lambda e:self.move(1),disabled=self.offset+PAGE >= len(rows))] if len(rows) > PAGE else []

    def detail_handler(self, row):
        async def handler(event):
            await self.ui.guard(lambda:self.detail(row))
        return handler

    def complex_handler(self, row):
        async def action():
            if not self.valid():
                return
            layers = [{'role':'receptor','path':row['receptor_path']},{'role':'pose','path':row['pose_path']}]
            name = ' · '.join(filter(None,[row['compound_id'],row['receptor_id'],(row['engine'] or '').upper()]))
            await self.ui.preview_complex(self.results.project_id,layers,name or row['pose_id'])

        async def handler(event):
            await self.ui.guard(action)
        return handler

    def pose_handler(self, row):
        async def handler(event):
            await self.ui.guard(lambda:self.ui.preview_artifact(self.results.project_id,row['pose_path']))
        return handler

    async def detail(self, row):
        if not self.valid():
            return
        approved = bool(row['pose_valid'])
        criteria = [ft.Row([ft.Icon(ft.Icons.CHECK if row.get(key) else ft.Icons.CLOSE if row.get(key) is False else ft.Icons.REMOVE,
                color=APPROVED if row.get(key) else REJECTED if row.get(key) is False else MUTED,size=18),
            ft.Text(label)],spacing=8) for key,label in self.report['criteria']]
        measurements = [ft.Row([ft.Text(label,expand=True,size=13),
            verbatim(ft.Text(quantity(row.get(key))+(' '+unit if unit and row.get(key) is not None else ''),size=13))],spacing=12)
            for key,label,unit in DETAILS if row.get(key) is not None]
        tolerances = self.report['summary'].get('tolerances') or {}
        contents = [ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE if approved else ft.Icons.ERROR_OUTLINE,
                color=APPROVED if approved else REJECTED),
            ft.Text('Pose aprovada' if approved else 'Pose reprovada',size=18,
                    color=APPROVED if approved else REJECTED)],spacing=8),
            verbatim(ft.Text(' · '.join(filter(None,[row['pose_id'],row['receptor_id'],row['compound_id'],
                (row['engine'] or '').upper()])),size=13,color=MUTED,selectable=True))]
        if row['rejection_details']:
            contents.append(verbatim(ft.Text(row['rejection_details'],color=REJECTED,selectable=True)))
        if row['error']:
            contents.append(verbatim(ft.Text(row['error'],color=REJECTED,size=12,selectable=True)))
        contents.extend([ft.Text('Critérios',size=15,weight=ft.FontWeight.W_600),*criteria,
            ft.Text('Medidas da pose',size=15,weight=ft.FontWeight.W_600),*measurements])
        if row['energy_method']:
            contents.append(verbatim(ft.Text('Campo de força da referência interna: '+row['energy_method'],size=12,color=MUTED)))
        if row['resolved_by']:
            contents.append(verbatim(ft.Text('Receptor associado por: '+row['resolved_by'],size=12,color=MUTED)))
        if tolerances:
            contents.append(ft.ExpansionTile(title=ft.Text('Tolerâncias aplicadas'),controls=[
                verbatim(ft.Text(TOLERANCES.get(key,key)+': '+quantity(value),size=13))
                for key,value in sorted(tolerances.items())]))
        dialog = ft.AlertDialog(title=ft.Text('Validação física e geométrica da pose'),
            content=ft.Container(ft.Column(contents,spacing=10,scroll=ft.ScrollMode.AUTO),
                width=min(760,max(280,(self.ui.page.width or 1000)-80)),height=min(560,max(300,(self.ui.page.height or 900)-260))),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)

    def summary_panel(self):
        summary = self.report['summary']
        counts = self.report['counts']
        published = summary.get('published_poses')
        strict = summary.get('strict_filter')
        cards = [card('Poses avaliadas',counts['total']),
                 card('Aprovadas',counts['valid'],APPROVED),
                 card('Reprovadas',counts['invalid'],REJECTED),
                 card('Entregues aos próximos blocos',published if published is not None else '—')]
        note = ('Somente as poses aprovadas seguiram para os próximos blocos.' if strict else
                'Todas as poses seguiram para os próximos blocos, com os metadados da validação.'
                if strict is False else '')
        rejections = summary.get('rejection_counts') or {}
        controls = [ft.ResponsiveRow(cards,run_spacing=4)]
        if note:
            controls.append(ft.Text(note,size=13,color=MUTED))
        if rejections:
            from biomolexplorer.pose_geometry import REASONS
            controls.append(ft.Text('Motivos de reprovação',size=15,weight=ft.FontWeight.W_600))
            controls.extend(verbatim(ft.Text(f'{REASONS.get(reason,reason)} · {total} pose(s)',size=13))
                            for reason,total in sorted(rejections.items(),key=lambda item:(-item[1],item[0])))
        unresolved = summary.get('unresolved_inputs') or 0
        if unresolved:
            controls.append(ft.Text(f'{unresolved} pose(s) sem receptor associado: o choque estérico não foi avaliado nelas.',
                                    size=13,color='#B45309'))
        if self.report['truncated']:
            controls.append(ft.Text('Relatório muito extenso: as primeiras poses são exibidas. Baixe o relatório completo nos arquivos deste bloco.',
                                    size=12,color='#B45309'))
        criteria = summary.get('criteria')
        if criteria:
            controls.append(ft.Text('Critérios inspirados em PoseBusters (Buttenschoen et al., Chem. Sci. 15, 3130, 2024) · '
                                    +str(criteria),size=12,color=MUTED))
        return ft.Column(controls,spacing=10)

    def build(self):
        controls = [self.summary_panel(),ft.Row([self.choice,ft.Container(self.count,expand=True)],
                                                vertical_alignment=ft.CrossAxisAlignment.CENTER),
                    ft.Row([self.table],scroll=ft.ScrollMode.AUTO),self.pager]
        if self.files is not None:
            controls.append(ft.ExpansionTile(title=ft.Text('Arquivos deste bloco'),
                subtitle=ft.Text('Relatório, poses aprovadas e receptores para download ou remoção',size=12,color=MUTED),
                controls=[ft.Container(self.files.build(),padding=20)]))
        return ft.Column(controls,spacing=16)
