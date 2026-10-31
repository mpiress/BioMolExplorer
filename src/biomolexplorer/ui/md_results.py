"""Paged MD status, failures and reuse with existing authorized file previews."""
import flet as ft
from biomolexplorer.md_results import MDResults
from .file_table import FileTable
from .localization import verbatim
from .feedback import close_dialog


class MDResultsTable:
    def __init__(self,results,rows):
        self.results=results;self.ui=results.ui;self.rows=rows;self.offset=0;self.limit=30;self.active=True
        self.service=MDResults(self.ui.store);self.children=[]
        self.count=ft.Text();self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(c)) for c in
            ('Receptor','Ligante','Estado','AMBER Score (kcal/mol)','Detalhes','Ações')])
        self.redraw()
    def valid(self):return self.active and self.results.valid()
    def redraw(self):
        self.table.rows=[]
        for row in self.rows[self.offset:self.offset+self.limit]:
            failed=row['status']=='failed'
            async def open(e,key=row['id']):await self.ui.guard(lambda:self.open(key))
            async def remove(e,key=row['id']):
                async def action():
                    if not self.valid():return
                    if getattr(self.ui,'dirty',False):raise ValueError('Salve as alterações do pipeline antes de remover um par.')
                    r=self.results
                    project=await self.ui.call(self.service.remove,r.token,r.project_id,r.run_id,r.stage['id'],key)
                    if not self.valid():return
                    self.ui.current=project
                    self.rows=await self.ui.call(self.service.pairs,r.token,r.project_id,r.run_id,r.stage['id'])
                    self.redraw();self.ui.page.update()
                await self.ui.guard(action)
            self.table.rows.append(ft.DataRow(color='#FEE2E2' if failed else None,cells=[
                ft.DataCell(verbatim(ft.Text(row['receptor']))),ft.DataCell(verbatim(ft.Text(row['ligand']))),
                ft.DataCell(ft.Text('Falha na dinâmica' if failed else 'Reutilizado' if row.get('reused') else 'Concluído',color='#B91C1C' if failed else '#166534')),
                ft.DataCell(ft.Text('—' if failed else f"{row['amber_score']:.4f}")),
                ft.DataCell(verbatim(ft.Text(row.get('error',''),width=280,max_lines=3,tooltip=row.get('error','')))),
                ft.DataCell(ft.Row([ft.IconButton(ft.Icons.FOLDER_OPEN,tooltip='Arquivos e gráficos do par',on_click=open),
                    ft.IconButton(ft.Icons.DELETE_OUTLINE,tooltip='Remover par da dinâmica',on_click=remove,disabled=not self.results.writable)]))]))
        self.count.value=str(len(self.rows))+' pares de dinâmica · '+str(sum(r['status']=='failed' for r in self.rows))+' falhas'
    def move(self,direction):
        self.offset=max(0,min(max(0,len(self.rows)-1),self.offset+direction*self.limit));self.redraw();self.ui.page.update()
    def build(self):
        async def compare(e):await self.ui.guard(self.compare)
        return ft.Column([self.count,ft.TextButton('Comparar pares por receptor e protocolo',on_click=compare),ft.Row([self.table],scroll=ft.ScrollMode.AUTO),
            ft.Row([ft.TextButton('Anterior',on_click=lambda e:self.move(-1)),ft.TextButton('Próxima',on_click=lambda e:self.move(1))])],spacing=16)
    def close(self):
        self.active=False
        for child in self.children:child.active=False
    async def compare(self):
        r=self.results
        rows=await self.ui.call(self.service.comparison,r.token,r.project_id,r.run_id,r.stage['id'])
        if not self.valid():return
        fields=('receptor','ligand','group','mean_protein_backbone_rmsd','mean_ligand_rmsd','mean_hydrogen_bonds','mean_protein_rg','mean_sasa','amber_score')
        table=ft.DataTable(columns=[ft.DataColumn(ft.Text(c)) for c in ('Receptor','Ligante','Grupo','RMSD backbone (Å)','RMSD ligante (Å)','Ligações H','Rg proteína (Å)','SASA (Å²)','AMBER Score')],
            rows=[ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text('—' if row[k] is None else f'{row[k]:.4g}' if isinstance(row[k],(int,float)) else row[k][:12] if k=='group' else row[k]))) for k in fields]) for row in rows])
        dialog=ft.AlertDialog(title=ft.Text('Comparação descritiva · sem ranking de afinidade'),
            content=ft.Container(ft.Column([ft.Text('Compare somente linhas do mesmo grupo: receptor e protocolo idênticos.'),ft.Row([table],scroll=ft.ScrollMode.AUTO)],scroll=ft.ScrollMode.AUTO),width=min(1000,max(280,(self.ui.page.width or 1000)-80)),height=450),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)
    async def open(self,key):
        r=self.results
        files=await self.ui.call(self.service.files,r.token,r.project_id,r.run_id,r.stage['id'],key)
        data=await self.ui.call(self.service.dashboard,r.token,r.project_id,r.run_id,r.stage['id'],key)
        if not self.valid():return
        async def preview(item):await self.ui.preview_artifact(r.project_id,item['path'])
        child=FileTable(self.ui,r.project_id,files,preview=preview);self.children.append(child)
        def finish(e):child.active=False;close_dialog(self.ui.page,dialog)
        summary=data.get('analysis',{}).get('summary',{})
        labels={'initial_complex_rg':('Raio de giro inicial do complexo','Å'),'final_complex_rg':('Raio de giro final do complexo','Å'),'hydrogen_atoms':('Átomos de hidrogênio',''),'frames':('Quadros',''),'simulation_duration_ps':('Duração simulada','ps'),
            'mean_protein_backbone_rmsd':('RMSD médio do backbone','Å'),
            'final_protein_backbone_rmsd':('RMSD final do backbone','Å'),
            'mean_ligand_rmsd':('RMSD médio do ligante','Å'),
            'final_ligand_rmsd':('RMSD final do ligante','Å'),
            'mean_ligand_site_distance':('Distância média ligante–sítio','Å'),
            'mean_complex_rg':('Raio de giro médio do complexo','Å'),
            'mean_protein_rg':('Raio de giro médio da proteína','Å'),
            'mean_hydrogen_bonds':('Ligações H médias na interface',''),
            'mean_sasa':('SASA média do complexo','Å²')}
        cards=[]
        for name,(label,unit) in labels.items():
            if name in summary:
                value=summary[name]
                cards.append(ft.Container(ft.Column([ft.Text(label,size=12),
                    ft.Text(f'{value:.4g} {unit}',size=20)]),padding=12,col={'xs':12,'sm':6,'lg':4}))
        definitions=data.get('analysis',{}).get('definitions',{})
        protocol=data.get('protocol',{})
        contents=[ft.Text('Refinamento curto em solvente implícito. Não equivale a produção equilibrada nem a energia livre de ligação.'),
            ft.Text(data.get('analysis',{}).get('unavailable','')),
            ft.ResponsiveRow(cards),
            ft.ExpansionTile(title=ft.Text('Definições e interpretação'),controls=[verbatim(ft.Text(name+': '+value)) for name,value in definitions.items()]),
            child.build()]
        for title,name in [('Resíduos mais flexíveis','most_flexible_residues'),('Contatos persistentes','most_persistent_contacts'),('Ligações H persistentes','most_persistent_hbonds')]:
            values=summary.get(name,[])
            if values:
                contents.insert(-1,ft.ExpansionTile(title=ft.Text(title),controls=[verbatim(ft.Text(
                    (str(value.get('name',''))+' '+str(value.get('residue','')) if 'residue' in value else str(value.get('donor',''))+' → '+str(value.get('acceptor','')))+
                    (f" · {value['occupancy']:.1%} dos quadros" if 'occupancy' in value else f" · RMSF {value['rmsf']:.3f} Å"+
                     (' · sítio inicial' if value.get('binding_site') else ' · fora do sítio inicial')))) for value in values]))
        energy=data.get('energies',{})
        if energy:
            contents.insert(-1,ft.ExpansionTile(title=ft.Text('Energias · kcal/mol · não são energia livre'),controls=[verbatim(ft.Text(f'{k}: {v:.5g}')) for k,v in energy.items() if isinstance(v,(int,float))]))
        sections={'Comparação estrutural inicial / final':('structural_rmsd','structural_rg'),'RMSD':('protein_backbone_rmsd','protein_rmsd','ligand_rmsd'),
            'RMSF':('rmsf',),'Ligações H · interface':('hydrogen_bonds',),
            'Raio de giro':('complex_rg','protein_rg'),'SASA':('sasa',),
            'Distância ligante–sítio':('ligand_site_distance',),'Energia potencial sander':('trajectory_energy',)}
        for title,names in sections.items():
            buttons=[]
            for name in names:
                item=next((f for f in files if f['name']==name+'.png'),None)
                if item:
                    async def show(e,item=item):await self.ui.guard(lambda:preview(item))
                    buttons.append(ft.TextButton(name,on_click=show))
            if buttons:contents.insert(-1,ft.ExpansionTile(title=ft.Text(title),controls=buttons))
        dialog=ft.AlertDialog(title=ft.Text('Análises, protocolo e arquivos do par'),
            content=ft.Container(ft.Column(contents,scroll=ft.ScrollMode.AUTO),width=min(1000,max(280,(self.ui.page.width or 1000)-80)),height=550),on_dismiss=lambda e:setattr(child,'active',False),
            actions=[ft.TextButton('Fechar',on_click=finish)])
        self.ui.page.show_dialog(dialog)
