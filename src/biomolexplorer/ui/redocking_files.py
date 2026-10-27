"""Bounded, target-grouped table for the redocking structure selection."""
from pathlib import Path
import flet as ft

from .artifact_picker import ArtifactPicker
from .localization import verbatim


class RedockingFilePicker(ArtifactPicker):
    def __init__(self,rows,origins,page=None,extras=None,writable=True):
        self.viewport_width=min(1180,max(360,(getattr(page,'width',None) or 1148)-148))
        self.target=ft.Dropdown(label='Alvo dos arquivos selecionados',value='all',
            options=[],enable_filter=True,enable_search=True,on_select=self.filter_changed)
        self.review_filter=ft.Dropdown(label='Revisão dos ligantes',value='all',options=[
            ft.DropdownOption(key='all',text='Todos os estados'),
            ft.DropdownOption(key='pending',text='A revisar'),
            ft.DropdownOption(key='reviewed',text='Revisados'),
            ft.DropdownOption(key='unavailable',text='Indisponíveis')],on_select=self.filter_changed)
        for field in (self.target,self.review_filter):
            field.content_padding=ft.Padding.symmetric(horizontal=16,vertical=18)
            field.expand=True
        self.overview=ft.Text(size=13,color='#475569')
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(label)) for label in
            ('PDB / arquivo','Alvo / origem','Ligantes selecionados','Revisão','Ações')],rows=[],
            column_spacing=20,data_row_min_height=64,data_row_max_height=64,
            heading_row_color='#EFF6FF',border=ft.Border.all(1,'#E2E8F0'),border_radius=12)
        self.table_host=ft.Container(content=ft.Row([self.table],scroll=ft.ScrollMode.AUTO),
            on_size_change=self.resize_table,size_change_interval=100)
        super().__init__(rows,origins,page,extras,writable)
        for control in self.control.controls:
            for child in getattr(control,'controls',[]) or []:
                if isinstance(child,ft.TextButton) and child.content=='Limpar seleção':child.content='Limpar toda a seleção'
        original=self.control.controls
        self.search.label='Pesquisar PDB, alvo ou ligante'
        self.search.disabled=False
        self.search.content_padding=ft.Padding.symmetric(horizontal=16,vertical=18)
        original[6].spacing=24
        original[6].alignment=ft.MainAxisAlignment.SPACE_BETWEEN
        self.add_panel=ft.ExpansionTile(title=ft.Text('Adicionar PDBs à seleção'),
            subtitle=ft.Text('Escolha uma origem e adicione arquivos individualmente ou em lote.',size=12),
            expanded=not any(flag.value for flag,_ in rows),maintain_state=True,
            tile_padding=16,controls_padding=ft.Padding.only(left=16,right=16,top=12,bottom=20),
            shape=ft.RoundedRectangleBorder(radius=12,side=ft.BorderSide(1,'#E2E8F0')),
            collapsed_shape=ft.RoundedRectangleBorder(radius=12,side=ft.BorderSide(1,'#E2E8F0')),
            controls=[ft.Column([original[i] for i in (0,2,3,4,5)],spacing=24,
                horizontal_alignment=ft.CrossAxisAlignment.STRETCH)])
        filters=ft.ResponsiveRow([
            ft.Container(ft.Row([self.target]),col={'xs':12,'md':8}),
            ft.Container(ft.Row([self.review_filter]),col={'xs':12,'md':4})],spacing=24,run_spacing=24)
        self.control.spacing=24
        self.control.horizontal_alignment=ft.CrossAxisAlignment.STRETCH
        self.selected.horizontal_alignment=ft.CrossAxisAlignment.STRETCH
        self.control.controls=[ft.Container(original[1],padding=ft.Padding.symmetric(vertical=8)),self.add_panel,
            ft.Column([original[6],
                ft.Text('Organize os PDBs por alvo e revise os ligantes sem perder a seleção. Verde indica uma revisão salva, não a validação da preparação.',size=12,color='#64748B')],spacing=12),
            ft.Container(self.overview,padding=16,bgcolor='#F8FAFC',border_radius=12),
            ft.Container(filters,padding=ft.Padding.symmetric(vertical=8)),
            ft.Row([ft.TextButton('Remover desta visualização',on_click=self.remove_view,disabled=not writable),
                self.count],alignment=ft.MainAxisAlignment.SPACE_BETWEEN,spacing=24,wrap=True),
            self.selected,original[9]]

    def resize_table(self,event):
        width=max(1,int(event.width))
        if abs(width-self.viewport_width)<2:return
        self.viewport_width=width
        self.redraw()
        if self.page:self.page.update()

    def table_layout(self):
        compact=self.viewport_width<900
        width=max(360,self.viewport_width)
        self.table.width=width
        self.table.horizontal_margin=16 if compact else 20
        self.table.column_spacing=16 if compact else 24
        self.table.data_row_min_height=self.table.data_row_max_height=112 if compact else 80
        titles=('PDB / alvo','Ligantes e revisão','Ações') if compact else (
            'PDB / arquivo','Alvo / origem','Ligantes selecionados','Revisão','Ações')
        self.table.columns=[ft.DataColumn(ft.Text(label)) for label in titles]
        if compact:
            available=width-32-32-88
            return True,(available*.46,available*.54,88)
        available=width-40-96-120-132-192
        return False,(120,available*.48,available*.52,132,192)

    def group(self,flag,ref):
        origin=self.origin(ref)
        label=getattr(flag,'target_label',None) or self.origins.get(origin,origin)
        return origin+'|'+label,label

    def origin_label(self,flag,ref):
        return getattr(flag,'origin_label',None) or self.origins.get(self.origin(ref),'')

    @staticmethod
    def status(flag):
        return 'unavailable' if not getattr(flag,'available',True) else 'reviewed' if getattr(flag,'reviewed',False) else 'pending'

    def available_rows(self):
        query=(self.search.value or '').casefold()
        return [(i,flag,ref) for i,(flag,ref) in enumerate(self.rows)
            if self.origin(ref)==self.source.value and getattr(flag,'available',True)
            and query in ' '.join((self.label(flag),flag.label,self.group(flag,ref)[1],getattr(flag,'ligand_summary',''))).casefold()]

    def selected_rows(self):
        query=(self.search.value or '').casefold()
        return [(i,flag,ref) for i,(flag,ref) in enumerate(self.rows) if flag.value
            and (self.target.value=='all' or self.group(flag,ref)[0]==self.target.value)
            and (self.review_filter.value=='all' or self.status(flag)==self.review_filter.value)
            and query in (' '.join((self.label(flag),flag.label,self.group(flag,ref)[1],getattr(flag,'ligand_summary','')))).casefold()]

    def redraw(self):
        super().redraw()
        groups={};selected=reviewed=unavailable=0
        for flag,ref in self.rows:
            if not flag.value:continue
            selected+=1;done=self.status(flag)=='reviewed';reviewed+=done
            unavailable+=self.status(flag)=='unavailable'
            key,label=self.group(flag,ref)
            entry=groups.setdefault(key,[label,0,0,self.origin_label(flag,ref)])
            entry[1]+=1;entry[2]+=done
        self.target.options=[ft.DropdownOption(key='all',text='Todos os alvos')]+[
            verbatim(ft.DropdownOption(key=key,text=f'{label} · {done}/{total} '+self.translate('revisados')+f' · {origin}'),'text')
            for key,(label,total,done,origin) in sorted(groups.items(),key=lambda item:(item[1][0].casefold(),item[0]))]
        if self.target.value not in {'all',*groups}:self.target.value='all'
        self.overview.value=f'{selected} PDBs selecionados · {len(groups)} alvos/origens · {reviewed} revisados · {selected-reviewed-unavailable} a revisar · {unavailable} indisponíveis'
        chosen=self.selected_rows()
        self.offset=min(self.offset,max(0,((len(chosen)-1)//self.PAGE_SIZE)*self.PAGE_SIZE))
        compact,widths=self.table_layout()
        table_rows=[]
        for i,flag,ref in chosen[self.offset:self.offset+self.PAGE_SIZE]:
            status=self.status(flag);done=status=='reviewed'
            color='#15803D' if done else '#B91C1C' if status=='unavailable' else '#B45309'
            label='Revisado' if done else 'Indisponível' if status=='unavailable' else 'A revisar'
            name=Path(getattr(flag,'physical_file',flag.label)).name
            def remove(e,index=i):self.set_selected([index],False)
            actions=list(self.extra_controls(flag))+[ft.IconButton(ft.Icons.REMOVE_CIRCLE_OUTLINE,
                tooltip='Remover arquivo',on_click=remove,disabled=not self.writable)]
            filename=verbatim(ft.Text(name,weight=ft.FontWeight.W_600,max_lines=2,
                overflow=ft.TextOverflow.ELLIPSIS,tooltip=flag.label),'value')
            origin=[verbatim(ft.Text(self.group(flag,ref)[1],max_lines=1,
                    overflow=ft.TextOverflow.ELLIPSIS,tooltip=self.group(flag,ref)[1]),'value'),
                    verbatim(ft.Text(self.origin_label(flag,ref),size=11,color='#64748B',max_lines=1,overflow=ft.TextOverflow.ELLIPSIS,tooltip=self.origin_label(flag,ref)),'value')]
            ligands=verbatim(ft.Text(getattr(flag,'ligand_summary','') or '—',max_lines=2,
                overflow=ft.TextOverflow.ELLIPSIS,tooltip=getattr(flag,'ligand_summary','')),'value')
            review=ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE_OUTLINE if done else
                    ft.Icons.ERROR_OUTLINE if status=='unavailable' else ft.Icons.PENDING_OUTLINED,color=color,size=18),
                    ft.Text(label,color=color,size=12)],spacing=8)
            if compact:
                controls=[ft.Column([filename,*origin],tight=True,spacing=6),
                    ft.Column([ligands,review],tight=True,spacing=12),
                    ft.Column([ft.Row(actions[index:index+2],spacing=4) for index in range(0,len(actions),2)],
                        tight=True,spacing=4)]
            else:
                controls=[filename,ft.Column(origin,tight=True,spacing=6),ligands,review,ft.Row(actions,spacing=4)]
            for control,width in zip(controls,widths):control.width=width
            table_rows.append(ft.DataRow(color='#F0FDF4' if done else '#FFFFFF',cells=[
                ft.DataCell(control) for control in controls]))
        self.table.rows=table_rows
        self.selected.controls=[self.table_host] if chosen else [
            ft.Text('Nenhum PDB selecionado corresponde aos filtros.' if selected else 'Nenhum arquivo selecionado.',size=12,color='#64748B')]
        start=self.offset+1 if chosen else 0
        self.count.value=f'{start}–{min(self.offset+self.PAGE_SIZE,len(chosen))} de {len(chosen)} PDBs nesta visualização'
        self.previous.disabled=self.offset==0;self.next.disabled=self.offset+self.PAGE_SIZE>=len(chosen)
        if not selected and hasattr(self,'add_panel'):self.add_panel.expanded=True

    def filter_changed(self,e=None):self.offset=0;self.update()
    def search_changed(self,e):self.offset=0;super().search_changed(e)
    def remove_view(self,e=None):self.set_selected([i for i,_,_ in self.selected_rows()],False)

    def set_selected(self,indices,value):
        if self.writable and hasattr(self,'add_panel'):
            if value:self.add_panel.expanded=False
            elif not any(flag.value for flag,_ in self.rows):self.add_panel.expanded=True
        super().set_selected(indices,value)

    def translate(self,value):
        translator=getattr(self.page,'translator',None)
        return translator(value) if translator else value

    def apply_reviews(self,statuses):
        for flag,_ in self.rows:
            review=statuses.get(getattr(flag,'physical_file',''))
            if review is None:continue
            flag.reviewed=bool(review.get('reviewed'))
            flag.ligand_summary='; '.join(f"{r['LIGAND']} · {r['RESNUM']} · {r['CHAIN']}" for r in review.get('ligands',[]))
        self.update()

    def view_state(self):
        return dict(target=self.target.value,review=self.review_filter.value,search=self.search.value or '',offset=self.offset,source=self.source.value)

    def restore_view(self,state):
        self.target.value=state.get('target','all')
        self.review_filter.value=state.get('review','all') if state.get('review','all') in ('all','pending','reviewed','unavailable') else 'all'
        self.search.value=state.get('search','')
        if state.get('source') in self.origins:self.source.value=state['source']
        self.offset=max(0,int(state.get('offset',0)))
        self.redraw()
