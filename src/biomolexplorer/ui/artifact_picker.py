"""Compact source/file picker with bounded rendering and incremental selection."""
import flet as ft
from .localization import verbatim


class ArtifactPicker:
    PAGE_SIZE=25
    def __init__(self, rows, origins, page=None, extras=None, writable=True):
        self.rows,self.origins,self.page,self.extras=rows,origins,page,extras or {}
        self.writable=writable;self.touched=False;self.offset=0;self.available_offset=0
        self.on_change=None
        keys=list(origins)
        self.source=ft.Dropdown(label='Bloco de origem',options=[
            verbatim(ft.DropdownOption(key=k,text=origins.get(k,k)+' · '+k.split(':')[-1][:8]),'text') for k in keys],
            value=keys[0] if keys else None,enable_filter=True,enable_search=True,menu_height=320,on_select=self.source_changed,disabled=not writable,expand=True)
        self.search=ft.TextField(label='Pesquisar arquivos',on_change=self.search_changed,disabled=not writable,expand=True)
        self.available=ft.Dropdown(label='Arquivos disponíveis',options=[],enable_filter=True,enable_search=True,menu_height=320,expand=True,disabled=not writable)
        self.count=ft.Text(size=12,color='#64748B');self.feedback=ft.Text(size=12,color='#B45309')
        self.selected=ft.Column(spacing=8)
        self.previous=ft.TextButton('Anterior',on_click=lambda e:self.move(-1))
        self.next=ft.TextButton('Próxima',on_click=lambda e:self.move(1))
        self.available_previous=ft.TextButton('Arquivos anteriores',on_click=lambda e:self.move_available(-1))
        self.available_next=ft.TextButton('Mais arquivos',on_click=lambda e:self.move_available(1))
        self.control=ft.Column([ft.Row([self.source]),ft.Row([self.search]),
            ft.Row([self.available],expand=False),
            ft.Row([ft.TextButton('Adicionar arquivo',on_click=self.add_file,disabled=not writable),
                ft.TextButton('Adicionar todos da origem',on_click=self.add_all,disabled=not writable),
                ft.TextButton('Adicionar todas as origens',on_click=self.add_all_sources,disabled=not writable)],wrap=True),
            ft.Row([self.available_previous,self.available_next],wrap=True),self.feedback,
            ft.Row([ft.Text('Arquivos selecionados',weight=ft.FontWeight.W_600),
                ft.TextButton('Limpar seleção',on_click=self.clear,disabled=not writable)],wrap=True),
            self.count,self.selected,ft.Row([self.previous,self.next],wrap=True)],spacing=12)
        self.control._artifact_picker=self
        self.redraw()

    @staticmethod
    def origin(ref):return ('stage:'+ref['stage']) if 'stage' in ref else 'asset:'+ref['asset']

    @staticmethod
    def label(flag):return getattr(flag,'display_label',flag.label)

    def available_rows(self):
        query=(self.search.value or '').casefold()
        return [(i,flag,ref) for i,(flag,ref) in enumerate(self.rows)
            if self.origin(ref)==self.source.value and getattr(flag,'available',True) and query in self.label(flag).casefold()]

    def redraw(self):
        self.source.tooltip=self.origins.get(self.source.value,'')+' · '+str(self.source.value or '')
        available=self.available_rows()
        self.available_offset=min(self.available_offset,max(0,((len(available)-1)//self.PAGE_SIZE)*self.PAGE_SIZE))
        self.available.options=[verbatim(ft.DropdownOption(key=str(i),text=self.label(flag)),'text')
            for i,flag,_ in available[self.available_offset:self.available_offset+self.PAGE_SIZE]]
        if self.available.value not in {o.key for o in self.available.options}:self.available.value=None
        self.available_previous.disabled=self.available_offset==0
        self.available_next.disabled=self.available_offset+self.PAGE_SIZE>=len(available)
        chosen=[(i,flag,ref) for i,(flag,ref) in enumerate(self.rows) if flag.value]
        self.offset=min(self.offset,max(0,((len(chosen)-1)//self.PAGE_SIZE)*self.PAGE_SIZE))
        self.count.value=f'{len(available)} disponíveis nesta busca · {len(chosen)} selecionados · Página {self.offset//self.PAGE_SIZE+1}'
        controls=[]
        for i,flag,ref in chosen[self.offset:self.offset+self.PAGE_SIZE]:
            def remove(e,i=i):self.set_selected([i],False)
            origin=self.origins.get(self.origin(ref),self.origin(ref))
            controls.append(ft.Row([ft.Column([
                verbatim(ft.Text(self.label(flag),max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,tooltip=self.label(flag)),'value'),
                verbatim(ft.Text(origin+' · '+self.origin(ref).split(':')[-1][:8],size=11,color='#64748B',max_lines=1,overflow=ft.TextOverflow.ELLIPSIS,tooltip=origin+' · '+self.origin(ref)),'value'),
                *([ft.Text('Arquivo indisponível — remova ou selecione novamente.',size=12,color='#B91C1C')] if not getattr(flag,'available',True) else [])],expand=True),
                ft.IconButton(ft.Icons.REMOVE_CIRCLE_OUTLINE,tooltip='Remover arquivo',on_click=remove,disabled=not self.writable)],vertical_alignment=ft.CrossAxisAlignment.START))
            controls.extend(self.extras.get(id(flag),[]))
        self.selected.controls=controls or [ft.Text('Nenhum arquivo selecionado.',size=12,color='#64748B')]
        self.previous.disabled=self.offset==0;self.next.disabled=self.offset+self.PAGE_SIZE>=len(chosen)
        if not self.rows:self.feedback.value='Nenhuma origem com arquivos disponível. Execute ou conecte uma origem compatível.'
        elif not available:self.feedback.value='Esta origem não possui arquivos disponíveis para esta busca.'

    def update(self):
        self.redraw()
        if self.page:self.page.update()

    def source_changed(self,e):self.available_offset=0;self.feedback.value='';self.update()
    def search_changed(self,e):self.available_offset=0;self.feedback.value='';self.update()
    def move(self,direction):self.offset=max(0,self.offset+direction*self.PAGE_SIZE);self.update()
    def move_available(self,direction):self.available_offset=max(0,self.available_offset+direction*self.PAGE_SIZE);self.update()

    def set_selected(self,indices,value):
        if not self.writable:return
        self.touched=True
        changed=0;callbacks={}
        for i in indices:
            flag,_=self.rows[i]
            if value and not getattr(flag,'available',True):continue
            if flag.value==value:continue
            flag.value=value;changed+=1
            if getattr(flag,'on_change',None):callbacks[id(flag.on_change)]=flag.on_change
        for callback in callbacks.values():callback(None)
        if self.on_change:self.on_change()
        self.feedback.value=f'{changed} arquivo(s) '+('adicionado(s).' if value else 'removido(s).') if changed else 'Arquivo já adicionado ou nenhuma alteração.'
        self.update()

    def add_file(self,e=None):
        if self.available.value is None:self.feedback.value='Escolha um arquivo disponível.';self.update();return
        self.set_selected([int(self.available.value)],True)
    def add_all(self,e=None):self.set_selected([i for i,(_,ref) in enumerate(self.rows) if self.origin(ref)==self.source.value],True)
    def add_all_sources(self,e=None):self.set_selected(range(len(self.rows)),True)
    def clear(self,e=None):self.set_selected(range(len(self.rows)),False)
