"""Compact selected-compound card and explicit, bounded scaffold exploration."""
import base64
import threading
import flet as ft
from biomolexplorer.molecular_network import numeric_properties
from biomolexplorer.network_scaffolds import individual,common
from biomolexplorer.network_analysis import METRICS
from .feedback import close_dialog
from .network_charts import card
from .localization import verbatim


class MolecularInspector:
    def __init__(self,explorer):
        self.explorer=explorer;self.view=explorer.view;self.properties={};self.structure={};self.dialog=None;self.relations_panel=None
        self.summary=ft.Column([ft.Text('Selecione um composto no grafo ou no ranking.',size=13)],spacing=8)
        self.scaffold=ft.Column(spacing=8);self.result_panel=ft.Column(spacing=8)
        self.scope=ft.Dropdown(label='Conjunto do MCS',value='mcc',expand=True,
            options=[ft.DropdownOption(key=k,text=v) for k,v in [('full','Full Graph'),('mcc','LCC'),('community','Comunidade selecionada'),('filtered','Subconjunto filtrado')]],on_select=self.change_scope)
        self.timeout=ft.TextField(label='Tempo máximo (s)',value='10',expand=True,on_change=self.change_parameters)
        self.limit=ft.TextField(label='Máximo de compostos',value='50',expand=True,on_change=self.change_parameters)
        self.timeout.tooltip='Limite da busca RDKit: 1–30 segundos. Preparação e desenho podem acrescentar tempo; um resultado parcial não confirma maximalidade.'
        self.limit.tooltip='Limite de segurança: 2–200 compostos. Restrinja conjuntos maiores; nenhum subconjunto é amostrado automaticamente.'
        self.running=False;self.result=None;self.cache={};self.structures={};self.cancelled=None;self.request=0;self.pending_scope=False
        self.compute_button=ft.TextButton('Calcular MCS do conjunto',on_click=self.calculate)
        self.cancel_button=ft.TextButton('Interromper MCS',on_click=self.cancel,disabled=True)
        self.notice=ft.Text(size=12,color='#64748B')
        saved=self.view.data.get('fragment')
        if saved and saved.get('smarts') and set(self.view.data.get('mcc',[]))==set(self.explorer.analysis.lcc) and saved.get('compounds')==len(self.explorer.analysis.lcc):
            self.result=dict(saved,ids=sorted(self.explorer.analysis.lcc),scope='LCC · resultado persistido',notice=('Tempo limite: maximalidade não confirmada. ' if saved.get('status')=='partial' else '')+'Resultado reutilizado do pipeline; os critérios originais estão preservados.')
            if not self.result.get('smiles'):
                self.result['smiles']=next((c.value for c in self.view.fragment_panel.controls if isinstance(c,ft.TextField) and c.label=='SMILES do fragmento'),'')
            self.show_result()

    def clear(self):
        self.summary.controls=[ft.Text('Selecione um composto no grafo ou no ranking.',size=13)]
        self.scaffold.controls=[];self.structure={};self.properties={}
        self.explorer.molecular.neighbor.options=[];self.explorer.molecular.neighbor.value=None;self.explorer.molecular.neighbor.disabled=True

    async def update(self,identifier,properties,version):
        view=self.view
        smarts=self.result.get('smarts','') if self.result and identifier in self.result['ids'] else ''
        smiles=properties.get('canonical_smiles')
        key=(smiles if isinstance(smiles,str) else None,smarts)
        if key not in self.structures:
            structure=await view.ui.call(individual,key[0],smarts)
            if len(self.structures)>=32:self.structures.pop(next(iter(self.structures)))
            self.structures[key]=structure
        if version!=view.selection_version or not view.active or view.ui.token!=view.token:return
        self.structure=self.structures[key];self.properties=properties
        neighbor=self.explorer.molecular.neighbor
        neighbors=sorted(view.graph.neighbors(identifier))
        neighbor.options=[ft.DropdownOption(key=n,text=n) for n in neighbors];neighbor.value=neighbors[0] if neighbors else None;neighbor.disabled=not bool(neighbors)
        graph=self.explorer.analysis.universe(view.mode)
        values=self.explorer.analysis.cache.get((view.mode,self.explorer.metric),({},''))[0]
        community=self.explorer.analysis.community_cache.get(view.mode,({},None,''))[0].get(identifier)
        async def locate(e):await view.center_on(identifier)
        self.summary.controls=[verbatim(ft.Text(identifier,size=18,weight=ft.FontWeight.W_600,selectable=True)),
            ft.Image(src=self.structure['image'],height=190,fit=ft.BoxFit.CONTAIN) if self.structure['image'] else ft.Text(self.structure['notice'],size=12),
            ft.Text(f'{view.mode} · Degree: {graph.degree(identifier) if identifier in graph else "fora do universo"} · C{community}' if community is not None else f'{view.mode} · Degree: {graph.degree(identifier) if identifier in graph else "fora do universo"}',size=12),
            *([ft.Text(METRICS[self.explorer.metric].split(' · ')[0]+': '+(f'{values[identifier]:.5g}' if identifier in values else 'indisponível'),size=12)] if self.explorer.metric!='degree' else []),
            ft.Row([ft.TextButton('Localizar',on_click=locate),ft.TextButton('Vizinhança',on_click=self.explorer.ego_selected)],wrap=True),
            ft.TextButton('View Molecular Details',on_click=self.details)]
        self.scaffold.controls=[ft.Text('Bemis–Murcko · scaffold individual',size=12),
            ft.Image(src=self.structure['scaffold'],height=125,fit=ft.BoxFit.CONTAIN) if self.structure['scaffold'] else ft.Text(self.structure['notice'],size=12),
            ft.Text(f'{self.structure["atoms"]} átomos · núcleo de anéis e conectores; não implica fragmento compartilhado.',size=11)]
        if smarts:self.scaffold.controls.append(ft.Text(f'MCS destacado: {self.structure.get("highlighted_atoms",0)} átomos no composto.',size=11))

    async def details(self,e=None):
        view=self.view
        async def action():
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if not view.active or view.ui.token!=view.token or view.selected is None:return
            page=view.ui.page
            width=max(260,min(1040,(getattr(page,'width',None) or 1200)-48))
            height=max(240,min(740,(getattr(page,'height',None) or 900)-140))
            structure=ft.Column([ft.Image(src=self.structure['image'],height=260,fit=ft.BoxFit.CONTAIN)] if self.structure.get('image') else [ft.Text('Estrutura indisponível')])
            if self.structure.get('scaffold'):structure.controls.extend([ft.Text('Bemis–Murcko'),ft.Image(src=self.structure['scaffold'],height=170,fit=ft.BoxFit.CONTAIN)])
            async def open_3d(e):await view.ui.guard(lambda:view.show_3d(view.selected,view.selection_version))
            structure.controls.append(ft.TextButton('Visualizar estrutura 3D',on_click=open_3d,disabled=not bool(self.structure.get('image'))))
            for key in ('molecule_chembl_id','canonical_smiles','standard_inchi','standard_inchi_key','inchi','inchikey'):
                if self.properties.get(key):structure.controls.extend([ft.Text(key,size=11),verbatim(ft.Text(str(self.properties[key]),size=12,selectable=True))])
            if self.result and self.result.get('smarts') and view.selected in self.result['ids']:
                structure.controls.extend([ft.Text('MCS · '+self.result['scope'],size=12),ft.Text(str(self.result['atoms'])+' átomos · '+str(self.result['compounds'])+' compostos · '+self.result['status'],size=11)])
                if self.result.get('image'):structure.controls.append(ft.Image(src=base64.b64decode(self.result['image']),height=150,fit=ft.BoxFit.CONTAIN))
                for key in ('smiles','smarts'):
                    if self.result.get(key):structure.controls.append(ft.TextField(label=key.upper()+' · fragmento',value=self.result[key],read_only=True,multiline=True,min_lines=1,max_lines=3))
            labels={'MW':('Molecular Weight','g/mol'),'molecule_properties.full_mwt':('Molecular Weight','g/mol'),
                'molecule_properties.mw_freebase':('Free-base Molecular Weight','g/mol'),'molecule_properties.alogp':('ALogP',''),'WLOGP':('WLOGP',''),'TPSA':('TPSA','Å²'),
                'molecule_properties.psa':('PSA','Å²'),'molecule_properties.hbd':('H-Bond Donors',''),
                'molecule_properties.hba':('H-Bond Acceptors',''),'molecule_properties.rtb':('Rotatable Bonds','')}
            rows=[]
            for key,value in self.properties.items():
                if value is None or value=='':continue
                if key in ('molecule_chembl_id','canonical_smiles','standard_inchi','standard_inchi_key','inchi','inchikey','degree','degree_density','component'):continue
                if key in ('molecule_properties','properties') and any(k.startswith(key+'.') for k in self.properties):continue
                name,unit=labels.get(key,(key,''))
                origin='grafo' if view.nodes[view.selected]['properties'].get(key.split('.')[0]) is not None else 'fonte selecionada'
                rows.append(ft.ResponsiveRow([ft.Text(name,col={'xs':12,'sm':5},size=12,tooltip=origin),
                    verbatim(ft.Text(str(value)+(' '+unit if unit else ''),col={'xs':12,'sm':7},size=12,selectable=True))],spacing=8,run_spacing=2))
            network=[ft.Text(f'{key} · registrado no Full Graph: {self.properties[key]}',size=12) for key in ('degree','degree_density','component') if self.properties.get(key) is not None]
            for (mode,metric),(values,notice) in self.explorer.analysis.cache.items():
                if mode==view.mode and view.selected in values:network.append(ft.Text(f'{metric}: {values[view.selected]:.6g}',tooltip=notice))
            component=next((i for i,g in view.components.items() if view.selected in g),None)
            mapping=self.explorer.analysis.community_cache.get(view.mode,({},None,''))[0]
            network.extend([ft.Text(f'{view.mode} · Componente: {component}'),ft.Text('Comunidade: '+str(mapping.get(view.selected,'não calculada')))])
            # Full relationship controls are mounted only in this internally scrolling modal.
            relation_index=next((i for i,c in enumerate(view.details.controls) if isinstance(c,ft.Text) and str(c.value).startswith('Relações deste composto')),len(view.details.controls))
            self.relations_panel=ft.Column(view.details.controls[relation_index:],spacing=6)
            source=(view.molecular_enrichment or {}).get('source',view.data.get('dashboard_context',{}))
            body=ft.Column([ft.Text('Origem: '+str(source.get('stage_name') or source.get('stage_id') or 'grafo persistido')+' · execução '+str(source.get('run_id','não registrada')),size=11),
                ft.ResponsiveRow([card('Molecular Structure',structure,{'xs':12,'lg':5}),card('Molecular Properties',ft.Column(rows,spacing=8),{'xs':12,'lg':7})]),
                card('Network Properties',ft.Row(network,wrap=True),12),self.relations_panel],scroll=ft.ScrollMode.AUTO,spacing=14)
            def dismissed(e=None):
                if self.dialog is dialog:self.dialog=None;self.relations_panel=None
            def close(e=None):dismissed();close_dialog(page,dialog)
            dialog=ft.AlertDialog(on_dismiss=dismissed,title=verbatim(ft.Text(view.selected+' · molécula 2D',size=20)),content=ft.Container(width=width,height=height,content=body),
                actions=[ft.TextButton('Fechar molécula',on_click=close)])
            self.dialog=dialog;page.show_dialog(dialog)
        await view.ui.guard(action)

    def refresh_relations(self,controls):
        if self.relations_panel is None:return
        start=next((i for i,c in enumerate(controls) if isinstance(c,ft.Text) and str(c.value).startswith('Relações deste composto')),len(controls))
        self.relations_panel.controls=controls[start:]

    def close_details(self):
        if self.dialog:close_dialog(self.view.ui.page,self.dialog)
        self.dialog=None;self.relations_panel=None

    def cancel(self,e=None):
        self.request+=1
        if self.cancelled:self.cancelled.set()
        self.cancel_button.disabled=True;self.compute_button.disabled=True
        self.notice.value='Interrupção solicitada; a busca termina no próximo callback do RDKit.'
        self.view.ui.page.update()

    def close(self):
        self.close_details()
        self.request+=1;self.pending_scope=False
        if self.cancelled:self.cancelled.set()
        self.cache.clear();self.structures.clear()

    def invalidate(self):
        self.request+=1
        if self.cancelled:self.cancelled.set()
        self.result=None;self.result_panel.controls=[]
        self.notice.value='Conjunto ou parâmetros alterados. Calcule o MCS para confirmar o novo resultado.'

    async def refresh_selected(self):
        view=self.view
        if view.active and view.ui.token==view.token and view.selected in view.displayed():
            await view.select(view.selected,show_popup=False)

    async def change_scope(self,e=None):
        if not self.view.active or self.view.ui.token!=self.view.token:return
        running=self.running
        self.invalidate();self.pending_scope=running
        await self.refresh_selected()
        if not running:await self.calculate()
        self.view.ui.page.update()

    async def change_parameters(self,e=None):
        if not self.view.active or self.view.ui.token!=self.view.token:return
        self.invalidate();self.pending_scope=False
        try:
            timeout,limit=int(self.timeout.value),int(self.limit.value)
            if not 1<=timeout<=30 or not 2<=limit<=200:
                self.notice.value='Use 1–30 segundos e limite de 2–200 compostos.'
        except (ValueError,TypeError):self.notice.value='Informe tempo e limite como números inteiros.'
        await self.refresh_selected();self.view.ui.page.update()

    def show_result(self):
        result=self.result
        if not result:return
        status={'complete':'Concluído','partial':'Parcial','unavailable':'Indisponível','cancelled':'Interrompido','no_common_fragment':'Sem fragmento comum'}.get(result['status'],result['status'])
        information=ft.Column([ft.Text(result['scope'],size=12,weight=ft.FontWeight.W_600),
            ft.Text(f'{result["compounds"]} compostos · {result["atoms"]} átomos',size=12),ft.Text(status,size=12),
            ft.Text(result.get('notice',''),size=12)],spacing=8)
        image=ft.Image(src=base64.b64decode(result['image']),height=150,fit=ft.BoxFit.CONTAIN) if result.get('image') else ft.Text('Estrutura indisponível',size=12)
        # The inspector is narrow even on desktop: the result uses its own
        # two-column row, while each parameter gets a full-width bounded row.
        self.result_panel.controls=[ft.Row([ft.Container(image,expand=True),ft.Container(information,expand=True)],spacing=12,vertical_alignment=ft.CrossAxisAlignment.START),
            ft.Text('Critérios: cobertura 100%',size=11),
            ft.Text('Anéis com anéis: '+('Sim' if result.get('ring_matches_ring_only') else 'Não' if result.get('ring_matches_ring_only') is False else 'não registrado'),size=11),
            ft.Text('Anéis completos: '+('Sim' if result.get('complete_rings_only') else 'Não' if result.get('complete_rings_only') is False else 'não registrado'),size=11),
            ft.Text('Tempo máximo: sem limite (resultado persistido)' if result.get('timeout')==0 else 'Tempo máximo: '+str(result.get('timeout','não registrado'))+' s',size=11)]
        if result.get('limit') is not None:self.result_panel.controls.append(ft.Text('Máximo de compostos: '+str(result['limit']),size=11))
        if result.get('reference_id'):self.result_panel.controls.append(ft.Text('Desenho extraído da referência: '+result['reference_id']+' · correspondência por elementos e ordem das ligações; valências não exigidas.',size=11))
        if result.get('smiles'):self.result_panel.controls.append(ft.TextField(label='SMILES do fragmento',value=result['smiles'],read_only=True,multiline=True,min_lines=1,max_lines=3))

    async def calculate(self,e=None):
        if self.running:return
        self.running=True;self.request+=1;request=self.request
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token or self.compute_button.disabled:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            try:timeout=int(self.timeout.value);limit=int(self.limit.value)
            except (ValueError,TypeError):self.notice.value='Informe tempo e limite como números inteiros.';view.ui.page.update();return
            scope=self.scope.value;mode=view.mode
            if scope=='community':
                mapping=self.explorer.analysis.community_cache.get(view.mode,({},None,''))[0]
                group=self.explorer.community_choice.value
                if group=='all':self.notice.value='Escolha uma comunidade no painel de configurações.';view.ui.page.update();return
                ids=sorted(n for n,c in mapping.items() if str(c)==group)
            elif scope=='filtered':ids=sorted(self.explorer.filtered_graph)
            else:ids=sorted(self.explorer.analysis.universe(scope))
            if not 1<=timeout<=30 or not 2<=limit<=200 or not 2<=len(ids)<=limit:
                self.notice.value=f'Conjunto com {len(ids)} compostos. Use 1–30 s e limite 2–200; restrinja o conjunto até o limite. Nenhum composto é amostrado.';view.ui.page.update();return
            properties,_=await view.ui.call(numeric_properties,[view.nodes[n] for n in ids],view.molecular_enrichment)
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            rows=tuple((n,properties[n].get('canonical_smiles') if isinstance(properties[n].get('canonical_smiles'),str) else None) for n in ids);key=(rows,timeout,limit)
            if self.cancelled:self.cancelled.set()
            cancelled=threading.Event();self.cancelled=cancelled
            self.compute_button.disabled=True;self.cancel_button.disabled=False;self.notice.value='Calculando MCS do conjunto explícito…';view.ui.page.update()
            try:
                result=self.cache.get(key)
                if result is None:result=await view.ui.call(common,rows,timeout,limit,cancelled)
                if request!=self.request or not view.active or view.ui.token!=view.token:return
                if len(self.cache)>=8:self.cache.pop(next(iter(self.cache)))
                self.cache[key]=result
                self.result=dict(result,scope=scope+' · '+mode+' · conjunto confirmado de '+str(len(ids))+' compostos',source=(view.molecular_enrichment or {}).get('source',view.data.get('dashboard_context',{})))
                self.show_result();self.notice.value=result.get('notice','')
                if view.selected in view.displayed():await view.select(view.selected,show_popup=False)
            finally:
                if self.cancelled is cancelled and view.active:self.compute_button.disabled=False;self.cancel_button.disabled=True;view.ui.page.update()
        try:await self.view.ui.guard(action)
        finally:
            self.running=False
            if self.view.active and self.view.ui.token==self.view.token:
                self.compute_button.disabled=False;self.cancel_button.disabled=True;self.view.ui.page.update()
                if self.pending_scope:
                    self.pending_scope=False
                    await self.calculate()

    def build(self):
        self.mcs_tile=ft.ExpansionTile(title=ft.Text('MCS · subestrutura compartilhada',size=13),controls=[
            ft.Container(padding=ft.Padding.only(top=16,bottom=12),content=ft.Column([ft.Row([self.scope]),
                ft.Row([self.timeout]),ft.Row([self.limit]),ft.Row([self.compute_button,self.cancel_button],wrap=True),self.notice,self.result_panel],
                spacing=16,horizontal_alignment=ft.CrossAxisAlignment.STRETCH))])
        return card('Selected Compound',ft.Column([self.summary,ft.Divider(),ft.Text('Molecular Scaffold',weight=ft.FontWeight.W_600),
            self.scaffold,self.mcs_tile],spacing=8),{'xs':12,'lg':3})
