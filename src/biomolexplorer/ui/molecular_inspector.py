"""Compact selected-compound card and explicit, bounded scaffold exploration."""
import base64
import threading
import flet as ft
from .processing import processing
from biomolexplorer.molecular_network import numeric_properties
from biomolexplorer.network_scaffolds import individual,common
from biomolexplorer.network_analysis import METRICS
from .feedback import close_dialog,readable_error
from .network_charts import card
from .localization import verbatim
from .molecular_cards import badge,metric,quick_metrics,image_panel,property_sections


class MolecularInspector:
    def __init__(self,explorer):
        self.explorer=explorer;self.view=explorer.view;self.properties={};self.structure={};self.dialog=None;self.relations_panel=None
        self.summary=ft.Column([ft.Text('Selecione um composto no grafo ou no ranking.',size=13)],spacing=8)
        self.scaffold=ft.Column(spacing=8);self.result_panel=ft.Column(spacing=10,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        self.scope=ft.Dropdown(label='Conjunto do MCS',value='mcc',expand=True,
            options=[ft.DropdownOption(key=k,text=v) for k,v in [('full','Full Graph'),('mcc','MCC / LCC · maior componente conectado'),('community','Comunidade selecionada'),('filtered','Subconjunto filtrado')]],on_select=self.change_scope)
        self.timeout=ft.TextField(label='Tempo máximo (s)',value='10',expand=True,on_change=self.change_parameters)
        self.limit=ft.TextField(label='Máximo de compostos',value=str(max(50,len(self.view.nodes))),expand=True,on_change=self.change_parameters)
        self.timeout.tooltip='Limite da busca RDKit: 1–30 segundos. Preparação e desenho podem acrescentar tempo; um resultado parcial não confirma maximalidade.'
        self.limit.tooltip='Inicia com capacidade para todos os compostos deste grafo. Você pode reduzir o limite; nenhum subconjunto é amostrado automaticamente.'
        self.community=ft.Dropdown(label='Comunidade do MCS',visible=False,expand=True,options=[],on_select=self.change_scope,enable_filter=True,enable_search=True,menu_height=280)
        self.community_row=ft.Row([self.community],visible=False)
        self.community_mode=None;self.global_community=self.explorer.community_choice.value
        self.running=False;self.result=None;self.cache={};self.structures={};self.cancelled=None;self.request=0;self.pending_scope=False
        self.compute_button=ft.TextButton('Calcular MCS do conjunto',on_click=self.calculate)
        self.cancel_button=ft.TextButton('Interromper MCS',on_click=self.cancel,disabled=True)
        self.notice=ft.Text(size=12,color='#64748B')
        self.progress=ft.Container(visible=False)
        self.context=ft.Text(size=12,color='#334155')
        self.scope_ids=None
        self.describe_scope()
        self.persisted_result=None;self.persisted_rows=None
        saved=self.view.data.get('fragment')
        if saved and saved.get('smarts') and set(self.view.data.get('mcc',[]))==set(self.explorer.analysis.lcc) and saved.get('compounds')==len(self.explorer.analysis.lcc):
            self.result=dict(saved,ids=sorted(self.explorer.analysis.lcc),scope='LCC · resultado persistido',notice=('Tempo limite: maximalidade não confirmada. ' if saved.get('status')=='partial' else '')+'Resultado reutilizado do pipeline; os critérios originais estão preservados.')
            if not self.result.get('smiles'):
                self.result['smiles']=next((c.value for c in self.view.fragment_panel.controls if isinstance(c,ft.TextField) and c.label=='SMILES do fragmento'),'')
            self.persisted_result=dict(self.result)
            self.persisted_rows=self.molecular_rows(sorted(self.explorer.analysis.lcc))
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
        membership='MCC / LCC' if identifier in self.explorer.analysis.lcc else 'Outro componente'
        tags=[badge('Grupo',membership,system=True),badge('Universo','MCC / LCC' if view.mode=='mcc' else 'Full Graph')]
        for key in ('BBB','HIA'):
            if properties.get(key):tags.append(badge(key,properties[key]))
        network=[metric('Conexões no universo',graph.degree(identifier) if identifier in graph else '—'),
            metric('Comunidade',community if community is not None else '—')]
        if self.explorer.metric!='degree':
            network.append(metric(METRICS[self.explorer.metric].split(' · ')[0],f'{values[identifier]:.5g}' if identifier in values else '—'))
        self.summary.controls=[ft.Row([ft.Icon(ft.Icons.SCIENCE_OUTLINED,color='#2552E8',size=22),
            verbatim(ft.Text(identifier,size=19,weight=ft.FontWeight.W_700,selectable=True,expand=True))],spacing=8),
            *([verbatim(ft.Text(str(properties['pref_name']),size=13,color='#475569',selectable=True))] if properties.get('pref_name') else []),
            ft.Row(tags,wrap=True,spacing=6,run_spacing=6),
            image_panel(self.structure['image'],190,self.structure['notice']),quick_metrics(properties),
            ft.ResponsiveRow(network,spacing=8,run_spacing=8),
            ft.Row([ft.TextButton('Localizar',icon=ft.Icons.MY_LOCATION,on_click=locate),
                ft.TextButton('Vizinhança',icon=ft.Icons.HUB_OUTLINED,on_click=self.explorer.ego_selected)],wrap=True),
            ft.FilledTonalButton('View Molecular Details',icon=ft.Icons.OPEN_IN_NEW,on_click=self.details)]
        self.scaffold.controls=[verbatim(ft.Text('Bemis–Murcko · '+identifier,size=12)),
            ft.Text('Scaffold apenas deste composto. Selecionar outro nó atualiza esta estrutura; o MCS do conjunto permanece separado.',size=12),
            ft.Image(src=self.structure['scaffold'],height=125,fit=ft.BoxFit.CONTAIN) if self.structure['scaffold'] else ft.Text(self.structure['notice'],size=12),
            ft.Text(f'{self.structure["atoms"]} átomos · núcleo de anéis e conectores; não implica fragmento compartilhado.',size=11)]
        if self.result:
            self.scaffold.controls.append(ft.Text('Este composto pertence ao conjunto do MCS exibido.' if identifier in self.result['ids'] else 'Este composto está fora do conjunto do MCS exibido.',size=11))
        if smarts:self.scaffold.controls.append(ft.Text(f'MCS destacado: {self.structure.get("highlighted_atoms",0)} átomos no composto.',size=11))

    async def details(self,e=None):
        view=self.view
        identifier,version=view.selected,view.selection_version
        properties=dict(self.properties);structure=dict(self.structure)
        def current():
            return (view.active and view.ui.token==view.token and view.selected==identifier
                and view.selection_version==version and identifier is not None)
        async def action():
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if not current():return
            page=view.ui.page
            width=max(260,min(1120,(getattr(page,'width',None) or 1200)-64))
            height=max(240,min(820,(getattr(page,'height',None) or 1000)-160))
            async def open_3d(e):await view.ui.guard(lambda:view.show_3d(identifier,version))
            molecule=ft.Column([image_panel(structure.get('image'),260,structure.get('notice','Estrutura indisponível')),
                ft.FilledTonalButton('Visualizar estrutura 3D',icon=ft.Icons.VIEW_IN_AR,on_click=open_3d,disabled=not bool(structure.get('image')))],spacing=12)
            left=[card('Estrutura molecular · 2D',molecule,12)]
            if structure.get('scaffold'):
                left.append(card('Scaffold individual · Bemis–Murcko',ft.Column([
                    image_panel(structure['scaffold'],170,framed=False),ft.Text('Scaffold individual do composto selecionado.',size=12,color='#64748B')],spacing=8),12))
            result=self.result
            if result and result.get('smarts') and identifier in result['ids']:
                shared=[ft.Text(result['scope'],size=12),ft.Row([badge('Compostos',result['compounds']),badge('Átomos',result['atoms']),badge('Status',result['status'])],wrap=True)]
                if result.get('image'):shared.append(image_panel(base64.b64decode(result['image']),160,framed=False))
                for key in ('smiles','smarts'):
                    if result.get(key):shared.append(ft.TextField(label=key.upper()+' · fragmento',value=result[key],read_only=True,multiline=True,min_lines=1,max_lines=3))
                left.append(card('MCS · scaffold comum ao conjunto',ft.Column(shared,spacing=10),12))
            graph=self.explorer.analysis.universe(view.mode)
            component=next((i for i,g in view.components.items() if identifier in g),None)
            mapping=self.explorer.analysis.community_cache.get(view.mode,({},None,''))[0]
            network=[metric('Conexões no universo',graph.degree(identifier) if identifier in graph else '—'),
                metric('Componente',component if component is not None else '—'),
                metric('Comunidade',mapping.get(identifier,'—'))]
            for (mode,key),(values,notice) in self.explorer.analysis.cache.items():
                if mode==view.mode and identifier in values and key!='degree':
                    item=metric(METRICS[key].split(' · ')[0],f'{values[identifier]:.6g}');item.tooltip=notice;network.append(item)
            for key in ('degree','degree_density','component'):
                if properties.get(key) is not None:network.append(metric(key+' · Full Graph',properties[key]))
            relation_index=next((i for i,c in enumerate(view.details.controls) if isinstance(c,ft.Text) and str(c.value).startswith('Relações deste composto')),len(view.details.controls))
            self.relations_panel=ft.Column(view.details.controls[relation_index:],spacing=6)
            source=(view.molecular_enrichment or {}).get('source',dict(view.data.get('dashboard_context',{}),origin=view.data.get('dataset_origin',{}).get('origin')))
            provenance=ft.Container(padding=14,bgcolor='#EFF6FF',border_radius=12,content=ft.Row([
                ft.Icon(ft.Icons.SOURCE_OUTLINED,color='#2552E8'),ft.Column([ft.Text('Origem dos dados',size=11,color='#64748B'),
                    verbatim(ft.Text(str(source.get('origin') or source.get('stage_name') or source.get('stage_id') or 'grafo persistido'),size=13,selectable=True)),
                    ft.Text('Execução: '+str(source.get('run_id','não registrada')),size=11,color='#64748B')],expand=True,spacing=4)],spacing=12))
            right=ft.Column([quick_metrics(properties,col={'xs':6,'md':4}),*property_sections(properties,view.nodes[identifier]['properties'])],spacing=14,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
            body=ft.Column([provenance,ft.ResponsiveRow([
                ft.Column(left,col={'xs':12,'lg':5},spacing=14),ft.Column([right],col={'xs':12,'lg':7})],spacing=18,run_spacing=18),
                card('Posição na rede',ft.Column([badge('Universo','MCC / LCC' if view.mode=='mcc' else 'Full Graph'),ft.ResponsiveRow(network,spacing=8,run_spacing=8)],spacing=12),12),
                card('Compostos relacionados',self.relations_panel,12)],scroll=ft.ScrollMode.AUTO,spacing=18)
            def dismissed(e=None):
                if self.dialog is dialog:self.dialog=None;self.relations_panel=None
            def close(e=None):dismissed();close_dialog(page,dialog)
            header=ft.Row([ft.Container(ft.Icon(ft.Icons.SCIENCE_OUTLINED,color='#2552E8',size=28),padding=12,bgcolor='#EFF6FF',border_radius=12),
                ft.Column([verbatim(ft.Text(identifier,size=22,weight=ft.FontWeight.W_700,selectable=True)),ft.Text('Detalhes moleculares · 2D',size=12,color='#64748B'),
                    *([verbatim(ft.Text(str(properties['pref_name']),size=13,color='#64748B'))] if properties.get('pref_name') else [])],expand=True,spacing=4)],spacing=14)
            dialog=ft.AlertDialog(on_dismiss=dismissed,title=header,bgcolor='#F8FAFC',content=ft.Container(width=width,height=height,content=body),
                actions=[ft.TextButton('Fechar molécula',on_click=close)])
            self.dialog=dialog;page.show_dialog(dialog)
        await view.ui.guard(action)

    def resize_details(self):
        if self.dialog is None:return
        page=self.view.ui.page
        self.dialog.content.width=max(260,min(1120,(getattr(page,'width',None) or 1200)-64))
        self.dialog.content.height=max(240,min(820,(getattr(page,'height',None) or 1000)-160))

    def refresh_relations(self,controls):
        if self.relations_panel is None:return
        start=next((i for i,c in enumerate(controls) if isinstance(c,ft.Text) and str(c.value).startswith('Relações deste composto')),len(controls))
        self.relations_panel.controls=controls[start:]

    def close_details(self):
        if self.dialog:close_dialog(self.view.ui.page,self.dialog)
        self.dialog=None;self.relations_panel=None

    def selected_set(self):
        scope=self.scope.value
        if scope=='community':
            mapping=self.explorer.analysis.community_cache.get(self.view.mode,({},None,''))[0]
            group=self.community.value or self.explorer.community_choice.value
            ids=sorted(n for n,c in mapping.items() if str(c)==group)
            label='Comunidade '+str(group)+' · '+('MCC / LCC' if self.view.mode=='mcc' else 'Full Graph')
        elif scope=='filtered':
            ids=sorted(self.explorer.graph())
            label='Subconjunto filtrado · '+('MCC / LCC' if self.view.mode=='mcc' else 'Full Graph')
        else:
            ids=sorted(self.explorer.analysis.universe(scope))
            label='MCC / LCC · maior componente conectado' if scope=='mcc' else 'Full Graph · todos os compostos'
        return ids,label

    def describe_scope(self):
        self.community.visible=self.community_row.visible=self.scope.value=='community'
        ids,label=self.selected_set()
        self.context.value=label+'\nConjunto do MCS: '+str(len(ids))+' compostos'
        return ids,label

    async def sync_filters(self):
        if self.scope.value not in ('community','filtered'):return
        if self.scope.value=='community':
            chosen=self.explorer.community_choice.value
            if chosen!=self.global_community and chosen!='all':self.community.value=chosen
            self.global_community=chosen
        ids,_=self.describe_scope()
        if self.scope_ids is not None and ids!=self.scope_ids:
            await self.change_scope()

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

    def molecular_rows(self,ids):
        matches=(self.view.molecular_enrichment or {}).get('matches',{})
        rows=[]
        for identifier in ids:
            smiles=self.view.nodes[identifier]['properties'].get('canonical_smiles')
            if smiles is None:smiles=matches.get(identifier,{}).get('properties',{}).get('canonical_smiles')
            rows.append((identifier,smiles if isinstance(smiles,str) else None))
        return tuple(rows)

    async def prepare_community(self):
        if self.scope.value!='community':return
        mode=self.view.mode
        mapping,_,_=await self.view.ui.call(self.explorer.analysis.communities,mode)
        if not self.view.active or self.view.ui.token!=self.view.token or mode!=self.view.mode or self.scope.value!='community':return
        from collections import Counter
        counts=Counter(mapping.values())
        groups={str(group):group for group in counts}
        self.community.options=[ft.DropdownOption(key=str(group),text=f'C{group} · {count} compostos')
            for group,count in sorted(counts.items(),key=lambda item:(-item[1],str(item[0])))]
        if self.community.value not in groups or self.community_mode!=mode:
            chosen=self.explorer.community_choice.value
            if chosen in groups:self.community.value=chosen
            elif self.view.selected in mapping:self.community.value=str(mapping[self.view.selected])
            else:self.community.value=str(max(counts,key=lambda group:(counts[group],str(group)))) if counts else None
        self.community_mode=mode;self.global_community=self.explorer.community_choice.value

    async def follow_mode(self):
        self.scope.value=self.view.mode
        await self.change_scope()

    async def change_scope(self,e=None):
        view=self.view
        if not view.active or view.ui.token!=view.token:return
        self.invalidate();self.pending_scope=False
        request=self.request
        self.community.visible=self.scope.value=='community'
        self.describe_scope()
        async def action():
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            await self.prepare_community()
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            ids,_=self.describe_scope()
            # Both graph controls and the MCS dropdown restore the same exact
            # saved MCC. Its original criteria remain visible, including for
            # sets larger than the interactive search limit.
            if (self.scope.value=='mcc' and self.persisted_result
                    and ids==self.persisted_result['ids']
                    and self.molecular_rows(ids)==self.persisted_rows):
                self.result=dict(self.persisted_result);self.scope_ids=ids
                self.show_result();self.notice.value=self.result.get('notice','')
                if hasattr(self,'mcs_tile'):self.mcs_tile.expanded=True
                await self.refresh_selected()
            else:
                self.pending_scope=self.running
                await self.refresh_selected()
                if not self.running and request==self.request:await self.calculate()
            view.ui.page.update()
        await view.ui.guard(action)

    async def change_parameters(self,e=None):
        if not self.view.active or self.view.ui.token!=self.view.token:return
        self.invalidate();self.pending_scope=False
        try:
            timeout,limit=int(self.timeout.value),int(self.limit.value)
            if not 1<=timeout<=30 or limit<2:
                self.notice.value='Use 1–30 segundos e limite de pelo menos 2 compostos.'
        except (ValueError,TypeError):self.notice.value='Informe tempo e limite como números inteiros.'
        await self.refresh_selected();self.view.ui.page.update()

    def show_result(self):
        result=self.result
        if not result:return
        status={'complete':'Concluído','partial':'Parcial','unavailable':'Indisponível','cancelled':'Interrompido','no_common_fragment':'Sem fragmento comum'}.get(result['status'],result['status'])
        information=ft.Column([ft.Text(result['scope'],size=12,weight=ft.FontWeight.W_600),
            ft.Text(f'{result["compounds"]} compostos · {result["atoms"]} átomos',size=12),ft.Text(status,size=12),
            ft.Text(result.get('notice',''),size=12)],spacing=8)
        image=image_panel(base64.b64decode(result['image']) if result.get('image') else None,180,framed=False)
        # Keep the shared scaffold legible across the full sidebar width.
        self.result_panel.controls=[image,information,
            ft.Text('Critérios: cobertura 100%',size=11),
            ft.Text('Anéis com anéis: '+('Sim' if result.get('ring_matches_ring_only') else 'Não' if result.get('ring_matches_ring_only') is False else 'não registrado'),size=11),
            ft.Text('Anéis completos: '+('Sim' if result.get('complete_rings_only') else 'Não' if result.get('complete_rings_only') is False else 'não registrado'),size=11),
            ft.Text('Tempo máximo: sem limite (resultado persistido)' if result.get('timeout')==0 else 'Tempo máximo: '+str(result.get('timeout','não registrado'))+' s',size=11)]
        if result.get('limit') is not None:self.result_panel.controls.append(ft.Text('Máximo de compostos: '+str(result['limit']),size=11))
        if result.get('reference_id'):self.result_panel.controls.append(ft.Text('Desenho extraído da referência: '+result['reference_id']+' · correspondência por elementos e ordem das ligações; valências não exigidas.',size=11))
        self.result_panel.controls.append(ft.TextField(label='Compostos utilizados no MCS',value=', '.join(result['ids']),read_only=True,multiline=True,min_lines=1,max_lines=3))
        if result.get('smiles'):self.result_panel.controls.append(ft.TextField(label='SMILES do fragmento',value=result['smiles'],read_only=True,multiline=True,min_lines=1,max_lines=3))

    async def calculate(self,e=None):
        if self.running:return
        if not self.view.active or self.view.ui.token!=self.view.token:return
        self.running=True;self.request+=1;request=self.request
        self.result=None;self.result_panel.controls=[]
        self.compute_button.disabled=True;self.cancel_button.disabled=False
        self.progress.visible=True
        self.notice.value='Preparando o conjunto para calcular o MCS…'
        if hasattr(self,'mcs_tile'):self.mcs_tile.expanded=True
        self.view.ui.page.update()
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            try:timeout=int(self.timeout.value);limit=int(self.limit.value)
            except (ValueError,TypeError):self.notice.value='Informe tempo e limite como números inteiros.';view.ui.page.update();return
            await self.prepare_community()
            if request!=self.request or not view.active or view.ui.token!=view.token:return
            scope=self.scope.value
            ids,label=self.describe_scope();self.scope_ids=ids
            if scope=='community' and not ids:
                self.notice.value='Nenhuma comunidade disponível neste universo. Selecione outro conjunto.';view.ui.page.update();return
            if not 1<=timeout<=30 or limit<2:
                self.notice.value='Use 1–30 segundos e limite de pelo menos 2 compostos.';view.ui.page.update();return
            if not 2<=len(ids)<=limit:
                self.notice.value=f'Conjunto com {len(ids)} compostos; limite configurado: {limit}. Selecione de 2 até o limite de compostos ou aumente o limite para incluir o conjunto completo. Nenhum composto é amostrado.';view.ui.page.update();return
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
                self.result=dict(result,scope=label,source=(view.molecular_enrichment or {}).get('source',view.data.get('dashboard_context',{})))
                self.show_result();self.notice.value=result.get('notice','')
            finally:
                if self.cancelled is cancelled and view.active:self.compute_button.disabled=False;self.cancel_button.disabled=True;view.ui.page.update()
        async def execute():
            try:
                async with processing(self.view.ui,'Calculando fragmento comum…',
                        'Analisando o conjunto de moléculas selecionado. Aguarde um instante.',on_cancel=self.cancel):
                    await action()
            except Exception as error:
                if request==self.request and self.view.active and self.view.ui.token==self.view.token:
                    self.notice.value='Falha ao calcular MCS: '+readable_error(error)
                raise
        try:await self.view.ui.guard(execute)
        finally:
            self.running=False
            self.progress.visible=False
            if self.view.active and self.view.ui.token==self.view.token:
                self.compute_button.disabled=False;self.cancel_button.disabled=True
                await self.refresh_selected()
                self.view.ui.page.update()
                if self.pending_scope:
                    self.pending_scope=False
                    await self.calculate()

    def build(self):
        self.mcs_tile=ft.ExpansionTile(title=ft.Text('MCS · scaffold comum ao conjunto',size=13),expanded=True,controls=[
            ft.Container(padding=ft.Padding.only(top=16,bottom=12),content=ft.Column([self.context,ft.Text('O MCS usa todos os compostos do conjunto indicado. Clicar em um nó altera somente o scaffold individual do composto selecionado.',size=12),ft.Row([self.scope]),self.community_row,
                ft.Row([self.timeout]),ft.Row([self.limit]),ft.Row([self.compute_button,self.cancel_button],wrap=True),self.progress,self.notice,self.result_panel],
                spacing=16,horizontal_alignment=ft.CrossAxisAlignment.STRETCH))])
        return card('Composto selecionado',ft.Column([self.summary,self.scaffold,ft.Divider(),self.mcs_tile],spacing=14),{'xs':12,'lg':3})
