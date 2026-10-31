"""Editable research hypotheses without asking users to write JSON."""
import flet as ft
from biomolexplorer.admet_config import ALERT_CATALOGS,PROFILES,DESCRIPTORS,PARAMETERS,resolve_config
from .settings_layout import grid,section

CHOICES={
 'admet_mode':('Como marcar a hipótese de pesquisa',{'annotate':'Anotar a hipótese','filter':'Marcar seleção pela hipótese'}),
 'prediction_backend':('Método de avaliação',{'rdkit':'RDKit local','admet_ai':'RDKit + ADMET-AI'}),
 'bbb_target':('Hipótese de penetração cerebral',{'any':'Não restringir BBB','BBB+':'Priorizar BBB+','BBB-':'Priorizar BBB−'}),
 'hia_target':('Hipótese de absorção intestinal',{'any':'Não restringir HIA','HIA+':'Priorizar HIA+','HIA-':'Priorizar HIA−'}),
 'rule_set':('Regras físico-químicas',{'none':'Não impor regra','lipinski':'Lipinski · regra de cinco','veber':'Veber · polaridade e flexibilidade','oral':'Lipinski + Veber','fragment':'Fragmentos · regra de três ampliada'}),
 'alert_policy':('Como tratar os alertas estruturais?',{'annotate':'Anotar para revisão','exclude':'Rejeitar na seleção da hipótese'}),
}
ML_PROPERTIES=['ML_HIA_Hou','ML_BBB_Martins','ML_hERG','ML_AMES','ML_DILI','ML_CYP3A4_Veith',
    'ML_Pgp_Broccatelli','ML_Solubility_AqSolDB','ML_PPBR_AZ','ML_Clearance_Hepatocyte_AZ','ML_Half_Life_Obach']


class AdmetSettings:
    def __init__(self,form):
        self.form=form;self.ui=form.ui;self.writable=form.writable
        config=resolve_config(form.stage['parameters']);self.fields={};self.limits=[];self.limit_box=ft.Column(spacing=12)
        self.profile=ft.Dropdown(label='Objetivo da pesquisa',value=config['admet_profile'],disabled=not self.writable,
            options=[ft.DropdownOption(key=k,text=v['label']) for k,v in PROFILES.items()],on_select=self.apply_profile)
        self.note=ft.Text(PROFILES[config['admet_profile']]['description'],size=13,color='#64748B')
        for key,(label,choices) in CHOICES.items():
            self.fields[key]=ft.Dropdown(label=label,value=config[key],disabled=not self.writable,
                options=[ft.DropdownOption(key=k,text=v) for k,v in choices.items()])
        self.violations=ft.TextField(label='Violações permitidas por regra (0–6)',value=str(config['max_rule_violations']),disabled=not self.writable)
        self.catalogs={name:ft.Checkbox(label=name,value=name in config['alert_catalogs'],disabled=not self.writable)
            for name in ALERT_CATALOGS}
        alert_controls=[]
        for name,checkbox in self.catalogs.items():
            alert_controls.append(ft.Container(col={'xs':12,'md':6},content=ft.Row([
                checkbox,ft.IconButton(icon=ft.Icons.INFO_OUTLINE,icon_size=18,
                    tooltip=ALERT_CATALOGS[name],on_click=lambda e,catalog=name:self.show_catalog(catalog))],spacing=4)))
        for key,interval in config['property_limits'].items():self.add_limit(key,interval)
        async def add(e):self.add_limit();self.ui.page.update()
        from biomolexplorer.admet_predictors import predictor_available
        self.available=predictor_available()
        self.readiness=ft.Text('ADMET-AI disponível para execução local.' if self.available else
            'ADMET-AI indisponível: utilize RDKit local ou instale o ambiente opcional antes de executar.',size=12,
            color='#15803D' if self.available else '#B45309')
        self.fields['prediction_backend'].on_select=self.backend_changed
        for key in ('admet_mode','bbb_target','hia_target','rule_set','alert_policy'):
            self.fields[key].on_select=self.update_summary
        self.criteria_summary=ft.Text(size=12,color='#475569')
        self.update_summary()
        def help(e):
            from .action_help import show_help
            paragraphs=[
                'BBB/HIA são heurísticas locais de permeação passiva. No perfil oral, Veber exige zero violações; o limite informado aplica-se a Lipinski.',
                'Alertas estruturais não demonstram toxicidade. No modo descritivo, até alertas configurados para rejeição continuam no relatório.',
                'RDKit não prediz CYP, hERG, AMES, DILI ou clearance. ADMET-AI requer instalação no worker; as saídas ML_ mantêm o endpoint original, sem limiar universal.',
                'Cada entrada gera BBB+, BBB−, HIA+ e HIA− independentemente dos demais critérios. A regra selecionada gera um arquivo com seu nome; os alertas geram uma única tabela true/false. A marca ADMET_selected preserva a hipótese combinada e seus motivos.']
            show_help(self.ui,'Métodos e limites ADMET',[ft.Text(p,size=13) for p in paragraphs],paragraphs)
        self.control=ft.Container(data='wide',content=ft.Column([
            section('Objetivo e saída',[grid([self.profile,self.fields['admet_mode']]),self.note,
                self.fields['prediction_backend'],self.readiness,ft.Text('Saídas por entrada: BBB+, BBB−, HIA+, HIA−; um arquivo da regra escolhida e uma tabela única de alertas com colunas por catálogo. As classificações são independentes e preservam a origem.',size=12,color='#475569')]),
            section('Critérios de seleção',[grid([self.fields['bbb_target'],self.fields['hia_target'],self.fields['rule_set'],self.violations]),self.criteria_summary]),
            ft.ExpansionTile(title=ft.Text('Alertas para revisão química'),expanded=bool(config['alert_policy']=='exclude'),controls=[
                ft.ResponsiveRow(alert_controls,spacing=12,run_spacing=8),self.fields['alert_policy']],controls_padding=16),
            ft.ExpansionTile(title=ft.Text('Intervalos de propriedades'),expanded=bool(self.limits),controls=[
                ft.Text('Mínimo ou máximo vazio não restringe esse lado. Todos os critérios são combinados (E). Use os nomes e as unidades do modelo para propriedades ML_.',size=12,color='#64748B'),
                self.limit_box,ft.TextButton('Adicionar intervalo',icon=ft.Icons.ADD,on_click=add,disabled=not self.writable)],controls_padding=16),
            ft.TextButton('Métodos e limites ADMET',icon=ft.Icons.INFO_OUTLINE,on_click=help)],spacing=18))
        form.field_controls.update(self.fields)
        for key in PARAMETERS:form.readers[key]=lambda key=key:self.read()[key]

    def add_limit(self,key='MW',interval=(None,None)):
        prop=ft.Dropdown(label='Propriedade',value=key,editable=True,enable_filter=True,expand=True,disabled=not self.writable,
            options=[ft.DropdownOption(key=k,text=v) for k,v in DESCRIPTORS.items()]+[ft.DropdownOption(key=k,text=k) for k in ML_PROPERTIES])
        low=ft.TextField(label='Mínimo',value='' if interval[0] is None else str(interval[0]),disabled=not self.writable)
        high=ft.TextField(label='Máximo',value='' if interval[1] is None else str(interval[1]),disabled=not self.writable)
        record=(prop,low,high)
        def remove(e):self.limits.remove(record);self.limit_box.controls.remove(row);self.ui.page.update()
        prop.data='wide'
        row=grid([prop,low,high,ft.IconButton(icon=ft.Icons.DELETE_OUTLINE,
            tooltip='Remover intervalo',on_click=remove,disabled=not self.writable)])
        self.limits.append(record);self.limit_box.controls.append(row)

    def apply_profile(self,e):
        profile=self.profile.value;config=resolve_config(dict(admet_profile=profile))
        self.note.value=PROFILES[profile]['description']
        # Applying a preset intentionally resets its criteria, while retaining model choice and custom ranges.
        for key in ('bbb_target','hia_target','rule_set','alert_policy'):self.fields[key].value=config[key]
        self.violations.value=str(config['max_rule_violations'])
        self.fields['admet_mode'].value='annotate' if profile in ('exploratory','natural_products','custom') else 'filter'
        self.update_summary();self.ui.page.update()

    def update_summary(self,e=None):
        chosen=[self.fields[k].value for k in ('bbb_target','hia_target','rule_set') if self.fields[k].value not in ('any','none')]
        self.criteria_summary.value='Critérios ativos: '+(' · '.join(chosen) if chosen else 'Sem restrição de barreiras ou regras.')
        self.violations.disabled=not self.writable or self.fields['rule_set'].value=='none'
        if e is not None:self.ui.page.update()

    def show_catalog(self,name):
        from .action_help import show_help
        description=ALERT_CATALOGS[name]
        show_help(self.ui,name,[ft.Text(description,size=14)],[description])

    def backend_changed(self,e):
        if self.fields['prediction_backend'].value=='admet_ai' and not self.available:
            self.readiness.color='#B45309'
        self.ui.page.update()

    def read(self):
        limits={}
        for prop,low,high in self.limits:
            key=prop.value
            if key in limits:raise ValueError('Defina apenas um intervalo para cada propriedade ADMET.')
            limits[key]=[float(low.value) if low.value.strip() else None,float(high.value) if high.value.strip() else None]
        return resolve_config(dict(admet_profile=self.profile.value,
            **{key:c.value for key,c in self.fields.items()},max_rule_violations=int(self.violations.value),
            alert_catalogs=[key for key,c in self.catalogs.items() if c.value],property_limits=limits))
