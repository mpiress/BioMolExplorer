"""Editable research hypotheses without asking users to write JSON."""
import flet as ft
from biomolexplorer.admet_config import PROFILES,DESCRIPTORS,PARAMETERS,resolve_config

CHOICES={
 'admet_mode':('Saída para o próximo bloco',{'annotate':'Anotar todos os compostos','filter':'Filtrar pela hipótese'}),
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
        self.catalogs={name:ft.Checkbox(label=name,value=name in config['alert_catalogs'],disabled=not self.writable) for name in ('PAINS','BRENK','NIH')}
        for key,interval in config['property_limits'].items():self.add_limit(key,interval)
        async def add(e):self.add_limit();self.ui.page.update()
        def grid(controls):
            for control in controls:control.expand=True
            return ft.ResponsiveRow([ft.Container(ft.Row([control]),col={'xs':12,'lg':6},padding=ft.Padding.only(top=8,bottom=4))
                for control in controls],spacing=20,run_spacing=12)
        def section(title,content):
            return ft.Container(padding=16,bgcolor='#F8FAFC',border_radius=12,border=ft.Border.all(1,'#E2E8F0'),
                content=ft.Column([ft.Text(title,size=15,weight=ft.FontWeight.W_600,color='#172B4D'),*content],spacing=14))
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
                'Um alerta PAINS, BRENK ou NIH não demonstra toxicidade. No modo descritivo, até alertas configurados para rejeição continuam no relatório.',
                'RDKit não prediz CYP, hERG, AMES, DILI ou clearance. ADMET-AI requer instalação no worker; as saídas ML_ mantêm o endpoint original, sem limiar universal.',
                'O CSV principal segue a opção de saída. A tabela evaluated preserva os compostos válidos avaliados e seus motivos; o relatório ADMET registra métodos, versões e limites.']
            show_help(self.ui,'Métodos e limites ADMET',[ft.Text(p,size=13) for p in paragraphs],paragraphs)
        self.control=ft.Container(data='wide',content=ft.Column([
            section('Objetivo e saída',[grid([self.profile,self.fields['admet_mode']]),self.note,
                self.fields['prediction_backend'],self.readiness]),
            section('Critérios de seleção',[grid([self.fields['bbb_target'],self.fields['hia_target'],self.fields['rule_set'],self.violations]),self.criteria_summary]),
            ft.ExpansionTile(title=ft.Text('Alertas para revisão química'),expanded=bool(config['alert_policy']=='exclude'),controls=[
                ft.Row(list(self.catalogs.values()),wrap=True),self.fields['alert_policy']],controls_padding=16),
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
        row=ft.Container(padding=ft.Padding.symmetric(vertical=8),content=ft.ResponsiveRow([
            ft.Container(ft.Row([prop]),col={'xs':12,'lg':6}),ft.Container(low,col={'xs':5,'lg':2}),ft.Container(high,col={'xs':5,'lg':2}),
            ft.Container(ft.IconButton(icon=ft.Icons.DELETE_OUTLINE,tooltip='Remover intervalo',on_click=remove,disabled=not self.writable),col=2)],spacing=12,run_spacing=12))
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
