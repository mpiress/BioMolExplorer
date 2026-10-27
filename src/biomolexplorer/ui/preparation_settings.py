"""Shared receptor and ligand preparation controls, independent of tool names."""
import copy
import re
import flet as ft
from biomolexplorer.redocking_config import DEFAULTS
from .settings_layout import section


def preparation_values(stage, config, role):
    values=dict(DEFAULTS, **config.get(role, {}))
    if role in config:return values
    values['charge_type']=stage['parameters'].get('charge_type','gas') if role=='ligand' else 'gas'
    templates=stage.get('templates',{})
    names=('prepare_ligand','prepare_better_conform') if role=='ligand' else ('prepare_receptor',)
    legacy=next((templates['chimera/'+n+'.template'] for n in names if 'chimera/'+n+'.template' in templates),None)
    if legacy is not None:
        lines=legacy.splitlines()
        values['add_hydrogens']='addh' in lines
        values['minimize']=any(line.startswith('minimize') for line in lines)
        match=re.search(r'method (gas|am1)',legacy)
        if match:values['charge_type']=match.group(1)
    complex_source=templates.get('chimera/prepare_complex.template')
    if complex_source is not None:
        values['remove_solvent']='delete solvent' in complex_source.splitlines()
        values['remove_hydrogens']='delete element.H' in complex_source.splitlines()
    return values


class PreparationSettings:
    def __init__(self,form,config=None):
        self.form=form;config=copy.deepcopy(config or {})
        self.extra={key:value for key,value in config.items() if key not in ('receptor','ligand','cofactors','preparation_backend','ambertools_path','ambertools_env','vina_charge_policy')}
        self.settings={};self.groups={}
        self.backend=ft.Dropdown(label='Ferramenta de preparação molecular',value='legacy_chimera',
            tooltip='Chimera aplica hidrogênios, cargas gas/am1 e minimização para Vina e DOCK6.',
            options=[ft.DropdownOption(key='legacy_chimera',text='UCSF Chimera')],disabled=True)
        self.ambertools=ft.TextField(label='Instalação AmberTools/Antechamber',value=config.get('ambertools_path',''),disabled=not form.writable)
        self.ambertools_env=ft.TextField(label='Ambiente Conda AmberTools',
            value=config.get('ambertools_env','BioMolExplorer-AmberTools'),disabled=not form.writable)
        self.vina_policy=ft.Dropdown(label='Política de cargas no Vina padrão',value=config.get('vina_charge_policy','fast'),
            tooltip='Rápido: Gasteiger/RDKit para vina/vinardo. Preservar método: respeita a seleção de cargas, incluindo AM1-BCC.',
            options=[ft.DropdownOption(key='fast',text='Rápido · Gasteiger/RDKit'),
                ft.DropdownOption(key='preserve',text='Preservar método selecionado')],disabled=not form.writable)
        self.receptor_inactive=False
        self.has_cofactors=ft.Switch(label='Considerar cofatores como parte do receptor',value=bool(config.get('cofactors')),disabled=not form.writable)
        self.cofactors=ft.TextField(label='Cofatores do receptor (códigos separados por vírgula)',value=', '.join(config.get('cofactors',[])))
        self.has_cofactors.on_change=lambda e:self.changed()
        groups=[]
        for role,title,notice in (
            ('receptor','Preparação do receptor','Separe a cadeia e o ligante de referência do complexo. Estas opções se aplicam ao receptor; os cofatores selecionados permanecem nele.'),
            ('ligand','Preparação e conformação do ligante','Estas opções se aplicam somente ao ligante de referência e à sua conformação. As opções do receptor são independentes.')):
            values=preparation_values(form.stage,config,role)
            options={key:ft.Switch(label=label,value=values[key],disabled=not form.writable) for key,label in (
                ('remove_solvent','Remover solvente'),('remove_hydrogens','Remover hidrogênios existentes'),
                ('add_hydrogens','Adicionar hidrogênios'),('minimize','Minimizar energia'))}
            options['charge_type']=ft.Dropdown(label='Método de cargas do receptor' if role=='receptor' else 'Método de cargas do ligante',
                tooltip='gas: cargas Gasteiger, geralmente mais rápidas. am1: cargas AM1-BCC, com custo maior de preparação. A minimização também pode aumentar esse tempo.',
                value={'gasteiger':'gas','am1bcc':'am1'}.get(values['charge_type'],values['charge_type']),options=[ft.DropdownOption(key=v,text=label) for v,label in [('gas','Gasteiger'),('am1','AM1-BCC')]],disabled=not form.writable)
            self.settings[role]=options
            if role=='ligand' and form.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6'):
                notice='Prepare os compostos selecionados de ChEMBL, PubChem, ZINC ou arquivos próprios. Todas as fontes usam estas opções.'
            controls=[ft.Text(notice,size=12)]
            if role=='ligand':
                controls += [self.backend,ft.Text('Chimera prepara receptores e ligantes. Open Babel converte os formatos para Vina e DOCK6, sem substituir o método de cargas selecionado.',size=12)]
            if role=='receptor':controls.append(ft.Text('Resíduos proteicos padrão usam templates AMBER. O método de cargas selecionado controla cofatores e resíduos não padrão, não a proteína inteira.',size=12))
            if role=='receptor':controls += [self.has_cofactors,self.cofactors]
            controls += list(options.values())
            groups.append(section(title,controls,expanded=True))
            self.groups[role]=groups[-1]
        self.control=ft.Column(groups,spacing=20,data='input')
        self.sync()

    def sync(self,inactive=False):
        inactive=inactive or not self.form.writable
        self.backend.disabled=True;self.ambertools.disabled=inactive;self.ambertools_env.disabled=inactive;self.vina_policy.disabled=inactive
        receptor_inactive=inactive or self.receptor_inactive
        self.groups['receptor'].visible=not self.receptor_inactive
        self.has_cofactors.disabled=receptor_inactive
        self.cofactors.disabled=receptor_inactive or not self.has_cofactors.value
        for role,options in self.settings.items():
            for control in options.values():control.disabled=receptor_inactive if role=='receptor' else inactive

    def changed(self):
        self.sync();self.form.ui.page.update()

    def read(self):
        result={role:{key:control.value for key,control in options.items()} for role,options in self.settings.items()}
        result['cofactors']=[v.strip().upper() for v in self.cofactors.value.split(',') if v.strip()] if self.has_cofactors.value else []
        result.update(self.extra)
        result['preparation_backend']=self.backend.value
        result['preparation_backend']='legacy_chimera'
        if self.has_cofactors.value and not result['cofactors'] and not self.receptor_inactive:raise ValueError('Informe quais cofatores devem ser mantidos no receptor.')
        return result
