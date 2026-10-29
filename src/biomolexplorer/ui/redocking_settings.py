"""Shared preparation and responsive, mutually exclusive redocking engine panels."""
import flet as ft
from .localization import Translator
from .runtime_messages import session_language

VINA={'pH','sizeof_box','exhaustiveness','num_modes'}
DOCK6={'dock6_app_path','density','radius','distance','conformer_search_type','max_orientations','num_scored_conformers','dock6_seed'}


class RedockingSettings:
    def __init__(self,form):
        self.form=form;self.engine=form.field_controls['redocking_engine']
        self.tr=Translator(session_language(form.ui))
        self.engine.on_select=self.changed
        self.panels={}
        self.changed()

    def active_field(self,key):
        return key not in (DOCK6 if self.engine.value=='vina' else VINA)

    def active_template(self,name):
        return not name.startswith('dock6/' if self.engine.value=='vina' else 'vina/')

    def changed(self,e=None):
        for engine,panel in self.panels.items():panel.visible=engine==self.engine.value
        for name,tile in self.form.template_tiles.items():tile.visible=self.active_template(name)
        if e is not None:self.form.ui.page.update()

    def fields(self,controls):
        cells=[]
        for control in controls:
            content=control
            if isinstance(control,(ft.TextField,ft.Dropdown)):
                control.width=None
                col={'xs':12,'md':6}
                if getattr(control,'read_only',False):col=12
                if isinstance(control,ft.TextField):control.helper_max_lines=3
                else:
                    control.expand=True
                    content=ft.Column([ft.Row([control])],spacing=10)
                    if getattr(control,'_redocking_helper',None):content.controls.append(ft.Text(control._redocking_helper,size=12,color='#475569'))
            else:col=12
            cells.append(ft.Container(content,col=col,padding=ft.Padding.only(top=12,bottom=16)))
        return ft.ResponsiveRow(cells,spacing=32,run_spacing=24)

    def layout(self):
        form=self.form;t=self.tr
        def card(title,description,content):
            return ft.Container(padding=24,bgcolor='#FFFFFF',border_radius=16,border=ft.Border.all(1,'#E2E8F0'),
                content=ft.Column([ft.Text(t(title),size=18,weight=ft.FontWeight.W_600),ft.Text(t(description),size=13,color='#64748B'),content],spacing=24))
        preparation=form.preparation_settings
        # Both role panels have exactly the same vertical structure. Optional
        # receptor cofactors live in a separate disclosure beneath the pair.
        for role,title in (('receptor','Preparação comum dos receptores'),('ligand','Preparação comum dos ligantes')):
            group=preparation.groups[role];group.title=ft.Text(t(title))
            charge=preparation.settings[role]['charge_type'];charge.width=None;charge.expand=True
            notice='Estas opções serão aplicadas a todos os receptores selecionados.' if role=='receptor' else 'Estas opções serão aplicadas a todos os ligantes selecionados.'
            group.controls=[ft.Container(ft.Text(t(notice),size=13),height=56),
                *[ft.Container(preparation.settings[role][key],height=56) for key in ('remove_solvent','remove_hydrogens','add_hydrogens','minimize')],
                ft.Container(ft.Row([charge]),padding=ft.Padding.only(top=24,bottom=20))]
            group.controls_padding=24;group.tile_padding=20
        roles=ft.ResponsiveRow([ft.Container(preparation.groups[role],col={'xs':12,'md':6},
            border=ft.Border.all(1,'#E2E8F0'),border_radius=12) for role in ('receptor','ligand')],spacing=24,run_spacing=24)
        preparation.cofactors.visible=bool(preparation.has_cofactors.value)
        cofactors=ft.ExpansionTile(title=ft.Text(t('Cofatores do receptor (opcional)')),
            expanded=bool(preparation.has_cofactors.value),tile_padding=16,controls_padding=24,
            controls=[preparation.has_cofactors,ft.Container(preparation.cofactors,padding=ft.Padding.only(top=16))])
        # Hiding a field must release its parent's space, too.
        cofactor_cell=cofactors.controls[-1];cofactor_cell.visible=preparation.cofactors.visible
        def toggle(e):
            preparation.sync();cofactor_cell.visible=preparation.cofactors.visible;form.ui.page.update()
        preparation.has_cofactors.on_change=toggle
        legacy=ft.Text(t('Este projeto tinha preparações diferentes por par. Ao salvar, a configuração comum exibida substituirá essas diferenças.'),color='#B45309',visible=form.legacy_preparation_conflict)
        helpers={
            'pH':'Aplicado na conversão PDBQT quando a adição de hidrogênios está ativada.',
            'exhaustiveness':'Inteiro positivo. Aumentar este valor amplia a busca e o tempo de processamento.',
            'num_modes':'Máximo de poses exportadas pelo Vina; use um inteiro positivo.',
            'density':'Densidade da superfície SES; padrão: 0,5.',
            'radius':'Raio da sonda em Å; padrão: 1,4.',
            'distance':'Distância máxima em Å entre as esferas e o ligante preparado.',
            'conformer_search_type':'Flexível explora ligações rotacionáveis; rígida mantém a geometria interna.',
            'max_orientations':'Inteiro positivo; padrão: 1000 orientações.',
            'num_scored_conformers':'Inteiro positivo. As poses são selecionadas pelo menor Grid Score.',
            'dock6_seed':'Inteiro positivo fixo para reprodutibilidade; padrão: 73.',
        }
        for key,helper in helpers.items():
            control=form.field_controls[key]
            if isinstance(control,ft.Dropdown):control._redocking_helper=t(helper)
            else:control.helper=t(helper)
        box=form.field_controls['sizeof_box']
        # Avoid placing an entire three-field row inside a half-width grid cell.
        if isinstance(box,ft.Row):
            for field in box.controls:field.width=None;field.helper=t('Dimensão total da caixa em Å.')
            box=ft.ResponsiveRow([ft.Container(c,col={'xs':12,'sm':4},padding=ft.Padding.symmetric(vertical=12)) for c in box.controls],spacing=24,run_spacing=24)
            form.field_controls['sizeof_box']=box
        for tile in form.template_tiles.values():
            if not getattr(tile,'_redocking_fields_wrapped',False):
                tile.controls=[self.fields(tile.controls)]
                tile.controls_padding=24;tile.tile_padding=20
                tile._redocking_fields_wrapped=True
        engine_cards=[]
        keys_by_engine=[('vina',('pH','exhaustiveness','num_modes')),('dock6',('dock6_app_path','conformer_search_type','density','radius','distance','max_orientations','num_scored_conformers','dock6_seed'))]
        for engine,keys in keys_by_engine:
            advanced=[tile for name,tile in form.template_tiles.items() if name.startswith(engine+'/')]
            controls=[self.fields([form.field_controls[k] for k in keys])]
            if engine=='vina':controls.insert(0,ft.Column([ft.Text(t('O centro da caixa é obtido do ligante de referência.'),size=13),box],spacing=16))
            panel=ft.Column([*controls,*advanced],spacing=28,visible=engine==self.engine.value)
            self.panels[engine]=panel;engine_cards.append(panel)
        identify=[form.title,form.enabled]
        if hasattr(form,'processing'):identify.append(form.processing)
        inputs=[editor.control for editor in form.input_editors.values()]
        if hasattr(form,'provided'):inputs.insert(0,form.provided.control)
        return ft.Column([
            card('Identificação e entradas','Escolha as estruturas e o modo de processamento.',ft.Column([self.fields(identify),*inputs],spacing=24)),
            card('Curadoria dos pares','Revise PDB, ligante, resíduo e cadeia. Não é necessário repetir a preparação em cada par.',form.redocking_pairs.control),
            card('Preparação única do conjunto','Chimera prepara as estruturas. Receptores preparados podem ser reutilizados sem novo preparo.',
                ft.Column([self.fields([form.field_controls[k] for k in ('target','prepare_complex')]),legacy,roles,cofactors],spacing=24)),
            card('Método de redocking','Vina utiliza PDBQT e caixa de busca. DOCK6 6.11 utiliza MOL2, superfície, esferas e grades.',ft.Column([self.fields([self.engine]),*engine_cards],spacing=28)),
        ],spacing=28,scroll=ft.ScrollMode.AUTO,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
