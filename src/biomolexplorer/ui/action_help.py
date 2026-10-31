"""Short contextual help, sized for its text and the available viewport."""
import math
import flet as ft
from .feedback import close_dialog

ACTION_HELP = {
    'order': ('Ordem de execução', 'Organize a prioridade dos blocos respeitando as dependências.', 'Mover um bloco leva seus dependentes. A ordem é salva para as próximas execuções.'),
    'projects': ('Meus projetos', 'Volte à lista dos projetos acessíveis à sua conta.',
        'Busque por nome ou descrição, consulte projetos arquivados e aceite ou recuse convites. Abra um projeto para continuar seu estudo.'),
    'import': ('Importar projeto', 'Continue um estudo exportado do BioMolExplorer.',
        'Selecione o ZIP do projeto e a pasta de destino. A importação recupera configurações, entradas, resultados e histórico compatíveis; não executa o pipeline.'),
    'new': ('Novo projeto', 'Crie um espaço independente para sua pesquisa.',
        'Escolha nome, descrição, cor e pasta principal. A aplicação cria a subpasta do projeto. Um destino existente exige confirmação antes da substituição.'),
    'workspace': ('Workspace', 'Volte à lista de projetos sem encerrar sua sessão.',
        'As alterações do pipeline são salvas automaticamente. Use Salvar para confirmar um rascunho antes de sair quando necessário.'),
    'project': ('Projeto', 'Edite os dados de apresentação deste projeto.',
        'Altere nome, descrição e cor. Essa ação não altera os parâmetros científicos dos blocos. Disponível para proprietários e editores.'),
    'share': ('Compartilhar', 'Gerencie quem pode acessar o projeto.',
        'Convide uma conta cadastrada como Editor ou Leitor e revogue acessos. O destinatário precisa aceitar o convite. Apenas o proprietário gerencia o compartilhamento.'),
    'history': ('Histórico', 'Consulte a autoria e as alterações do projeto.',
        'Explore versões anteriores da configuração, dos arquivos e dos resultados. O proprietário pode restaurar uma versão; a restauração exige confirmação.'),
    'save': ('Salvar', 'Confirme as alterações atuais do pipeline.',
        'Posições, configurações e conexões também são salvas automaticamente após uma breve pausa. Salvar não inicia cálculos. Disponível para proprietários e editores.'),
    'templates': ('Modelos prontos', 'Use um pipeline de exemplo como ponto de partida.',
        'Escolha um modelo para configurar o fluxo. O modelo adiciona etapas ao pipeline atual; revise as conexões e os parâmetros depois de aplicá-lo. Disponível para proprietários e editores.'),
    'run': ('Executar pipeline', 'Execute as etapas ativas na ordem definida pelas conexões.',
        'A configuração é salva antes de iniciar. Se houver resultados compatíveis, escolha se deseja reaproveitá-los. Entradas ainda não configuradas podem pausar o fluxo para seleção. Disponível para proprietários e editores.'),
    'selected': ('Executar até seleção', 'Execute o bloco selecionado e as dependências necessárias para ele.',
        'Selecione um bloco no canvas ou no menu Bloco selecionado. Etapas posteriores não fazem parte dessa execução. A escolha de reaproveitamento também se aplica. Disponível para proprietários e editores.'),
    'canvas': ('Controles do canvas', 'Organize e explore o pipeline sem alterar os cálculos.',
        'Arraste o cabeçalho para mover um bloco; duplo clique abre a configuração. Organizar blocos distribui as etapas pelas dependências. Duplicar cria uma cópia em espaço livre. Use zoom e Redefinir visão para navegar.',
        'Ctrl+Z desfaz, Ctrl+Y refaz e Delete remove o bloco selecionado. Escape cancela uma conexão em andamento. A lateral pode ser ocultada para ampliar a área de desenho.',
        'A câmera exporta todo o pipeline como PNG, incluindo blocos fora da área visível. A imagem preserva as posições e conexões e pode ser salva para artigos e apresentações.'),
}


def help_dimensions(page, paragraphs):
    width=min(580,max(180,(getattr(page,'width',None) or 1200)-64))
    characters=max(18,int((width-16)/7))
    lines=sum(max(1,sum(math.ceil(len(line)/characters) or 1 for line in value.splitlines())) for value in paragraphs)
    wanted=lines*22+max(0,len(paragraphs)-1)*12+12
    available=max(80,(getattr(page,'height',None) or 900)-210)
    return width,min(wanted,available)


def resize_help(page,dialog):
    if dialog is not None and getattr(dialog,'open',False):
        dialog.content.width,dialog.content.height=help_dimensions(page,dialog._help_paragraphs)


def show_help(ui,title,controls,paragraphs):
    translated=[ui.tr(value) if hasattr(ui,'tr') else value for value in paragraphs]
    width,height=help_dimensions(ui.page,translated)
    dialog=ft.AlertDialog(title=ft.Text(title,size=21,weight=ft.FontWeight.W_600),
        content=ft.Container(width=width,height=height,content=ft.Column(controls,tight=True,spacing=12,scroll=ft.ScrollMode.AUTO)),
        inset_padding=ft.Padding.symmetric(horizontal=24,vertical=24),
        actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(ui.page,dialog))])
    dialog._help_paragraphs=translated
    ui.help_dialog=dialog
    ui.page.show_dialog(dialog)
    return dialog


def show_action_help(ui,key):
    title,*paragraphs=ACTION_HELP[key]
    return show_help(ui,title,[ft.Text(value,size=14) for value in paragraphs],paragraphs)


def action_row(ui,key,action):
    action.expand=True
    action.style=ft.ButtonStyle(padding=ft.Padding.symmetric(horizontal=12,vertical=12),
        alignment=ft.Alignment.CENTER_LEFT,shape=ft.RoundedRectangleBorder(radius=10),
        bgcolor='#2563EB' if key=='run' else '#FFFFFF',color='#FFFFFF' if key=='run' else '#172B4D')
    # Permit long translated action names to wrap in the narrow sidebar.
    if isinstance(action.content,str):action.content=ft.Text(action.content,size=13,max_lines=2)
    return ft.Row([action,ft.IconButton(icon=ft.Icons.INFO_OUTLINE,icon_size=18,
        icon_color='#64748B',tooltip='Informações · '+ui.tr(ACTION_HELP[key][0]),
        on_click=lambda e:show_action_help(ui,key))],spacing=2)


def action_group(ui,key,title,icon,actions):
    if not hasattr(ui,'sidebar_expanded'):ui.sidebar_expanded={}
    def changed(e):
        ui.sidebar_expanded[key]=str(e.data).lower()=='true'
    return ft.ExpansionTile(title=ft.Text(title,size=13,weight=ft.FontWeight.W_600,color='#172B4D'),
        leading=ft.Icon(icon,size=18,color='#64748B'),controls=actions,
        expanded=ui.sidebar_expanded.get(key,True),on_change=changed,maintain_state=True,dense=True,
        min_tile_height=44,tile_padding=ft.Padding.symmetric(horizontal=4),
        controls_padding=ft.Padding.only(bottom=12),expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
        shape=ft.RoundedRectangleBorder(radius=10,side=ft.BorderSide(0,ft.Colors.TRANSPARENT)),
        collapsed_shape=ft.RoundedRectangleBorder(radius=10,side=ft.BorderSide(0,ft.Colors.TRANSPARENT)))
