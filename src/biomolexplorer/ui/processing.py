"""One branded splash for foreground work, shared across nested UI actions."""
import asyncio
from contextlib import asynccontextmanager
from contextvars import ContextVar
from functools import wraps

from .configuration_loading import ProcessingSplash

_CURRENT=ContextVar('biomol_processing',default=None)
_BACKGROUND=ContextVar('biomol_background',default=False)

MESSAGES={
    'set_pdb_ligands':('Salvando ligantes PDB…','Validando os ligantes e atualizando os pares de redocking. Aguarde um instante.'),
    'pdb_ligands':('Analisando ligantes PDB…','Identificando os resíduos disponíveis na estrutura. Aguarde um instante.'),
    'save_pipeline':('Salvando pipeline…','Gravando as configurações e o histórico do projeto. Aguarde um instante.'),
    'project':('Carregando projeto…','Lendo as informações e permissões do projeto. Aguarde um instante.'),
    'list_projects':('Carregando projetos…','Preparando a lista de projetos. Aguarde um instante.'),
    'list_runs':('Carregando execuções…','Consultando as etapas e os resultados disponíveis. Aguarde um instante.'),
    'assets':('Carregando arquivos…','Preparando os arquivos disponíveis para seleção. Aguarde um instante.'),
    'export_project':('Exportando projeto…','Preparando o pacote com os arquivos e o histórico. Aguarde um instante.'),
    'import_project':('Importando projeto…','Verificando e organizando o projeto recebido. Aguarde um instante.'),
    'rollback':('Restaurando projeto…','Restaurando os arquivos e as configurações selecionadas. Aguarde um instante.'),
    'history':('Carregando histórico…','Consultando as alterações do projeto. Aguarde um instante.'),
    'submit':('Preparando execução…','Validando o pipeline e organizando os dados de entrada. Aguarde um instante.'),
    'resume':('Continuando execução…','Validando a seleção e preparando o próximo bloco. Aguarde um instante.'),
    'cancel':('Cancelando execução…','Interrompendo os processos e atualizando o estado da execução. Aguarde um instante.'),
    'available_results':('Verificando resultados preservados…','Validando os arquivos concluídos e suas origens. Aguarde um instante.'),
    'existing_results':('Verificando resultados…','Identificando os cálculos que podem ser reaproveitados. Aguarde um instante.'),
    'read_file':('Carregando arquivo…','Lendo os dados solicitados. Aguarde um instante.'),
    'read_bytes':('Preparando download…','Lendo o arquivo para download. Aguarde um instante.'),
    'import_local_file':('Importando arquivo…','Verificando e registrando o arquivo no projeto. Aguarde um instante.'),
    'finish_upload':('Concluindo envio…','Verificando e registrando o arquivo recebido. Aguarde um instante.'),
    'create_project':('Criando projeto…','Preparando o diretório e as configurações do projeto. Aguarde um instante.'),
    'create_project_in_parent':('Criando projeto…','Preparando o diretório e as configurações do projeto. Aguarde um instante.'),
    'update_project':('Salvando projeto…','Gravando as informações do projeto. Aguarde um instante.'),
    'delete_project':('Excluindo projeto…','Removendo os arquivos e registros do projeto. Aguarde um instante.'),
    'remove_project_results':('Excluindo todos os resultados…','Removendo os arquivos calculados e atualizando as etapas. Aguarde um instante.'),
    'remove_all':('Excluindo resultados da etapa…','Removendo os arquivos calculados e suas referências. Aguarde um instante.'),
    'remove':('Removendo dados…','Gravando a alteração e atualizando os resultados. Aguarde um instante.'),
    'remove_asset':('Excluindo arquivo…','Removendo o arquivo e atualizando suas referências. Aguarde um instante.'),
    'page':('Carregando tabela…','Preparando os registros solicitados. Aguarde um instante.'),
    'tables':('Carregando resultados…','Organizando as tabelas disponíveis. Aguarde um instante.'),
    'metadata_records':('Preparando pares de redocking…','Lendo os ligantes, resíduos e cadeias disponíveis. Aguarde um instante.'),
    'pdb_review_statuses':('Carregando revisões PDB…','Organizando as estruturas e os ligantes revisados. Aguarde um instante.'),
    'compound_properties':('Calculando propriedades moleculares…','Analisando a estrutura selecionada. Aguarde um instante.'),
    'molecule_image':('Preparando visualização molecular…','Gerando a representação da molécula. Aguarde um instante.'),
    'individual':('Analisando molécula…','Calculando as propriedades e o scaffold da estrutura. Aguarde um instante.'),
    'common':('Calculando fragmento comum…','Analisando o conjunto de moléculas selecionado. Aguarde um instante.'),
    'layout':('Organizando grafo…','Calculando a disposição dos nós e conexões. Aguarde um instante.'),
    'metrics':('Analisando grafo…','Calculando as métricas da rede. Aguarde um instante.'),
    'communities':('Identificando comunidades…','Analisando os grupos de moléculas do grafo. Aguarde um instante.'),
    'footprint_pdf':('Carregando footprint…','Lendo as energias de interação por resíduo. Aguarde um instante.'),
    'pdf_preview':('Preparando gráfico footprint…','Convertendo o gráfico para visualização. Aguarde um instante.'),
    'login':('Entrando no BioMolExplorer…','Verificando suas credenciais. Aguarde um instante.'),
    'register':('Criando conta…','Registrando seu acesso ao BioMolExplorer. Aguarde um instante.'),
}
DEFAULT=('Processando dados…','O BioMolExplorer está trabalhando. Aguarde um instante.')

for names,message in [
    (('admet_datasets','admet_model'),('Preparando resultados ADMET…','Carregando as propriedades, os critérios e os gráficos ADMET. Aguarde um instante.')),
    (('graph_datasets','graph_enrichment','graph_model','graph_sources','load_view','ResultsViewer'),
        ('Preparando análise de grafos…','Lendo a rede e organizando os gráficos e seus dados. Aguarde um instante.')),
    (('statistics','dashboard','aggregate','comparison'),('Calculando estatísticas…','Analisando os resultados e preparando as comparações. Aguarde um instante.')),
    (('simulations','simulation'),('Carregando resultados de redocking…','Organizando os complexos, as poses e os valores de RMSD. Aguarde um instante.')),
    (('pairs',),('Carregando pares receptor–ligante…','Organizando as simulações e os resultados disponíveis. Aguarde um instante.')),
    (('pdb_view_url','compound_view_url','issue_document','issue_result','scene_payload'),
        ('Preparando visualização 3D…','Organizando a estrutura e os dados da cena molecular. Aguarde um instante.')),
    (('archive','export_csv','csv','svg','image_png','graph_exchange','copyfile'),
        ('Preparando exportação…','Gerando os arquivos solicitados. Aguarde um instante.')),
    (('render_pipeline',),('Gerando imagem do pipeline…','Desenhando os blocos e as conexões. Aguarde um instante.')),
    (('reference_choices',),('Organizando compostos…','Lendo as moléculas disponíveis para seleção. Aguarde um instante.')),
    (('validate_bundle','validate_file'),('Validando arquivos…','Conferindo o formato e o conteúdo dos arquivos selecionados. Aguarde um instante.')),
    (('directories','project_destination'),('Carregando pastas…','Consultando os diretórios disponíveis. Aguarde um instante.')),
    (('mkdir',),('Criando pasta…','Preparando o diretório selecionado. Aguarde um instante.')),
    (('prepare_upload',),('Preparando envio…','Registrando os arquivos que serão enviados ao projeto. Aguarde um instante.')),
    (('members','invitations'),('Carregando compartilhamento…','Consultando os participantes e os convites do projeto. Aguarde um instante.')),
    (('invite','revoke','accept_invitation'),('Atualizando compartilhamento…','Gravando os convites e as permissões do projeto. Aguarde um instante.')),
    (('logout',),('Saindo do BioMolExplorer…','Encerrando sua sessão. Aguarde um instante.')),
    (('get_run',),('Carregando execução…','Consultando o estado das etapas e seus resultados. Aguarde um instante.')),
    (('numeric_properties',),MESSAGES['compound_properties']),
    (('remove_failure',),MESSAGES['remove']),
]:
    for name in names:MESSAGES[name]=message


def splash_for(ui):
    splash=getattr(ui,'_configuration_loading',None)
    if splash is None:
        splash=ProcessingSplash(ui.page)
        ui._configuration_loading=splash
    return splash


class ProcessingScope:
    def __init__(self,ui,title,description,modal,on_cancel=None):
        self.ui=ui;self.title=title;self.description=description;self.modal=modal
        self.request=object();self.shown=False;self.closed=False
        self.message=None
        self.on_cancel=on_cancel

    async def start(self,function=None,args=()):
        if self.closed or _BACKGROUND.get():return
        title,description=MESSAGES.get(getattr(function,'__name__',''),DEFAULT)
        title=self.title or title;description=self.description or description
        cancel=self.on_cancel
        if cancel is None and getattr(function,'__name__','')=='common' and args and hasattr(args[-1],'set'):
            cancel=lambda e:args[-1].set()
        message=(title,description,cancel)
        if not self.shown or message!=self.message:
            splash_for(self.ui).show(self.request,title,description,modal=self.modal,on_cancel=cancel)
            self.shown=True;self.message=message
            # Deliver the splash before validation, disk access or form rebuilds.
            await asyncio.sleep(0)

    def close(self):
        self.closed=True
        if self.shown:splash_for(self.ui).hide(self.request)


@asynccontextmanager
async def processing(ui,title=None,description=None,*,immediate=True,modal=True,on_cancel=None):
    if _BACKGROUND.get():
        yield None
        return
    parent=_CURRENT.get()
    if parent is not None and parent.ui is ui and not parent.closed:
        previous=(parent.title,parent.description,parent.on_cancel)
        if title:parent.title=title;parent.description=description
        if on_cancel:parent.on_cancel=on_cancel
        try:
            if immediate:await parent.start()
            yield parent
        finally:
            parent.title,parent.description,parent.on_cancel=previous
        return
    scope=ProcessingScope(ui,title,description,modal,on_cancel)
    context=_CURRENT.set(scope)
    try:
        if immediate:await scope.start()
        yield scope
    finally:
        scope.close()
        _CURRENT.reset(context)


def background(function):
    """Polling and automatic draft writes must never steal foreground focus."""
    @wraps(function)
    async def wrapped(*args,**kwargs):
        context=_CURRENT.set(None);quiet=_BACKGROUND.set(True)
        try:return await function(*args,**kwargs)
        finally:_BACKGROUND.reset(quiet);_CURRENT.reset(context)
    return wrapped


def with_processing(title,description=None,*,foreground=False):
    def decorate(function):
        @wraps(function)
        async def wrapped(owner,*args,**kwargs):
            ui=getattr(owner,'ui',None) or getattr(getattr(owner,'view',None),'ui',None) or owner
            quiet=_BACKGROUND.set(False) if foreground else None
            try:
                async with processing(ui,title,description,immediate=False):return await function(owner,*args,**kwargs)
            finally:
                if quiet is not None:_BACKGROUND.reset(quiet)
        return wrapped
    return decorate


@asynccontextmanager
async def pause_processing(ui):
    """A native file dialog waits for the user, rather than processing data."""
    scope=_CURRENT.get()
    if scope is not None and scope.ui is ui and not scope.closed and scope.shown:
        splash_for(ui).hide(scope.request)
        scope.shown=False
    quiet=_BACKGROUND.set(True)
    try:yield
    finally:_BACKGROUND.reset(quiet)


async def save_file(ui,*args,**kwargs):
    async with pause_processing(ui):return await ui.picker.save_file(*args,**kwargs)


async def pick_files(ui,*args,**kwargs):
    async with pause_processing(ui):return await ui.picker.pick_files(*args,**kwargs)
