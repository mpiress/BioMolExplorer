"""Packaged WebGL viewer with project-authorized, bounded browser access."""
import html
import hashlib
import json
import secrets
import threading
import time
from collections import OrderedDict
from functools import lru_cache
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit

from .visualizations import MAX_VIEW_BYTES
from .workspace import AccessDenied

RESOURCE=Path(__file__).parent/'resources/viewer'
PREFIX='/molecular-viewer'
ASSETS={'3Dmol-min.js','pdb-viewer.js','pdb-viewer.css','docking-viewer.js'}
PRIVATE_HEADERS={'Cache-Control':'no-store','Referrer-Policy':'no-referrer','X-Content-Type-Options':'nosniff',
    'Content-Security-Policy':"default-src 'none'; script-src 'self'; style-src 'self'; connect-src 'self'; img-src 'self' data: blob:; worker-src 'self' blob:; base-uri 'none'; frame-ancestors 'none'"}
TEXT={
 'pt':{'viewer':'Estrutura 3D','representation':'Representação','cartoon':'Fitas + ligantes','sticks':'Bastões','spheres':'Esferas','lines':'Linhas',
       'color':'Colorir por','chain':'Cadeia','element':'Elemento','spectrum':'Sequência','model':'Modelo','all':'Todas as cadeias',
       'ligands':'Ligantes','water':'Água','center':'Recentrar','fullscreen':'Tela inteira','snapshot':'Salvar imagem',
       'hint':'Arraste para girar · Roda do mouse para zoom · Botão direito para deslocar',
       'loading':'Carregando estrutura…','error':'Não foi possível abrir a estrutura. Volte ao projeto e abra o visualizador novamente.',
       'webgl':'Seu navegador não conseguiu iniciar o WebGL. Ative a aceleração gráfica ou use outro navegador.',
       'atoms':'átomos','chains':'cadeias','picked':'Selecione um átomo para ver seus detalhes.',
       'conformer':'Conformero gerado localmente a partir do SMILES; não representa uma pose de docking.'},
 'en':{'viewer':'3D structure','representation':'Representation','cartoon':'Ribbons + ligands','sticks':'Sticks','spheres':'Spheres','lines':'Lines',
       'color':'Color by','chain':'Chain','element':'Element','spectrum':'Sequence','model':'Model','all':'All chains',
       'ligands':'Ligands','water':'Water','center':'Recenter','fullscreen':'Fullscreen','snapshot':'Save image',
       'hint':'Drag to rotate · Mouse wheel to zoom · Right drag to pan',
       'loading':'Loading structure…','error':'Unable to open the structure. Return to the project and open the viewer again.',
       'webgl':'Your browser could not start WebGL. Enable graphics acceleration or use another browser.',
       'atoms':'atoms','chains':'chains','picked':'Select an atom to see its details.',
       'conformer':'Conformer generated locally from SMILES; this is not a docking pose.'}}

TEXT['pt'].update(receptor_style='Estilo do receptor',ligand_style='Estilo do ligante',hydrogens='Hidrogênios do ligante',
    hydrogen_hint='Exibe apenas os hidrogênios presentes no arquivo.',interaction_view='Interações 2D',close='Voltar ao 3D')
TEXT['en'].update(receptor_style='Receptor style',ligand_style='Ligand style',hydrogens='Ligand hydrogens',
    hydrogen_hint='Shows only hydrogen atoms present in the file.',interaction_view='2D interactions',close='Back to 3D')

@lru_cache(maxsize=32)
def _asset_digest(name,mtime_ns,size):
    """Version cache URLs whenever a packaged asset changes, including in development."""
    return hashlib.sha256((RESOURCE/name).read_bytes()).hexdigest()[:16]


def viewer_document(name,language='pt',prefix=PREFIX,generated=False,scene=False):
    language=language if language in TEXT else 'pt'
    labels=TEXT[language]
    file_format=Path(name).suffix.lower().lstrip('.')
    ligand='.lig.' in name.lower() or file_format=='sdf'
    config=json.dumps({'name':name,'labels':labels,'format':file_format if file_format in ('pdb','pdbqt','mol2','sdf') else 'pdb',
        'ligand':ligand,'generated':generated,'scene':scene,'language':language,
        'representation':'sticks' if ligand or file_format=='mol2' else 'cartoon'},ensure_ascii=False).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    source=(RESOURCE/'pdb-viewer.html').read_text(encoding='utf-8')
    for asset in ASSETS:
        stat=(RESOURCE/asset).stat()
        digest=_asset_digest(asset,stat.st_mtime_ns,stat.st_size)
        source=source.replace('/assets/'+asset+'"','/assets/'+asset+'?v='+digest+'"')
    for key,value in dict(labels,name=name,language=language,prefix=prefix,config=config).items():
        source=source.replace('{{'+key+'}}',value if key=='config' else html.escape(value,quote=True))
    return source


class StructureViewers:
    """Opaque capabilities never contain login tokens or expose project folders."""
    TTL=24*60*60
    LIMIT=4096

    def __init__(self,store,clock=time.monotonic):
        self.store=store;self.clock=clock;self.tickets=OrderedDict();self.lock=threading.RLock()
        self.server=None;self.thread=None

    def issue(self,token,pid,filename,language='pt'):
        self.store.project(token,pid)
        path=self.store.scoped_path(pid,filename)
        if path.suffix.lower() not in ('.pdb','.pdbqt','.mol2','.sdf') or not path.is_file():raise ValueError('Estrutura molecular indisponível.')
        if path.stat().st_size>MAX_VIEW_BYTES:raise ValueError('Visualização maior que o limite de 32 MB. Baixe o arquivo para consultá-lo localmente.')
        return self._issue({'token':token,'pid':pid,'path':str(path),'language':language})

    def issue_compound(self,token,pid,smiles,name,language='pt'):
        """Keep generated conformers ephemeral; authorize every browser fetch."""
        from .visualizations import molecule_sdf
        self.store.project(token,pid)
        data=molecule_sdf(smiles).encode('utf-8')
        if len(data)>MAX_VIEW_BYTES:raise ValueError('Estrutura molecular indisponível.')
        self.store.project(token,pid)
        return self._issue({'token':token,'pid':pid,'name':str(name)+'.sdf',
                            'data':data,'language':language,'generated':True})

    def issue_result(self,token,pid,rid,sid,kind,selection,language='pt'):
        record=dict(token=token,pid=pid,rid=rid,sid=sid,kind=kind,selection=selection,language=language)
        self._result_spec(record)
        return self._issue(record)

    def _result_spec(self,record):
        if record['kind']=='redocking':
            from .redocking_results import RedockingResults
            return RedockingResults(self.store).scene(record['token'],record['pid'],record['rid'],record['sid'],record['selection'])
        if record['kind']=='docking':
            from .docking_results import DockingResults
            return DockingResults(self.store).scene(record['token'],record['pid'],record['rid'],record['sid'],**record['selection'])
        raise AccessDenied('Conformação não autorizada.')

    def issue_document(self,token,pid,rid,sid,selection,language='pt'):
        from .docking_results import DockingResults
        DockingResults(self.store).footprint(token,pid,rid,sid,**selection)
        return self._issue(dict(token=token,pid=pid,rid=rid,sid=sid,selection=selection,language=language,document=True))

    def _issue(self,record):
        with self.lock:
            now=self.clock()
            for key in [k for k,v in self.tickets.items() if v['expires']<=now]:self.tickets.pop(key)
            while len(self.tickets)>=self.LIMIT:self.tickets.popitem(last=False)
            key=secrets.token_urlsafe(32)
            self.tickets[key]=dict(record,expires=now+self.TTL)
            return key

    def ticket(self,key):
        with self.lock:
            record=self.tickets.get(key)
            if record is None or record['expires']<=self.clock():
                self.tickets.pop(key,None);raise AccessDenied('Visualização expirada ou indisponível.')
            record=dict(record)
        self.store.project(record['token'],record['pid'])
        return record

    def response(self,path,prefix=PREFIX):
        """Shared HTTP adapter for same-origin web and loopback desktop hosts."""
        language = 'en'
        try:
            if path.startswith(PREFIX+'/assets/'):
                name=path.removeprefix(PREFIX+'/assets/')
                if name not in ASSETS:return 404,'text/plain',b'Not found',PRIVATE_HEADERS
                mime={'js':'text/javascript','css':'text/css'}[name.rsplit('.',1)[1]]
                return 200,mime,(RESOURCE/name).read_bytes(),{'Cache-Control':'no-cache','X-Content-Type-Options':'nosniff'}
            parts=path.removeprefix(PREFIX+'/').split('/')
            if not path.startswith(PREFIX+'/') or len(parts)>2 or (len(parts)==2 and parts[1] not in ('structure','scene')):
                return 404,'text/plain',b'Not found',PRIVATE_HEADERS
            with self.lock:
                language = self.tickets.get(parts[0], {}).get('language', 'en')
            record=self.ticket(parts[0])
            if record.get('document'):
                from .docking_results import DockingResults
                filename=DockingResults(self.store).footprint(record['token'],record['pid'],record['rid'],record['sid'],**record['selection'])
                data=self.store.read_file(record['token'],record['pid'],filename,MAX_VIEW_BYTES+1)
                if len(data)>MAX_VIEW_BYTES or not data.startswith(b'%PDF-'):raise ValueError('PDF indisponível.')
                return 200,'application/pdf',data,dict(PRIVATE_HEADERS,**{'Content-Disposition':'inline'})
            if 'kind' in record:
                spec=self._result_spec(record)
                if len(parts)==2:
                    if parts[1]!='scene':raise AccessDenied('Conformação não autorizada.')
                    from .docking_scene import scene_payload
                    data=scene_payload(self.store,record['token'],record['pid'],spec,cutoff=8.,include_interactions=True)
                    return 200,'application/json',json.dumps(data,ensure_ascii=False).encode(),PRIVATE_HEADERS
                document=viewer_document(spec['name'],record['language'],prefix,scene=True)
                return 200,'text/html; charset=utf-8',document.encode(),PRIVATE_HEADERS
            if len(parts)==2:
                if parts[1]!='structure':raise AccessDenied('Conformação não autorizada.')
                data=record['data'] if 'data' in record else self.store.read_file(record['token'],record['pid'],record['path'],MAX_VIEW_BYTES+1)
                if len(data)>MAX_VIEW_BYTES:raise ValueError('PDB exceeds preview limit.')
                return 200,'text/plain; charset=utf-8',data,PRIVATE_HEADERS
            document=viewer_document(record['name'] if 'name' in record else Path(record['path']).name,
                                     record['language'],prefix,generated=record.get('generated',False))
            return 200,'text/html; charset=utf-8',document.encode(),PRIVATE_HEADERS
        except (AccessDenied,FileNotFoundError,ValueError,OSError):
            message = ('Visualização expirada ou indisponível. Retorne ao projeto e abra novamente.' if language == 'pt' else
                       'Visualization expired or unavailable. Return to the project and reopen it.')
            return 403,'text/plain; charset=utf-8',message.encode('utf-8'),PRIVATE_HEADERS

    def url(self,key,web=False,page_url=None):
        if web:
            if not page_url:raise ValueError('Endereço da interface web indisponível.')
            parsed=urlsplit(page_url)
            scheme={'ws':'http','wss':'https','http':'http','https':'https'}.get(parsed.scheme)
            if scheme is None or not parsed.netloc:
                raise ValueError('Endereço HTTP/HTTPS da interface web inválido.')
            return urlunsplit((scheme,parsed.netloc,PREFIX+'/'+key,'',''))
        self.start_desktop()
        return f'http://127.0.0.1:{self.server.server_port}{PREFIX}/{key}'

    def start_desktop(self):
        with self.lock:
            if self.server is not None:return
            owner=self
            class Handler(BaseHTTPRequestHandler):
                def do_GET(self):
                    status,mime,data,headers=owner.response(urlsplit(self.path).path)
                    self.send_response(status);self.send_header('Content-Type',mime)
                    self.send_header('Content-Length',str(len(data)))
                    for name,value in headers.items():self.send_header(name,value)
                    self.end_headers()
                    try:self.wfile.write(data)
                    except (BrokenPipeError,ConnectionResetError):pass
                def log_message(self,*args):pass  # Never log capability URLs.
            self.server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
            self.server.daemon_threads=True
            self.thread=threading.Thread(target=self.server.serve_forever,name='biomol-viewer',daemon=True);self.thread.start()

    def close(self):
        with self.lock:
            server=self.server;self.server=None;self.tickets.clear()
        if server:server.shutdown();server.server_close()
        if self.thread:self.thread.join(timeout=2)
