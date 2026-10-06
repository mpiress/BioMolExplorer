"""Build the project's static bilingual documentation using only the standard library.

Supports the Markdown constructs used by these guides: headings, paragraphs,
lists, tables, fenced code, emphasis, links and the architecture Mermaid graph.
Markdown remains the source; do not edit generated HTML pages directly.
"""
import html
import os
from pathlib import Path
import re
import unicodedata

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
GUIDES = ('installation', 'user_manual', 'frontend', 'backend_usage', 'architecture', 'projects', 'pipeline_validation', 'chimera_migration', 'dms_migration', 'retrieval', 'redocking_configuration', 'logging')
LABELS = {
    'en': {'language':'English','other':'Português','home':'Platform overview','docs':'Documentation',
           'source':'Markdown source','contents':'On this page','diagram':'Application architecture',
           'diagram_source':'Mermaid source','intro':'Choose your path into BioMolExplorer.',
           'description':'Guides for using the workspace, running scientific operations and extending the platform.',
           'titles':['Installation and configuration','User manual','Flet workspace','Backend usage','Architecture and technical review','Projects and versions','Pipeline validation','Chimera replacement assessment','Native DMS port and validation','Flexible information retrieval','Redocking configuration','Logs and diagnostics'],
           'summaries':['Download from GitHub, install Chimera 1.17 and DOCK6 6.11, and check your environment before starting.',
                        'Detailed walkthrough from language selection and your first study to results and migration.',
                        'Accounts, private projects, sharing, flexible pipelines and your own files.',
                        'Scientific environment, CLI, job services, operation parameters and outputs.',
                        'Code structure, contracts, execution lifecycle, verification and limitations.',
                        'Project folders, migration, history, rollback and multiple validated inputs.',
                        'Stage contracts, connection matrix, regression tests and verification limits.',
                        'Current Chimera functions, Python alternatives and requirements for a validated migration.',
                        'Native SES generation and comparisons with the original UCSF C algorithm.',
                        'PDB search criteria, direct ChEMBL searches, optional filters and query reports.',
                        'Explicit pairs, shared chain, cofactors, preparation, RMSD tables, downloads and 3D inspection.',
                        'Correlated events, failure codes, per-job summaries and filtered log reports.'],
           'read':'Read guide','start':'First installation? Download the project and configure the required tools before starting.',
           'footer':'Documentation generated from the project Markdown sources. MIT license.'},
    'pt': {'language':'Português','other':'English','home':'Apresentação da plataforma','docs':'Documentação',
           'source':'Fonte Markdown','contents':'Nesta página','diagram':'Arquitetura da aplicação',
           'diagram_source':'Fonte Mermaid','intro':'Escolha seu caminho no BioMolExplorer.',
           'description':'Guias para usar o workspace, executar operações científicas e estender a plataforma.',
           'titles':['Instalação e configuração','Manual do usuário','Workspace Flet','Uso do backend','Arquitetura e revisão técnica','Projetos e versões','Validação do pipeline','Avaliação da substituição do Chimera','Porte nativo DMS e validação','Recuperação flexível da informação','Configuração do redocking','Logs e diagnóstico'],
           'summaries':['Download pelo GitHub, instalação de Chimera 1.17 e DOCK6 6.11 e verificação antes de iniciar.',
                        'Passo a passo do idioma e primeiro estudo aos resultados e à migração.',
                        'Contas, projetos privados, compartilhamento, pipelines flexíveis e arquivos próprios.',
                        'Ambiente científico, CLI, serviços de tarefas, parâmetros e resultados das operações.',
                        'Organização do código, contratos, ciclo de execução, verificações e limites.',
                        'Pastas, migração, histórico, restauração e várias entradas validadas.',
                        'Contratos das etapas, matriz de conexões, regressões e limites de verificação.',
                        'Funções atuais do Chimera, alternativas Python e condições para uma migração validada.',
                        'Geração SES nativa e comparação com o algoritmo C original da UCSF.',
                        'Critérios PDB, busca direta ChEMBL, filtros opcionais e relatórios das consultas.',
                        'Pares explícitos, cadeia comum, cofatores, preparação, tabela de RMSD, downloads e visualização 3D.',
                        'Eventos correlacionados, códigos de falha, resumos por job e consultas filtradas.'],
           'read':'Ler guia','start':'Primeira instalação? Baixe o projeto e configure as ferramentas obrigatórias antes de iniciar.',
           'footer':'Documentação gerada a partir dos arquivos Markdown do projeto. Licença MIT.'},
}

CSS = '''
:root{color-scheme:light;--ink:#172b4d;--muted:#526580;--teal:#087f74;--line:#dce5ef;--bg:#f4f7fb}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.7 system-ui,sans-serif}
a{color:var(--teal);text-underline-offset:3px}a:hover{color:#065b53}a:focus-visible,summary:focus-visible{outline:3px solid #0ea5e9;outline-offset:4px}
.top{background:#102b46;color:white;padding:20px max(24px,calc((100vw - 1240px)/2));display:flex;justify-content:space-between;gap:20px;flex-wrap:wrap}
.top a{color:#99f6e4}.brand{font-weight:750;text-decoration:none}.top nav{display:flex;gap:24px;flex-wrap:wrap}
.layout{max-width:1240px;margin:36px auto;display:grid;grid-template-columns:250px minmax(0,1fr);gap:32px;padding:0 24px}
aside{align-self:start;position:sticky;top:24px}aside a{display:block;margin:8px 0;text-decoration:none;padding:6px 12px;border-radius:8px}
aside a[aria-current=page]{background:#d9f2ed;font-weight:650}aside h2{font-size:13px;text-transform:uppercase;color:var(--muted);letter-spacing:.06em;margin-top:28px}
.toc a{font-size:14px;margin:2px 0}article{background:white;border:1px solid var(--line);border-radius:20px;padding:32px 40px;min-width:0}
h1,h2,h3{line-height:1.25;scroll-margin-top:24px}h1{font-size:clamp(28px,4vw,40px);margin-top:0}h2{font-size:25px;margin-top:40px}h3{font-size:20px}
article p,article li{overflow-wrap:anywhere}code{font: .88em/1.5 ui-monospace,monospace;background:#edf3f8;border-radius:4px;padding:2px 5px}
pre{background:#102b46;color:#e0f2fe;border-radius:12px;padding:20px;overflow:auto;font-size:14px;line-height:1.6}pre code{background:none;padding:0;color:inherit}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:12px;margin:24px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:12px 14px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#edf3f8}tr:last-child td{border-bottom:0}
.source{display:inline-block;font-size:14px;margin-bottom:20px}.footer{max-width:1240px;margin:36px auto;padding:0 24px;color:var(--muted);font-size:13px}
.portal{max-width:1140px;margin:64px auto;padding:0 24px}.eyebrow{font-size:13px;text-transform:uppercase;letter-spacing:.12em;font-weight:700;color:var(--teal)}
.portal h1{max-width:750px;font-size:clamp(32px,5vw,54px)}.lead{font-size:20px;color:var(--muted);max-width:780px}.cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:22px;margin:36px 0}
.card{border:1px solid var(--line);border-radius:18px;background:white;padding:28px;display:flex;flex-direction:column;gap:16px}.card h2{margin:0;font-size:23px}.card p{margin:0;flex:1;color:var(--muted)}.card .number{color:var(--teal);font-weight:750;font-size:14px}
.note{padding:22px 28px;background:#e0f4ef;border-radius:14px}.architecture-edges{padding:0;list-style:none;display:grid;gap:10px}.architecture-edges li{display:flex;align-items:center;gap:12px}.architecture-edges span{padding:8px 12px;background:#edf3f8;border:1px solid var(--line);border-radius:8px;font-size:14px}
details{margin:16px 0}summary{cursor:pointer;color:var(--teal)}
@media(max-width:900px){.layout{grid-template-columns:1fr}aside{position:static}.toc{display:none}aside nav{display:flex;gap:12px;flex-wrap:wrap}article{padding:24px}.cards{grid-template-columns:1fr}}
@media(max-width:480px){.layout{padding:0 12px;margin:20px auto}article{padding:20px 16px}.top{padding:16px}.portal{margin:36px auto}.architecture-edges li{flex-wrap:wrap}}
'''


def source_path(language, guide):
    return DOCS / ('en' if language == 'en' else '') / f'{guide}.md'


def output_path(language, guide):
    return DOCS / language / f'{guide}.html'


def href_for(url, source, output):
    if re.match(r'^(https?://|#|mailto:)', url):
        return url
    if ':' in url.split('/')[0]:
        raise ValueError(f'Unsupported link scheme: {url}')
    path, _, fragment = url.partition('#')
    destination = (source.parent / path).resolve()
    for language in LABELS:
        for guide in GUIDES:
            if destination == source_path(language, guide).resolve():
                destination = output_path(language, guide)
    if destination == DOCS / 'README.md':
        destination = DOCS / ('pt/index.html' if output.parent.name == 'pt' else 'index.html')
    return Path(os.path.relpath(destination, output.parent)).as_posix() + ('#' + fragment if fragment else '')


def inline(value, source, output):
    tokens = []

    def keep(fragment):
        tokens.append(fragment)
        return f'\x01{len(tokens)-1}\x02'

    value = re.sub(r'`([^`]+)`', lambda m: keep('<code>' + html.escape(m[1]) + '</code>'), value)
    value = re.sub(r'\[([^\]]+)\]\(([^\s)]+)\)',
                   lambda m: keep(f'<a href="{html.escape(href_for(m[2],source,output),quote=True)}">{html.escape(m[1])}</a>'), value)
    value = html.escape(value)
    value = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', value)
    value = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<em>\1</em>', value)
    return re.sub(r'\x01(\d+)\x02', lambda m: tokens[int(m[1])], value)


def render_markdown(source, output, language):
    lines = source.read_text(encoding='utf-8').splitlines()
    fragments, toc, ids = [], [], set()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if line.startswith('```'):
            syntax = line[3:].strip()
            block = []
            i += 1
            while i < len(lines) and not lines[i].startswith('```'):
                block.append(lines[i])
                i += 1
            if i == len(lines):
                raise ValueError(f'Unclosed code fence in {source}')
            escaped = html.escape('\n'.join(block))
            if syntax == 'mermaid':
                labels = dict(re.findall(r'(\w+)\[(.*?)\]', '\n'.join(block)))
                edges = re.findall(r'^\s*(\w+)(?:\[.*?\])?\s*-->\s*(\w+)', '\n'.join(block), re.M)
                graph = ''.join(f'<li><span>{html.escape(labels.get(a,a).strip("()"))}</span><b aria-hidden="true">→</b><span>{html.escape(labels.get(b,b).strip("()"))}</span></li>' for a,b in edges)
                fragments.append(f'<figure><figcaption>{LABELS[language]["diagram"]}</figcaption><ul class="architecture-edges">{graph}</ul></figure><details><summary>{LABELS[language]["diagram_source"]}</summary><pre><code>{escaped}</code></pre></details>')
            else:
                fragments.append(f'<pre><code class="language-{html.escape(syntax)}">{escaped}</code></pre>')
            i += 1
            continue
        heading = re.match(r'^(#{1,6})\s+(.+)$', line)
        if heading:
            level, value = len(heading[1]), heading[2]
            identifier = re.sub(r'[^a-z0-9]+','-',unicodedata.normalize('NFKD',value).encode('ascii','ignore').decode().lower()).strip('-') or 'section'
            base, counter = identifier, 2
            while identifier in ids:
                identifier = f'{base}-{counter}'
                counter += 1
            ids.add(identifier)
            fragments.append(f'<h{level} id="{identifier}">{inline(value,source,output)}</h{level}>')
            if level == 2:
                toc.append((identifier,value))
            i += 1
            continue
        if line.startswith('|') and i+1 < len(lines) and re.fullmatch(r'[| :\-]+',lines[i+1]):
            header = [cell.strip() for cell in line.strip('|').split('|')]
            rows = []
            i += 2
            while i < len(lines) and lines[i].startswith('|'):
                cells = [cell.strip() for cell in lines[i].strip('|').split('|')]
                if len(cells)!=len(header):
                    raise ValueError(f'Invalid table in {source}: {lines[i]}')
                rows.append('<tr>' + ''.join(f'<td>{inline(cell,source,output)}</td>' for cell in cells) + '</tr>')
                i += 1
            fragments.append('<div class="table-wrap"><table><thead><tr>' + ''.join(f'<th scope="col">{inline(cell,source,output)}</th>' for cell in header) + '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table></div>')
            continue
        item = re.match(r'^(?:- |\d+\. )(.+)$',line)
        if item:
            ordered = bool(re.match(r'^\d+\.',line))
            pattern = r'^\d+\. (.+)$' if ordered else r'^- (.+)$'
            items = []
            while i<len(lines) and (match:=re.match(pattern,lines[i])):
                items.append('<li>' + inline(match[1],source,output) + '</li>')
                i += 1
            tag = 'ol' if ordered else 'ul'
            fragments.append(f'<{tag}>' + ''.join(items) + f'</{tag}>')
            continue
        paragraph = [line.strip()]
        i += 1
        while i<len(lines) and lines[i].strip() and not re.match(r'^(#|```|\||- |\d+\. )',lines[i]):
            paragraph.append(lines[i].strip())
            i += 1
        fragments.append('<p>' + inline(' '.join(paragraph),source,output) + '</p>')
    return '\n'.join(fragments), toc


def page(title, language, content, home, portal, other):
    labels = LABELS[language]
    return f'''<!DOCTYPE html>
<!-- Generated by scripts/build_docs.py. Edit the Markdown sources instead. -->
<html lang="{'pt-BR' if language=='pt' else 'en'}"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1"><meta name="description" content="{html.escape(title)} — BioMolExplorer documentation"><title>{html.escape(title)} | BioMolExplorer</title><link rel="stylesheet" href="{'../' if home=='../../index.html' else ''}assets/docs.css"></head>
<body><header class="top"><a class="brand" href="{home}">BioMolExplorer</a><nav aria-label="{'Navegação' if language=='pt' else 'Navigation'}"><a href="{home}">{labels['home']}</a><a href="{portal}">{labels['docs']}</a><a href="{other}" hreflang="{'en' if language=='pt' else 'pt-BR'}" lang="{'en' if language=='pt' else 'pt-BR'}">{labels['other']}</a></nav></header>{content}<footer class="footer">{labels['footer']}</footer></body></html>'''


def build():
    (DOCS/'assets').mkdir(exist_ok=True)
    (DOCS/'assets/docs.css').write_text(CSS.strip()+'\n',encoding='utf-8')
    for language, labels in LABELS.items():
        (DOCS/language).mkdir(exist_ok=True)
        for index,guide in enumerate(GUIDES):
            source, output = source_path(language,guide), output_path(language,guide)
            body, toc = render_markdown(source,output,language)
            navigation = ''.join(f'<a href="{name}.html"'+(' aria-current="page"' if name==guide else '')+f'>{labels["titles"][n]}</a>' for n,name in enumerate(GUIDES))
            contents = ''.join(f'<a href="#{identifier}">{html.escape(title)}</a>' for identifier,title in toc)
            source_link = os.path.relpath(source,output.parent)
            content = f'<main class="layout"><aside><nav aria-label="{labels["docs"]}">{navigation}</nav><div class="toc"><h2>{labels["contents"]}</h2>{contents}</div></aside><article><a class="source" href="{source_link}">{labels["source"]}</a>{body}</article></main>'
            output.write_text(page(labels['titles'][index],language,content,'../../index.html','../index.html' if language=='en' else 'index.html',f'../{"pt" if language=="en" else "en"}/{guide}.html'),encoding='utf-8')
        portal = DOCS/'index.html' if language=='en' else DOCS/'pt/index.html'
        prefix = 'en/' if language=='en' else ''
        cards = ''.join(f'<section class="card"><span class="number">0{n+1}</span><h2>{labels["titles"][n]}</h2><p>{labels["summaries"][n]}</p><a href="{prefix}{guide}.html">{labels["read"]} →</a></section>' for n,guide in enumerate(GUIDES))
        content = f'<main class="portal"><p class="eyebrow">BioMolExplorer · {labels["docs"]}</p><h1>{labels["intro"]}</h1><p class="lead">{labels["description"]}</p><div class="cards">{cards}</div><p class="note">{labels["start"]} <a href="{prefix}{GUIDES[0]}.html">{labels["titles"][0]} →</a></p></main>'
        portal.write_text(page(labels['docs'],language,content,'../index.html' if language=='en' else '../../index.html','./index.html' if language=='en' else 'index.html','pt/index.html' if language=='en' else '../index.html'),encoding='utf-8')
    print(f'Built {len(GUIDES)*len(LABELS)} documentation pages and {len(LABELS)} language portals.')


if __name__=='__main__':
    build()
