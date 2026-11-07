"""Portable, path-free experimental methods report."""
import copy
import csv
import html
import io
import json
import math
import time
import zipfile
from base64 import b64decode,b64encode
from datetime import datetime,timezone
from pathlib import Path
from uuid import uuid4
from xml.sax.saxutils import escape

from .artifact_contracts import operation_contract
from .ui.localization import Translator,default_stage_names

MAX_RESULT_TABLES=12
MAX_RESULT_ROWS=5000
MAX_FIGURES=8
MAX_FIGURE_BYTES=4*1024*1024
MAX_DOCUMENT_RUNS=50
MAX_DOCUMENT_ARTIFACTS=250
MAX_PARAMETER_TEXT=700


def _semantic_stage(stage):
    result={key:copy.deepcopy(stage.get(key)) for key in
            ('id','block_id','name','operation','enabled','parameters','bindings','input_processing') if key in stage}
    # Physical paths are runtime adapter details and must not appear as method
    # configuration. Uploaded assets remain represented by stable IDs.
    for key in ('base_input_path','base_selected_mols','base_vina_path','base_dock6_path','similarity_path','dock6_app_path'):
        result.get('parameters',{}).pop(key,None)
    result['contract']=operation_contract(stage['operation'])
    return result


def _result_evidence(store,project_id,run):
    """Bounded, read-only summaries from the latest run for publication use."""
    tables=[];figures=[]
    for stage in run.get('stages',[]):
        for value in stage.get('artifacts',[]):
            if len(tables)>=MAX_RESULT_TABLES and len(figures)>=MAX_FIGURES:break
            try:path=store.scoped_path(project_id,value)
            except Exception:continue
            if not path.is_file():continue
            suffix=path.suffix.lower()
            if suffix=='.csv' and len(tables)<MAX_RESULT_TABLES and path.stat().st_size<=8*1024*1024:
                try:
                    with path.open(encoding='utf-8-sig',errors='replace',newline='') as stream:
                        reader=csv.DictReader(stream);fields=(reader.fieldnames or [])[:40];rows=[]
                        for index,row in enumerate(reader):
                            if index>=MAX_RESULT_ROWS:break
                            rows.append(row)
                    numeric={}
                    for field in fields:
                        values=[]
                        for row in rows:
                            try:value=float(row.get(field,''))
                            except (TypeError,ValueError):continue
                            if math.isfinite(value):values.append(value)
                        if values and len(values)>=max(2,len(rows)//2):
                            numeric[field]={'n':len(values),'min':min(values),'max':max(values),'mean':sum(values)/len(values)}
                    tables.append({'stage':stage.get('name'),'operation':stage.get('operation'),'name':path.name,
                        'rows_sampled':len(rows),'columns':fields,'numeric_summary':numeric})
                except (OSError,csv.Error):pass
            elif suffix in ('.png','.jpg','.jpeg') and len(figures)<MAX_FIGURES and path.stat().st_size<=MAX_FIGURE_BYTES:
                try:
                    from PIL import Image as PILImage
                    image=PILImage.open(path);image.thumbnail((1400,900));buffer=io.BytesIO()
                    if image.mode not in ('RGB','L'):image=image.convert('RGB')
                    image.save(buffer,'JPEG',quality=82,optimize=True)
                    figures.append({'stage':stage.get('name'),'operation':stage.get('operation'),'name':path.name,
                        'width':image.width,'height':image.height,'mime':'image/jpeg','data':b64encode(buffer.getvalue()).decode()})
                except (OSError,ValueError):pass
    return {'tables':tables,'figures':figures,'limits':{'rows_per_table':MAX_RESULT_ROWS,
        'tables':MAX_RESULT_TABLES,'figures':MAX_FIGURES}}


def build_report_preview(store,token,project_id):
    """Read the dialog overview without opening results or rendering figures."""
    project=store.project(token,project_id)
    with store.connect() as db:
        runs=db.execute('SELECT COUNT(*) FROM runs WHERE project_id=?',(project_id,)).fetchone()[0]
        artifacts=db.execute('SELECT COUNT(*) FROM artifact_instances WHERE project_id=?',(project_id,)).fetchone()[0]
    return {'generated_at':datetime.now(timezone.utc).isoformat(),
        'project':{key:project[key] for key in ('id','name','description','revision')},
        'pipeline':[{key:stage.get(key) for key in ('id','block_id','name','operation')} for stage in project['pipeline']],
        'document_scope':{'execution_count':runs,'artifact_count':artifacts}}


def build_report(store,token,project_id,language='en'):
    project=store.project(token,project_id)
    with store.connect() as db:
        run_count=db.execute('SELECT COUNT(*) FROM runs WHERE project_id=?',(project_id,)).fetchone()[0]
        runs=[store._run_row(row) for row in db.execute(
            'SELECT * FROM runs WHERE project_id=? ORDER BY created DESC LIMIT ?',
            (project_id,MAX_DOCUMENT_RUNS)).fetchall()]
        artifact_count=db.execute('SELECT COUNT(*) FROM artifact_instances WHERE project_id=?',(project_id,)).fetchone()[0]
        artifacts=[dict(row) for row in db.execute('''SELECT artifact_id,occurrence_id,indexed_run_id,producer_run_id,
            stage_id,block_id,output_port,artifact_type,schema_version,content_sha256,state,metadata,dependencies,provenance
            FROM artifact_instances WHERE project_id=? ORDER BY created DESC LIMIT ?''',
            (project_id,MAX_DOCUMENT_ARTIFACTS)).fetchall()]
    for artifact in artifacts:
        for key in ('metadata','dependencies','provenance'):artifact[key]=json.loads(artifact[key])
    executions=[]
    for run in runs:
        stages=[]
        for stage in run['stages']:
            stages.append({key:copy.deepcopy(stage.get(key)) for key in
                ('id','name','operation','status','processing_version','producer_software','processing_runtime','input_manifest','artifact_manifest','error') if key in stage})
        executions.append({'id':run['id'],'status':run['status'],'created':run['created'],'updated':run['updated'],'stages':stages})
    evidence=_result_evidence(store,project_id,runs[0]) if runs else {'tables':[],'figures':[],'limits':{}}
    from .scientific_report import analyze_run,pipeline_diagram
    analyses=analyze_run(store,project_id,runs[0],project['pipeline'],language) if runs else []
    from . import __version__
    return {'schema':'biomolexplorer.experimental-report','schema_version':1,
        'generated_at':datetime.now(timezone.utc).isoformat(),
        'application':{'name':'BioMolExplorer','version':__version__},
        'project':{'id':project_id,'name':project['name'],'description':project['description'],'revision':project['revision']},
        'pipeline':[_semantic_stage(stage) for stage in project['pipeline']],
        'executions':executions,'artifacts':artifacts,'results':evidence,'analyses':analyses,
        'document_scope':{'execution_count':run_count,'included_executions':len(runs),
            'artifact_count':artifact_count,'included_artifacts':len(artifacts),
            'bounded':run_count>len(runs) or artifact_count>len(artifacts)},
        'pipeline_figure':pipeline_diagram(project['pipeline'],language),
        'notes':{'pt':'Configurações, versões, hashes e proveniência para reprodução do experimento.',
                 'en':'Configurations, versions, hashes and provenance for experiment reproduction.'}}


def public_report(report):
    """JSON remains machine-readable; embedded document images stay out of it."""
    result=copy.deepcopy(report)
    for figure in result.get('results',{}).get('figures',[]):figure.pop('data',None)
    for analysis in result.get('analyses',[]):
        for chart in analysis.get('charts',[]):chart.pop('data',None)
    if result.get('pipeline_figure'):result['pipeline_figure'].pop('data',None)
    return result


def report_bytes(store,token,project_id):
    return json.dumps(public_report(build_report(store,token,project_id)),ensure_ascii=False,indent=2,allow_nan=False).encode()


def printable_file(store,report,language='en'):
    directory=Path(store.staging)/'report-previews';directory.mkdir(parents=True,exist_ok=True)
    now=time.time()
    for old in directory.glob('*.html'):
        try:
            if now-old.stat().st_mtime>3600:old.unlink()
        except OSError:pass
    path=directory/(uuid4().hex+'.html');path.write_bytes(report_html(report,language,True));return path


def report_content(report,language='en'):
    """Human-readable view derived from the portable scientific manifest."""
    tr=Translator(language);project=report['project'];pipeline=report.get('pipeline',[])
    executions=report.get('executions',[]);artifacts=report.get('artifacts',[])
    latest=executions[0] if executions else None
    labels={
        'title':tr('Relatório Experimental Reprodutível'),
        'generated':tr('Gerado em'),
        'project':tr('Projeto'),'description':tr('Descrição'),
        'revision':tr('Revisão do projeto'),'summary':tr('Resumo experimental'),
        'pipeline':tr('Pipeline e parâmetros'),
        'runs':tr('Histórico de execuções'),
        'artifacts':tr('Artefatos científicos'),
        'blocks':tr('blocos configurados'),'executions':tr('execuções'),
        'outputs':tr('artefatos indexados'),'status':'Status','operation':tr('Operação'),
        'parameters':tr('Parâmetros'),'bindings':tr('Referências de entrada'),
        'latest':tr('Execução mais recente'),'none':tr('Não disponível'),
        'repro':tr('Reprodutibilidade'),
        'results':tr('Resultados processados e métricas descritivas'),
        'methods':tr('Materiais e métodos'),
        'software':tr('Software e ambiente de execução'),
    }
    def metric_value(metric):
        item=metric.get('value')
        if isinstance(item,dict):
            parts=[f"n={item.get('n')}",f"mean={item.get('mean'):.5g}",f"median={item.get('median'):.5g}",f"min={item.get('min'):.5g}",f"max={item.get('max'):.5g}"]
            if item.get('sd') is not None:parts.append(f"SD={item['sd']:.5g}")
            if item.get('unit'):parts.append(str(item['unit']))
            return '; '.join(parts)
        if item is None:return labels['none']
        suffix='%' if metric.get('kind')=='percent' else ''
        return (f'{item:.5g}' if isinstance(item,float) else str(item))+suffix
    def display_parameter(item):
        text=json.dumps(item,ensure_ascii=False,sort_keys=True) if isinstance(item,(dict,list,tuple)) else str(item)
        if len(text)>MAX_PARAMETER_TEXT:
            suffix=tr(' … [valor completo no manifesto JSON]')
            text=text[:MAX_PARAMETER_TEXT-len(suffix)]+suffix
        return text
    metric_names={'retrieved_records': 'Retrieved records', 'unique_compounds': 'Unique compounds', 'retrieved_structures': 'Retrieved structures', 'evaluated_compounds': 'Evaluated compounds', 'accepted_compounds': 'Accepted compounds', 'rejected_compounds': 'Rejected compounds', 'retention': 'Retention', 'molecules_processed': 'Processed molecules', 'valid_fingerprints': 'Valid fingerprints', 'invalid_fingerprints': 'Invalid fingerprints', 'bit_density': 'Fingerprint bit density', 'molecules': 'Molecules', 'relationships': 'Relationships', 'similarity_coefficients': 'Similarity coefficients', 'connected_components': 'Connected components', 'largest_component': 'Largest connected component', 'network_density': 'Network density', 'evaluated_complexes': 'Evaluated complexes', 'rmsd': 'RMSD', 'success_rate': 'Success rate', 'reported_results': 'Reported results', 'prepared_receptor_files': 'Prepared receptors', 'prepared_ligand_files': 'Prepared ligands'}
    def metric_label(name):return tr(metric_names.get(name,str(name).replace('_',' ').title()))
    def scientific_title(text):
        if not text:return text
        prefix,separator,payload=str(text).partition(' · ')
        return prefix+separator+tr(payload) if separator else tr(text)
    state=(tr('configurado, ainda não executado')) if not latest else tr(latest.get('status','unknown'))
    targets=[];sources=[];source_names={'retrieve_compounds':'ChEMBL','retrieve_pubchem':'PubChem','retrieve_zinc':'ZINC','retrieve_structures':'PDB'}
    for stage in pipeline:
        parameters=stage.get('parameters',{})
        for key in ('target','search_term','pdb_ids','pdb_codes'):
            target_value=parameters.get(key)
            if target_value not in (None,'',[]):targets.append(str(target_value))
        if stage.get('operation') in source_names:sources.append(source_names[stage['operation']])
    identity=[project.get('description') or (tr('Nenhum objetivo foi informado.')),
        {'type':'table','columns':[labels['project'],labels['status'],'BioMolExplorer',labels['generated']],
         'rows':[[project.get('name'),state,report.get('application',{}).get('version',''),report.get('generated_at')]]}]
    if targets:identity.append((tr('Alvos moleculares: '))+', '.join(dict.fromkeys(targets)))
    if sources:identity.append((tr('Fontes de dados: '))+', '.join(dict.fromkeys(sources)))
    sections=[(labels['project'],[
        *identity])]
    latest_stages={stage.get('id'):stage for stage in (latest or {}).get('stages',[])}
    from collections import Counter
    overview=[];by_stage={stage.get('id'):stage for stage in pipeline};stage_indices={stage.get('id'):index for index,stage in enumerate(pipeline,1)}
    block_counts=Counter(stage.get('block_id') for stage in pipeline if stage.get('block_id'))
    from .bindings import sources as binding_sources
    from .flow import input_ports
    def stage_label(stage,index=None):
        if not stage:return labels['none']
        stored_block=stage.get('block_id')
        block=stored_block if stored_block and block_counts[stored_block]==1 else (f'B{index:02d}' if index is not None else '')
        return (block+' · ' if block else '')+(tr(stage.get('name')) if stage.get('name') in default_stage_names() else str(stage.get('name') or stage.get('operation')))
    for index,stage in enumerate(pipeline,1):
        execution=latest_stages.get(stage.get('id'),{})
        incoming=[]
        for field,binding in stage.get('bindings',{}).items():
            for reference in binding_sources(binding):
                source=by_stage.get(reference.get('stage'))
                incoming.append(stage_label(source,stage_indices.get(source.get('id')) if source else None)+' → '+field)
        inputs='; '.join(incoming) or (tr('Fonte configurada'))
        outputs=', '.join((stage.get('contract') or {}).get('outputs',{})) or labels['none']
        overview.append([index,stage_label(stage,index),inputs,outputs,tr(execution.get('status',tr('não executado')))])
    pipeline_items=[]
    if report.get('pipeline_figure'):
        figure=dict(report['pipeline_figure']);figure['caption']=(
            tr('Captura do canvas: as posições dos blocos são preservadas; as setas coloridas terminam nas portas de entrada correspondentes.'))
        pipeline_items.append(figure)
    connections=[];connection_index=0
    for target_index,target in enumerate(pipeline,1):
        fields=[port['field'] for port in input_ports(target)]
        for field in fields:
            binding=target.get('bindings',{}).get(field)
            for reference in binding_sources(binding) if binding else []:
                source=by_stage.get(reference.get('stage'))
                if source is None:continue
                connection_index+=1
                connections.append([f'C{connection_index:02d}',stage_label(source,stage_indices.get(source.get('id'))),reference.get('output_port','auto'),
                    stage_label(target,target_index),field])
    if connections:
        pipeline_items.append({'type':'table','title':tr('Mapa de conexões'),
            'columns':['ID',tr('Bloco de origem'),tr('Porta de saída'),
                tr('Bloco de destino'),tr('Porta de entrada')],
            'rows':connections})
    pipeline_items.append({'type':'table','columns':['#',labels['operation'],tr('Entrada'),tr('Resultado'),labels['status']],'rows':overview})
    sections.append(((tr('Visão geral do pipeline')),pipeline_items))
    analyses=report.get('analyses',[])
    if analyses:
        for analysis in analyses:
            method=analysis.get('method','');items=[method.get(language,method.get('en','')) if isinstance(method,dict) else tr(method)]
            parameters=analysis.get('parameters',{})
            if parameters:items.append({'type':'table','columns':[labels['parameters'],tr('Valor')],
                'rows':[[key,display_parameter(parameter)] for key,parameter in sorted(parameters.items())]})
            metrics=analysis.get('metrics',[])
            if metrics:items.append({'type':'table','columns':[tr('Conjunto'),tr('Métrica'),tr('Valor'),tr('Amostra/critério')],
                'rows':[[metric.get('dataset','—'),metric_label(metric['name']),metric_value(metric),metric.get('criterion') or (('n='+str(metric['denominator'])) if metric.get('denominator') is not None else '')] for metric in metrics]})
            for table in analysis.get('tables',[]):
                rows=table.get('rows',[])[:25];columns=table.get('columns') or list(dict.fromkeys(key for row in rows if isinstance(row,dict) for key in row))
                items.append({'type':'table','title':scientific_title(table.get('title')),'columns':columns,'rows':[[row.get(column,'') for column in columns] for row in rows if isinstance(row,dict)]})
            for chart in analysis.get('charts',[]):items.append({'type':'figure',**chart,'title':scientific_title(chart.get('title'))})
            for warning in analysis.get('warnings',[]):
                warning=warning.get(language,warning.get('en','')) if isinstance(warning,dict) else warning
                items.append(tr('Atenção: ')+tr(str(warning)))
            sections.append((tr(analysis.get('name')) if analysis.get('name') in default_stage_names() else (analysis.get('name') or analysis.get('operation')),items))
    sections.append((labels['repro'],[
        (tr('A configuração completa, valores padrão, referências de artefatos, hashes e proveniência são preservados no manifesto JSON complementar.')),
        tr('Este documento inclui as {0} execuções mais recentes e os {1} registros de artefatos mais recentes; as contagens totais permanecem registradas no manifesto.').format(
            report.get('document_scope',{}).get('included_executions',len(executions)),
            report.get('document_scope',{}).get('included_artifacts',len(artifacts))),
        'Schema: '+str(report.get('schema'))+' v'+str(report.get('schema_version'))]))
    return labels,sections


def report_html(report,language='en',auto_print=False):
    tr=Translator(language)
    from .ui.branding import image_bytes
    labels,sections=report_content(report,language);project=report['project']
    logo=b64encode(image_bytes('logo.png')).decode()
    def body(section):
        title,items=section;rows=[]
        for item in items:
            if isinstance(item,dict) and item.get('type')=='table':
                caption='<h3>'+html.escape(str(item.get('title','')))+'</h3>' if item.get('title') else ''
                head=''.join('<th>'+html.escape(str(value))+'</th>' for value in item.get('columns',[]))
                values=''.join('<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in row)+'</tr>' for row in item.get('rows',[]))
                rows.append(caption+'<div class="table"><table><thead><tr>'+head+'</tr></thead><tbody>'+values+'</tbody></table></div>')
            elif isinstance(item,dict) and item.get('type')=='figure':
                rows.append('<figure><img src="data:'+item.get('mime','image/png')+';base64,'+item['data']+'"><figcaption>'+html.escape(item.get('caption') or item.get('title',''))+'</figcaption></figure>')
            elif isinstance(item,tuple):
                rows.append('<article><h3>'+html.escape(str(item[0]))+'</h3>'+''.join('<p>'+html.escape(str(v))+'</p>' for v in item[1])+'</article>')
            else:rows.append('<p>'+html.escape(str(item))+'</p>')
        pipeline_class=' class="pipeline-section"' if title==tr('Visão geral do pipeline') else ''
        return '<section'+pipeline_class+'><h2>'+html.escape(title)+'</h2>'+''.join(rows)+'</section>'
    script='<script>window.addEventListener("load",()=>window.print())</script>' if auto_print else ''
    figures=''.join('<figure><img src="data:'+item['mime']+';base64,'+item['data']+'"><figcaption>'+html.escape(item.get('stage','')+' · '+item['name'])+'</figcaption></figure>' for item in report.get('results',{}).get('figures',[]))
    figures=('<section><h2>'+html.escape(labels['results'])+' · '+(tr('Figuras'))+'</h2><div class="figures">'+figures+'</div></section>') if figures else ''
    title=project['name']+' — '+labels['title']
    return ('<!doctype html><html><head><meta charset="utf-8"><title>'+html.escape(title)+'</title>'+script+
        '<style>@page{size:A4;margin:20mm 18mm 18mm}@page pipeline{size:A4 landscape;margin:15mm}*{box-sizing:border-box}body{font:10.5pt Georgia,"Times New Roman",serif;color:#172b4d;margin:0;line-height:1.52}header{display:flex;align-items:center;gap:18px;border-bottom:3px solid #2552e8;padding-bottom:14px;margin-bottom:26px}header img{width:58px;height:auto}h1{font:700 20pt Arial,sans-serif;margin:0;color:#091d48}header p{font:9.5pt Arial,sans-serif;margin:5px 0 0;color:#52637a}h2{font:700 13.5pt Arial,sans-serif;color:#2552e8;border-bottom:1px solid #dce4f3;padding-bottom:5px;margin-top:24px}h3{font:700 10.5pt Arial,sans-serif;margin:13px 0 3px;color:#091d48}p{margin:4px 0;text-align:justify;overflow-wrap:anywhere}article{break-inside:avoid;padding:5px 0}.table{overflow:hidden;margin:9px 0 14px}table{width:100%;border-collapse:collapse;font:8.5pt Arial,sans-serif}th{background:#eaf0ff;color:#091d48;text-align:left}th,td{border-bottom:1px solid #dce4f3;padding:5px 7px;vertical-align:top}td:nth-child(n+2){font-variant-numeric:tabular-nums}.figures{display:grid;grid-template-columns:1fr 1fr;gap:14px}figure{margin:10px 0;break-inside:avoid;border:1px solid #dce4f3;padding:8px;background:#fff}figure img{display:block;width:100%;height:auto;max-height:95mm;object-fit:contain}figcaption{font:8.5pt Arial,sans-serif;color:#52637a;margin-top:6px;text-align:center}.pipeline-section{page:pipeline;break-before:page;break-after:page}.pipeline-section figure img{max-height:125mm}@media print{header{break-after:avoid}.figures{display:block}figure{margin:0 0 12mm}}</style></head><body><header><img src="data:image/png;base64,'+logo+'"><div><h1>'+html.escape(title)+'</h1><p>'+html.escape(labels['generated'])+': '+html.escape(report['generated_at'])+'</p></div></header>'+''.join(body(s) for s in sections)+figures+'</body></html>').encode()


def report_docx(report,language='en'):
    """Create an interoperable DOCX with no optional runtime dependency."""
    tr=Translator(language)
    from .ui.branding import image_bytes
    labels,sections=report_content(report,language);project=report['project']
    def paragraph(text,style=None):
        prop=f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>' if style else ''
        return '<w:p>'+prop+'<w:r><w:t xml:space="preserve">'+escape(str(text))+'</w:t></w:r></w:p>'
    def drawing(rid,index,width=548640,height=513360,name='logo.png'):
        return (f'''<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing><wp:inline xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"><wp:extent cx="{width}" cy="{height}"/><wp:docPr id="{index}" name="{escape(name)}"/><a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:nvPicPr><pic:cNvPr id="{index}" name="{escape(name)}"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed="{rid}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{width}" cy="{height}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>''')
    def table(item):
        rows=[item.get('columns',[])]+item.get('rows',[])
        return '<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/></w:tblPr>'+''.join('<w:tr>'+''.join('<w:tc><w:tcPr><w:shd w:fill="EAF0FF"/></w:tcPr>'+(paragraph(value,'TableText'))+'</w:tc>' for value in row)+'</w:tr>' for row in rows)+'</w:tbl>'
    def section_break(landscape=False):
        size='<w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>' if landscape else '<w:pgSz w:w="11906" w:h="16838"/>'
        return '<w:p><w:pPr><w:sectPr>'+size+'<w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360"/></w:sectPr></w:pPr></w:p>'
    title=project['name']+' — '+labels['title']
    content=[paragraph(title,'Title'),paragraph(labels['generated']+': '+report['generated_at'])];embedded=[]
    for section_title,items in sections:
        pipeline_section=section_title==tr('Visão geral do pipeline')
        if pipeline_section:content.append(section_break(False))
        content.append(paragraph(section_title,'Heading1'))
        for item in items:
            if isinstance(item,dict) and item.get('type')=='table':
                if item.get('title'):content.append(paragraph(item['title'],'Heading2'))
                content.append(table(item))
            elif isinstance(item,dict) and item.get('type')=='figure':
                embedded.append(item);index=len(embedded)
                width,height=(8500000,3570000) if pipeline_section else (4572000,2940000)
                content.extend([drawing('rId'+str(index+2),index+1,width,height,item.get('title','figure.png')),paragraph(item.get('caption') or item.get('title',''),'Caption')])
            elif isinstance(item,tuple):content.extend([paragraph(item[0],'Heading2'),*(paragraph(v) for v in item[1])])
            else:content.append(paragraph(item))
        if pipeline_section:content.append(section_break(True))
    figures=report.get('results',{}).get('figures',[])
    if figures:
        content.append(paragraph(labels['results']+' · '+(tr('Figuras')),'Heading1'))
        for index,figure in enumerate(figures,1):
            embedded.append(figure);number=len(embedded);content.extend([drawing('rId'+str(number+2),number+1,4572000,2940000,figure['name']),
                paragraph(figure.get('stage','')+' · '+figure['name'],'Caption')])
    header=drawing('rId1',1)+paragraph(title)
    document=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>'+''.join(content)+'<w:sectPr><w:headerReference w:type="default" r:id="rId2"/><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1440" w:right="1020" w:bottom="1020" w:left="1020" w:header="360"/></w:sectPr></w:body></w:document>')
    output=io.BytesIO()
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('[Content_Types].xml','''<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/><Default Extension="jpg" ContentType="image/jpeg"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/header1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.header+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>''')
        archive.writestr('_rels/.rels','''<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>''')
        relations=['<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>','<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/header" Target="header1.xml"/>']
        for index,figure in enumerate(embedded,1):
            extension='jpg' if figure.get('mime')=='image/jpeg' else 'png'
            relations.append(f'<Relationship Id="rId{index+2}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/figure{index}.{extension}"/>')
        archive.writestr('word/_rels/document.xml.rels','<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+''.join(relations)+'</Relationships>')
        archive.writestr('word/_rels/header1.xml.rels','''<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/logo.png"/></Relationships>''')
        archive.writestr('word/header1.xml','<?xml version="1.0" encoding="UTF-8"?><w:hdr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'+header+'</w:hdr>')
        archive.writestr('word/styles.xml','''<?xml version="1.0"?><w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:pPr><w:jc w:val="both"/><w:spacing w:after="120" w:line="276" w:lineRule="auto"/></w:pPr><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="22"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:rPr><w:b/><w:color w:val="091D48"/><w:sz w:val="36"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:rPr><w:b/><w:color w:val="2552E8"/><w:sz w:val="28"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:rPr><w:b/><w:color w:val="091D48"/><w:sz w:val="23"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="Caption"/><w:pPr><w:jc w:val="center"/></w:pPr><w:rPr><w:i/><w:color w:val="52637A"/><w:sz w:val="18"/></w:rPr></w:style></w:styles>''')
        archive.writestr('word/document.xml',document);archive.writestr('word/media/logo.png',image_bytes('logo.png'))
        for index,figure in enumerate(embedded,1):
            extension='jpg' if figure.get('mime')=='image/jpeg' else 'png'
            archive.writestr(f'word/media/figure{index}.{extension}',b64decode(figure['data']))
    return output.getvalue()


def report_pdf(report,language='en'):
    tr=Translator(language)
    from reportlab.lib.colors import HexColor
    from reportlab.lib.enums import TA_JUSTIFY,TA_LEFT
    from reportlab.lib.pagesizes import A4,landscape
    from reportlab.lib.styles import ParagraphStyle,getSampleStyleSheet
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import BaseDocTemplate,Frame,Image,NextPageTemplate,PageBreak,PageTemplate,Paragraph,Spacer,Table,TableStyle
    tr=Translator(language)
    from .ui.branding import image_bytes
    labels,sections=report_content(report,language);output=io.BytesIO();styles=getSampleStyleSheet()
    title=report['project']['name']+' — '+labels['title']
    styles.add(ParagraphStyle(name='BioTitle',parent=styles['Title'],textColor=HexColor('#091D48'),alignment=TA_LEFT,fontName='Helvetica-Bold'))
    styles.add(ParagraphStyle(name='Caption',parent=styles['Normal'],textColor=HexColor('#52637A'),fontSize=8,alignment=1,spaceAfter=8))
    styles['BodyText'].alignment=TA_JUSTIFY;styles['BodyText'].leading=15
    styles['Heading1'].textColor=HexColor('#2552E8')
    story=[Paragraph(html.escape(title),styles['BioTitle']),Paragraph(html.escape(labels['generated']+': '+report['generated_at']),styles['Normal']),Spacer(1,16)]
    for title,items in sections:
        pipeline_section=title==tr('Visão geral do pipeline')
        if pipeline_section:story.extend([NextPageTemplate('landscape'),PageBreak()])
        story.append(Paragraph(html.escape(title),styles['Heading1']))
        for item in items:
            if isinstance(item,dict) and item.get('type')=='table':
                if item.get('title'):story.append(Paragraph(html.escape(str(item['title'])),styles['Heading2']))
                data=([[Paragraph(html.escape(str(value)),styles['Normal']) for value in item.get('columns',[])]]+
                    [[Paragraph(html.escape(str(value)),styles['Normal']) for value in row] for row in item.get('rows',[])])
                if data and data[0]:
                    table=Table(data,repeatRows=1,hAlign='LEFT')
                    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#EAF0FF')),('TEXTCOLOR',(0,0),(-1,0),HexColor('#091D48')),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),7.5),('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),.25,HexColor('#DCE4F3')),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4)]));story.extend([table,Spacer(1,8)])
            elif isinstance(item,dict) and item.get('type')=='figure':
                width,height=(735,310) if pipeline_section else (430,275)
                story.extend([Image(io.BytesIO(b64decode(item['data'])),width=width,height=height,kind='proportional'),Paragraph(html.escape(item.get('caption') or item.get('title','')),styles['Caption']),Spacer(1,10)])
            elif isinstance(item,tuple):
                story.append(Paragraph(html.escape(str(item[0])),styles['Heading2']))
                story.extend(Paragraph(html.escape(str(value)),styles['BodyText']) for value in item[1])
            else:story.append(Paragraph(html.escape(str(item)),styles['BodyText']))
        story.append(Spacer(1,8))
        if pipeline_section:story.extend([NextPageTemplate('portrait'),PageBreak()])
    figures=report.get('results',{}).get('figures',[])
    if figures:
        story.append(Paragraph(html.escape(labels['results']+' · '+(tr('Figuras'))),styles['Heading1']))
        for figure in figures:
            story.extend([Image(io.BytesIO(b64decode(figure['data'])),width=430,height=275,kind='proportional'),
                Paragraph(html.escape(figure.get('stage','')+' · '+figure['name']),styles['Caption']),Spacer(1,12)])
    def header_footer(canvas,document):
        page_width,page_height=canvas._pagesize
        canvas.saveState();canvas.drawImage(ImageReader(io.BytesIO(image_bytes('logo.png'))),45,page_height-48,width=24,height=22,mask='auto')
        canvas.setFont('Helvetica-Bold',8);canvas.setFillColor(HexColor('#091D48'));canvas.drawString(76,page_height-40,title[:100])
        canvas.setStrokeColor(HexColor('#2552E8'));canvas.line(45,page_height-53,page_width-45,page_height-53)
        canvas.setFont('Helvetica',8);canvas.setFillColor(HexColor('#52637A'));canvas.drawCentredString(page_width/2,24,str(document.page));canvas.restoreState()
    document=BaseDocTemplate(output,pagesize=A4,rightMargin=45,leftMargin=45,topMargin=70,bottomMargin=42,title=title,author='BioMolExplorer')
    portrait_frame=Frame(45,42,A4[0]-90,A4[1]-112,id='portrait-frame')
    wide=landscape(A4);landscape_frame=Frame(45,42,wide[0]-90,wide[1]-112,id='landscape-frame')
    document.addPageTemplates([PageTemplate(id='portrait',pagesize=A4,frames=[portrait_frame],onPage=header_footer),
        PageTemplate(id='landscape',pagesize=wide,frames=[landscape_frame],onPage=header_footer)])
    document.build(story)
    return output.getvalue()
