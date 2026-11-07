"""Bounded scientific interpretation of persisted pipeline results."""
import ast
import csv
import io
import json
import math
from collections import Counter
from pathlib import Path

MAX_ROWS=5000
MAX_ANALYZED_FILES_PER_STAGE=32

from .i18n import locales,localized_dictionary
from .ui.localization import Translator

ESSENTIAL={'target','search_term','pdb_ids','reference','reference_type','metric','threshold','similarity_threshold',
 'pH','charge_type','preparation_tool','redocking_engine','rmsd_threshold','exhaustiveness','num_modes',
 'scoring_protocol','consensus_method','admet_profile','selection_rule','morgan','maccs','pharmacophore','radius','morgan_n_bits'}


def _number(value):
    try:value=float(value)
    except (TypeError,ValueError):return None
    return value if math.isfinite(value) else None


def _stats(values,unit=None):
    values=sorted(value for value in values if value is not None)
    if not values:return None
    middle=len(values)//2;median=values[middle] if len(values)%2 else (values[middle-1]+values[middle])/2
    mean=sum(values)/len(values);variance=sum((v-mean)**2 for v in values)/(len(values)-1) if len(values)>1 else None
    return {'n':len(values),'mean':mean,'median':median,'min':values[0],'max':values[-1],
        'sd':math.sqrt(variance) if variance is not None else None,'unit':unit}


def _rows(path):
    with path.open(encoding='utf-8-sig',errors='replace',newline='') as stream:
        reader=csv.DictReader(stream);rows=[]
        for index,row in enumerate(reader):
            if index>=MAX_ROWS:break
            rows.append(row)
    return reader.fieldnames or [],rows


def _first(fields,*candidates):
    lookup={field.casefold():field for field in fields}
    for candidate in candidates:
        if candidate.casefold() in lookup:return lookup[candidate.casefold()]
    return None


def _histogram(values,title,xlabel,language='en'):
    if len(values)<3 or min(values)==max(values):return None
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    figure,axis=plt.subplots(figsize=(6.4,3.6));axis.hist(values,bins=min(20,max(5,int(math.sqrt(len(values))))),color='#4878CF',edgecolor='white')
    axis.set_title(title);axis.set_xlabel(xlabel);axis.set_ylabel(Translator(language)('Count'));axis.grid(axis='y',alpha=.2);figure.tight_layout()
    output=io.BytesIO();figure.savefig(output,format='png',dpi=150);plt.close(figure);return output.getvalue()


def relevant_parameters(stage):
    from .catalog import new_stage
    configured=stage.get('configuration',stage);operation=configured.get('operation',stage.get('operation'))
    current=configured.get('parameters',{});defaults=new_stage(operation).get('parameters',{}) if operation else {}
    return {key:value for key,value in current.items() if key in ESSENTIAL or value!=defaults.get(key)}


def pipeline_diagram(pipeline,language='en'):
    if not pipeline:return None
    import base64
    from .ui.pipeline_png import render_pipeline
    data=render_pipeline(pipeline,Translator(language))
    return {'type':'figure','title':'Pipeline overview',
        'caption':'Canvas snapshot preserving block positions, input ports and data-flow connections.',
        'mime':'image/png','data':base64.b64encode(data).decode()}


def _analyze_table(operation,fields,rows,parameters,*,max_charts=8,language='en'):
    tr=Translator(language)
    def histogram(values,title,xlabel):
        if len(charts)>=max_charts:return None
        return _histogram(values,tr(title),tr(xlabel),language)
    metrics=[];tables=[];charts=[]
    identifier=_first(fields,'molecule_chembl_id','compound_id','name','PDB_CODE','pdb_id')
    if operation in ('retrieve_compounds','retrieve_pubchem','retrieve_zinc'):
        unique=len({row.get(identifier) for row in rows if identifier and row.get(identifier)}) if identifier else None
        metrics=[{'name':'retrieved_records','value':len(rows),'kind':'count'}]
        if unique is not None:metrics.append({'name':'unique_compounds','value':unique,'kind':'count'})
        activity=_first(fields,'standard_value','pchembl_value','activity_value')
        measure=_first(fields,'standard_type','activity_type');unit=_first(fields,'standard_units','units')
        if activity:
            groups={}
            for row in rows:
                value=_number(row.get(activity));key=(row.get(measure) or activity,row.get(unit) or None)
                if value is not None:groups.setdefault(key,[]).append(value)
            for (name,u),values in groups.items():
                summary=_stats(values,u);tables.append({'title':'Activity '+str(name),'rows':[summary]})
                image=histogram(values,'Activity distribution: '+str(name),str(name)+(f' ({u})' if u else ''))
                if image:charts.append({'title':'Activity distribution: '+str(name),'caption':'Observed activity values; incompatible measures are shown separately.','data':image})
    elif operation=='retrieve_structures':
        pdb=_first(fields,'PDB_CODE','pdb_id');resolution=_first(fields,'RESOLUTION','resolution')
        ids=sorted({row.get(pdb) for row in rows if pdb and row.get(pdb)})
        metrics=[{'name':'retrieved_structures','value':len(ids) or len(rows),'kind':'count'}]
        columns=[c for c in (pdb,_first(fields,'TARGET','target'),_first(fields,'ORGANISM','organism'),_first(fields,'METHOD','experimental_method'),resolution,_first(fields,'LIGAND','ligand')) if c]
        tables.append({'title':'PDB structures','columns':columns,'rows':[{c:r.get(c) for c in columns} for r in rows[:50]]})
        if resolution:
            values=[_number(r.get(resolution)) for r in rows];summary=_stats(values,'Å')
            if summary:metrics.append({'name':'resolution','value':summary,'kind':'statistics'})
    elif operation=='admet':
        selected=_first(fields,'ADMET_selected','selected','accepted');total=len(rows);accepted=sum(str(r.get(selected,'')).casefold() in ('true','1','yes') for r in rows) if selected else None
        metrics=[{'name':'evaluated_compounds','value':total,'kind':'count'}]
        if accepted is not None:metrics.extend([{'name':'accepted_compounds','value':accepted,'kind':'count'},
            {'name':'retention','value':100*accepted/total if total else None,'kind':'percent','denominator':total},
            {'name':'rejected_compounds','value':total-accepted,'kind':'count'}])
        for field in fields:
            values=[_number(r.get(field)) for r in rows];values=[v for v in values if v is not None]
            if len(values)>=max(3,len(rows)//2):
                summary=_stats(values);tables.append({'title':field,'rows':[summary]})
                image=histogram(values,'Distribution of '+field,field)
                if image and len(charts)<3:charts.append({'title':'Distribution of '+field,'caption':'Descriptive distribution of model-derived molecular properties.','data':image})
    elif operation=='fingerprints':
        fp=_first(fields,'fingerprint','morgan','maccs','pharmacophore');valid=0;densities=[]
        for row in rows:
            try:bits=ast.literal_eval(row.get(fp,'')) if fp else []
            except (ValueError,SyntaxError):bits=[]
            if isinstance(bits,(list,tuple)) and bits:valid+=1;densities.append(sum(bool(v) for v in bits)/len(bits))
        metrics=[{'name':'molecules_processed','value':len(rows),'kind':'count'},{'name':'valid_fingerprints','value':valid,'kind':'count'},
            {'name':'invalid_fingerprints','value':len(rows)-valid,'kind':'count'}]
        if densities:metrics.append({'name':'bit_density','value':_stats(densities),'kind':'statistics'})
    elif operation in ('similarity','graphs'):
        source=_first(fields,'source');target=_first(fields,'target');score=_first(fields,'value','similarity','score')
        values=[_number(r.get(score)) for r in rows] if score else [];values=[v for v in values if v is not None]
        nodes={r.get(key) for r in rows for key in (source,target) if key and r.get(key)}
        metrics=[{'name':'molecules','value':len(nodes),'kind':'count'},{'name':'relationships','value':len(rows),'kind':'count'}]
        if values:metrics.append({'name':'similarity_coefficients','value':_stats(values),'kind':'statistics'})
        if operation=='graphs' and source and target:
            import networkx as nx
            graph=nx.Graph();graph.add_edges_from((r.get(source),r.get(target)) for r in rows if r.get(source) and r.get(target))
            components=list(nx.connected_components(graph));metrics.extend([{'name':'connected_components','value':len(components),'kind':'count'},
                {'name':'largest_component','value':max(map(len,components),default=0),'kind':'count'},
                {'name':'network_density','value':nx.density(graph) if graph.number_of_nodes()>1 else None,'kind':'ratio'}])
            degrees=[degree for _,degree in graph.degree()];image=histogram(degrees,'Degree distribution','Degree')
        else:image=histogram(values,'Similarity distribution','Similarity coefficient')
        if image:charts.append({'title':'Degree distribution' if operation=='graphs' else 'Similarity distribution','caption':'Descriptive distribution calculated from persisted relationships.','data':image})
        if score:tables.append({'title':'Highest recorded relationships','columns':[source,target,score],
            'rows':sorted(({source:r.get(source),target:r.get(target),score:r.get(score)} for r in rows),key=lambda r:_number(r.get(score)) or -math.inf,reverse=True)[:20]})
    elif operation=='redocking':
        rmsd=_first(fields,'RMSD','rmsd');values=[_number(r.get(rmsd)) for r in rows] if rmsd else [];values=[v for v in values if v is not None]
        threshold=_number(parameters.get('rmsd_threshold')) or 2.0;success=sum(v<=threshold for v in values)
        metrics=[{'name':'evaluated_complexes','value':len(values),'kind':'count'},{'name':'rmsd','value':_stats(values,'Å'),'kind':'statistics'},
            {'name':'success_rate','value':100*success/len(values) if values else None,'kind':'percent','denominator':len(values),'criterion':f'RMSD ≤ {threshold} Å'}]
        if rmsd:tables.append({'title':'Redocking complexes','columns':[c for c in ('PDB_CODE','LIGAND','CHAIN',rmsd) if c in fields],'rows':rows[:50]})
        image=histogram(values,'Redocking RMSD distribution','RMSD (Å)')
        if image:charts.append({'title':'Redocking RMSD distribution','caption':f'Success criterion: RMSD ≤ {threshold} Å.','data':image})
    elif operation in ('docking_vina','docking_dock6','consensus','dynamics_dock6'):
        candidates=('affinity','vina_score','Grid_Score','grid_score','score','amber_score','total_energy')
        score=_first(fields,*candidates);values=[_number(r.get(score)) for r in rows] if score else [];values=[v for v in values if v is not None]
        unit='kcal/mol' if score and any(term in score.casefold() for term in ('affinity','amber','energy')) else None
        metrics=[{'name':'reported_results','value':len(rows),'kind':'count'}]
        if values:metrics.append({'name':score,'value':_stats(values,unit),'kind':'statistics'})
        if score:
            identity=[c for c in (_first(fields,'receptor','PDB_CODE','target'),_first(fields,'ligand','compound','molecule_chembl_id'),score) if c]
            tables.append({'title':'Top-ranked persisted results','columns':identity,'rows':sorted(rows,key=lambda r:_number(r.get(score)) if _number(r.get(score)) is not None else math.inf)[:20]})
            image=histogram(values,'Distribution of '+score,score+(f' ({unit})' if unit else ''))
            if image:charts.append({'title':'Distribution of '+score,'caption':'Computational scores are ranking quantities and are not experimental binding affinities.','data':image})
    else:
        metrics=[{'name':'persisted_records','value':len(rows),'kind':'count'}]
    return metrics,tables,charts


def analyze_run(store,project_id,run,pipeline,language='en'):
    configured={stage['id']:stage for stage in pipeline};analyses=[];chart_budget=8
    for stage in run.get('stages',[]):
        operation=stage.get('operation');configuration=stage.get('configuration') or configured.get(stage.get('id'),stage)
        analysis={'stage_id':stage.get('id'),'name':stage.get('name'),'operation':operation,'status':stage.get('status'),
            'method':{code:localized_dictionary(code,'methods').get(operation,Translator(code)('Operação configurada do pipeline BioMolExplorer.')) for code in locales()},'parameters':relevant_parameters(dict(stage,configuration=configuration)),
            'metrics':[],'tables':[],'charts':[],'warnings':[]}
        if stage.get('status')!='succeeded':
            if stage.get('error'):analysis['warnings'].append(str(stage['error']))
            analyses.append(analysis);continue
        paths=[]
        inspected=0;structured_skipped=False
        for artifact in stage.get('artifacts',[]):
            try:path=store.scoped_path(project_id,artifact)
            except Exception:continue
            if not path.is_file():continue
            paths.append(path)
            if path.suffix.lower() not in ('.csv','.json'):continue
            if inspected>=MAX_ANALYZED_FILES_PER_STAGE:
                structured_skipped=True;continue
            inspected+=1
            if path.suffix.lower()=='.json' and path.stat().st_size<=2*1024*1024:
                try:data=json.loads(path.read_text(encoding='utf-8'))
                except (OSError,ValueError):data=None
                if isinstance(data,dict):
                    keys=('compounds','duplicates','conformations','selected','evaluated','excluded_records','failed','failures','accepted','rejected')
                    for key in keys:
                        value=data.get(key)
                        if isinstance(value,(int,float)) and not isinstance(value,bool):analysis['metrics'].append({'name':key,'value':value,'kind':'count'})
            if path.suffix.lower()!='.csv' or path.stat().st_size>8*1024*1024:continue
            try:fields,rows=_rows(path)
            except (OSError,csv.Error):continue
            if not rows:continue
            metrics,tables,charts=_analyze_table(operation,fields,rows,configuration.get('parameters',{}),max_charts=chart_budget,language=language)
            for metric in metrics:metric['dataset']=path.name
            for table in tables:table['title']=path.stem+' · '+table.get('title','Results')
            analysis['metrics'].extend(metrics);analysis['tables'].extend(tables)
            for chart in charts:
                if chart_budget<=0:break
                import base64
                chart['mime']='image/png';chart['data']=base64.b64encode(chart['data']).decode();analysis['charts'].append(chart);chart_budget-=1
        if operation=='prepare_structures':
            receptors={path.name for path in paths if path.suffix.lower() in ('.pdb','.pdbqt') and not any(word in path.stem.casefold() for word in ('ligand','lig.','pose'))}
            ligands={path.name for path in paths if path.suffix.lower() in ('.mol2','.sdf','.pdbqt') and any(word in path.stem.casefold() for word in ('ligand','lig.','mol','compound'))}
            analysis['metrics'].extend([{'name':'prepared_receptor_files','value':len(receptors),'kind':'count'},
                {'name':'prepared_ligand_files','value':len(ligands),'kind':'count'}])
        if structured_skipped:
            analysis['warnings'].append({'en':
                f'Analysis was bounded to {MAX_ANALYZED_FILES_PER_STAGE} structured result files for this stage; all artifacts remain available in the project.',
                'pt':f'A análise foi limitada a {MAX_ANALYZED_FILES_PER_STAGE} arquivos estruturados nesta etapa; todos os artefatos continuam disponíveis no projeto.'})
        if stage.get('error') and str(stage['error']) not in analysis['warnings']:analysis['warnings'].append(str(stage['error']))
        analyses.append(analysis)
    return analyses
