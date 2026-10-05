"""Independent molecular graph experiments, common fragments and export figures."""
import ast
import base64
import io
import math
import re
import textwrap
from collections import Counter
from pathlib import Path

import networkx as nx
import pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import Draw,rdFMCS

from biomolexplorer.input_validation import ALIASES,columns,validate_file
from biomolexplorer.visualizations import SUFFIX,graph_view,write_view
from biomolexplorer.progress import report_progress
from biomolexplorer.molecule_quality import CURRENT_REPORT, QualityReport, clean_dataframe, filter_edges


def fragment_molecule(smarts, reference_smiles):
    """Extract actual atoms/bonds for the query from its reference structure."""
    reference = Chem.MolFromSmiles(reference_smiles) if isinstance(reference_smiles, str) else None
    query = Chem.MolFromSmarts(smarts) if isinstance(smarts, str) and smarts else None
    if reference is None or query is None:
        return None
    match = reference.GetSubstructMatch(query)
    if not match:
        return None
    matched_atoms = set(match)
    matched_bonds = {frozenset((match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()]))
                     for bond in query.GetBonds()}

    def extract(source):
        # Retain the reference's atom/neighbor order and stereo bond metadata.
        # Rebuilding atoms in query order can invert a tetrahedral center and
        # drops the direction/stereo atoms that define double-bond geometry.
        fragment = Chem.RWMol(source)
        for bond in source.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in matched_atoms and b in matched_atoms and frozenset((a, b)) not in matched_bonds:
                fragment.RemoveBond(a, b)
        for index in reversed(range(source.GetNumAtoms())):
            if index not in matched_atoms:
                fragment.RemoveAtom(index)
        molecule = fragment.GetMol()
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        return molecule

    molecule = extract(reference)
    # A partial aromatic ring cannot be written as valid aromatic SMILES. Use
    # the reference's Kekulé bond orders in that case, then sanitize the cut.
    molecule.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(molecule)
    partial_aromatic = any(atom.GetIsAromatic() and not atom.IsInRing() for atom in molecule.GetAtoms())
    if partial_aromatic or Chem.SanitizeMol(molecule, catchErrors=True) != Chem.SanitizeFlags.SANITIZE_NONE:
        reference = Chem.Mol(reference)
        Chem.Kekulize(reference, clearAromaticFlags=True)
        molecule = extract(reference)
        Chem.SanitizeMol(molecule)
    # A cut may remove a substituent required to define a stereocenter. Clear
    # only those assignments; preserve stereo supported by the retained graph.
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule


def common_fragment(smiles,timeout=30,ring_matches_ring_only=True,complete_rings_only=False):
    result={'status':'empty','smarts':'','smiles':'','atoms':0,'bonds':0,'image':'','compounds':len(smiles),
            'timeout':timeout,'ring_matches_ring_only':ring_matches_ring_only,'complete_rings_only':complete_rings_only}
    if not smiles:return result
    molecules=[];seen=set()
    for value in smiles:
        mol=Chem.MolFromSmiles(value) if isinstance(value,str) and value else None
        if mol is None:
            return dict(result,status='missing_structures')
        canonical=Chem.MolToSmiles(mol)
        if canonical not in seen:seen.add(canonical);molecules.append(mol)
    if len(molecules)==1:
        query=molecules[0];smarts=Chem.MolToSmarts(query);partial=False
    else:
        match=rdFMCS.FindMCS(molecules,timeout=timeout,threshold=1.0,
            ringMatchesRingOnly=ring_matches_ring_only,completeRingsOnly=complete_rings_only)
        smarts=match.smartsString;partial=bool(match.canceled)
        query=Chem.MolFromSmarts(smarts) if smarts else None
    if query is None or query.GetNumAtoms()==0:return dict(result,status='partial' if partial else 'no_common_fragment')
    # SMARTS remains the matching query; SMILES describes the actual fragment
    # rendered from the reference molecule, preserving atom and bond identities.
    depiction=fragment_molecule(smarts,Chem.MolToSmiles(molecules[0]))
    image=io.BytesIO();Draw.MolToImage(depiction,size=(360,250),kekulize=False).save(image,format='PNG')
    return dict(result,status='partial' if partial else 'complete',smarts=smarts,smiles=Chem.MolToSmiles(depiction),atoms=query.GetNumAtoms(),
                bonds=query.GetNumBonds(),image=base64.b64encode(image.getvalue()).decode('ascii'))


def read_compounds(paths,report=None):
    frames=[]
    for filename in paths:
        frame=pd.read_csv(filename,dtype={'molecule_chembl_id':str,'name':str}).rename(columns=ALIASES)
        if not {'molecule_chembl_id','canonical_smiles'}<=set(frame):continue
        validate_file(filename,'compounds',validate_rows=False)
        frame=clean_dataframe(frame,source=filename,report=report)
        if any(not isinstance(code,str) or not code for code in frame['molecule_chembl_id']):
            raise ValueError('A tabela de compostos desta análise precisa de códigos preenchidos em molecule_chembl_id.')
        frame=frame.drop(columns=['fingerprint'],errors='ignore')
        frames.append(frame)
    if not frames:return pd.DataFrame(columns=['molecule_chembl_id','canonical_smiles'])
    data=pd.concat(frames,ignore_index=True)
    conflicted=[]
    for identifier,group in data.groupby('molecule_chembl_id'):
        structures={Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in group['canonical_smiles'] if isinstance(s,str) and Chem.MolFromSmiles(s)}
        if len(structures)>1:
            (report or CURRENT_REPORT.get() or QualityReport()).exclude('tabelas de compostos',None,[identifier],'código associado a estruturas conflitantes')
            conflicted.append(identifier)
    data=data[~data['molecule_chembl_id'].isin(conflicted)]
    return data.drop_duplicates('molecule_chembl_id').reset_index(drop=True)


def exact_edges(filename,metric,threshold,report=None):
    """Compare each pair once; no approximate candidates are discarded."""
    validate_file(filename,'fingerprints',validate_rows=False)
    data=pd.read_csv(filename,dtype={'molecule_chembl_id':str,'name':str}).rename(columns=ALIASES)
    data=clean_dataframe(data,'fingerprints',filename,report)
    for identifier,group in data.groupby('molecule_chembl_id'):
        if len({tuple(ast.literal_eval(v)) for v in group['fingerprint']})>1:
            raise ValueError('O código '+identifier+' possui fingerprints diferentes no mesmo arquivo.')
    data=data.drop_duplicates('molecule_chembl_id')
    identifiers=list(data['molecule_chembl_id'])
    vectors=[]
    for bits in data['fingerprint']:
        values=ast.literal_eval(bits)
        vector=DataStructs.ExplicitBitVect(len(values))
        for i,value in enumerate(values):
            if value:vector.SetBit(i)
        vectors.append(vector)
    name=getattr(metric,'value',metric)
    function=getattr(DataStructs,name+'Similarity')
    rows=[]
    for index,left in enumerate(vectors):
        if index%100==0:report_progress(f'Calculando similaridade exata: {index}/{len(vectors)} compostos…')
        for j in range(index+1,len(vectors)):
            # Directed metrics are represented as an undirected molecular graph
            # only when both directions satisfy the chosen threshold.
            score=min(float(function(left,vectors[j])),float(function(vectors[j],left)))
            lower=-1 if name=='McConnaughey' else 0
            if not math.isfinite(score) or not lower<=score<=1:
                raise ValueError('Esta métrica produziu valores fora de 0–1. Use uma métrica compatível com grafos de similaridade.')
            if score>=threshold/100:
                rows.append({'source':identifiers[index],'target':identifiers[j],'value':score})
    return pd.DataFrame(rows,columns=['source','target','value']),identifiers


def build_graph(dataset,edges):
    graph=nx.Graph();graph.add_nodes_from(map(str,dataset['molecule_chembl_id']))
    known=set(graph)
    for row in edges.itertuples(index=False):
        source,target=str(row.source),str(row.target)
        if source not in known or target not in known:continue
        if source==target:continue
        weight=float(row.value)
        if graph.has_edge(source,target):weight=max(weight,graph[source][target]['value'])
        graph.add_edge(source,target,value=weight)
    return graph


def report_png(model):
    """MCC, query depiction and three degree panels, using isolated Agg figures."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.image import imread
    from matplotlib.ticker import MaxNLocator
    graph=nx.Graph();graph.add_nodes_from(n['id'] for n in model['nodes'])
    graph.add_edges_from((e['source'],e['target']) for e in model['edges'])
    mcc=graph.subgraph(model['mcc']);degrees=dict(mcc.degree())
    fig=Figure(figsize=(13,11),facecolor='#F8FAFD');FigureCanvasAgg(fig)
    fig.subplots_adjust(left=.07,right=.96,bottom=.07,top=.88)
    grid=fig.add_gridspec(3,6,height_ratios=[1.65,1,1],hspace=.55,wspace=1.05)
    ax=fig.add_subplot(grid[0,:4]);fragment_ax=fig.add_subplot(grid[0,4:])
    rank_ax=fig.add_subplot(grid[1,:3]);hist_ax=fig.add_subplot(grid[1,3:]);dist_ax=fig.add_subplot(grid[2,:])
    for item in (ax,fragment_ax,rank_ax,hist_ax,dist_ax):
        item.set_facecolor('white');item.spines[['top','right']].set_visible(False)
    positions=model.get('mcc_positions') or {n['id']:(n['x'],n['y']) for n in model['nodes'] if n['id'] in mcc}
    values=[degrees[n] for n in mcc];norm=Normalize(vmin=0,vmax=max([1]+values))
    nx.draw_networkx_edges(mcc,positions,ax=ax,edge_color='#94A3B8',alpha=.35,width=.8)
    if mcc:
        nx.draw_networkx_nodes(mcc,positions,ax=ax,node_color=values,cmap='viridis',vmin=norm.vmin,vmax=norm.vmax,
            node_size=20,edgecolors='white',linewidths=.4)
        if len(mcc)<=20:
            offset=max(.1,max((p[1] for p in positions.values()),default=0)-min((p[1] for p in positions.values()),default=0))*.04
            nx.draw_networkx_labels(mcc,{n:(p[0],p[1]+offset) for n,p in positions.items()},ax=ax,font_size=8,font_color='#172B4D')
    else:ax.text(.5,.5,'No compounds in MCC',ha='center',transform=ax.transAxes)
    colorbar=fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),ax=ax,label='Node degree',shrink=.8)
    colorbar.locator=MaxNLocator(integer=True);colorbar.update_ticks()
    ax.set_title(f'MCC · {len(mcc)} compounds · {mcc.number_of_edges()} relationships',loc='left',color='#172B4D',fontsize=13);ax.set_axis_off()
    fragment=model.get('fragment',{});fragment_ax.set_axis_off()
    if fragment.get('image'):
        fragment_ax.imshow(imread(io.BytesIO(base64.b64decode(fragment['image'])),format='png'))
        qualifier='Best match · time limit reached' if fragment['status']=='partial' else 'Maximum common fragment'
        fragment_ax.set_title(qualifier+'\n'+str(fragment['atoms'])+' atoms · '+str(fragment['bonds'])+' bonds',fontsize=10,color='#172B4D')
    else:
        fragment_ax.text(.5,.5,'Common fragment unavailable\n'+fragment.get('status','not calculated'),ha='center',va='center',transform=fragment_ax.transAxes,color='#64748B')
    sequence=sorted(values,reverse=True)
    rank_ax.plot(range(1,len(sequence)+1),sequence,color='#2552E8',lw=2);rank_ax.fill_between(range(1,len(sequence)+1),sequence,color='#2552E8',alpha=.1)
    rank_ax.set(title='MCC degree rank',xlabel='Rank',ylabel='Degree')
    counts=Counter(sequence);hist_ax.bar(list(counts),list(counts.values()),color='#21918C',width=.8)
    hist_ax.set(title='MCC degree histogram',xlabel='Degree',ylabel='Compounds')
    for current,label,color in ((graph,'Full graph','#94A3B8'),(mcc,'MCC','#2552E8')):
        distribution=Counter(dict(current.degree()).values());x=list(range(max(distribution,default=0)+1)) if distribution else []
        dist_ax.plot(x,[distribution[d] for d in x],marker='o',color=color,label=label)
    dist_ax.set(title='Full graph and MCC · degree distribution',xlabel='Degree',ylabel='Compounds');dist_ax.legend(frameon=False)
    for item in (rank_ax,hist_ax,dist_ax):
        item.grid(axis='y',alpha=.15);item.set_axisbelow(True)
        item.xaxis.set_major_locator(MaxNLocator(integer=True));item.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.suptitle(textwrap.fill(model['title'],90),fontsize=13,color='#172B4D',y=.97)
    out=io.BytesIO();fig.savefig(out,format='png',dpi=180)
    return out.getvalue()


def export_analysis(dataset,edges,output,identifier,title,origin=None,mcs_timeout=30,ring_matches_ring_only=True,complete_rings_only=False):
    output=Path(output);graph=build_graph(dataset,edges)
    model=graph_view(graph,dataset,title)
    mcc=graph.subgraph(model['mcc']);selected=dataset[dataset['molecule_chembl_id'].isin(model['mcc'])].copy()
    report_progress('Buscando o fragmento comum do MCC: '+title+'…')
    model['fragment']=common_fragment(selected['canonical_smiles'].tolist(),mcs_timeout,ring_matches_ring_only,complete_rings_only)
    model['origin']=origin or {};model['analysis_id']=identifier
    model['statistics']={'nodes':len(graph),'edges':graph.number_of_edges(),'components':nx.number_connected_components(graph),
        'density':nx.density(graph),'mcc_nodes':len(mcc),'mcc_edges':mcc.number_of_edges(),'mcc_density':nx.density(mcc)}
    selected['degree']=selected['molecule_chembl_id'].map(dict(mcc.degree()))
    for folder in ('plots','Molecules','data/maxcomp','centroids'):(output/folder).mkdir(parents=True,exist_ok=True)
    write_view(output/'plots'/(identifier+SUFFIX),model)
    selected.to_csv(output/'Molecules'/(identifier+'.csv'),index=False)
    pd.DataFrame([{'source':a,'target':b,'value':attrs['value']} for a,b,attrs in mcc.edges(data=True)],
        columns=['source','target','value']).to_csv(output/'data/maxcomp'/(identifier+'.csv'),index=False)
    (output/'plots'/(identifier+'.png')).write_bytes(report_png(model))
    if model['fragment']['image']:(output/'centroids'/(identifier+'.png')).write_bytes(base64.b64decode(model['fragment']['image']))
    return model,selected


def folder_inputs(compounds_path,fingerprints_path,similarity_path,metric,fingerprint):
    compound_files=[str(p) for p in sorted(Path(compounds_path).glob('*.csv')) if {'canonical_smiles'}<=columns(p)] if compounds_path else []
    result=[]
    for kind,folder in [('fingerprints',fingerprints_path),('similarity',similarity_path or (str(Path(compounds_path)/'Similarity') if compounds_path else None))]:
        if folder is None:continue
        required={'fingerprint','molecule_chembl_id'} if kind=='fingerprints' else {'source','target','value'}
        for path in sorted(Path(folder).glob('*.csv')):
            if required<=columns(path):
                # Preserve legacy per-dataset matching when a collection folder
                # contains independent datasets with isolated compounds.
                stem=path.stem
                prefixes=[getattr(metric,'value',metric)+'_'+getattr(fingerprint,'value',fingerprint)+'_',getattr(fingerprint,'value',fingerprint)+'_']
                name=next((stem.removeprefix(prefix) for prefix in prefixes if stem.startswith(prefix)),stem)
                matching=[p for p in compound_files if Path(p).stem==name]
                result.append({'kind':kind,'file':str(path),'compound_files':matching or compound_files,'label':path.name})
    if not result:raise ValueError('Selecione uma entrada de fingerprints ou de similaridade pronta para gerar os grafos.')
    return result


def analyze_inputs(entries,output,metric='Tanimoto',fingerprint='morgan',threshold=70,mcs_timeout=30,ring_matches_ring_only=True,complete_rings_only=False):
    output=Path(output);selected_frames=[]
    quality=CURRENT_REPORT.get() or QualityReport('graphs')
    for index,entry in enumerate(entries):
        path=Path(entry['file']);kind=entry['kind'];report_progress(f'Análise de grafos {index+1}/{len(entries)} · {path.name}')
        validate_file(path,kind,validate_rows=False)
        dataset=read_compounds(entry.get('compound_files',[])+([str(path)] if 'canonical_smiles' in columns(path) else []),quality)
        if kind=='fingerprints':
            edges,ids=exact_edges(path,metric,threshold,quality)
            dataset=dataset[dataset['molecule_chembl_id'].isin(ids)]
            absent=set(ids)-set(dataset['molecule_chembl_id'])
            for code in sorted(absent):quality.exclude(path,None,[code],'fingerprint sem molécula correspondente na entrada')
        else:
            edges=pd.read_csv(path,dtype={'source':str,'target':str})
            edges=clean_dataframe(edges,'similarity',path,quality)
            if not entry.get('compound_files') and dataset.empty:
                # External edge tables are sufficient for topology. Structures
                # and isolated nodes require the optional compound table.
                identifiers=sorted(set(edges['source'])|set(edges['target']))
                dataset=pd.DataFrame({'molecule_chembl_id':identifiers,'canonical_smiles':[None]*len(identifiers)})
        before=len(quality.records)
        edges=filter_edges(edges,dataset['molecule_chembl_id'],path,quality)
        base=re.sub('[^A-Za-z0-9_.+-]+','_',path.stem)[:90]
        identifier=f'{index+1:03d}_{kind}_{base}'
        origin={'kind':kind,'file':path.name,'metric':getattr(metric,'value',metric) if kind=='fingerprints' else entry.get('metric','pronta'),
            'fingerprint':entry.get('fingerprint',getattr(fingerprint,'value',fingerprint)),
            'threshold':threshold if kind=='fingerprints' else None,'label':entry.get('label',path.name),
            'excluded_relationships':len(quality.records)-before}
        title=entry.get('label',path.name)+' · '+('similaridade '+str(origin['metric']) if kind=='fingerprints' else 'similaridade pronta')
        model,selected=export_analysis(dataset,edges,output,identifier,title,origin,mcs_timeout,ring_matches_ring_only,complete_rings_only)
        if kind=='fingerprints':
            (output/'data/similarity').mkdir(parents=True,exist_ok=True)
            edges.to_csv(output/'data/similarity'/(identifier+'.csv'),index=False)
        selected_frames.append(selected)
    # Keep the legacy convenience union for existing downstream pipelines. Each
    # analysis also has its own CSV, view and figure and can be selected directly.
    combined=pd.concat(selected_frames,ignore_index=True) if selected_frames else pd.DataFrame(columns=['molecule_chembl_id','canonical_smiles'])
    (output/'Molecules').mkdir(parents=True,exist_ok=True)
    # Independent uploads can reuse an identifier for different molecules. The
    # convenience union must preserve those structures without ambiguous IDs.
    if not combined.empty:
        combined['_structure']=combined['canonical_smiles'].map(
            lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s)) if isinstance(s,str) and s else '')
        conflicts=set(combined.groupby('molecule_chembl_id')['_structure'].nunique().loc[lambda n:n>1].index)
        if conflicts:
            import hashlib
            combined['original_code']=combined['molecule_chembl_id']
            for row in combined.index:
                if combined.at[row,'molecule_chembl_id'] in conflicts:
                    code=combined.at[row,'molecule_chembl_id'][:75]
                    combined.at[row,'molecule_chembl_id']=code+'_'+hashlib.sha256(combined.at[row,'_structure'].encode()).hexdigest()[:12]
        combined=combined.drop(columns='_structure')
    combined.drop_duplicates('molecule_chembl_id').to_csv(output/'Molecules/molecules.csv',index=False)
    quality.write(output)
