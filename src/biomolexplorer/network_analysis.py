"""Read-only topological analysis of persisted molecular relationship networks.

Relationship scores are deliberately not used as shortest-path distances.
All metrics and communities describe unweighted, undirected topology.
"""
import csv
import io
import math
from collections import Counter
import networkx as nx

METRICS={
    'degree':'Degree · número de conexões no universo analisado.',
    'density':'Degree Centrality · grau / (n − 1); depende do universo analisado.',
    'betweenness':'Betweenness · participação nos menores caminhos topológicos; aproximação informada em redes grandes.',
    'closeness':'Closeness · proximidade topológica, com correção de Wasserman–Faust para componentes desconectados.',
    'eigenvector':'Eigenvector · conexões a nós centrais; disponível apenas em universo conectado.',
    'pagerank':'PageRank · importância topológica com amortecimento 0,85; sem pesos de escore.',
    'core':'K-core · maior núcleo de grau mínimo k ao qual o composto pertence.',
    'clustering':'Clustering · fração de conexões possíveis entre os vizinhos.',
    'triangles':'Triangles · número de triângulos que contêm o composto.',
}
PALETTE=['#2563EB','#D97706','#0D9488','#7C3AED','#DB2777','#475569','#65A30D','#0891B2']


class NetworkAnalysis:
    DRAW_NODES=3000
    DRAW_EDGES=12000

    def __init__(self,data):
        self.data=data
        self.graph=nx.Graph()
        self.graph.add_nodes_from(n['id'] for n in data['nodes'])
        for edge in data['edges']:
            self.graph.add_edge(edge['source'],edge['target'],**{k:v for k,v in edge.items() if k not in ('source','target')})
        self.components=sorted(nx.connected_components(self.graph),key=lambda c:(-len(c),min(c)))
        self.lcc=set(self.components[0]) if self.components else set()
        self.component_edges={i+1:self.graph.subgraph(group).number_of_edges() for i,group in enumerate(self.components)}
        self.cache={};self.layouts={};self.community_cache={};self.statistics_cache={}

    def universe(self,mode):
        return self.graph.subgraph(self.lcc) if mode=='mcc' else self.graph

    def metrics(self,mode,metric):
        key=(mode,metric)
        if key in self.cache:return self.cache[key]
        graph=self.universe(mode);n=len(graph);m=graph.number_of_edges();notice=''
        if not n:values={}
        elif metric=='degree':values=dict(graph.degree())
        elif metric=='density':values=nx.degree_centrality(graph) if n>1 else {next(iter(graph)):0.0}
        elif metric=='core':
            simple=graph.copy();simple.remove_edges_from(nx.selfloop_edges(simple))
            values=nx.core_number(simple)
        elif metric in ('clustering','triangles'):
            from .network_centrality import triangles
            values,notice=triangles(graph)
            if metric=='clustering' and values:
                values={node:2*value/(degree*(degree-1)) if (degree:=len(set(graph[node])-{node}))>1 else 0. for node,value in values.items()}
        elif metric=='betweenness':
            from .network_centrality import betweenness
            values,notice=betweenness(graph,list(nx.connected_components(graph)))
        elif metric=='closeness':
            from .network_centrality import closeness
            values,notice=closeness(graph,list(nx.connected_components(graph)))
        elif metric=='eigenvector':
            from .network_centrality import eigenvector
            values,notice=eigenvector(graph)
        elif metric=='pagerank':
            try:
                values=nx.pagerank(graph,alpha=.85,weight=None,max_iter=1000,tol=1e-8)
                notice='PageRank · solver esparso · amortecimento 0,85 · tolerância 1e-8.'
            except nx.PowerIterationFailedConvergence:values={};notice='PageRank não convergiu em 1.000 iterações; selecione um componente menor.'
        else:raise ValueError('Métrica topológica desconhecida.')
        self.cache[key]=(values,notice)
        return values,notice

    def communities(self,mode):
        if mode in self.community_cache:return self.community_cache[mode]
        graph=self.universe(mode)
        if len(graph)>100000 or graph.number_of_edges()>500000:
            result=({},None,'Comunidades indisponíveis acima de 100.000 nós ou 500.000 arestas.')
        else:
            groups=nx.community.louvain_communities(graph,weight=None,seed=42) if graph.number_of_edges() else [{n} for n in graph]
            groups=sorted(groups,key=lambda c:(-len(c),min(c)))
            mapping={node:i+1 for i,c in enumerate(groups) for node in c}
            modularity=nx.community.modularity(graph,groups,weight=None) if graph.number_of_edges() else None
            result=(mapping,modularity,'Louvain sem pesos · resolução 1 · semente 42. Comunidades são grupos topológicos.')
        self.community_cache[mode]=result
        return result

    def statistics(self,mode,calculate=False):
        if mode in self.statistics_cache and (not calculate or self.statistics_cache[mode].get('clustering_calculated')):return self.statistics_cache[mode]
        graph=self.universe(mode);n=len(graph);m=graph.number_of_edges()
        groups=list(nx.connected_components(graph))
        clustering=nx.average_clustering(graph) if calculate and n and m<=20000 and n<=5000 else None
        result=dict(nodes=n,edges=m,density=nx.density(graph),components=len(groups),
            lcc_nodes=max(map(len,groups),default=0),isolates=nx.number_of_isolates(graph),
            mean_degree=2*m/n if n else 0,clustering=clustering,
            degrees=dict(sorted(Counter(dict(graph.degree()).values()).items())),
            component_sizes=dict(sorted(Counter(map(len,groups)).items())))
        result['clustering_calculated']=calculate
        result.update(diameter_lcc=None,mean_path_lcc=None,assortativity=None)
        if calculate and n:
            if len(self.lcc)<=1500 and self.graph.subgraph(self.lcc).number_of_edges()<=20000:
                lcc=graph.subgraph(self.lcc)
                if len(lcc):
                    diameter=total=pairs=0
                    for _,lengths in nx.all_pairs_shortest_path_length(lcc):
                        diameter=max(diameter,max(lengths.values(),default=0))
                        total+=sum(lengths.values());pairs+=len(lengths)-1
                    result.update(diameter_lcc=diameter,mean_path_lcc=total/pairs if pairs else 0.)
            if m and len({d for _,d in graph.degree() if d})>1:
                result['assortativity']=float(nx.degree_assortativity_coefficient(graph,weight=None))
        self.statistics_cache[mode]=result
        return result

    def filtered(self,mode,component=None,community=None,query='',degree=(None,None),
                 metric='degree',centrality=(None,None),score=(None,None),ego=None,hops=1):
        graph=self.universe(mode);allowed=set(graph)
        if component is None and community is None and not query and degree==(None,None) and centrality==(None,None) and score==(None,None) and ego is None:
            return graph
        if component is not None:
            allowed &= self.components[component-1] if 1<=component<=len(self.components) else set()
        if community is not None:
            mapping=self.community_cache.get(mode,({},None,''))[0]
            allowed &= {n for n,c in mapping.items() if c==community}
        if query:allowed={n for n in allowed if query.casefold() in n.casefold()}
        values=self.metrics(mode,metric)[0]
        def inside(value,bounds):return value is not None and all(limit is None or (value>=limit if i==0 else value<=limit) for i,limit in enumerate(bounds))
        allowed={n for n in allowed if inside(graph.degree(n),degree) and
            (centrality==(None,None) or inside(values.get(n),centrality))}
        if ego is not None:
            if ego not in graph:return graph.subgraph([]).copy()
            allowed &= set(nx.single_source_shortest_path_length(graph,ego,cutoff=max(1,min(3,hops))))
        view=graph.subgraph(allowed)
        return nx.subgraph_view(view,filter_edge=lambda a,b:inside(graph[a][b].get('value'),score)) if score!=(None,None) else view

    def layout(self,mode,kind):
        key=(mode,kind)
        if key in self.layouts:return self.layouts[key]
        graph=self.universe(mode);n=len(graph);notice=''
        if kind=='kamada' and (n>400 or graph.number_of_edges()>8000):
            return None,'Kamada–Kawai indisponível acima de 400 nós ou 8.000 arestas.'
        if kind=='force' and (n>1500 or graph.number_of_edges()>20000):
            return None,'Force-directed indisponível acima de 1.500 nós ou 20.000 arestas.'
        if not n:positions={}
        elif kind in ('fr_components','fa2','fa2_linlog'):
            from .graph_layouts import packed_force
            positions,notice=packed_force(graph,kind)
        elif kind=='circular':positions=nx.circular_layout(graph)
        elif kind=='random':
            from .visualizations import _separate_positions
            positions=_separate_positions({k:(float(p[0])*math.sqrt(n)*2,float(p[1])*math.sqrt(n)*2) for k,p in nx.random_layout(graph,seed=42).items()},graph)
        elif kind=='radial':
            from .visualizations import separated_graph_layout
            positions=separated_graph_layout(graph,circular=True)
        elif kind=='force':positions=nx.spring_layout(graph,seed=42,weight=None,iterations=60)
        elif kind=='kamada':positions=nx.kamada_kawai_layout(graph,weight=None)
        elif kind=='hierarchical':
            # BFS shells are an explicit exploratory rule, never a causal hierarchy.
            positions={};offset=0
            for component in sorted(nx.connected_components(graph),key=lambda c:(-len(c),min(c))):
                root=min(component,key=lambda node:(-graph.degree(node),node))
                distances=nx.single_source_shortest_path_length(graph,root)
                counts=Counter()
                for node in sorted(component,key=lambda n:(distances[n],n)):
                    depth=distances[node];positions[node]=(depth*2,offset+counts[depth]);counts[depth]+=1
                offset+=max(counts.values(),default=0)+3
            notice='Hierarquia exploratória BFS: raiz de maior grau por componente; sem direção causal.'
        elif kind=='community':
            mapping,_,notice=self.communities(mode)
            if not mapping:return None,notice
            centers=nx.circular_layout(nx.empty_graph(sorted(set(mapping.values()))))
            positions={}
            for community in sorted(set(mapping.values())):
                members=sorted(n for n in graph if mapping[n]==community)
                cx,cy=centers[community]
                for i,node in enumerate(members):
                    angle=i*2*math.pi/len(members)
                    positions[node]=(cx*5+.8*math.cos(angle),cy*5+.8*math.sin(angle))
        else:raise ValueError('Layout topológico desconhecido.')
        if kind in ('force','kamada','community') and positions:
            from .visualizations import _separate_positions
            scale=max(1,math.sqrt(n))
            positions=_separate_positions({node:(float(point[0])*scale,float(point[1])*scale) for node,point in positions.items()},graph)
        result=({n:(float(p[0]),float(p[1])) for n,p in positions.items()},notice)
        self.layouts[key]=result
        return result

    def aggregate(self,mode):
        mapping,_,notice=self.communities(mode)
        if not mapping:return nx.Graph(),notice
        graph=self.universe(mode);aggregate=nx.Graph()
        counts=Counter(mapping.values())
        aggregate.add_nodes_from((c,dict(size=count,internal_edges=0)) for c,count in counts.items())
        for a,b in graph.edges():
            left,right=mapping[a],mapping[b]
            if left==right:aggregate.nodes[left]['internal_edges']+=1
            else:
                if not aggregate.has_edge(left,right):aggregate.add_edge(left,right,count=0)
                aggregate[left][right]['count']+=1
        return aggregate,notice

    def export_csv(self,mode,metric,scope=None):
        values,notice=self.metrics(mode,metric);mapping=self.community_cache.get(mode,({},None,''))[0]
        stream=io.StringIO();writer=csv.writer(stream)
        writer.writerow(['compound_id','universe','metric','value','community','method_notice'])
        for node in sorted(self.universe(mode)):writer.writerow([node,scope or mode,metric,values.get(node,''),mapping.get(node,''),METRICS[metric]+' '+notice])
        return stream.getvalue().encode('utf-8')
