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
        elif metric=='betweenness':
            if n>5000 or m>50000:values={};notice='Betweenness indisponível acima de 5.000 nós ou 50.000 arestas.'
            else:
                k=None if n<=400 else min(64,n)
                values=nx.betweenness_centrality(graph,k=k,seed=42,weight=None)
                if k:notice=f'Betweenness aproximada: {k} pivôs, semente 42.'
        elif metric=='closeness':
            if n>1500 or m>20000:values={};notice='Closeness indisponível acima de 1.500 nós ou 20.000 arestas.'
            else:values=nx.closeness_centrality(graph,wf_improved=True)
        elif metric=='eigenvector':
            if not nx.is_connected(graph):values={};notice='Eigenvector exige um universo conectado; escolha LCC.'
            elif m>50000:values={};notice='Eigenvector indisponível acima de 50.000 arestas.'
            else:
                try:values=nx.eigenvector_centrality(graph,max_iter=300,weight=None)
                except nx.PowerIterationFailedConvergence:values={};notice='Eigenvector não convergiu em 300 iterações.'
        elif metric=='pagerank':
            if m>50000 or n>20000:values={};notice='PageRank indisponível acima de 20.000 nós ou 50.000 arestas.'
            else:
                try:values=nx.pagerank(graph,alpha=.85,weight=None,max_iter=200)
                except nx.PowerIterationFailedConvergence:values={};notice='PageRank não convergiu em 200 iterações.'
        else:raise ValueError('Métrica topológica desconhecida.')
        self.cache[key]=(values,notice)
        return values,notice

    def communities(self,mode):
        if mode in self.community_cache:return self.community_cache[mode]
        graph=self.universe(mode)
        if len(graph)>10000 or graph.number_of_edges()>50000:
            result=({},None,'Comunidades indisponíveis acima de 10.000 nós ou 50.000 arestas.')
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
        elif kind=='circular':positions=nx.circular_layout(graph)
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

    def export_csv(self,mode,metric):
        values,notice=self.metrics(mode,metric);mapping=self.community_cache.get(mode,({},None,''))[0]
        stream=io.StringIO();writer=csv.writer(stream)
        writer.writerow(['compound_id','universe','metric','value','community','method_notice'])
        for node in sorted(self.universe(mode)):writer.writerow([node,mode,metric,values.get(node,''),mapping.get(node,''),METRICS[metric]+' '+notice])
        return stream.getvalue().encode('utf-8')
