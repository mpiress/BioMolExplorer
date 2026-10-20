"""Component-aware unweighted centralities, with explicit work budgets.

Local path metrics are normalized against the requested universe, never the
component size by accident. SciPy batches bound shortest-path matrix memory.
"""
import networkx as nx

PATH_WORK=200_000_000
BRANDES_WORK=20_000_000
TRIANGLE_WORK=20_000_000


def _tree(graph):
    root=next(iter(graph));parent={root:None};order=[root];depth={root:0}
    for node in order:
        for neighbor in graph[node]:
            if neighbor not in parent:
                parent[neighbor]=node;depth[neighbor]=depth[node]+1;order.append(neighbor)
    sizes=dict.fromkeys(graph,1)
    for node in reversed(order[1:]):sizes[parent[node]]+=sizes[node]
    return parent,order,depth,sizes


def betweenness(graph,components):
    n=len(graph);values=dict.fromkeys(graph,0.);budget=BRANDES_WORK
    plans=[];pivots=approximate=0
    for members in components:
        sub=graph.subgraph(members);size=len(sub);edges=sub.number_of_edges()
        if size<=2:continue
        if edges==size*(size-1)//2 and not nx.number_of_selfloops(sub):continue
        if edges==size-1 and not nx.number_of_selfloops(sub):
            plans.append((sub,'tree',None));continue
        step=size+2*edges
        if size<=400 and size*step<=min(2_000_000,budget):k=None;budget-=size*step
        else:
            k=min(64,size,budget//max(1,step))
            if k<min(8,size):
                return {},f'Betweenness: orçamento de {BRANDES_WORK} visitas excedido; selecione um componente ou comunidade para reduzir o custo.'
            budget-=k*step;pivots+=k;approximate+=1
        plans.append((sub,'brandes',k))
    scale=2/((n-1)*(n-2)) if n>2 else 0.
    for sub,method,k in plans:
        if method=='tree':
            parent,order,_,sizes=_tree(sub);size=len(sub)
            for node in order:
                branches=[sizes[other] if parent.get(other)==node else size-sizes[node] for other in sub[node]]
                values[node]=((size-1)**2-sum(s*s for s in branches))/2*scale
        else:
            raw=nx.betweenness_centrality(sub,k=k,normalized=False,weight=None,seed=42)
            values.update({node:value*scale for node,value in raw.items()})
    notice=(f'Betweenness aproximada: {pivots} pivôs em {approximate} componentes, até 64 por componente; semente 42. Demais componentes exatos.'
            if approximate else 'Betweenness exata por componente.')
    return values,notice


def closeness(graph,components):
    n=len(graph);values=dict.fromkeys(graph,0.)
    work=0
    for members in components:
        size=len(members);edges=graph.subgraph(members).number_of_edges()
        if edges not in (size-1,size*(size-1)//2):work+=size*(size+2*edges)
    if work>PATH_WORK:
        return {},f'Closeness: custo estimado {work} acima do orçamento {PATH_WORK}; selecione um componente ou comunidade. Nenhuma aproximação foi aplicada.'
    from scipy.sparse.csgraph import shortest_path
    import numpy as np
    for members in components:
        sub=graph.subgraph(members);size=len(sub)
        if size<2:continue
        correction=(size-1)/max(1,n-1)
        if sub.number_of_edges()==size*(size-1)//2 and not nx.number_of_selfloops(sub):
            values.update(dict.fromkeys(sub,correction))
        elif sub.number_of_edges()==size-1:
            parent,order,depth,sizes=_tree(sub);totals={order[0]:sum(depth.values())}
            for node in order[1:]:totals[node]=totals[parent[node]]+size-2*sizes[node]
            values.update({node:(size-1)/total*correction for node,total in totals.items()})
        else:
            nodes=sorted(sub)
            matrix=nx.to_scipy_sparse_array(sub,nodelist=nodes,weight=None,dtype=float,format='csr')
            # SciPy csgraph in the supported environment requires int32 indices.
            matrix.indices=matrix.indices.astype(np.int32);matrix.indptr=matrix.indptr.astype(np.int32)
            for start in range(0,size,128):
                distances=shortest_path(matrix,directed=False,unweighted=True,indices=np.arange(start,min(start+128,size)),method='D')
                totals=distances.sum(axis=1)
                for node,total in zip(nodes[start:start+128],totals):values[node]=(size-1)/float(total)*correction
    return values,'Closeness exata por componente · correção Wasserman–Faust · caminhos em lotes de até 128 origens.'


def triangles(graph):
    neighbors={node:set(graph[node])-{node} for node in graph}
    degree={node:len(adj) for node,adj in neighbors.items()}
    forward={node:{other for other in adj if (degree[node],node)<(degree[other],other)} for node,adj in neighbors.items()}
    work=sum(min(len(adj),len(forward[other])) for node,adj in forward.items() for other in adj)
    if work>TRIANGLE_WORK:
        return {},f'Clustering / triângulos: custo estimado {work} acima do orçamento {TRIANGLE_WORK}; selecione um componente ou comunidade.'
    values=dict.fromkeys(graph,0)
    for node,adj in forward.items():
        for other in adj:
            for third in adj & forward[other]:
                values[node]+=1;values[other]+=1;values[third]+=1
    return values,''


def eigenvector(graph):
    if not nx.is_connected(graph):
        return {},'Eigenvector exige um universo conectado; escolha LCC ou um componente no escopo do ranking.'
    if len(graph)<=2:return nx.eigenvector_centrality(graph,max_iter=1000,weight=None),''
    import numpy as np
    from scipy.sparse.linalg import eigsh,ArpackNoConvergence
    nodes=list(graph);matrix=nx.to_scipy_sparse_array(graph,nodelist=nodes,weight=None,dtype=float,format='csr')
    try:
        _,vectors=eigsh(matrix,k=1,which='LA',v0=np.ones(len(nodes)),tol=1e-8,maxiter=2000)
    except ArpackNoConvergence:
        return {},'Eigenvector: o solver esparso não convergiu em 2.000 iterações; selecione um componente menor.'
    vector=np.abs(vectors[:,0]);vector/=np.linalg.norm(vector)
    return dict(zip(nodes,map(float,vector))),'Eigenvector · solver esparso · tolerância 1e-8.'
