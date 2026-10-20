"""Deterministic CPU layouts with bounded work and no node sampling.

ForceAtlas2 force laws and adaptive speed follow Jacomy et al. (2014),
https://doi.org/10.1371/journal.pone.0098679 . This implementation uses exact
repulsion, not Gephi's Barnes–Hut acceleration or its interactive engine.
"""
import math
import numpy as np
import networkx as nx


def forceatlas2(graph, linlog=False, iterations=80):
    nodes=sorted(graph);n=len(nodes)
    if n<2:return {node:(0.,0.) for node in nodes}
    if n>1500:raise ValueError('ForceAtlas2 CPU: limite de 1.500 nós por componente.')
    index={node:i for i,node in enumerate(nodes)}
    mass=np.array([graph.degree(node)+1. for node in nodes])
    points=np.random.default_rng(42).uniform(-1,1,(n,2))*math.sqrt(n)
    edges=np.array([(index[a],index[b]) for a,b in graph.edges() if a!=b],dtype=int).reshape(-1,2)
    previous=np.zeros_like(points);speed=1.;efficiency=1.
    for _ in range(iterations):
        delta=points[:,None,:]-points[None,:,:]
        distance2=np.einsum('ijk,ijk->ij',delta,delta)
        np.fill_diagonal(distance2,np.inf)
        repulsion=2*mass[:,None]*mass[None,:]/np.maximum(distance2,1e-9)
        force=np.einsum('ij,ijk->ik',repulsion,delta)
        distance=np.linalg.norm(points,axis=1)
        force-=points*(mass/np.maximum(distance,1e-9))[:,None]
        if len(edges):
            a,b=edges.T;attraction=points[a]-points[b]
            if linlog:
                length=np.linalg.norm(attraction,axis=1)
                attraction*= (np.log1p(length)/np.maximum(length,1e-9))[:,None]
            np.add.at(force,a,-attraction);np.add.at(force,b,attraction)
        swinging=np.linalg.norm(force-previous,axis=1)
        traction=np.linalg.norm(force+previous,axis=1)/2
        total_swing=float(mass@swinging);total_traction=float(mass@traction)
        jitter=max(.05,min(10.,.05*math.sqrt(n)))
        if total_swing>jitter*total_traction:efficiency=max(.05,efficiency*.7)
        elif speed<1000:efficiency=min(1.,efficiency*1.3)
        target=jitter*jitter*efficiency*total_traction/max(total_swing,1e-9)
        speed=max(1e-6,speed+min(target-speed,.5*speed))
        step=force*(speed/(1+np.sqrt(speed*swinging)))[:,None]
        # Bound a single integration step, including nearly coincident nodes.
        step*=np.minimum(1.,10/np.maximum(np.linalg.norm(step,axis=1),1e-9))[:,None]
        points+=step;previous=force
    points-=points.mean(axis=0)
    points/=max(float(np.abs(points).max()),1e-9)
    return {node:tuple(points[i]) for i,node in enumerate(nodes)}


def packed_force(graph,kind):
    """Run forces per component; pack all boxes including cost-limited rings."""
    from .visualizations import separated_graph_layout,_separate_positions
    groups=sorted(nx.connected_components(graph),key=lambda c:(-len(c),min(c)))
    boxes=[];budget=8_000_000;fallback=0
    for members in groups:
        sub=graph.subgraph(sorted(members));n=len(sub)
        cost=n*n
        if n<2:positions={next(iter(sub)):(0.,0.)}
        elif n<=1500 and sub.number_of_edges()<=20000 and cost<=budget:
            budget-=cost
            raw=(nx.spring_layout(sub,seed=42,weight=None,iterations=60) if kind=='fr_components'
                 else forceatlas2(sub,linlog=kind=='fa2_linlog'))
            positions=_separate_positions({k:(float(p[0])*math.sqrt(n),float(p[1])*math.sqrt(n)) for k,p in raw.items()},sub)
        else:
            positions=separated_graph_layout(sub,circular=True);fallback+=1
        xmin=min(p[0] for p in positions.values());ymin=min(p[1] for p in positions.values())
        width=max(p[0] for p in positions.values())-xmin+3
        height=max(p[1] for p in positions.values())-ymin+3
        boxes.append((positions,xmin,ymin,width,height))
    limit=max([b[3] for b in boxes]+[math.sqrt(sum(b[3]*b[4] for b in boxes)*1.6)])
    result={};x=y=row=0.
    for positions,xmin,ymin,width,height in boxes:
        if x and x+width>limit:x=0.;y+=row;row=0.
        result.update({k:(p[0]-xmin+x+1.5,p[1]-ymin+y+1.5) for k,p in positions.items()})
        x+=width;row=max(row,height)
    method='Fruchterman–Reingold · 60 iterações' if kind=='fr_components' else 'ForceAtlas2 CPU · 80 iterações · repulsão exata'+(' · LinLog' if kind=='fa2_linlog' else '')
    notice=method+' · por componente · semente 42 · sem pesos de similaridade.'
    if fallback:notice+=f' {fallback} componentes em anéis por limite de custo; todos os nós preservados.'
    return result,notice
