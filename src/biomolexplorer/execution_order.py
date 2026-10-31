"""Durable scheduling priorities and dependency-preserving group movement."""
import json
import time
from .bindings import dependencies


def ordered(stages,preference=()):
    by_id={s['id']:s for s in stages}
    priority=list(dict.fromkeys([identifier for identifier in preference if identifier in by_id]+list(by_id)))
    pending={identifier:dependencies(by_id[identifier]) for identifier in priority};result=[];completed=set()
    while pending:
        identifier=next((identifier for identifier,parents in pending.items() if parents<=completed),None)
        if identifier is None:raise ValueError('Há um ciclo ou dependência ausente na ordem de execução.')
        result.append(identifier);completed.add(identifier);pending.pop(identifier)
    return result


def move(stages,order,identifier,direction):
    if direction not in (-1,1):raise ValueError('Movimento inválido na ordem de execução.')
    order=ordered(stages,order)
    if identifier not in order:raise ValueError('O bloco não pertence à ordem de execução.')
    carried={identifier}
    while True:
        grown=carried|{s['id'] for s in stages if dependencies(s)&carried}
        if grown==carried:break
        carried=grown
    group=[key for key in order if key in carried];remaining=[key for key in order if key not in carried]
    position=sum(key not in carried for key in order[:order.index(identifier)])
    destination=position+direction
    if not 0<=destination<=len(remaining):return order
    proposed=remaining[:destination]+group+remaining[destination:]
    indexes={key:index for index,key in enumerate(proposed)}
    if any(indexes[parent]>=indexes[s['id']] for s in stages for parent in dependencies(s)):
        return order
    return proposed


def read(store,token,pid):
    project=store.project(token,pid)
    with store.connect() as db:
        row=db.execute('SELECT stage_ids,confirmed FROM execution_preferences WHERE project_id=?',(pid,)).fetchone()
        last=db.execute('SELECT stages FROM runs WHERE project_id=? ORDER BY created DESC LIMIT 1',(pid,)).fetchone()
    preference=json.loads(row[0]) if row else [s['id'] for s in json.loads(last[0])] if last else []
    return dict(order=ordered(project['pipeline'],preference),confirmed=bool(row[1]) if row else bool(last))


def save(store,token,pid,order):
    project=store.project(token,pid,'editor');expected={s['id'] for s in project['pipeline']}
    if not isinstance(order,list) or len(order)!=len(set(order)) or set(order)!=expected:
        raise ValueError('A ordem deve incluir exatamente os blocos atuais do pipeline.')
    if ordered(project['pipeline'],order)!=order:raise ValueError('As dependências devem executar antes dos blocos que as utilizam.')
    with store.connect() as db:
        db.execute('BEGIN IMMEDIATE')
        current=store.project(token,pid,'editor')
        if current['revision']!=project['revision']:raise ValueError('O pipeline foi alterado. Reabra a organização da execução.')
        previous=db.execute('SELECT stage_ids,confirmed FROM execution_preferences WHERE project_id=?',(pid,)).fetchone()
        if previous and json.loads(previous[0])==order and previous[1]:return dict(order=order,confirmed=True)
        from .project_state import snapshot,record
        before=snapshot(store,db,pid)
        db.execute('INSERT INTO execution_preferences VALUES(?,?,1,?) ON CONFLICT(project_id) DO UPDATE SET stage_ids=excluded.stage_ids,confirmed=1,updated=excluded.updated',
            (pid,json.dumps(order),time.time()))
        record(store,db,pid,store.user(token)['id'],'execution_order','Ordem de execução atualizada.',before)
    return dict(order=order,confirmed=True)
