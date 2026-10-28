"""Human-readable dataset origins from exact, transitive artifact provenance."""
from pathlib import Path

from .artifacts import records


class ResultOrigins:
    def __init__(self, stages, run_id=None, load_run=None):
        self.load_run = load_run
        self.loaded = {run_id}
        self.catalog = {}
        self.by_path = {}
        self.memo = {}
        self.stages = {}
        self.add(stages, run_id)

    def add(self, stages, run_id):
        for stage in stages:
            self.stages[(run_id, stage['id'])] = stage
            for record in records(stage, run_id):
                key = (record['run_id'], record['stage_id'], record['path'])
                self.catalog[key] = (record, stage)
                self.by_path.setdefault(record['path'], []).append(key)

    def find(self, reference):
        path = reference.get('path')
        def matches():
            return [key for key in self.by_path.get(path, [])
                    if not reference.get('stage_id') or key[1] == reference['stage_id']]
        candidates = matches()
        run_id = reference.get('run_id')
        exact = [key for key in candidates if key[0] == run_id]
        if exact:
            return exact[0]
        if self.load_run and run_id and run_id not in self.loaded:
            self.loaded.add(run_id)
            run = self.load_run(run_id)
            if run:
                self.add(run['stages'], run_id)
            candidates = matches()
            exact = [key for key in candidates if key[0] == run_id]
            if exact:
                return exact[0]
        return candidates[0] if len(candidates) == 1 else None

    def parents(self, record, stage):
        inputs = record.get('inputs') or {}
        if not inputs:
            batch = next((b for b in stage.get('batches', []) if record['path'] in b.get('artifacts', [])), {})
            inputs = {field: [dict(path=str(p), name=Path(p).name) for p in paths]
                      for field, paths in (batch or stage).get('input_files', {}).items()}
        if record.get('operation') == 'graphs' and inputs.get('similarity_path'):
            # Molecular annotation files do not change the origin of the edges.
            inputs = {'similarity_path': inputs['similarity_path']}
        return [ref for rows in inputs.values() for ref in rows]

    def routes(self, key, seen=frozenset()):
        if key in seen:
            return []
        if key in self.memo:
            return self.memo[key]
        record, stage = self.catalog[key]
        config = stage.get('configuration', stage)
        params = config.get('parameters', {})
        step = dict(stage_id=record['stage_id'], stage_name=record['stage_name'],
                    operation=record['operation'], name=record['name'], path=record['path'],
                    run_id=record['run_id'], identifier=params.get('search_term') or
                    record.get('context', {}).get('search_term') or '',
                    target=record.get('context', {}).get('target') or params.get('target') or '')
        parents = self.parents(record, stage)
        paths = []
        for parent in parents:
            parent_key = self.find(parent)
            if parent_key:
                paths.extend(route + [step] for route in self.routes(parent_key, seen | {key}))
            else:
                leaf = dict(stage_id=parent.get('stage_id', parent.get('asset_id', '')),
                            stage_name=parent.get('stage_name', parent.get('name', '')),
                            operation='external', name=parent.get('name', Path(parent.get('path', '')).name),
                            path=parent.get('path', ''), run_id=parent.get('run_id'), identifier='')
                paths.append([leaf, step])
        if not paths:
            paths = [[step]]
        unique = {}
        for route in paths:
            unique[tuple((s['stage_id'], s['path'], s['run_id']) for s in route)] = route
        self.memo[key] = list(unique.values())
        return self.memo[key]

    def describe(self, path, stage_id=None, run_id=None):
        key = self.find(dict(path=str(path), stage_id=stage_id, run_id=run_id))
        if key is None:
            return dict(source_label=Path(path).name, origin=Path(path).name,
                        display_name=Path(path).name, source_trace='', source_routes=[])
        routes = self.routes(key)
        roots = []; branches = []; traces = []
        for route in routes:
            root = route[0]
            label = root['stage_name'] or root['name']
            if root['identifier'] and root['identifier'] != label:
                label += ' · ' + str(root['identifier'])
            if label not in roots:
                roots.append(label)
            branch = ' → '.join(s['stage_name'] + (' · ' + s['name'] if s['operation'] == 'admet' else '')
                                for s in route[1:-1])
            if branch and branch not in branches:
                branches.append(branch)
            traces.append(' → '.join(s['stage_name'] + ' / ' + s['name'] for s in route))
        source = ' + '.join(roots)
        origin = source + (' — ' + ' + '.join(branches) if branches else '')
        scope=''
        if self.catalog[key][0]['operation']=='graphs':
            from .graph_exports import subset_label
            scope=subset_label(path)
        return dict(source_label=source, origin=origin, display_name=source + ' — ' + (scope+' — ' if scope else '') + Path(path).name,
                    source_trace='\n'.join(dict.fromkeys(traces)), source_routes=routes)
