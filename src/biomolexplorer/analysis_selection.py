"""Follow growing generated datasets without replacing intentional subsets."""
import copy
from pathlib import Path

from .artifact_choices import matches_selector
from .artifacts import records, reference
from .bindings import sources, pack
from .fingerprint_selection import filename_kind

FIELDS = {'similarity': ('base_input_path', 'fingerprints'), 'graphs': ('similarity_path', 'similarity')}


def candidates(producer, consumer):
    operation = consumer['operation']
    paths = [p for p in producer.get('artifacts', []) if Path(p).suffix == '.csv']
    if operation == 'similarity':
        kind = consumer['parameters'].get('fingerprint', 'morgan')
        kind = getattr(kind, 'value', kind)
        return [p for p in paths if filename_kind(p) == kind]
    # Similarity workers publish edges in Similarity; metadata/report files
    # are not additional graph experiments.
    return [p for p in paths if Path(p).parent.name == 'Similarity']


def covers(refs, producer, consumer):
    files = candidates(producer, consumer)
    if not files:return False
    catalog = {r['path']: r for r in records(producer)}
    return all(any(ref.get('artifact_id') == catalog[path]['artifact_id'] if ref.get('artifact_id') else
                   ref.get('selector', 'auto') != 'auto' and matches_selector(path, ref['selector'])
                   for ref in refs) for path in files)


def expand_selection(stage, completed, histories=()):
    """Migrate an old complete selection, or honor an explicit follow choice."""
    if stage['operation'] not in FIELDS or stage.get('provided_results') or stage.get('follow_inputs') is False:
        return stage, []
    field, operation = FIELDS[stage['operation']]
    group = stage.get('bindings', {}).get(field)
    refs = sources(group) if group else []
    lookup = {s['id']: s for s in completed}
    replacements = {}; changes = []
    for producer_id in dict.fromkeys(r['stage'] for r in refs if 'stage' in r):
        producer = lookup.get(producer_id)
        if not producer or producer['operation'] != operation or producer.get('provided'):
            continue
        selected = [r for r in refs if r.get('stage') == producer_id]
        # Unconfigured sources still go through the normal input confirmation.
        if any(r.get('selector', 'auto') == 'auto' for r in selected):continue
        follow = stage.get('follow_inputs') is True
        if not follow:
            for history in histories:
                old = next((s for s in history if s['id'] == stage['id']), None)
                parent = next((s for s in history if s['id'] == producer_id and s['status'] == 'succeeded'), None)
                config = old.get('configuration', {}) if old else {}
                prior_group = config.get('bindings', {}).get(field)
                prior = [r for r in sources(prior_group) if r.get('stage') == producer_id] if prior_group else []
                if config.get('follow_inputs') is False:break
                if parent and prior == selected and covers(prior, parent, stage):
                    follow = True;break
        if not follow:continue
        paths = candidates(producer, stage)
        if not paths:continue
        chosen_paths=set(paths)
        new = [reference(r) for r in records(producer) if r['path'] in chosen_paths]
        if new != selected:
            replacements[producer_id] = new
            changes.append(dict(stage_id=producer_id, previous=len(selected), current=len(new)))
    if not changes:return stage, []
    result = copy.deepcopy(stage);expanded = [];seen = set()
    for ref in refs:
        identity = ref.get('stage')
        if identity in replacements:
            if identity not in seen:expanded.extend(replacements[identity]);seen.add(identity)
        else:expanded.append(ref)
    result['bindings'][field] = pack(expanded)
    # A single generated source is the common analysis chain. Persist its
    # inferred intent so future expansions need no historical inference.
    if len({r.get('stage') for r in refs}) == 1:result['follow_inputs'] = True
    return result, changes
