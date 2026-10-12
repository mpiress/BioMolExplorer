"""Source identities for read-only dashboards, honoring persisted provenance."""
import hashlib
import json
from pathlib import Path
from .artifact_choices import logical_path


def records(stage,run_id):
    previous={record['path']:record for record in stage.get('artifact_records') or []}
    params=stage.get('configuration',{}).get('parameters',stage.get('parameters',{}))
    for path in stage.get('artifacts',[]):
        prior=previous.get(path,{})
        parts=Path(logical_path(path)).parts
        if stage['id'] in parts:
            index=len(parts)-1-list(reversed(parts)).index(stage['id'])
            locator='/'.join(parts[index+1:])
        else:locator=logical_path(path)
        native_run=None
        if stage['id'] in parts and 'runs' in parts[:index]:
            run_index=len(parts[:index])-1-list(reversed(parts[:index])).index('runs')
            candidate=parts[run_index+1]
            if len(candidate)==32 and all(c in '0123456789abcdef' for c in candidate):native_run=candidate
        occurrence=prior.get('run_id') or native_run or stage.get('reused_from_run') or run_id
        identity=prior.get('artifact_id') or hashlib.sha256(json.dumps([stage['id'],occurrence,locator]).encode()).hexdigest()
        yield dict(artifact_id=identity,path=path,name=Path(path).name,stage_id=stage['id'],
            stage_name=stage.get('name',stage['id']),run_id=occurrence,
            context=prior.get('context') or stage.get('artifact_context') or {'configured_target':params.get('target'),'search_term':params.get('search_term')})
