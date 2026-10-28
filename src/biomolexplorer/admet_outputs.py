"""Explicit ADMET output roles and a narrowly scoped legacy-name migration."""
import copy
import json
from pathlib import Path
import re

SCHEMA_VERSION=4


def export_standard(evaluated, configuration, destination, stem):
    """Partition each endpoint independently; filtering one must not erase another."""
    from .admet_config import selection_reasons
    destination=Path(destination);outputs={}
    for endpoint in ('BBB','HIA'):
        for sign in ('+','-'):
            role=endpoint+sign;name=stem+'_'+role+'.csv'
            evaluated[evaluated[endpoint]==role].to_csv(destination/name,index=False)
            outputs[role]=name
    rules=configuration['rule_set']
    if rules!='none':
        label={'oral':'lipinski-veber'}.get(rules,rules)
        config=copy.deepcopy(configuration)
        config.update(bbb_target='any',hia_target='any',alert_policy='annotate',property_limits={})
        selected=[not selection_reasons(row,config) for row in evaluated.to_dict('records')]
        role='rules-'+label;name=stem+'_'+role+'.csv'
        evaluated.loc[selected].to_csv(destination/name,index=False);outputs[role]=name
    alerts=evaluated.copy()
    for catalog in configuration['alert_catalogs']:
        alerts['Alert_'+catalog]=alerts['StructuralAlerts'].fillna('').map(lambda value:catalog+': ' in str(value))
    name=stem+'_alerts.csv';alerts.to_csv(destination/name,index=False);outputs['alerts']=name
    return outputs


def legacy_reference(ref, files):
    """Map only a proven old ADMET slot to one declared role in the same batch.

    No filename guessing across batches, producers, datasets or rule sets. An
    unknown/forged identity and an ambiguous/missing output remain errors.
    """
    from .artifacts import locator,slot_id
    selector=Path(str(ref.get('selector','')))
    match=re.fullmatch(r'B\d+__admet__(bbb-pos|bbb-neg|hia-pos|hia-neg|rules-(?:lipinski|veber|lipinski-veber|fragment)|alerts)\.csv',selector.name)
    if not match or not ref.get('stage'):return None
    parent=selector.parent.as_posix()
    if parent=='.':return None  # The old filename alone cannot identify its dataset.
    qualified=parent if parent.startswith('individual/') else 'individual/'+parent
    if len(Path(qualified).parts)!=2:return None
    old_locator=qualified+'/'+selector.name
    if ref.get('artifact_id') and ref['artifact_id']!=slot_id(ref['stage'],old_locator):return None
    role={'bbb-pos':'BBB+','bbb-neg':'BBB-','hia-pos':'HIA+','hia-neg':'HIA-'}.get(match[1],match[1])
    allowed={str(Path(p)) for p in files};matches=[]
    for file in files:
        path=Path(file)
        if not path.name.endswith('_admet_report.json'):continue
        if Path(locator(ref['stage'],path)).parent.as_posix()!=qualified:continue
        try:report=json.loads(path.read_text())
        except (OSError,ValueError):continue
        name=report.get('outputs',{}).get(role)
        # Schema 3 already declares its tagged stem and endpoint partitions.
        if not name and report.get('schema_version')==3 and role in ('BBB+','BBB-','HIA+','HIA-'):
            stem=report.get('output_stem','')
            if stem and Path(stem).name==stem:name=stem+'_'+role+'.csv'
        if not name or Path(name).name!=name:continue
        candidate=path.parent/name
        if str(candidate) in allowed and candidate.is_file():matches.append(str(candidate))
    return matches[0] if len(matches)==1 else None
