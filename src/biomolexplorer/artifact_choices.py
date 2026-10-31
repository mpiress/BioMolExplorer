"""Qualified output selectors prevent ambiguous filenames across data sources."""
from pathlib import Path
from .input_validation import columns
import re


def logical_path(path):
    """Remove the worker job identity while retaining the batch and output tree."""
    return re.sub(r'\.biomolexplorer/jobs/[^/]+/artifacts/', '', Path(path).as_posix())


def matches_selector(path, selector):
    path=logical_path(path); selector=logical_path(selector)
    return path==selector or path.endswith('/'+selector)


def selector_for(path, files):
    """Find the shortest unique suffix after removing transient job directories."""
    parts=Path(logical_path(path)).parts
    for length in range(1,len(parts)+1):
        suffix='/'.join(parts[-length:])
        if sum(matches_selector(p,suffix) for p in files)==1:
            return suffix
    raise ValueError('Os resultados contêm caminhos de arquivo duplicados.')


def produced_files(item):
    """Public outputs; ADMET inputs are lineage, not downstream datasets."""
    files=list(dict.fromkeys(Path(p) for p in item.get('artifacts',[])))
    config=item.get('configuration',item)
    if item.get('operation',config.get('operation'))!='admet' or config.get('provided_results'):return files
    inputs={str(Path(p)) for owner in [item,*item.get('batches',item.get('_execution_batches',[]))]
            for paths in owner.get('input_files',{}).values() for p in paths}
    result=[]
    for path in files:
        if str(path) in inputs or '.quality-inputs' in path.parts or '.inputs' in path.parts:continue
        if path.suffix.lower()=='.csv':
            canonical=re.fullmatch(r'B\d+__admet__(?:compounds|evaluated|selected|rejected|bbb-pos|bbb-neg|hia-pos|hia-neg|alerts|rules-lipinski-veber|rules-lipinski|rules-veber|rules-fragment)(?:-[a-f0-9]{8})?',path.stem)
            if not canonical:
                fields=columns(path)
                if not ({'BBB','HIA'}<=fields or 'ADMET_selected' in fields):continue
        result.append(path)
    return result


def choices(item):
    files=produced_files(item)
    if item['operation'] in ('retrieve_compounds','retrieve_pubchem','expand_similar_compounds') and not item.get('configuration',{}).get('provided_results'):
        files=[p for p in files if p.suffix=='.csv' and (any('__'+role in p.stem for role in ('compounds','mol','sim','merged')) or p.name=='compounds.csv' and p.parent.parent.name=='compounds'
                or p.stem.endswith(('_FULL','_MOLS','_SIMS')))]
    if item['operation']=='retrieve_structures':
        files=[p for p in files if p.suffix.lower()=='.pdb']
    if item['operation'] in ('docking_vina','docking_dock6','consensus') and not item.get('configuration',{}).get('provided_results'):
        files=consolidated_docking_tables(files,item['operation'])
    if item['operation']=='graphs' and not item.get('configuration',{}).get('provided_results'):
        from .graph_exports import molecular_files
        files=molecular_files(files)
    if item['operation']=='redocking':files=redocking_receptors(files)
    result=set()
    for p in files:
        parts=p.parts
        if item.get('batches') and item['id'] in parts:
            result.add(logical_path('/'.join(parts[parts.index(item['id'])+1:])))
        elif 'artifacts' in parts:
            index=len(parts)-1-list(reversed(parts)).index('artifacts')
            result.add('/'.join(parts[index+1:]))
        elif p.suffix.lower() in ('.csv','.pdb','.pdbqt','.mol2','.txt'):
            result.add(p.name)
    return result


def consolidated_docking_tables(paths,operation=None):
    """Final full compound score tables, excluding contexts and pose intermediates."""
    tables=[]
    for path in dict.fromkeys(Path(p) for p in paths):
        if path.suffix.lower()!='.csv' or path.parent.name in ('context','poses','Prepared','inputs','.inputs'):continue
        fields=columns(path)
        docking={'molecule_chembl_id','engine','score','conformer_file'}<=fields
        consensus={'molecule_chembl_id','vina','dock6'}<=fields
        if (operation!='consensus' and docking) or (operation not in ('docking_vina','docking_dock6') and consensus):tables.append(path)
    return tables


def redocking_receptors(paths):
    """Public redocking handoffs expose receptors, never crystal ligands."""
    from .docking_inputs import preferred_receptors,prepared_receptor
    from .md_inputs import receptor_input_files
    paths=[Path(p) for p in paths]
    return list(dict.fromkeys([*preferred_receptors(p for p in paths if prepared_receptor(p)),*receptor_input_files(paths)]))
