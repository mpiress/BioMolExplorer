"""Qualified output selectors prevent ambiguous filenames across data sources."""
from pathlib import Path
from .input_validation import columns


def choices(item):
    files=[Path(p) for p in item.get('artifacts',[])]
    if item['operation'] in ('retrieve_compounds','expand_similar_compounds') and not item.get('configuration',{}).get('provided_results'):
        files=[p for p in files if p.suffix=='.csv' and (p.name=='compounds.csv' and p.parent.parent.name=='compounds'
                or p.stem.endswith(('_FULL','_MOLS','_SIMS')))]
    if item['operation']=='retrieve_structures':
        files=[p for p in files if p.suffix.lower()=='.pdb']
    result=set()
    for p in files:
        parts=p.parts
        if item.get('batches') and item['id'] in parts:
            result.add('/'.join(parts[parts.index(item['id'])+1:]))
        elif 'artifacts' in parts:
            index=len(parts)-1-list(reversed(parts)).index('artifacts')
            result.add('/'.join(parts[index+1:]))
        elif p.suffix.lower() in ('.csv','.pdb','.pdbqt','.mol2','.txt'):
            result.add(p.name)
    return result
