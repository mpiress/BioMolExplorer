"""Names and labels for molecular subsets exported by graph analyses."""
from pathlib import Path
import re


def subset_label(path, translate=None):
    """Describe a molecular CSV without reading it or changing its selector."""
    translate=translate or (lambda value:value)
    path=Path(path)
    match=re.fullmatch(r'(.+)__full_graph',path.stem)
    if match:return translate('Grafo completo')+' · '+match[1]
    match=re.fullmatch(r'(.+)__component_(\d+)(_mcc)?',path.stem)
    if match:
        scope=('MCC / LCC · ' if match[3] else '')+translate('Componente')+' '+str(int(match[2]))
        return scope+' · '+match[1]
    if path.parent.name=='Molecules':
        if path.name=='molecules.csv':return translate('União dos MCCs')+' · '+translate('todas as análises')
        return 'MCC / LCC · '+translate('compatibilidade')+' · '+path.stem
    return ''


def molecular_files(files):
    """Graph ports supply compounds, while edge CSVs remain downloadable results.

    Keep legacy/imported layouts that do not contain the known Molecules folder.
    """
    files=[Path(p) for p in files]
    compounds=[p for p in files if p.parent.name=='Molecules' and p.suffix.lower()=='.csv']
    return compounds if compounds else files
