"""Entity-level selections materialized only at external-tool boundaries."""
import csv
import hashlib
import json
import math
from pathlib import Path

from .storage import write_dataframe


IDENTIFIER_COLUMNS=('molecule_chembl_id','PubChem_CID','pubchem_cid','zinc_id','ZINC_ID','compound_id','id')


def _ids(selection):
    values=list(selection.get('entity_ids',[]))
    for key in ('pdb_id','ligand_id','receptor_id','pose_id','community_id'):
        if selection.get(key) is not None:values.append(str(selection[key]))
    return {str(value).strip().casefold() for value in values}


def path_matches(path,selection):
    """Cheap structure/pose filtering; metadata tables are filtered separately."""
    if not selection:return True
    path=Path(path);stem=path.name.split('.',1)[0].casefold();parts={p.casefold() for p in path.parts}
    pdb=str(selection.get('pdb_id','')).casefold()
    chain=str(selection.get('chain_id','')).casefold()
    receptor=str(selection.get('receptor_id','')).casefold()
    ligand=str(selection.get('ligand_id','')).casefold()
    if pdb and not (stem==pdb or stem.startswith(pdb+'_') or pdb in parts):return False
    if chain and not (stem.endswith('_'+chain) or '_'+chain+'_' in stem):
        # A retrieved PDB commonly contains several chains in one immutable
        # artifact.  Validate membership from the structure without creating
        # a derivative file; the preparation adapter applies the configured
        # chain/residue choice later.
        available=set()
        if path.suffix.lower() in ('.pdb','.pdbqt') and path.is_file():
            try:
                with path.open(encoding='utf-8',errors='replace') as stream:
                    available={line[21].strip().casefold() for line in stream
                               if line.startswith(('ATOM  ','HETATM')) and len(line)>21 and line[21].strip()}
            except OSError:return False
        if chain not in available:return False
    assembly=str(selection.get('assembly_id','')).casefold()
    if assembly and not (f'assembly_{assembly}' in stem or f'assembly-{assembly}' in stem):
        # Assembly selection is supported only when the producer publishes a
        # distinct assembly artifact; silently using the asymmetric unit would
        # be scientifically incorrect.
        return False
    if receptor and receptor not in stem and receptor not in parts:return False
    if ligand and ligand not in stem and ligand not in parts:return False
    return True


def _property_match(value,rule):
    if isinstance(rule,dict):
        try:number=float(value)
        except (TypeError,ValueError):return False
        return (rule.get('min') is None or number>=float(rule['min'])) and (rule.get('max') is None or number<=float(rule['max']))
    if isinstance(rule,(list,tuple,set)):return str(value) in {str(v) for v in rule}
    return str(value).casefold()==str(rule).casefold()


def row_matches(row,selection):
    wanted=_ids(selection)
    if wanted:
        values={str(row.get(column,'')).strip().casefold() for column in IDENTIFIER_COLUMNS}
        # Docking and structure tables use these identifiers as well.
        values|={str(row.get(column,'')).strip().casefold() for column in ('PDB_CODE','pdb_id','receptor_id','pose_id','LIGAND','CHAIN')}
        if not wanted & values:return False
    mappings={'pdb_id':('PDB_CODE','pdb_id'),'chain_id':('CHAIN','chain_id'),
              'ligand_id':('LIGAND','ligand_id'),'receptor_id':('receptor_id',),
              'pose_id':('pose_id',),'community_id':('community','community_id','graph_component')}
    for key,columns in mappings.items():
        if selection.get(key) is not None and str(selection[key]).casefold() not in {
                str(row.get(column,'')).casefold() for column in columns}:return False
    return all(_property_match(row.get(name),rule) for name,rule in selection.get('properties',{}).items())


def csv_view(source,selection,root):
    """Create a deterministic cached adapter view without altering the artifact."""
    source=Path(source)
    signature=hashlib.sha256(json.dumps([str(source.resolve()),source.stat().st_size,source.stat().st_mtime_ns,selection],
        sort_keys=True,separators=(',',':')).encode()).hexdigest()
    destination=Path(root)/signature/source.name
    if destination.is_file():return destination
    import pandas as pd
    frame=pd.read_csv(source)
    mask=[row_matches({key:value for key,value in row.items() if not (isinstance(value,float) and math.isnan(value))},selection)
          for row in frame.to_dict('records')]
    selected=frame.loc[mask]
    if selected.empty:raise ValueError('A seleção não corresponde a nenhuma entidade disponível.')
    write_dataframe(selected,destination)
    return destination
