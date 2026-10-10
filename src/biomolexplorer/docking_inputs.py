"""Pair docking poses with their reference receptor and selected compounds."""
from .artifact_choices import matches_selector
import csv
from pathlib import Path


def target_input(stage, field):
    return field == 'base_input_path' and (stage['operation'] in ('docking_vina', 'docking_dock6') and stage['parameters'].get('receptor_prepared',True) or
        stage['operation']=='prepare_structures' and stage['parameters'].get('receptor_prepared',False))


def prepared_receptor(path):
    """The preparation stage gives ligand-free receptors this dedicated suffix."""
    return Path(path).name.endswith('.dockprep.pdbqt')


def ligand_input_file(path):
    path=Path(path)
    return path.suffix.lower() in ('.csv','.sdf','.mol2','.pdbqt') and '.dockprep.' not in path.name and '.noH.' not in path.name


def compound_codes(path):
    from .molecule_quality import _validated_row
    codes=set()
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        for row in csv.DictReader(stream):
            try:_,code,_=_validated_row(row,'compounds')
            except (ValueError,TypeError,KeyError):continue
            if code:codes.add(code)
    return codes


def matching_poses(poses, records, codes):
    identities={identifier for record in records for code in codes for identifier in (
        f'{record[0]}_{record[1]}_{record[2]}{record[3]}_{code}',
        f'{record[0]}_{record[3]}_{code}', f'{record[0]}_{code}')}
    # Standalone externally prepared poses can use only the compound code.
    identities.update(codes)
    return [Path(p) for p in poses if Path(p).name.removesuffix('.pdbqt').removesuffix('.lig') in identities]


def dock6_variant_matches(stage, refs, results, store, project_id):
    def files(reference):
        if 'asset' in reference:
            with store.connect() as db:
                row=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(project_id,reference['asset'])).fetchone()
            return [store.scoped_path(project_id,row[0])] if row else []
        paths=[store.scoped_path(project_id,p) for p in results.get(reference['stage'],[])]
        selector=reference.get('selector','auto')
        return paths if selector=='auto' else [p for p in paths if matches_selector(p,selector)]
    receptor_ref=refs['base_input_path']
    receptor_files=files(receptor_ref)
    receptors={p.name.split('.',1)[0] for p in receptor_files if p.name.endswith('.dockprep.pdbqt')}
    all_paths=[store.scoped_path(project_id,p) for p in results.get(receptor_ref.get('stage'),[])]
    records=[]
    for path in all_paths:
        if path.name!='pdb_codes.csv':continue
        with path.open(encoding='utf-8-sig',newline='') as stream:
            for record in csv.DictReader(stream):
                if f"{record.get('PDB_CODE')}_{record.get('CHAIN')}" in receptors:
                    records.append([record[k] for k in ('PDB_CODE','LIGAND','RESNUM','CHAIN')])
    if not records and stage['parameters'].get('pdb_code'):
        record=stage['parameters']['pdb_code']
        if record and isinstance(record[0],(list,tuple)):record=record[0]
        records=[record]
    codes=set()
    for path in files(refs['base_selected_mols']):
        if path.suffix.lower()=='.csv':codes.update(compound_codes(path))
    return bool(matching_poses(files(refs['base_vina_path']),records,codes))
