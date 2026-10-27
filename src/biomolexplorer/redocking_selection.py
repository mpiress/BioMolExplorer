"""Keep a curated structure's ligands in connected redocking execution plans."""
from copy import deepcopy
from pathlib import Path

from .artifact_choices import matches_selector, selector_for
from .bindings import pack, sources
from .redocking_config import pair_key, validate_pairs


def sync_curated_pairs(pipeline, producer_id, path, artifacts, previous, curated):
    """Apply one explicit ligand edit; retain selections from other structures."""
    code = path.stem.upper()
    def record(row):
        result = [code, row['LIGAND'], int(row['RESNUM']), row['CHAIN']]
        if row.get('RESOLUTION'):
            result.append(float(row['RESOLUTION']))
        return result
    selected = [record(row) for row in curated]
    previous_keys = set()
    for row in previous:
        # Retrievals can retain placeholder rows for a PDB with no ligand.
        try:
            if row.get('LIGAND') and row.get('CHAIN'):
                previous_keys.add(pair_key([code, row['LIGAND'], int(row['RESNUM']), row['CHAIN']]))
        except (ValueError, TypeError, KeyError):
            continue
    removed = previous_keys - {pair_key(r) for r in selected}
    updated = []; added = 0
    for stage in pipeline:
        if stage['operation'] != 'redocking' or stage.get('provided_results'):
            continue
        binding = stage.get('bindings', {}).get('base_input_path')
        refs = sources(binding) if binding else []
        connected = [ref for ref in refs if ref.get('stage') == producer_id]
        if not connected:
            continue
        before = deepcopy(stage)
        params = stage['parameters']
        original = params.get('pdb_codes') or []
        settings = deepcopy(params.get('preparation_pairs') or {})
        # Newly curated receptors join the same general preparation policy.
        first=next((settings[pair_key(row)] for row in original if pair_key(row) in settings),{})
        general={key:deepcopy(first[key]) for key in
            ('receptor','ligand','preparation_backend','ambertools_path','ambertools_env','vina_charge_policy') if key in first}
        receptors = {}
        for row in original:
            config = settings.get(pair_key(row), {})
            receptors.setdefault((str(row[0]).upper(), row[3]),
                {key: deepcopy(config[key]) for key in ('cofactors','receptor','ligand','preparation_backend','ambertools_path','ambertools_env','vina_charge_policy') if key in config})
        records = [row for row in original if pair_key(row) not in removed]
        keys = {pair_key(row) for row in records}
        for row in selected:
            key = pair_key(row)
            if key in keys:
                continue
            records.append(row); keys.add(key); added += 1
            config = dict(deepcopy(general),**deepcopy(receptors.get((code,row[3]),{})))
            if config:
                settings[key] = deepcopy(config)
        settings = {key: value for key, value in settings.items() if key in keys}
        if records:
            validate_pairs(records, settings)
        params['pdb_codes'] = records
        if settings or 'preparation_pairs' in params:
            params['preparation_pairs'] = settings
        # Saved explicit file choices can otherwise exclude a newly curated PDB
        # even after adding its ligand tuple to pdb_codes.
        if selected and not any(ref.get('selector', 'auto') == 'auto' or
                Path(ref['selector']).name == 'pdb_codes.csv' or matches_selector(path, ref['selector'])
                for ref in connected):
            stage['bindings']['base_input_path'] = pack(refs + [
                dict(stage=producer_id, selector=selector_for(path, artifacts))])
        if stage != before:
            updated.append(stage['id'])
    return dict(updated_redocking_stages=updated, added_redocking_pairs=added)
