"""Chemical component definitions available to the Chimera preparation.

Only the cache identity of the local definition lives here.  The IDATM
command repair that used to accompany it was removed on 2026-11-28, by
decision C8 of the triage of 2026-11-27: the two-phase preparation reads
the formal charge from the hydrogenated graph and reports an open-shell
one, instead of rewriting atom types to get past the charge assignment.
"""


def local_definition(folder,code,config):
    """Fingerprint available chemistry without fetching it during cache lookup."""
    import os
    from pathlib import Path
    from .stage_cache import file_digest
    folder=Path(folder)
    explicit=config.get('chemistry_templates',{}).get(code)
    if explicit:
        path=Path(explicit)
        if not path.is_absolute():path=folder/path
        return {str(path):file_digest(path) if path.is_file() else None}
    roots=[folder/'chemistry',folder]
    if os.environ.get('BIOMOL_CCD_PATH'):roots.insert(0,Path(os.environ['BIOMOL_CCD_PATH']))
    preparation=os.environ.get('BIOMOL_PREPARATION_CACHE')
    cache=Path(preparation).parent/'calculations' if preparation else folder/'.calculation-cache'
    roots.append(cache/'ccd')
    for root in roots:
        for suffix in ('.cif','.sdf','.mol2'):
            path=root/(code+suffix)
            if path.is_file():return {code:file_digest(path)}
    return {}
