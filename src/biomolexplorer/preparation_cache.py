"""Project-scoped, content-addressed native preparation shared by all blocks.

Immutable entries are published by atomic rename only after every required file
has been validated. Readers verify every digest; failures are never cache hits.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from uuid import uuid4

from .stage_cache import file_digest
from .redocking_config import DEFAULTS


def receptor_key(source,record,options,ph):
    import ast
    import inspect
    import textwrap
    from .paths import resolve_path
    from . import chimera_preparation
    from .redocking_config import configure_template
    from caad.docking import Docking
    protocol=hashlib.sha256()
    for function in (configure_template,Docking.generate_docking_script,Docking.prepare_on_chimera,Docking.validate_prepared_file):
        protocol.update(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(function))),include_attributes=False).encode())
    tools=tool_manifest()
    templates={name:file_digest(resolve_path('src/scripts/chimera/'+name+'.template'))
               for name in ('prepare_complex','prepare_receptor')}
    payload=dict(schema=1,source=file_digest(source),chain=str(record[3]),
        receptor=dict(DEFAULTS,**options.get('receptor',{})),cofactors=sorted(set(options.get('cofactors',[]))),
        ph=ph,templates=templates,tools=tools,implementation=file_digest(chimera_preparation.__file__),protocol=protocol.hexdigest())
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()


def tool_manifest(search_path=None):
    tools={}
    for tool in ('chimera','obabel'):
        path=shutil.which(tool,path=search_path)
        if path:
            executable=Path(path).resolve();tools[tool]=file_digest(executable)
            if tool=='chimera':
                native=executable.parent.parent
                for relative in ('share/AddCharge/__init__.py','share/AddH/__init__.py','share/DockPrep/__init__.py',
                                 'share/MMMD/MMTKinter.py','bin/amber18/bin/sqm','bin/amber18/bin/antechamber'):
                    dependency=native/relative
                    if dependency.is_file():tools[relative]=file_digest(dependency)
    return tools


def root():
    if os.environ.get('BIOMOL_REUSE_RESULTS')=='0':return None
    value=os.environ.get('BIOMOL_CACHE_DIR')
    return Path(value)/'preparation'/'receptors' if value else None


def restore(key,prepared,prefix):
    base=root()
    if base is None:return False
    entry=base/key;prepared=Path(prepared)
    try:
        manifest=json.loads((entry/'manifest.json').read_text())
        if manifest.get('key')!=key:return False
        required={'.dockprep.mol2','.noH.pdb'}
        if not required<=set(manifest['files']):return False
        for suffix,digest in manifest['files'].items():
            if not suffix.startswith('.') or '/' in suffix or '\\' in suffix:return False
            if file_digest(entry/('receptor'+suffix))!=digest:return False
        for suffix in manifest['files']:shutil.copy2(entry/('receptor'+suffix),prepared/(prefix+suffix))
        return True
    except (OSError,ValueError,KeyError):return False


def publish(key,prepared,prefix):
    base=root()
    if base is None:return
    prepared=Path(prepared);required=('.dockprep.mol2','.noH.pdb')
    if not all((prepared/(prefix+s)).is_file() and (prepared/(prefix+s)).stat().st_size for s in required):return
    base.mkdir(parents=True,exist_ok=True)
    if restore_valid(base/key):return
    # Scientific caller validates native exports before publication.
    with tempfile.TemporaryDirectory(dir=base,prefix='.building-') as temp:
        work=Path(temp);files={}
        for suffix in (*required,'.dockprep.pdbqt','.dockprep.preparation.json'):
            source=prepared/(prefix+suffix)
            if source.is_file():
                shutil.copy2(source,work/('receptor'+suffix));files[suffix]=file_digest(source)
        (work/'manifest.json').write_text(json.dumps(dict(key=key,files=files),sort_keys=True)+'\n')
        target=base/key
        # An invalid entry is quarantined atomically, never overwritten in place.
        if target.exists() and not restore_valid(target):
            quarantine=base/('.invalid-'+uuid4().hex)
            try:target.rename(quarantine)
            except FileNotFoundError:pass
            else:shutil.rmtree(quarantine,ignore_errors=True)
        try:work.rename(target)
        except OSError:
            if not restore_valid(target):raise


def restore_valid(entry):
    try:
        meta=json.loads((entry/'manifest.json').read_text())
        return meta.get('key')==entry.name and {'.dockprep.mol2','.noH.pdb'}<=set(meta['files']) and all(
            s.startswith('.') and '/' not in s and '\\' not in s and file_digest(entry/('receptor'+s))==digest for s,digest in meta['files'].items())
    except (OSError,ValueError,KeyError):return False
