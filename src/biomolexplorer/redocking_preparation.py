"""Isolate structure preparation failures without hiding installation failures."""
import csv
import io
from pathlib import Path
import subprocess

from .redocking_config import pair_key, validate_structure_pairs
from .storage import write_text

REPORT_NAME='redocking_failures.csv'
FIELDS=['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION','STATUS','FAILED_STEP','ERROR']


def recoverable(error):
    if getattr(error,'error_code',None) in ('TOOL_NOT_FOUND','TOOL_START_FAILED','TOOL_VERSION_UNSUPPORTED'):return False
    if isinstance(error,FileNotFoundError) and Path(str(error.filename)).name in ('chimera','obabel','vina'):return False
    return isinstance(error,(ValueError,RuntimeError,FileNotFoundError,subprocess.SubprocessError))


def prepare_pairs(dock,records,settings,charge_type,pH,ready=None,engine='vina'):
    """Reuse receptors, but isolate ligand failures and clean partial outputs."""
    ready=set(ready or ());failed_receptors={};prepared=[];failures=[]
    for record in records:
        receptor=(record[0],record[3]);key=pair_key(record)
        try:
            if receptor in failed_receptors:raise ValueError(failed_receptors[receptor])
            validate_structure_pairs(dock.complexpath,[record],{key:settings[key]} if key in settings else {},engine=engine)
            completed=dock.prepare_for_docking([record],charge_type,pH,True,
                preparation_pairs=settings,_prepared_receptors=ready,**({'docking_engine':engine} if engine!='vina' else {}))
            validate_structure_pairs(dock.complexpath,[record],{key:settings[key]} if key in settings else {},prepared=True,engine=engine)
            prepared.extend(completed)
        except Exception as error:
            if not recoverable(error):raise
            reason=f'{type(error).__name__}: {error}'
            if receptor not in ready:failed_receptors[receptor]=reason
            prefix=f'{record[0]}_{record[1]}_{record[2]}{record[3]}.'
            receptor_prefix=f'{record[0]}_{record[3]}.'
            for path in Path(dock.outputpath).iterdir():
                if path.is_file() and (path.name.startswith(prefix) or
                        receptor not in ready and path.name.startswith(receptor_prefix)):
                    path.unlink()
            failures.append(dict(zip(FIELDS[:5],list(record[:4])+[record[4] if len(record)>4 else None]),
                STATUS='preparation_failed',FAILED_STEP='Preparação do complexo',ERROR=reason))
            dock.logger.warning('Redocking: complexo ignorado %s: %s',key,reason,exc_info=True)
    return prepared,failures


def write_failures(records,directory):
    if not records:return None
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=FIELDS)
    writer.writeheader();writer.writerows(records)
    path=Path(directory)/REPORT_NAME;write_text(stream.getvalue(),path)
    return path
