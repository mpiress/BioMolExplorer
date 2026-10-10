"""Keep compound preparation failures separate from usable docking inputs."""
from contextlib import contextmanager
from contextvars import ContextVar
import csv
from pathlib import Path
import subprocess

REPORT_NAME='preparation_failures.csv'
FIELDS=['compound_identifier','submitted_smiles','input_file','input_row','preparation_step','preparation_context','reason']
CURRENT_REPORT=ContextVar('preparation_failure_report',default=None)

class PreparationFailures:
    def __init__(self):self.records=[]

    @contextmanager
    def compound(self,dataset,row_number,row,context,cleanup,step='Preparação do ligante'):
        try:
            yield
        except (ValueError,RuntimeError,FileNotFoundError,subprocess.SubprocessError) as error:
            # Missing engines/configuration and filesystem failures are not bad compounds.
            if getattr(error,'error_code',None) in ('TOOL_NOT_FOUND','TOOL_START_FAILED'):raise
            if isinstance(error,FileNotFoundError) and Path(str(error.filename)).name in ('chimera','obabel'):raise
            for path in cleanup:Path(path).unlink(missing_ok=True)
            reason=f'{type(error).__name__}: {error}'
            self.records.append(dict(compound_identifier=row.get('molecule_chembl_id',''),
                submitted_smiles=row.get('canonical_smiles',''),input_file=str(dataset),input_row=row_number,
                preparation_step=step,preparation_context=str(context),reason=reason))
            from .molecule_quality import CURRENT_REPORT as QUALITY
            report=QUALITY.get()
            if report:report.exclude(dataset,row_number,[row.get('molecule_chembl_id')],'Falha em '+step+': '+reason)
            from .diagnostics import get_logger
            get_logger('scientific_tools').warning('Composto %s não utilizado (%s): %s',row.get('molecule_chembl_id',''),step,reason,exc_info=True)

    def write(self,directory):
        if not self.records:return None
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        path=directory/REPORT_NAME
        from .storage import write_text
        import io
        stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(self.records);write_text(stream.getvalue(),path)
        return path
