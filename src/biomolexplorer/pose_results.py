"""Authorized reading of a pose validation report for the results view.

The block already publishes everything the view needs: ``pose_validation.csv``
with one row per evaluated pose and ``pose_validation.json`` with the counts,
the tolerances and the criteria that produced them.  This reader turns both
into the table the interface renders, and resolves each pose to the copies the
block published, so the complex can be opened without reaching into the run
directory of the docking block that produced it.
"""
import csv
import json
from pathlib import Path

from .result_files import ResultFiles
from .workspace import AccessDenied

REPORT = 'pose_validation.csv'
SUMMARY = 'pose_validation.json'
MAX_REPORT_BYTES = 64 * 1024 * 1024
MAX_ROWS = 20000
# Criteria in the order the wrapper evaluates them.
CRITERIA = (('valid_chemistry','Química'),('valid_connectivity','Conectividade'),
            ('valid_geometry','Geometria'),('valid_planarity','Planaridade'),
            ('valid_energy','Energia interna'),('valid_clash','Choque estérico'))
NUMERIC = ('score','min_vdw_ratio','min_vdw_overlap','min_heavy_distance','delta_energy',
           'max_bond_deviation','max_angle_deviation','max_planarity_deviation')
INTEGER = ('clashing_contacts','planar_groups','atoms')
TEXT = ('pose_id','receptor_id','ligand_id','compound_id','engine','energy_method','bond_orders',
        'resolved_by','rejection_details','error','source_table')


def number(value):
    try:
        result = float(value)
    except (TypeError,ValueError):
        return None
    return result if result == result and abs(result) != float('inf') else None


def integer(value):
    value = number(value)
    return int(value) if value is not None else None


def flag(value):
    """The report writes the verdicts as the booleans the evaluation produced."""
    text = str(value).strip().casefold()
    return True if text in ('true','1','yes') else False if text in ('false','0','no') else None


class PoseValidationResults:
    """Report of one pose validation block, read with the project's authorization."""

    def __init__(self, store):
        self.store = store
        self.results = ResultFiles(store)

    def stage(self, token, pid, rid, sid):
        stage = self.results.stage(token,pid,rid,sid)
        if stage['operation'] != 'pose_validation':
            raise AccessDenied('Resultados de validação de poses indisponíveis.')
        return stage

    def _published(self, stage):
        """Poses and receptors this block copied next to its report."""
        published = {}
        for name in stage.get('artifacts',[]):
            path = Path(name)
            if path.parent.name in ('poses','receptors'):
                published.setdefault((path.parent.name,path.name),str(path))
        return published

    def _locate(self, pid, published, kind, value):
        """Prefer the copy this block published; keep the original otherwise."""
        if not value:
            return None
        for candidate in (published.get((kind,Path(value).name)),value):
            if not candidate:
                continue
            try:
                path = self.store.scoped_path(pid,candidate)
            except Exception:                           # noqa: BLE001 - outside the project
                continue
            if path.is_file():
                return str(path)
        return None

    def report(self, token, pid, rid, sid):
        """Summary, evaluated poses and the tolerances that judged them."""
        stage = self.stage(token,pid,rid,sid)
        published = self._published(stage)
        summary, rows, truncated = {}, [], False
        for name in stage.get('artifacts',[]):
            path = self.store.scoped_path(pid,name)
            if path.name == SUMMARY:
                try:
                    summary = json.loads(self.store.read_file(token,pid,str(path),MAX_REPORT_BYTES).decode('utf-8'))
                except (ValueError,UnicodeDecodeError):
                    summary = {}
            elif path.name == REPORT:
                data = self.store.read_file(token,pid,str(path),MAX_REPORT_BYTES+1)
                truncated = len(data) > MAX_REPORT_BYTES
                for record in csv.DictReader(data[:MAX_REPORT_BYTES].decode('utf-8-sig',errors='replace').splitlines()):
                    if len(rows) >= MAX_ROWS:
                        truncated = True
                        break
                    rows.append(self._row(pid,published,record))
        if not isinstance(summary,dict):
            summary = {}
        for row in rows:
            row['failed'] = [label for key,label in CRITERIA if row.get(key) is False]
        # Rejected poses first: they are the ones that need a decision.
        rows.sort(key=lambda row:(row['pose_valid'] is not False,str(row['receptor_id']),
                                  str(row['compound_id']),str(row['pose_id'])))
        return dict(summary=summary,rows=rows,truncated=truncated,criteria=list(CRITERIA),
                    counts=dict(total=len(rows),valid=sum(1 for row in rows if row['pose_valid']),
                                invalid=sum(1 for row in rows if row['pose_valid'] is False)))

    def _row(self, pid, published, record):
        row = {key:(record.get(key) or '').strip() for key in TEXT}
        row.update({key:number(record.get(key)) for key in NUMERIC})
        row.update({key:integer(record.get(key)) for key in INTEGER})
        row['pose_valid'] = flag(record.get('pose_valid'))
        for key,_ in CRITERIA:
            row[key] = flag(record.get(key))
        row['rejection_reasons'] = [item for item in (record.get('rejection_reasons') or '').split(';') if item]
        row['compound_id'] = row['compound_id'] or row['ligand_id']
        row['pose_path'] = self._locate(pid,published,'poses',record.get('pose_file'))
        row['receptor_path'] = self._locate(pid,published,'receptors',record.get('receptor_file'))
        return row
