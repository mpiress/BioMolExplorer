"""Remove bad molecular records from execution copies, with an auditable report."""
import ast
import csv
import hashlib
import math
from collections import Counter
from contextvars import ContextVar
from pathlib import Path

from .input_validation import ALIASES, _code, columns, validate_file
from .storage import write_json

REPORT_NAME = 'molecule_exclusions.json'
CURRENT_REPORT = ContextVar('molecule_quality_report', default=None)


class QualityReport:
    def __init__(self, operation='analysis'):
        self.operation = operation
        self.records = []
        self._seen = set()
        self._written = {}

    def exclude(self, source, row, identifiers, reason):
        identifiers = [str(i) for i in identifiers if i is not None and str(i)]
        key = (str(source), row, tuple(identifiers), reason)
        if key in self._seen:return
        self._seen.add(key)
        self.records.append(dict(source=str(source), row=row, identifiers=identifiers, reason=reason))

    def write(self, output):
        if not self.records:return None
        path = Path(output) / REPORT_NAME
        if self._written.get(path)==len(self.records) and path.exists():return path
        write_json(dict(version=1, operation=self.operation, excluded_records=len(self.records),
                        records=self.records), path)
        from .diagnostics import get_logger
        get_logger('backend').warning('%s: %s registros excluídos; relatório=%s', self.operation, len(self.records), path)
        self._written[path] = len(self.records)
        return path


def table_kind(path):
    fields = columns(path)
    if {'source','target','value'} <= fields:return 'similarity'
    if {'molecule_chembl_id','fingerprint'} <= fields:return 'fingerprints'
    if 'canonical_smiles' in fields:return 'compounds'
    if fields & {'molecule_chembl_id','molecule','id','Unnamed: 0'} and fields & {'vina','dock6','score'}:return 'scores'
    return None


def _validated_row(row, kind):
    """Return normalized data and its molecular identity, or a record error."""
    from rdkit import Chem
    if None in row or any(v is None for v in row.values()):raise ValueError('quantidade de campos diferente do cabeçalho')
    row = {ALIASES.get(k,k):v for k,v in row.items()}
    code = row.get('molecule_chembl_id')
    signature = None
    if 'canonical_smiles' in row:
        smiles = row['canonical_smiles']
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles,str) and smiles.strip() else None
        if mol is None or not mol.GetNumAtoms():raise ValueError('SMILES ausente ou inválido')
        signature = Chem.MolToSmiles(mol)
        if kind=='compounds' and not code:
            code = row['molecule_chembl_id'] = 'USER_' + hashlib.sha256(signature.encode()).hexdigest()[:12]
    if kind in ('compounds','fingerprints') and not _code(code):raise ValueError('código do composto ausente ou inválido')
    if kind=='fingerprints':
        bits = ast.literal_eval(row['fingerprint']) if isinstance(row['fingerprint'],str) else row['fingerprint']
        if not isinstance(bits,(list,tuple)) or not bits or any(type(v) is not int or v not in (0,1) for v in bits):
            raise ValueError('fingerprint inválido')
        signature = (signature, tuple(bits))
    if kind=='similarity':
        if not _code(row['source']) or not _code(row['target']):raise ValueError('código na relação de similaridade inválido')
        score = float(row['value'])
        if not math.isfinite(score) or not 0<=score<=1:raise ValueError('valor de similaridade inválido')
    if kind=='scores':
        code = row.get('molecule_chembl_id') or row.get('molecule') or row.get('id') or row.get('Unnamed: 0')
        if not _code(code):raise ValueError('código do composto ausente ou inválido')
    numeric = ('TPSA','WLOGP') if kind=='compounds' else ('vina','dock6','score','z-score','min-max') if kind=='scores' else ()
    for field in numeric:
        if field in row and (kind=='scores' or row[field] not in ('',None)) and not math.isfinite(float(row[field])):
            raise ValueError(field+' inválido')
    return row, code, signature


def clean_rows(rows, kind, report, source):
    valid = []; signatures = {}; conflicts = set(); widths = Counter()
    for number, original in enumerate(rows, 2):
        ids = [original.get('source'),original.get('target')] if kind=='similarity' else [original.get('molecule_chembl_id') or original.get('name')]
        try:
            row, code, signature = _validated_row(original, kind)
        except (ValueError, TypeError, SyntaxError, KeyError, OverflowError) as exc:
            report.exclude(source, number, ids, str(exc));continue
        if code and signature is not None:
            if code in signatures and signatures[code] != signature:conflicts.add(code)
            signatures[code] = signature
        width = len(signature[1]) if kind=='fingerprints' else None
        if width is not None:widths[width] += 1
        valid.append((number,row,code,width))
    width = widths.most_common(1)[0][0] if widths else None
    for number,row,code,current_width in valid:
        reason = ('código associado a dados moleculares conflitantes' if code in conflicts else
                  'tamanho de fingerprint incompatível' if current_width is not None and current_width!=width else None)
        if reason:report.exclude(source,number,[code],reason)
        else:yield row


def sanitize_csv(source, destination, kind=None, report=None):
    source, destination = Path(source), Path(destination)
    kind = kind or table_kind(source)
    if kind is None:raise ValueError('Tabela molecular sem as colunas esperadas.')
    # Broken schemas, paths and permissions remain errors; only records are excluded.
    validate_file(source,kind,validate_rows=False)
    report = report or CURRENT_REPORT.get() or QualityReport()
    with source.open(encoding='utf-8-sig',newline='') as stream:
        reader = csv.DictReader(stream)
        fields = [ALIASES.get(k,k) for k in reader.fieldnames]
        if kind=='compounds' and 'molecule_chembl_id' not in fields:fields.insert(0,'molecule_chembl_id')
        rows = list(clean_rows(reader,kind,report,source))
    from tempfile import NamedTemporaryFile
    destination.parent.mkdir(parents=True,exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode='w',encoding='utf-8',newline='',dir=destination.parent,delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        temporary.replace(destination)
    finally:
        if temporary:temporary.unlink(missing_ok=True)
    return destination


def clean_dataframe(frame, kind='compounds', source='dataset', report=None):
    import pandas as pd
    report = report or CURRENT_REPORT.get() or QualityReport()
    frame = frame.rename(columns=ALIASES)
    fields = list(frame.columns)
    if kind=='compounds' and 'molecule_chembl_id' not in fields:fields.insert(0,'molecule_chembl_id')
    # CSV blanks are represented as NaN by pandas; let the same record validator
    # handle blanks in memory and on disk, including optional missing metadata.
    rows = frame.astype(object).where(pd.notna(frame),'').to_dict('records')
    return pd.DataFrame(clean_rows(rows,kind,report,source),columns=fields)


def filter_edges(edges, identifiers, source, report=None):
    report = report or CURRENT_REPORT.get() or QualityReport()
    known = set(map(str,identifiers))
    keep = edges['source'].isin(known) & edges['target'].isin(known)
    for number, row in enumerate(edges.itertuples(index=False), 2):
        if str(row.source) not in known or str(row.target) not in known:
            report.exclude(source,number,[row.source,row.target],'relação sem molécula correspondente na entrada')
    return edges.loc[keep].reset_index(drop=True)


def merge_clean_csv(paths, destination, kind, report):
    import pandas as pd
    frames=[]
    for index,path in enumerate(paths):
        cleaned=Path(destination).parent/'.quality-inputs'/str(index)/Path(path).name
        sanitize_csv(path,cleaned,kind,report)
        frames.append(pd.read_csv(cleaned,dtype={'molecule_chembl_id':str,'source':str,'target':str}))
    combined=clean_dataframe(pd.concat(frames,ignore_index=True),kind,'entradas combinadas',report)
    keys=['source','target'] if kind=='similarity' else ['molecule_chembl_id']
    from .storage import write_dataframe
    write_dataframe(combined.drop_duplicates(keys),destination)
    return destination


def prepare_inputs(operation, parameters, output, report):
    """Sanitize immutable copies of molecular CSVs for every scientific adapter."""
    import copy
    import shutil
    from .paths import resolve_path
    parameters = copy.deepcopy(parameters)
    cache = {}; root = Path(output)/'.quality-inputs'
    def clean(path):
        path = Path(path)
        if path not in cache:
            kind = table_kind(path)
            destination = root/str(len(cache))/path.name
            cache[path] = sanitize_csv(path,destination,kind,report) if kind else path
        return str(cache[path])
    for entry in parameters.get('graph_inputs',[]):
        entry['file'] = clean(resolve_path(entry['file']))
        entry['compound_files'] = [clean(resolve_path(p)) for p in entry['compound_files']]
    # Structural folders and ZINC address lists use their existing validators.
    # Their compound tables still pass through this same boundary at docking.
    if operation=='graphs':return parameters
    for field in ('base_input_path','base_selected_mols'):
        if not parameters.get(field):continue
        folder = resolve_path(parameters[field])
        if not folder.is_dir():continue
        files = [p for p in sorted(folder.glob('*.csv')) if table_kind(p)]
        filename = (parameters.get('mol_filename','')+'.csv' if field=='base_selected_mols' else
                    parameters.get('input_file') if operation in ('admet','expand_similar_compounds') else
                    parameters.get('filename') if operation=='similarity' else None)
        selected = parameters.get('files') if operation=='fingerprints' else [filename] if filename else None
        if selected:files = [p for p in files if p.name in selected]
        if not files:continue
        destination = root/field;destination.mkdir(parents=True,exist_ok=True)
        for p in folder.iterdir():
            if p.is_file() and not p.is_symlink():shutil.copy2(p,destination/p.name)
        for path in files:sanitize_csv(path,destination/path.name,report=report)
        parameters[field] = str(destination)
    return parameters
