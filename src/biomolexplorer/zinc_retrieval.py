"""Read ZINC tranche download lists as data and normalize 2D/3D candidates."""
import bz2
import csv
import gzip
import hashlib
import re
import tempfile
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from collections import deque
from pathlib import Path
from urllib.parse import urljoin, urlsplit, urlunsplit

HOSTS={'zinc.docking.org','zinc15.docking.org','zinc20.docking.org',
       'files.docking.org','files2.docking.org','files3.docking.org','cache.docking.org'}
LIST_SUFFIXES={'.txt','.uri','.urls','.sh','.wget','.curl','.ps1'}


def tranche_url(value):
    url=urlsplit(value)
    if (url.scheme not in ('http','https') or url.hostname not in HOSTS or url.username or url.password
            or url.port not in (None,80 if url.scheme=='http' else 443) or url.fragment):
        raise ValueError('Endereço ZINC não autorizado: '+value)
    name=url.path.lower()
    if not name.endswith(('.smi','.smi.gz','.smi.bz2','.mol2','.mol2.gz','.mol2.bz2')):
        raise ValueError('Formato de tranche ZINC não suportado: '+value)
    return urlunsplit(('https',url.hostname,url.path,url.query,''))


def read_download_list(path):
    """Accept URI lists and URL-bearing cURL/wget/PowerShell exports; never execute."""
    urls=[]
    text=re.sub(r'\\\r?\n',' ',Path(path).read_text(encoding='utf-8-sig'))
    for number,line in enumerate(text.splitlines(),1):
        line=line.strip()
        if not line or line.startswith('#'):continue
        matches=re.findall(r'https?://[^\s\"\'<>`]+',line)
        if not matches:
            if line.startswith(('mkdir ','set ','cd ','echo ','$','REM ','rem ')):continue
            raise ValueError(f'Linha {number}: informe um link de tranche ZINC SMI ou MOL2.')
        for match in matches:
            url=tranche_url(match.rstrip(';'))
            if url not in urls:urls.append(url)
    if not urls:raise ValueError('O arquivo não contém links de tranches ZINC SMI ou MOL2.')
    return urls


def _download(session,url,destination):
    for _ in range(6):
        with session.get(url,stream=True,allow_redirects=False,timeout=(20,120)) as response:
            if response.status_code in (301,302,303,307,308):
                url=tranche_url(urljoin(url,response.headers.get('Location','')))
                continue
            response.raise_for_status()
            with destination.open('wb') as output:
                for chunk in response.iter_content(chunk_size=1024*1024):
                    if chunk:output.write(chunk)
            if not destination.stat().st_size:raise ValueError('O download ZINC está vazio: '+url)
            return url
    raise ValueError('Redirecionamentos excessivos no download ZINC: '+url)


@contextmanager
def _text(path):
    opener=gzip.open if path.suffix=='.gz' else bz2.open if path.suffix=='.bz2' else Path.open
    with opener(path,'rt',encoding='utf-8-sig') as stream:yield stream


def zinc_code(value):
    value=value.strip()
    if value.isdigit():value='ZINC'+value.zfill(12)
    match=re.fullmatch(r'ZINC\d+',value,re.I)
    if not match:raise ValueError('Identificador ZINC inválido: '+value)
    return value.upper()


def _smiles_rows(path):
    order=(0,1)
    with _text(path) as stream:
        for number,line in enumerate(stream,1):
            if not line.strip() or line.lstrip().startswith('#'):continue
            fields=line.split();lower=[f.lower() for f in fields]
            if any(f in ('smiles','smile','canonical_smiles') for f in lower):
                smiles=next(i for i,f in enumerate(lower) if f in ('smiles','smile','canonical_smiles'))
                identifier=next((i for i,f in enumerate(lower) if f in ('zinc_id','zincid','id','name','molecule_chembl_id')),None)
                if identifier is None:raise ValueError('O cabeçalho SMI está sem o identificador ZINC.')
                order=(smiles,identifier);continue
            if len(fields)<=max(order):raise ValueError(f'Linha {number}: tabela SMI ZINC incompleta.')
            yield dict(molecule_chembl_id=zinc_code(fields[order[1]]),canonical_smiles=fields[order[0]])


def _mol2_blocks(path):
    block=[]
    with _text(path) as stream:
        for line in stream:
            if line.strip()=='@<TRIPOS>MOLECULE':
                if block:yield ''.join(block)
                block=[]
            if block or line.strip()=='@<TRIPOS>MOLECULE':block.append(line)
        if block:yield ''.join(block)


def validate_download_workers(value):
    if type(value) is not int or not 1 <= value <= 16:
        raise ValueError('Downloads simultâneos deve ser um inteiro entre 1 e 16.')


def _download_job(index,url,temporary):
    """Each download owns its HTTP session; sessions are never shared by threads."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    destination=Path(temporary)/(str(index)+'_'+Path(urlsplit(url).path).name)
    with requests.Session() as session:
        session.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
        session.headers.update({'User-Agent':'BioMolExplorer/0.2 (ZINC tranche retrieval)'})
        resolved=_download(session,url,destination)
    return destination,resolved


@contextmanager
def _downloads(urls,temporary,workers):
    """Bound prefetch to the worker count and consume in deterministic list order."""
    pool=ThreadPoolExecutor(max_workers=workers,thread_name_prefix='zinc-download')
    pending=deque()
    remaining=iter(enumerate(urls))
    def submit():
        item=next(remaining,None)
        if item is not None:
            index,url=item
            pending.append((url,pool.submit(_download_job,index,url,temporary)))
    def results():
        while pending:
            url,future=pending.popleft()
            downloaded,resolved=future.result()
            yield url,downloaded,resolved
            downloaded.unlink()
            submit()
    try:
        for _ in range(workers):submit()
        yield results()
    finally:
        # Running requests finish before the temporary directory is removed.
        pool.shutdown(wait=True,cancel_futures=True)


def retrieve_tranches(list_file,output,download_workers=4):
    """Export compounds.csv and per-molecule MOL2 conformations in one bundle."""
    from rdkit import Chem
    from .input_validation import validate_file
    from .storage import write_json
    validate_download_workers(download_workers)
    urls=read_download_list(list_file)
    workers=min(download_workers,len(urls))
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    report={'provider':'ZINC','download_workers':workers,'downloads':[],'compounds':0,'duplicates':0,'conformations':0}
    seen={}
    fields=['molecule_chembl_id','canonical_smiles','source','source_url','structure_format','conformer_file','conformer_origin']
    with tempfile.TemporaryDirectory(prefix='bme-zinc-') as temporary, _downloads(urls,temporary,workers) as downloads:
        table=Path(temporary)/'compounds.csv'
        with table.open('w',encoding='utf-8',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
            def write(row,url,format,conformation=None):
                code=row['molecule_chembl_id'];molecule=Chem.MolFromSmiles(row['canonical_smiles'])
                if molecule is None or not molecule.GetNumAtoms():raise ValueError('SMILES ZINC inválido: '+code)
                row['canonical_smiles']=Chem.MolToSmiles(molecule)
                previous=seen.get(code)
                if previous is not None:
                    if previous[0]!=row['canonical_smiles']:raise ValueError('Estruturas ZINC conflitantes para '+code)
                    report['duplicates']+=1
                    # An existing 2D record can acquire its 3D file later; the CSV uses
                    # a deterministic path filled during final consolidation below.
                    if conformation is not None and not previous[1]:
                        target=output/'Conformers'/(code+'.mol2');target.parent.mkdir(exist_ok=True)
                        target.write_text(conformation,encoding='utf-8');seen[code]=(previous[0],True)
                        report['conformations']+=1
                    return
                if conformation is not None:
                    target=output/'Conformers'/(code+'.mol2');target.parent.mkdir(exist_ok=True)
                    target.write_text(conformation,encoding='utf-8');report['conformations']+=1
                seen[code]=(row['canonical_smiles'],conformation is not None)
                writer.writerow(dict(row,source='ZINC',source_url=url,structure_format=format,conformer_file='',conformer_origin=''))
                report['compounds']+=1
            for url,downloaded,resolved in downloads:
                format='mol2' if '.mol2' in downloaded.name.lower() else 'smi'
                before=report['compounds'];count=0
                if format=='smi':
                    for row in _smiles_rows(downloaded):write(row,url,format);count+=1
                else:
                    for block in _mol2_blocks(downloaded):
                        name=block.splitlines()[1].strip()
                        match=re.search(r'ZINC\d+',name,re.I)
                        code=zinc_code(match.group() if match else name)
                        structure=Path(temporary)/'molecule.mol2';structure.write_text(block,encoding='utf-8')
                        validate_file(structure,'prepared_structures')
                        molecule=Chem.MolFromMol2File(str(structure),removeHs=False)
                        if molecule is None:
                            from .docking_data import structure_smiles
                            smiles=structure_smiles(structure)
                        else:smiles=Chem.MolToSmiles(Chem.RemoveHs(molecule))
                        write(dict(molecule_chembl_id=code,canonical_smiles=smiles),url,format,block);count+=1
                if not count:raise ValueError('A tranche ZINC não contém moléculas: '+url)
                digest=hashlib.sha256()
                with downloaded.open('rb') as downloaded_stream:
                    for chunk in iter(lambda:downloaded_stream.read(1024*1024),b''):digest.update(chunk)
                report['downloads'].append(dict(url=url,resolved_url=resolved,format=format,records=count,
                    added=report['compounds']-before,sha256=digest.hexdigest()))
        final=output/'compounds.csv'
        with table.open(encoding='utf-8',newline='') as source,final.open('w',encoding='utf-8',newline='') as destination:
            writer=csv.DictWriter(destination,fieldnames=fields);writer.writeheader()
            for row in csv.DictReader(source):
                if seen[row['molecule_chembl_id']][1]:
                    row['conformer_file']=str((output/'Conformers'/(row['molecule_chembl_id']+'.mol2')).resolve())
                    row['structure_format']='mol2'
                    row['conformer_origin']='library'
                writer.writerow(row)
        validate_file(final,'compounds')
    write_json(report,output/'retrieval_report.json')
    return report
