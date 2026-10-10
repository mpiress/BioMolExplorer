"""Human-readable scientific file contracts, without executing uploaded code."""
import ast
import csv
import hashlib
import math
import re
from pathlib import Path

ALIASES={'smiles':'canonical_smiles','smile':'canonical_smiles','Canonical_SMILES':'canonical_smiles','name':'molecule_chembl_id','zinc_id':'molecule_chembl_id'}
EXPECTED={
    'compounds':'CSV UTF-8 com canonical_smiles e molecule_chembl_id. Também aceitamos smiles e name; códigos ausentes são gerados. Exemplo: molecule_chembl_id,canonical_smiles\nMOL1,CCO',
    'graph_compounds':'CSV UTF-8 com molecule_chembl_id e canonical_smiles preenchidos. Também aceitamos name e smiles. Use os mesmos códigos de source/target; nós isolados podem ter códigos adicionais. Exemplo: molecule_chembl_id,canonical_smiles\nMOL1,CCO',
    'admet':'CSV de resultados com molecule_chembl_id, canonical_smiles, TPSA e WLOGP. Os dois últimos campos devem ser números finitos.',
    'fingerprints':'CSV com molecule_chembl_id e fingerprint. Cada fingerprint é uma lista de bits 0/1, por exemplo "[0, 1, 0]", com o mesmo tamanho em todas as linhas.',
    'similarity':'CSV com source,target,value. source e target são códigos; value é um número entre 0 e 1. Exemplo: MOL1,MOL2,0.85',
    'structures':'Arquivo .pdb com registros ATOM/HETATM e coordenadas válidas. Metadados opcionais: pdb_codes.csv com PDB_CODE,LIGAND,RESNUM,CHAIN.',
    'prepared_structures':'Receptor .pdbqt com ATOM/HETATM, pdb_codes.csv (PDB_CODE,LIGAND,RESNUM,CHAIN) e centers.csv (uma coluna por complexo, com três coordenadas). Selecione o conjunto completo.',
    'vina':'Pose .pdbqt com coordenadas ATOM/HETATM e a linha REMARK VINA RESULT com o score.',
    'dock6':'Resultado *_scored.mol2 com as seções @<TRIPOS>MOLECULE e @<TRIPOS>ATOM e um Grid_Score numérico.',
    'scores':'CSV com um código de composto e colunas numéricas de score (vina/dock6 ou score).',
    'visualization':'PNG/JPEG ou arquivo *.biomol-view.json exportado pela plataforma.',
    'other':'Arquivo de dados gerais. Para listas de download ZINC, escolha o tipo Lista de downloads ZINC.',
    'zinc_urls':'Lista de downloads ZINC (TXT/URI ou script exportado), com links SMI/MOL2, inclusive .gz ou .bz2. Links HTTP oficiais são convertidos para HTTPS; comandos não são executados.',
}


def output_kind(operation):
    return {'retrieve_structures':'structures','prepare_structures':'prepared_structures','redocking':'prepared_structures',
            'fingerprints':'fingerprints','similarity':'similarity','docking_vina':'vina','docking_dock6':'dock6',
            'consensus':'scores'}.get(operation,'compounds')


def contract(kind,operation=None):
    if kind=='compounds' and operation=='graphs':return EXPECTED['graph_compounds']
    return EXPECTED['admet' if kind=='compounds' and operation=='admet' else kind]


def columns(path):
    try:
        with Path(path).open(encoding='utf-8-sig',newline='') as stream:
            return {ALIASES.get(name,name) for name in next(csv.reader(stream),[])}
    except (OSError,UnicodeError):
        return set()


def _code(value):
    return isinstance(value,str) and bool(re.fullmatch('[A-Za-z0-9_.+-]{1,100}',value)) and '..' not in value


def validate_file(path,kind,operation=None,validate_rows=True):
    path=Path(path)
    try:
        if not path.is_file() or path.is_symlink():
            raise ValueError('O arquivo está indisponível.')
        if kind=='other':
            return
        if kind=='zinc_urls':
            from .zinc_retrieval import read_download_list
            read_download_list(path)
            return
        if kind=='visualization':
            if path.name.endswith('.biomol-view.json'):
                from .visualizations import load_view
                load_view(path.read_bytes())
            else:
                from PIL import Image
                with Image.open(path) as image:image.verify()
            return
        if kind in ('vina','dock6') and path.suffix=='.csv':
            if not {'molecule_chembl_id','score'}<=columns(path):
                raise ValueError('Tabela de docking incompleta.')
            validate_file(path,'scores',operation,validate_rows)
            return
        if kind in ('compounds','fingerprints','similarity','scores'):
            if path.suffix.lower()!='.csv':raise ValueError('Selecione um arquivo CSV.')
            with path.open(encoding='utf-8-sig',newline='') as stream:
                reader=csv.DictReader(stream)
                fields=[ALIASES.get(k,k) for k in reader.fieldnames or []]
                if len(fields)!=len(set(fields)):raise ValueError('O cabeçalho possui colunas duplicadas.')
                required={'canonical_smiles'} if kind=='compounds' else {'molecule_chembl_id','fingerprint'} if kind=='fingerprints' else {'source','target','value'} if kind=='similarity' else set()
                if kind=='compounds' and operation=='admet':required|={'molecule_chembl_id','TPSA','WLOGP'}
                if not required.issubset(fields):raise ValueError('Faltam as colunas: '+', '.join(sorted(required-set(fields))))
                if kind=='scores' and (not set(fields)&{'molecule_chembl_id','molecule','id','Unnamed: 0'} or not set(fields)&{'vina','dock6','score'}):
                    raise ValueError('Inclua a coluna de código e ao menos uma coluna de score.')
                if not validate_rows:return
                seen={};width=None
                if kind=='compounds':
                    from rdkit import Chem
                for number,row in enumerate(reader,2):
                    if None in row or any(value is None for value in row.values()):raise ValueError(f'Linha {number}: quantidade de campos diferente do cabeçalho.')
                    row={ALIASES.get(k,k):v for k,v in row.items()}
                    identifier=row.get('molecule_chembl_id')
                    if identifier and not _code(identifier):raise ValueError(f'Linha {number}: código do composto inválido.')
                    if kind=='compounds':
                        smiles=row['canonical_smiles']
                        molecule=Chem.MolFromSmiles(smiles)
                        if molecule is None or not molecule.GetNumAtoms():raise ValueError(f'Linha {number}: SMILES inválido.')
                        canonical=Chem.MolToSmiles(molecule)
                        if identifier in seen and seen[identifier]!=canonical:raise ValueError(f'Linha {number}: o código {identifier} representa moléculas diferentes.')
                        if identifier:seen[identifier]=canonical
                        if operation=='admet' and any(not math.isfinite(float(row[key])) for key in ('TPSA','WLOGP')):
                            raise ValueError(f'Linha {number}: TPSA/WLOGP inválido.')
                    elif kind=='fingerprints':
                        if not identifier:raise ValueError(f'Linha {number}: código vazio.')
                        bits=ast.literal_eval(row['fingerprint'])
                        if not isinstance(bits,(list,tuple)) or not bits or any(type(v) is not int or v not in (0,1) for v in bits):raise ValueError(f'Linha {number}: fingerprint inválido.')
                        if width is not None and width!=len(bits):raise ValueError(f'Linha {number}: os fingerprints têm tamanhos diferentes.')
                        width=len(bits)
                    elif kind=='similarity':
                        value=float(row['value'])
                        if not _code(row['source']) or not _code(row['target']) or not math.isfinite(value) or not 0<=value<=1:raise ValueError(f'Linha {number}: relação de similaridade inválida.')
                    else:
                        for field in set(fields)&{'vina','dock6','score','z-score','min-max'}:
                            if not math.isfinite(float(row[field])):raise ValueError(f'Linha {number}: score inválido.')
            return
        if path.name=='pdb_codes.csv':
            if not {'PDB_CODE','LIGAND','RESNUM','CHAIN'}.issubset(columns(path)):
                raise ValueError('pdb_codes.csv está sem as colunas de metadados.')
            with path.open(encoding='utf-8-sig',newline='') as stream:
                for row in csv.DictReader(stream):
                    if any(not _code(row[k]) for k in ('PDB_CODE','LIGAND','CHAIN')):raise ValueError('Código PDB, ligante ou cadeia inválido.')
                    int(row['RESNUM'])
            return
        if path.name=='centers.csv':
            with path.open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.reader(stream))
            if len(rows)!=4 or not rows[0] or any(len(row)!=len(rows[0]) for row in rows):raise ValueError('centers.csv precisa de cabeçalho e três linhas de coordenadas.')
            if any(not math.isfinite(float(value)) for row in rows[1:] for value in row):raise ValueError('Centro de docking inválido.')
            return
        if kind=='prepared_structures' and path.suffix.lower()=='.mol2':
            text=path.read_text(encoding='utf-8')
            if '@<TRIPOS>MOLECULE' not in text or '@<TRIPOS>ATOM' not in text:
                raise ValueError('Estrutura MOL2 preparada incompleta.')
            atoms=text.split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0].strip().splitlines()
            if not atoms or any(len(line.split())<6 or any(not math.isfinite(float(value)) for value in line.split()[2:5]) for line in atoms):
                raise ValueError('Coordenadas MOL2 inválidas.')
            return
        if kind in ('structures','prepared_structures','vina'):
            allowed={'.pdb'} if kind=='structures' else {'.pdb','.pdbqt'} if kind=='prepared_structures' else {'.pdbqt'}
            if path.suffix.lower() not in allowed:raise ValueError('Extensão de estrutura inválida.')
            count=0;score=False
            with path.open(encoding='utf-8') as stream:
                for line in stream:
                    if line.startswith(('ATOM  ','HETATM')):
                        if len(line)<54 or any(not math.isfinite(float(line[a:b])) for a,b in ((30,38),(38,46),(46,54))):raise ValueError('Coordenadas ATOM/HETATM inválidas.')
                        count+=1
                    if line.startswith('REMARK VINA RESULT:'):
                        score=math.isfinite(float(line.split(':',1)[1].split()[0]))
            if not count:raise ValueError('Não foram encontrados átomos com coordenadas.')
            if kind=='vina' and not score:raise ValueError('A pose está sem o score Vina.')
            return
        if kind=='dock6':
            text=path.read_text()
            if not path.name.endswith('_scored.mol2') or any(marker not in text for marker in ('@<TRIPOS>MOLECULE','@<TRIPOS>ATOM')) or not re.search(r'Grid_Score:\s*-?\d',text):raise ValueError('Resultado DOCK6 incompleto.')
    except (ValueError,SyntaxError,UnicodeError,ImportError,KeyError,TypeError,IndexError,OSError) as exc:
        raise ValueError(f'{path.name}: {exc}\nPadrão esperado: {contract(kind,operation)}') from None


def validate_bundle(paths,kind,operation=None):
    if not paths:raise ValueError('Selecione seus arquivos. Padrão esperado: '+contract(kind,operation))
    data=[]
    for path in paths:
        if Path(path).suffix.lower() in ('.png','.jpg','.jpeg') or str(path).endswith('.biomol-view.json'):
            validate_file(path,'visualization')
        else:
            validate_file(path,kind,operation);data.append(Path(path))
    if not data:raise ValueError('Inclua os dados além do gráfico. Padrão esperado: '+contract(kind,operation))
    if kind=='prepared_structures' and (not {'pdb_codes.csv','centers.csv'}.issubset({p.name for p in data}) or not any(p.suffix=='.pdbqt' for p in data)):
        raise ValueError('O conjunto preparado está incompleto. Padrão esperado: '+contract(kind,operation))


def merge_csv(paths,destination,kind='compounds'):
    """Union selected datasets, preserving metadata and rejecting ID conflicts."""
    paths=list(dict.fromkeys(map(Path,paths)))
    fields=[];seen={};fingerprint_width=None
    for path in paths:
        validate_file(path,kind)
        with path.open(encoding='utf-8-sig',newline='') as stream:
            reader=csv.DictReader(stream)
            for name in reader.fieldnames:
                name=ALIASES.get(name,name)
                if name not in fields:fields.append(name)
    if kind=='compounds' and 'molecule_chembl_id' not in fields:fields.insert(0,'molecule_chembl_id')
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    from tempfile import NamedTemporaryFile
    temporary=None
    try:
      with NamedTemporaryFile(mode='w',encoding='utf-8',newline='',dir=destination.parent,delete=False) as out:
        temporary=Path(out.name)
        writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader()
        for path in paths:
            with path.open(encoding='utf-8-sig',newline='') as stream:
                for original in csv.DictReader(stream):
                    row={ALIASES.get(k,k):v for k,v in original.items()}
                    if kind=='compounds':
                        from rdkit import Chem
                        signature=Chem.MolToSmiles(Chem.MolFromSmiles(row['canonical_smiles']))
                        if not row.get('molecule_chembl_id'):
                            row['molecule_chembl_id']='USER_'+hashlib.sha256(signature.encode()).hexdigest()[:12]
                    else:
                        signature=tuple(ast.literal_eval(row['fingerprint'])) if kind=='fingerprints' else float(row['value']) if kind=='similarity' else None
                        if kind=='fingerprints':
                            if fingerprint_width is not None and len(signature)!=fingerprint_width:
                                raise ValueError('Os arquivos selecionados possuem fingerprints com tamanhos diferentes. Padrão esperado: '+contract(kind))
                            fingerprint_width=len(signature)
                    key=row['molecule_chembl_id'] if kind!='similarity' else (row['source'],row['target'])
                    if key in seen:
                        if seen[key]!=signature:raise ValueError(f'O código {key} possui dados conflitantes entre os arquivos selecionados. Use códigos únicos.')
                        continue
                    seen[key]=signature;writer.writerow(row)
      temporary.replace(destination)
    finally:
      if temporary is not None:temporary.unlink(missing_ok=True)
    return destination
