"""Validate the MOL2 graph and numeric data before DOCK6's native routines."""
import math
from pathlib import Path


def compound_code(filename):
    name=Path(filename).name
    for suffix in ('.lig.min_scored.mol2','.lig.mol2','_scored.mol2','.mol2'):
        if name.endswith(suffix):return name.removesuffix(suffix)
    raise ValueError('Nome de arquivo MOL2 não reconhecido: '+name)


def validate_mol2(path, *, connected=True, scored=False, sanitize=True):
    from rdkit import Chem
    path=Path(path)
    text=path.read_text()
    if text.count('@<TRIPOS>MOLECULE')!=1:
        raise ValueError('O arquivo precisa conter exatamente uma molécula MOL2.')
    sections={};current=None
    for line in text.splitlines():
        if line.startswith('@<TRIPOS>'):
            current=line;sections[current]=[]
        elif current is not None and line.strip():sections[current].append(line.split())
    atoms=sections.get('@<TRIPOS>ATOM',[])
    bonds=sections.get('@<TRIPOS>BOND',[])
    header=sections.get('@<TRIPOS>MOLECULE',[])
    try:
        counts=list(map(int,header[1][:2]))
        if counts!=[len(atoms),len(bonds)] or not atoms:raise ValueError()
        identities=[int(a[0]) for a in atoms]
        if len(set(identities))!=len(identities):raise ValueError()
        for a in atoms:
            if len(a)<9 or not all(math.isfinite(float(a[i])) for i in (2,3,4,8)):raise ValueError()
        for b in bonds:
            if int(b[1]) not in identities or int(b[2]) not in identities or b[1]==b[2]:raise ValueError()
    except (IndexError,TypeError,ValueError) as error:
        raise ValueError('MOL2 inválido: contagens, identificadores, ligações, coordenadas ou cargas inconsistentes.') from error
    mol=Chem.MolFromMol2Block(text,sanitize=sanitize,removeHs=False)
    if mol is None:raise ValueError('Topologia ou tipos atômicos MOL2 inválidos para o docking.')
    # SYBYL P.3 is also used for neutral pentavalent phosphorus. RDKit's
    # MOL2 reader may assign +2 although the explicit five bond valences
    # describe the neutral phosphoryl group (phosphonium valence 4 stays +1).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==15 and atom.GetFormalCharge()==2 and sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())==5:
            atom.SetFormalCharge(0)
    Chem.SanitizeMol(mol)
    if connected and len(Chem.GetMolFrags(mol))!=1:
        raise ValueError('O ligante MOL2 contém fragmentos desconectados; revise a topologia da pose e os contraíons antes do docking flexível.')
    if scored:
        from .docking_data import score
        score(path,'dock6')
    return mol
