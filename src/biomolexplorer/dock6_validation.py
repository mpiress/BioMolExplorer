"""Validate the MOL2 graph and numeric data before DOCK6's native routines."""
import math
from pathlib import Path


def compound_code(filename):
    name=Path(filename).name
    for suffix in ('.lig.min_scored.mol2','.lig.mol2','_scored.mol2','.mol2'):
        if name.endswith(suffix):return name.removesuffix(suffix)
    raise ValueError('Nome de arquivo MOL2 não reconhecido: '+name)


def validate_mol2(path, *, connected=True, scored=False):
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
    if any(a[5] in ('N.4','O.co2') for a in atoms) or '@<TRIPOS>UNITY_ATOM_ATTR' in sections:
        # RDKit's MOL2 reader does not reliably preserve these standard SYBYL
        # valence/charge conventions (including UNITY_ATOM_ATTR). Decode the
        # graph with Open Babel instead of rejecting valid cations/carboxylates.
        from openbabel import openbabel as ob
        allowed={'C.3','C.2','C.1','C.ar','C.cat','N.3','N.2','N.1','N.ar','N.am','N.pl3','N.4',
            'O.3','O.2','O.co2','S.3','S.2','S.O','S.o','S.O2','S.o2','P.3','F','Cl','Br','I','H','H.spc','H.t3p'}
        if any(a[5] not in allowed for a in atoms) or any(b[3] not in ('1','2','3','ar','am') for b in bonds):
            raise ValueError('Topologia ou tipos atômicos MOL2 inválidos para o docking.')
        mol=Chem.MolFromMol2Block(text,sanitize=True,removeHs=False)
        if mol is None:
            converter=ob.OBConversion();converter.SetInAndOutFormats('mol2','sdf');native=ob.OBMol()
            mol=Chem.MolFromMolBlock(converter.WriteString(native),removeHs=False) if converter.ReadString(native,text) else None
        if mol is not None and (mol.GetNumAtoms()!=len(atoms) or mol.GetNumBonds()!=len(bonds)):mol=None
    else:
        mol=Chem.MolFromMol2Block(text,sanitize=True,removeHs=False)
        # RDKit can reject otherwise valid charged SYBYL records whose formal
        # charge is encoded in a UNITY block (for example a protonated N.2).
        # Open Babel understands that extension and is already the canonical
        # graph bridge used by the native MOL2 loader.
        if mol is None:
            from openbabel import openbabel as ob
            converter=ob.OBConversion();converter.SetInAndOutFormats('mol2','sdf');native=ob.OBMol()
            mol=Chem.MolFromMolBlock(converter.WriteString(native),removeHs=False) if converter.ReadString(native,text) else None
            if mol is not None and (mol.GetNumAtoms()!=len(atoms) or mol.GetNumBonds()!=len(bonds)):mol=None
    if mol is None:
        # MOL2 has no portable formal-charge field for an iminium/nitrone N.2
        # or an anionic O.2. Resolve their closed-shell valences explicitly;
        # bond orders and coordinates are never inferred or changed here.
        candidate=Chem.MolFromMol2Block(text,sanitize=False,removeHs=False)
        if candidate is not None:
            for atom in candidate.GetAtoms():
                valence=sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
                if atom.GetAtomicNum()==7 and valence==4:atom.SetFormalCharge(1)
                elif atom.GetAtomicNum()==8 and valence==1:atom.SetFormalCharge(-1)
            try:
                Chem.SanitizeMol(candidate)
                if candidate.GetNumAtoms()==len(atoms) and candidate.GetNumBonds()==len(bonds):mol=candidate
            except ValueError:pass
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
