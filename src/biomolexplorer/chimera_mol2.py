"""Canonicalize Chimera's over-valent MOL2 carboxylate resonance export."""
from pathlib import Path


def normalize_carboxylates(path):
    """Keep coordinates/types/partial charges; choose one C=O resonance bond.

    Classic Chimera WriteMol2 infers bond orders from atom geometry and can
    emit two double bonds for a C.2 bonded to two terminal O.co2 atoms. This
    is the same carboxylate graph, but its integer bond orders must sum to
    four at carbon for DOCK6 and molecular readers.
    """
    path=Path(path);lines=path.read_text().splitlines(keepends=True)
    section=None;atoms={};bonds=[];neighbors={}
    for index,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):section=line.strip();continue
        fields=line.split()
        if section=='@<TRIPOS>ATOM' and len(fields)>=9:atoms[fields[0]]=fields[5]
        elif section=='@<TRIPOS>BOND' and len(fields)>=4:
            bonds.append((index,fields))
            for atom,other in ((fields[1],fields[2]),(fields[2],fields[1])):
                neighbors.setdefault(atom,[]).append((other,index,fields))
    changed=0
    for atom,kind in atoms.items():
        if kind!='C.2':continue
        attached=neighbors.get(atom,[])
        oxygen=[entry for entry in attached if atoms.get(entry[0])=='O.co2' and len(neighbors.get(entry[0],[]))==1]
        if len(attached)!=3 or len(oxygen)!=2 or any(entry[2][3]!='2' for entry in oxygen):continue
        other=[entry for entry in attached if entry not in oxygen]
        if other[0][2][3]!='1':continue
        # Stable choice among chemically equivalent oxygens, without moving atoms.
        _,index,fields=sorted(oxygen,key=lambda entry:int(entry[0]))[-1]
        rewritten=list(fields);rewritten[3]='1'
        lines[index]=' '.join(rewritten)+'\n';changed+=1
    if changed:path.write_text(''.join(lines))
    return changed


def canonicalize_component(path,folder,record,config):
    """Repair exported heavy-bond orders from CCD without changing atom data.

    Chimera assigns GAFF types/charges for minimization but WriteMol2 infers
    SYBYL bond orders again from geometry. Keep every prepared atom, coordinate
    and partial charge; use the component definition for heavy-bond orders.
    """
    from .structure_selection import component_path,chemical_template
    from rdkit import Chem
    definition=component_path(record[1],folder,config)
    molecule,_=chemical_template(definition,record[1])
    Chem.Kekulize(molecule,clearAromaticFlags=True)
    path=Path(path);lines=path.read_text().splitlines();section=None;atoms={};hydrogens=[];bonds=[]
    for line in lines:
        if line.startswith('@<TRIPOS>'):section=line;continue
        fields=line.split()
        if section=='@<TRIPOS>ATOM' and len(fields)>=9:atoms[fields[0]]=fields
        elif section=='@<TRIPOS>BOND' and len(fields)>=4:bonds.append(fields)
    heavy={}
    for identity,fields in atoms.items():
        if fields[5].split('.')[0]=='H':continue
        if fields[1] in heavy:raise ValueError('Átomos alternativos duplicados no ligante preparado; revise a seleção estrutural.')
        heavy[fields[1]]=identity
    expected={a.GetProp('_TriposAtomName') for a in molecule.GetAtoms()}
    if set(heavy)-expected:raise ValueError('Os átomos pesados do ligante não correspondem à definição química do componente.')
    missing=sorted(expected-set(heavy));edits={}
    for atom in molecule.GetAtoms():
        name=atom.GetProp('_TriposAtomName')
        if name in heavy and atom.GetAtomicNum()==8 and atom.GetFormalCharge()==-1 and any(n.GetAtomicNum()==7 for n in atom.GetNeighbors()):
            identity=heavy[name]
            if atoms[identity][5]=='O.co2':
                edits[identity]=list(atoms[identity]);edits[identity][5]='O.3'
    for fields in bonds:
        if any(atoms[a][5].split('.')[0]=='H' for a in fields[1:3]):hydrogens.append(fields[1:4])
    original={frozenset(fields[1:3]):fields[3] for fields in bonds}
    phosphorus={identity for identity,fields in atoms.items() if fields[5].split('.')[0]=='P'}
    rebuilt=[]
    for bond in molecule.GetBonds():
        names=[bond.GetBeginAtom().GetProp('_TriposAtomName'),bond.GetEndAtom().GetProp('_TriposAtomName')]
        if any(name not in heavy for name in names):continue
        left,right=[heavy[name] for name in names]
        order=str(int(bond.GetBondTypeAsDouble()))
        # Keep the phosphorus tautomer/protonation actually prepared by
        # Chimera, especially in truncated organophosphate adduct references.
        # Imposing the full free component's P=O on an observed P-OH would
        # create an over-valent oxygen and alter the chosen chemical state.
        if left in phosphorus or right in phosphorus:
            order=original.get(frozenset((left,right)),order)
        if order=='1' and any(atoms[a][5]=='N.am' for a in (left,right)) and any(atoms[a][5]=='C.2' for a in (left,right)):order='am'
        rebuilt.append([left,right,order])
    rebuilt+=hydrogens
    output=[];section=None;written=False;header_row=0
    for line in lines:
        if line.startswith('@<TRIPOS>'):
            if section=='@<TRIPOS>BOND' and not written:
                output += [str(i)+' '+' '.join(b) for i,b in enumerate(rebuilt,1)];written=True
            section=line;header_row=0;output.append(line);continue
        if section=='@<TRIPOS>BOND':continue
        if section=='@<TRIPOS>ATOM' and line.split() and line.split()[0] in edits:
            line=' '.join(edits[line.split()[0]])
        if section=='@<TRIPOS>MOLECULE' and line.strip():
            header_row+=1
            if header_row==2:
                values=line.split();values[1]=str(len(rebuilt));line=' '.join(values)
        output.append(line)
    if section=='@<TRIPOS>BOND' and not written:output += [str(i)+' '+' '.join(b) for i,b in enumerate(rebuilt,1)]
    path.write_text('\n'.join(output)+'\n')
    from .storage import write_json
    from .stage_cache import file_digest
    write_json(dict(component=record[1],definition_sha256=file_digest(definition),
        missing_heavy_atoms=missing,scope='observed_fragment' if missing else 'complete_component',
        prepared_phosphorus_bonds_preserved=bool(phosphorus),
        repair='CCD heavy-bond export; prepared coordinates and partial charges preserved'),
        path.with_suffix('.export.json'))
