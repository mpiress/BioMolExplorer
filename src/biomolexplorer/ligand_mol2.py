"""Repair Chimera MOL2 heavy-bond orders from named PDB CCD definitions.

Coordinates, atom inventory, hydrogen attachments and partial charges are kept.
Definitions are cached; unknown atoms or an invalid repaired graph are rejected.
"""
import hashlib
import io
import json
import os
from pathlib import Path
import re

from .dock6_validation import validate_mol2

PREPARATION_SCHEMA_VERSION=3


def component_definition(component):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from .paths import resolve_path
    from .storage import write_text
    if not re.fullmatch(r'[A-Z0-9]{1,5}',component):raise ValueError('Código CCD do ligante inválido.')
    root=Path(os.environ.get('BIOMOL_CACHE_DIR',resolve_path('.biomolexplorer/cache')))
    path=root/'calculations/ccd'/(component+'.cif')
    if path.is_file():text=path.read_text()
    else:
        import requests
        try:
            response=requests.get('https://files.rcsb.org/ligands/download/'+component+'.cif',timeout=30)
            response.raise_for_status();text=response.text
        except requests.RequestException as error:
            raise ValueError('Não foi possível obter a definição CCD do ligante '+component+'.') from error
    data=MMCIF2Dict(io.StringIO(text))
    if data.get('_chem_comp.id')!=[component]:raise ValueError('A definição CCD não corresponde ao ligante selecionado.')
    if not path.is_file():write_text(text,path)
    return data,hashlib.sha256(text.encode()).hexdigest()


def ensure_ligand_mol2(path):
    """Validate ordinary exports; repair a failed export only with exact CCD names."""
    path=Path(path)
    try:
        validate_mol2(path)
        return False
    except ValueError:
        # Retain numeric/count/connectivity validation; never repair malformed bytes.
        validate_mol2(path,sanitize=False)
    match=re.fullmatch(r'.+_([A-Za-z0-9]{1,5})_[+-]?\d+[^.]*\.lig\.mol2',path.name)
    if not match:raise ValueError('Não foi possível identificar o componente CCD do ligante MOL2.')
    component=match[1].upper();definition,sha=component_definition(component)
    symbols=dict(zip(definition['_chem_comp_atom.atom_id'],definition['_chem_comp_atom.type_symbol']))
    lines=path.read_text().splitlines();atoms={};names={};bonds=[];section=None
    for line in lines:
        if line.startswith('@<TRIPOS>'):section=line
        elif line.strip() and section=='@<TRIPOS>ATOM':
            a=line.split();atoms[a[0]]=a
            if a[1] in names:raise ValueError('Nomes atômicos duplicados no ligante MOL2.')
            names[a[1]]=a[0]
        elif line.strip() and section=='@<TRIPOS>BOND':bonds.append(line.split())
    heavy={a[1] for a in atoms.values() if a[5].split('.')[0] not in ('H','D')}
    if any(symbols.get(name,'').upper()!=atoms[names[name]][5].split('.')[0].upper() for name in heavy):
        raise ValueError('Os nomes e elementos do ligante não correspondem à definição CCD.')
    original={frozenset((b[1],b[2])):b[3] for b in bonds}
    corrected=[]
    for first,second,order in zip(definition['_chem_comp_bond.atom_id_1'],definition['_chem_comp_bond.atom_id_2'],definition['_chem_comp_bond.value_order']):
        if first not in heavy or second not in heavy:continue
        kind={'SING':'1','DOUB':'2','TRIP':'3','AROM':'ar'}.get(order)
        if kind is None:raise ValueError('Ordem de ligação CCD não suportada no ligante.')
        ids=(names[first],names[second])
        if kind=='1' and original.get(frozenset(ids))=='am':kind='am'
        corrected.append((*ids,kind))
    # Chimera supplies protonation and hydrogen geometry. Do not rebuild either.
    corrected.extend((b[1],b[2],b[3]) for b in bonds if any(atoms[i][5].split('.')[0] in ('H','D') for i in b[1:3]))
    output=[];section=None
    for index,line in enumerate(lines):
        if line=='@<TRIPOS>BOND':
            output.append(line);output.extend(f'{i} {a} {b} {kind}' for i,(a,b,kind) in enumerate(corrected,1));section=line;continue
        if line.startswith('@<TRIPOS>'):section=line
        elif section=='@<TRIPOS>BOND':continue
        if index==2:
            fields=line.split();fields[1]=str(len(corrected));line=' '.join(fields)
        output.append(line)
    temporary=path.with_suffix('.checked.mol2')
    try:
        temporary.write_text('\n'.join(output)+'\n');mol=validate_mol2(temporary)
        if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            raise ValueError('O preparo do ligante perdeu hidrogênios; refaça o preparo com a topologia CCD antes de aplicar cargas.')
        temporary.replace(path)
    finally:temporary.unlink(missing_ok=True)
    missing=sorted(name for name,element in symbols.items() if element not in ('H','D') and name not in heavy)
    path.with_suffix('.export.json').write_text(json.dumps(dict(component=component,definition_sha256=sha,
        missing_heavy_atoms=missing,scope='observed_fragment' if missing else 'complete_component',
        repair='CCD heavy-bond export; prepared coordinates and partial charges preserved'),indent=2)+'\n')
    return True


def prepare_chimera_ligand(script):
    """Supply CCD topology before Chimera adds H/charges, not after chemistry.

    PDB intermediates carry coordinates and connectivity but no bond orders.
    Reopening one can misclassify saturated carbon and lose required hydrogens.
    """
    from rdkit import Chem
    from .dock6_validation import validate_mol2
    script=Path(script);text=script.read_text()
    prefix=re.match(r'open (.+\.complex\.pdb)\nselect :([+-]?\d+)\.([A-Za-z0-9_+-]+)\nselect invert\ndelete selected\n',text)
    if not prefix:return False
    output=re.search(r'^write format mol2 #0 (.+\.lig\.mol2)$',text,re.M)
    if not output:return False
    component=re.fullmatch(r'.+_([A-Za-z0-9]{1,5})_[+-]?\d+[^.]*\.lig\.mol2',Path(output[1]).name)
    if not component:return False
    component=component[1].upper();definition,sha=component_definition(component)
    symbols=dict(zip(definition['_chem_comp_atom.atom_id'],definition['_chem_comp_atom.type_symbol']))
    charges=dict(zip(definition['_chem_comp_atom.atom_id'],definition['_chem_comp_atom.charge']))
    positions={}
    for line in Path(prefix[1]).read_text().splitlines():
        if line.startswith('ENDMDL'):break
        if not line.startswith(('ATOM  ','HETATM')):continue
        if (line[17:20].strip(),line[22:26].strip(),line[21:22].strip())!=(component,prefix[2],prefix[3]):continue
        if line[16:17] not in (' ','A'):continue
        name=line[12:16].strip();element=line[76:78].strip().upper()
        if element in ('H','D'):continue
        if name not in symbols or name in positions:raise ValueError('Os átomos do ligante não correspondem à definição CCD.')
        if element and element!=symbols[name].upper():raise ValueError('Elemento atômico incompatível com a definição CCD.')
        positions[name]=tuple(float(line[i:i+8]) for i in (30,38,46))
    if not positions:raise ValueError('Ligante cristalográfico ausente na referência de RMSD.')
    graph=Chem.RWMol();names=list(positions);ids={}
    for name in names:
        atom=Chem.Atom(symbols[name].capitalize());atom.SetFormalCharge(int(charges[name]));ids[name]=graph.AddAtom(atom)
    bondtypes={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for first,second,order in zip(definition['_chem_comp_bond.atom_id_1'],definition['_chem_comp_bond.atom_id_2'],definition['_chem_comp_bond.value_order']):
        if first in positions and second in positions:graph.AddBond(ids[first],ids[second],bondtypes[order])
    mol=graph.GetMol();Chem.SanitizeMol(mol)
    conformer=Chem.Conformer(len(names))
    for name,xyz in positions.items():conformer.SetAtomPosition(ids[name],xyz)
    mol.AddConformer(conformer);Chem.AssignAtomChiralTagsFromStructure(mol)
    def atom_type(atom):
        symbol=atom.GetSymbol()
        if symbol in ('C','N') and atom.GetIsAromatic():return symbol+'.ar'
        if symbol=='N' and atom.GetFormalCharge()>0 and atom.GetTotalDegree()==4:return 'N.4'
        if symbol=='N' and any(n.GetSymbol()=='C' and any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol()=='O' for b in n.GetBonds()) for n in atom.GetNeighbors()):return 'N.am'
        hybrid=str(atom.GetHybridization())
        if symbol in ('C','N'):return symbol+('.1' if hybrid=='SP' else '.2' if hybrid=='SP2' else '.3')
        if symbol=='O':
            if atom.GetFormalCharge()<0:
                return 'O.co2' if any(n.GetSymbol()=='C' and sum(x.GetSymbol()=='O' for x in n.GetNeighbors())>=2 for n in atom.GetNeighbors()) else 'O.2' if hybrid=='SP2' else 'O.3'
            return 'O.2' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in atom.GetBonds()) else 'O.3'
        if symbol=='S':
            oxy=sum(n.GetSymbol()=='O' and mol.GetBondBetweenAtoms(n.GetIdx(),atom.GetIdx()).GetBondType()==Chem.BondType.DOUBLE for n in atom.GetNeighbors())
            return 'S.o2' if oxy>=2 else 'S.o' if oxy else 'S.2' if hybrid=='SP2' else 'S.3'
        return 'P.3' if symbol=='P' else symbol
    lines=['@<TRIPOS>MOLECULE',component,f'{len(names)} {mol.GetNumBonds()} 1 0 0','SMALL','USER_CHARGES','','@<TRIPOS>ATOM']
    for index,name in enumerate(names):
        xyz=positions[name];atom=mol.GetAtomWithIdx(index)
        lines.append(f'{index+1} {name} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f} {atom_type(atom)} 1 {component} {atom.GetFormalCharge():.4f}')
    lines.append('@<TRIPOS>BOND')
    for index,bond in enumerate(mol.GetBonds()):
        kind='ar' if bond.GetIsAromatic() else str(int(bond.GetBondTypeAsDouble()))
        lines.append(f'{index+1} {bond.GetBeginAtomIdx()+1} {bond.GetEndAtomIdx()+1} {kind}')
    lines.extend(['@<TRIPOS>SUBSTRUCTURE',f'1 {component} 1 RESIDUE 4 {prefix[3]} {component} 0 ROOT'])
    destination=Path(output[1]).with_suffix('.input.mol2')
    # The CCD graph has already been sanitized above. MOL2 cannot encode an
    # implicit pyrrole H reliably; validate its syntax here and the full graph
    # strictly after Chimera applies the requested hydrogen preparation.
    destination.write_text('\n'.join(lines)+'\n');validate_mol2(destination,sanitize=False)
    # Keep explicit source hydrogens when requested, using a PDB overlay would
    # lose topology again. Current normal preparation removes them explicitly.
    if 'delete element.H' not in text:
        hydrogen_lines=[l for l in Path(prefix[1]).read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and
            (l[17:20].strip(),l[22:26].strip(),l[21:22].strip())==(component,prefix[2],prefix[3]) and l[76:78].strip().upper() in ('H','D')]
        if hydrogen_lines:
            destination.unlink();return False
    # Classic Chimera re-infers IDATM types from 3D geometry even for MOL2.
    # Pin only types derived from the validated CCD graph before addh/addcharge.
    types={}
    for name,index in ids.items():
        atom=mol.GetAtomWithIdx(index);symbol=atom.GetSymbol();hybrid=str(atom.GetHybridization())
        if symbol=='C':types[name]='Car' if atom.GetIsAromatic() else 'C3' if hybrid=='SP3' else 'C1' if hybrid=='SP' else 'C2'
        elif symbol=='N':
            if atom.GetFormalCharge()>0:types[name]='N3+' if atom.GetTotalDegree()==4 else 'N2+'
            elif hybrid=='SP3':types[name]='N3'
            elif hybrid=='SP':types[name]='N1'
            elif atom.GetIsAromatic():types[name]='Npl' if atom.GetTotalNumHs() or atom.GetDegree()==3 else 'N2'
            else:types[name]='N2' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in atom.GetBonds()) else 'Npl'
        elif symbol=='O':
            types[name]='O2-' if atom.GetFormalCharge()<0 else 'O2' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in atom.GetBonds()) else 'O3'
    typing=destination.with_suffix('.types.py')
    typing.write_text('import chimera\natom_types = '+repr(types)+'\nfor model in chimera.openModels.list(modelTypes=[chimera.Molecule]):\n    for atom in model.atoms:\n        if atom.name in atom_types:\n            atom.idatmType = atom_types[atom.name]\n')
    script.write_text('open '+str(destination)+'\nopen '+str(typing)+'\n'+text[prefix.end():])
    missing=sorted(name for name,symbol in symbols.items() if symbol not in ('H','D') and name not in positions)
    destination.with_suffix('.json').write_text(json.dumps(dict(component=component,definition_sha256=sha,
        missing_heavy_atoms=missing,scope='observed_fragment' if missing else 'complete_component',
        topology='CCD before Chimera hydrogen and charge assignment'),indent=2)+'\n')
    return True


def restore_pdbqt_names(source,destination):
    """Restore atom names after Open Babel pH conversion using unchanged positions.

    This runs before docking, where coordinates prove the one-to-one mapping;
    it must never be used to infer identities in a calculated docking pose.
    """
    import math
    source,destination=Path(source),Path(destination)
    section=None;heavy=[]
    for line in source.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):section=line
        elif section=='@<TRIPOS>ATOM' and line.strip():
            a=line.split()
            if a[5].split('.')[0] not in ('H','D'):heavy.append((a[1],a[5].split('.')[0].upper(),tuple(map(float,a[2:5]))))
    used=set();lines=[]
    elements={'A':'C','NA':'N','OA':'O','SA':'S'}
    for line in destination.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            kind=line.split()[-1].upper()
            if kind not in ('H','HD','HS'):
                xyz=tuple(float(line[i:i+8]) for i in (30,38,46));element=elements.get(kind,kind)
                matches=[(index,name) for index,(name,symbol,point) in enumerate(heavy) if symbol==element and math.dist(xyz,point)<.002]
                if len(matches)!=1 or matches[0][0] in used:raise ValueError('Mapeamento dos átomos pesados da pose incompleto.')
                index,name=matches[0];used.add(index)
                if len(name)>4:raise ValueError('Nome atômico não suportado no formato PDBQT.')
                line=line[:12]+f'{name:>4}'+line[16:]
        lines.append(line)
    if len(used)!=len(heavy):raise ValueError('Mapeamento dos átomos pesados da pose incompleto.')
    destination.write_text('\n'.join(lines)+'\n')
