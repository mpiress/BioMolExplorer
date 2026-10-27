"""Targeted Chimera IDATM repair from a trusted chemical component graph."""
import re


def local_definition(folder,code,config):
    """Fingerprint available chemistry without fetching it during cache lookup."""
    import os
    from pathlib import Path
    from .stage_cache import file_digest
    folder=Path(folder)
    explicit=config.get('chemistry_templates',{}).get(code)
    if explicit:
        path=Path(explicit)
        if not path.is_absolute():path=folder/path
        return {str(path):file_digest(path) if path.is_file() else None}
    roots=[folder/'chemistry',folder]
    if os.environ.get('BIOMOL_CCD_PATH'):roots.insert(0,Path(os.environ['BIOMOL_CCD_PATH']))
    preparation=os.environ.get('BIOMOL_PREPARATION_CACHE')
    cache=Path(preparation).parent/'calculations' if preparation else folder/'.calculation-cache'
    roots.append(cache/'ccd')
    for root in roots:
        for suffix in ('.cif','.sdf','.mol2'):
            path=root/(code+suffix)
            if path.is_file():return {code:file_digest(path)}
    return {}


def charge_commands(folder,record,config):
    """Restore CCD carbon hybridization, oxime hydroxyls and N+(=C)-O- geometry.

    Use an available local definition before AddH, or retrieve one after a
    failed charge assignment. Repair only CCD-confirmed motifs; never switch
    the requested charge method to bypass an AM1 failure.
    """
    from .structure_selection import component_path,chemical_template,source_structure,read_selected
    molecule,_=chemical_template(component_path(record[1],folder,config),record[1])
    result=[]
    chain=config.get('ligand_chain',record[3]);icode=config.get('ligand_insertion_code','')
    present=None
    try:raw=source_structure(folder,record[0])
    except ValueError:pass
    else:
        structure,_=read_selected(raw,config.get('structure',{}))
        present={atom.name for c in structure[0] for residue in c for atom in residue
            if c.name==chain and residue.name==record[1] and residue.seqid.num==int(record[2])
            and residue.seqid.icode.strip()==icode}
    def append(atom,kind):
        name=atom.GetProp('_TriposAtomName')
        if not re.fullmatch(r"[A-Za-z0-9_']+",name):raise ValueError('Nome atômico inválido na definição química do ligante.')
        if present is not None and name not in present:return
        result.append(f'setattr a idatmType {kind} :{record[2]}{icode}.{chain}@{name}')
    # Crystallographic geometry can make IDATM misclassify saturated carbons,
    # yielding too few hydrogens and an odd electron count for closed-shell AM1.
    from rdkit import Chem
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum()==6 and atom.GetFormalCharge()==0:
            kind='Car' if atom.GetIsAromatic() else {Chem.HybridizationType.SP:'C1',Chem.HybridizationType.SP2:'C2',Chem.HybridizationType.SP3:'C3'}.get(atom.GetHybridization())
            if kind is None:continue
            append(atom,kind)
    for oxygen in molecule.GetAtoms():
        if oxygen.GetAtomicNum()==8 and oxygen.GetFormalCharge()==0 and oxygen.GetDegree()==1 and oxygen.GetBonds()[0].GetBondTypeAsDouble()==2 and oxygen.GetNeighbors()[0].GetAtomicNum()!=15:
            append(oxygen,'O2')
            continue
        if oxygen.GetAtomicNum()!=8 or oxygen.GetFormalCharge()!=0 or oxygen.GetDegree()!=1:continue
        nitrogen=oxygen.GetNeighbors()[0]
        if nitrogen.GetAtomicNum()!=7 or molecule.GetBondBetweenAtoms(oxygen.GetIdx(),nitrogen.GetIdx()).GetBondTypeAsDouble()!=1:continue
        if not any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(nitrogen).GetAtomicNum()==6 for b in nitrogen.GetBonds()):continue
        append(oxygen,'O3')
    for nitrogen in molecule.GetAtoms():
        if nitrogen.GetAtomicNum()!=7 or nitrogen.GetFormalCharge()!=1 or nitrogen.GetDegree()!=3:continue
        oxygen=[a for a in nitrogen.GetNeighbors() if a.GetAtomicNum()==8 and a.GetFormalCharge()==-1 and a.GetDegree()==1]
        imines=[b for b in nitrogen.GetBonds() if b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(nitrogen).GetAtomicNum()==6]
        if len(oxygen)!=1 or len(imines)!=1:continue
        for atom,kind in ((nitrogen,'N2+'),(oxygen[0],'O3-')):
            append(atom,kind)
    return result
