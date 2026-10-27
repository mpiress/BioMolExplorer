"""In-process PDB/mmCIF selection and CCD-backed crystallographic chemistry."""
from pathlib import Path
import os
import re

import numpy as np
from rdkit import Chem

from .calculation_cache import CalculationCache
from .storage import write_text
from .md_analysis import PROTEIN

WATER={'HOH','WAT','DOD','SOL'}
CAPS={'ACE','NME','NHE'}
METALS={'Li','Na','K','Rb','Cs','Mg','Ca','Sr','Ba','Mn','Fe','Co','Ni','Cu','Zn','Cd','Hg'}


def protein_residue(residue):
    return residue.name in PROTEIN|CAPS or gemmi_module().find_tabulated_residue(residue.name).is_amino_acid()


def gemmi_module():
    try:import gemmi
    except ImportError:raise ValueError('Gemmi is required for native structural preparation; update requirements.yml.') from None
    return gemmi


def source_structure(folder,code):
    for suffix in ('.pdb','.cif','.mmcif'):
        path=Path(folder)/(code+suffix)
        if path.is_file():return path
    raise ValueError('Structure not found (PDB/mmCIF): '+code)


def atom_identity(chain,residue,atom,model=1):
    return dict(model=model,chain=chain.name,residue_name=residue.name,
                residue_number=residue.seqid.num,insertion_code=residue.seqid.icode.strip(),
                atom_name=atom.name,altloc=atom.altloc.replace('\x00',''),serial=atom.serial,
                element=atom.element.name,coordinates=[atom.pos.x,atom.pos.y,atom.pos.z])


def read_selected(path,settings=None):
    """One explicit model; choose a coherent residue alternate, never a mosaic."""
    settings=settings or {};gemmi=gemmi_module();structure=gemmi.read_structure(str(path))
    model=int(settings.get('model',1))
    if not 1<=model<=len(structure):raise ValueError('Requested structural model does not exist.')
    selected=structure.clone()
    for index in reversed(range(len(selected))):
        if index!=model-1:del selected[index]
    records=[]
    for chain in selected[0]:
        for residue in chain:
            occupancy={}
            for atom in residue:
                if atom.altloc!='\x00':occupancy[atom.altloc]=occupancy.get(atom.altloc,0)+atom.occ
            choice=settings.get('altloc','occupancy')
            choice=sorted(occupancy,key=lambda key:(-occupancy[key],key))[0] if occupancy and choice=='occupancy' else choice
            if occupancy and choice not in occupancy:raise ValueError('Requested alternate location absent in '+str(residue.seqid))
            for index in reversed(range(len(residue))):
                atom=residue[index]
                if atom.altloc not in ('\x00',choice):del residue[index]
            names=[atom.name for atom in residue]
            if len(names)!=len(set(names)):raise ValueError('Ambiguous duplicate atoms after alternate-location selection.')
            for atom in residue:
                if not np.isfinite([atom.pos.x,atom.pos.y,atom.pos.z]).all():raise ValueError('Non-finite structural coordinates.')
                records.append(atom_identity(chain,residue,atom,model))
    return selected,records


def subset(structure,predicate):
    result=structure.clone()
    for chain_index in reversed(range(len(result[0]))):
        chain=result[0][chain_index]
        for residue_index in reversed(range(len(chain))):
            if not predicate(chain,chain[residue_index]):del chain[residue_index]
        if not len(chain):del result[0][chain_index]
    if not any(len(residue) for chain in result[0] for residue in chain):raise ValueError('Structural selection contains no atoms.')
    return result


def write_structure(structure,path,remove_hydrogens=False):
    """PDB contracts must never silently truncate mmCIF identifiers."""
    structure=structure.clone()
    if remove_hydrogens:structure.remove_hydrogens()
    path=Path(path)
    if path.suffix in ('.cif','.mmcif'):
        write_text(structure.make_mmcif_document().as_string(),path)
    else:
        for chain in structure[0]:
            if len(chain.name)>1:raise ValueError('Chain ID exceeds the PDB docking contract; use an explicit chain mapping, not silent truncation.')
            for residue in chain:
                if len(residue.name)>3 or not -999<=residue.seqid.num<=9999:raise ValueError('Residue identity cannot be represented in the PDB docking contract.')
        write_text(structure.make_pdb_string(),path)
    return path


def component_path(code,folder,config):
    """Prefer explicit/local topology; otherwise retrieve and retain wwPDB CCD."""
    if not re.fullmatch('[A-Za-z0-9]{1,5}',code):raise ValueError('Invalid chemical component identifier.')
    explicit=config.get('chemistry_templates',{}).get(code)
    if explicit:
        path=Path(explicit)
        if not path.is_absolute():path=Path(folder)/path
        if not path.is_file():raise ValueError('Chemical topology is missing: '+str(path))
        return path
    roots=[Path(folder)/'chemistry',Path(folder)]
    if os.environ.get('BIOMOL_CCD_PATH'):roots.insert(0,Path(os.environ['BIOMOL_CCD_PATH']))
    for root in roots:
        for suffix in ('.cif','.sdf','.mol2'):
            path=root/(code+suffix)
            if path.is_file():return path
    cache=CalculationCache('chemical-components',folder);path=cache.root/'ccd'/(code+'.cif')
    if not path.is_file():
        import requests
        with cache.locked('ccd-'+code):
            if not path.is_file():
                try:
                    response=requests.get('https://files.rcsb.org/ligands/download/'+code+'.cif',timeout=30)
                    response.raise_for_status()
                    block=gemmi_module().cif.read_string(response.text).sole_block()
                    if not len(block.find_values('_chem_comp_atom.atom_id')):raise ValueError('CCD contains no chemical graph.')
                    write_text(response.text,path)
                except (requests.RequestException,ValueError) as error:
                    raise ValueError('No reliable chemical definition for '+code+'. Provide chemistry/'+code+'.cif, SDF or MOL2 with atom names; topology will not be inferred from distances.') from error
    return path


def chemical_template(path,code):
    path=Path(path)
    if path.suffix.lower()!='.cif':
        from .molecular_core import MolecularPreparationCore
        molecule,_=MolecularPreparationCore(path.parent).load(path)
        return molecule,{}
    blocks=gemmi_module().cif.read(str(path))
    block=next((b for b in blocks if b.find_value('_chem_comp.id')==code),blocks[0])
    columns=['atom_id','type_symbol','charge','pdbx_stereo_config']
    table=block.find('_chem_comp_atom.',columns);rw=Chem.RWMol();names={};stereo={}
    for row in table:
        name=row[0]
        if row[1].upper() in ('H','D'):continue
        atom=Chem.Atom(row[1].title());atom.SetFormalCharge(int(row[2]) if row[2] not in ('.','?') else 0)
        atom.SetProp('_TriposAtomName',name);names[name]=rw.AddAtom(atom)
        if row[3] in ('R','S'):stereo[name]=row[3]
    orders={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for row in block.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order']):
        if row[0] not in names or row[1] not in names:continue
        if row[2] not in orders:raise ValueError('Unsupported CCD bond order: '+row[2])
        rw.AddBond(names[row[0]],names[row[1]],orders[row[2]])
    molecule=rw.GetMol()
    if not molecule.GetNumAtoms():raise ValueError('Chemical component has no heavy atoms.')
    Chem.SanitizeMol(molecule)
    declared=block.find_value('_chem_comp.pdbx_formal_charge')
    if declared and declared not in ('.','?') and Chem.GetFormalCharge(molecule)!=int(declared):
        raise ValueError('CCD molecular and per-atom formal charges disagree.')
    return molecule,stereo


def crystal_molecule(structure,record,config,folder):
    """Join trusted bond orders to coordinates by atom names, not distances."""
    chain_id=config.get('ligand_chain',record[3]);icode=config.get('ligand_insertion_code','')
    residues=[(chain,residue) for chain in structure[0] for residue in chain
              if chain.name==chain_id and residue.name==record[1] and residue.seqid.num==int(record[2])
              and residue.seqid.icode.strip()==icode]
    if len(residues)!=1:raise ValueError('Crystallographic ligand selection is missing or ambiguous; specify chain and insertion code.')
    chain,residue=residues[0]
    ligand_atoms=[atom for atom in residue if not atom.element.is_hydrogen]
    metal_atoms=[]
    for other_chain in structure[0]:
        for other_residue in other_chain:
            if other_chain.name==chain_id and other_residue.seqid==residue.seqid:continue
            for atom in other_residue:
                if atom.element.is_hydrogen:continue
                if atom.element.name in METALS:metal_atoms.append(atom)
    from scipy.spatial import cKDTree
    ligand_xyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in ligand_atoms])
    if metal_atoms and cKDTree([[a.pos.x,a.pos.y,a.pos.z] for a in metal_atoms]).query(ligand_xyz)[0].min()<3.:
        raise ValueError('Metal within ligand coordination distance requires a dedicated parameterized protocol.')
    from .complex_classification import classify_complex,COVALENT,AMBIGUOUS,CovalentComplexError
    classification=classify_complex(structure,record,config)
    if classification['interaction_type']==COVALENT:raise CovalentComplexError(classification)
    if classification['interaction_type']==AMBIGUOUS:
        raise ValueError('Ambiguous receptor-ligand interaction requires explicit chemistry: '+str(classification['ambiguities'][0]))
    path=component_path(residue.name,folder,config);template,stereo=chemical_template(path,residue.name)
    template=Chem.RemoveHs(template);atoms={a.name:a for a in residue if not a.element.is_hydrogen}
    names=[a.GetProp('_TriposAtomName') if a.HasProp('_TriposAtomName') else '' for a in template.GetAtoms()]
    if not all(names) or len(set(names))!=len(names) or set(names)!=set(atoms):
        raise ValueError('Chemical template heavy-atom names/count do not match the crystal ligand.')
    molecule=Chem.Mol(template);molecule.RemoveAllConformers();conformer=Chem.Conformer(len(names));mapping=[]
    for index,atom in enumerate(molecule.GetAtoms()):
        source=atoms[names[index]]
        if source.element.atomic_number!=atom.GetAtomicNum():raise ValueError('Crystal/template element mismatch.')
        conformer.SetAtomPosition(index,(source.pos.x,source.pos.y,source.pos.z))
        info=Chem.AtomPDBResidueInfo();info.SetName(source.name.ljust(4));info.SetResidueName(residue.name)
        info.SetResidueNumber(residue.seqid.num);info.SetChainId(chain.name);info.SetInsertionCode(icode)
        info.SetIsHeteroAtom(True);atom.SetMonomerInfo(info)
        mapping.append(dict(prepared_index=index,original=atom_identity(chain,residue,source,config.get('structure',{}).get('model',1))))
    conformer.Set3D(True);molecule.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(molecule);Chem.AssignStereochemistry(molecule,cleanIt=True,force=True)
    for atom,name in zip(molecule.GetAtoms(),names):
        expected=stereo.get(name)
        if not expected and template.GetAtomWithIdx(atom.GetIdx()).HasProp('_CIPCode'):
            expected=template.GetAtomWithIdx(atom.GetIdx()).GetProp('_CIPCode')
        if expected and (not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=expected):
            raise ValueError('Crystallographic stereochemistry conflicts with the chemical definition at '+name)
    if any(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()):raise ValueError('Radical crystal ligand requires a dedicated spin protocol.')
    return molecule,mapping,path
