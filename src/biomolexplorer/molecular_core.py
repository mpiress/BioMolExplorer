"""Engine-independent small-molecule preparation; exporters never assign charges."""
from dataclasses import dataclass, field
from contextlib import contextmanager
import hashlib
import json
import math
import os
import re
from pathlib import Path
import tempfile

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdPartialCharges
from Bio.Data.PDBData import protein_letters_3to1

from .calculation_cache import CalculationCache, identity
from .storage import write_json

CHARGE_ALIASES={'gas':'gasteiger','am1':'am1bcc','gasteiger':'gasteiger','am1bcc':'am1bcc'}


@contextmanager
def preparation_workspace(root):
    if os.environ.get('BIOMOL_PREPARATION_DEBUG','0')=='1':
        debug=Path(root)/'debug';debug.mkdir(parents=True,exist_ok=True)
        yield Path(tempfile.mkdtemp(prefix='charge-',dir=debug))
    else:
        with tempfile.TemporaryDirectory(prefix='bme-charge-') as folder:yield Path(folder)


def atom_name(atom,index):
    symbol=atom.GetSymbol();width=4-len(symbol)
    value=np.base_repr(index+1,base=36).zfill(width)
    if len(value)>width:raise ValueError('Ligand atom count exceeds the stable Amber atom-name mapping.')
    return symbol+value


def charge_method(value):
    try:return CHARGE_ALIASES[value]
    except (KeyError,TypeError):raise ValueError('Unsupported molecular charge method: '+str(value)) from None


def charge_strategy(method,engines,policy='fast',minimize=False):
    requested=charge_method(method);effective=requested;scoring=None
    if engines=='vina':
        from .paths import resolve_path
        text=resolve_path('src/scripts/vina/config.template').read_text()
        match=re.search(r'^\s*scoring\s*=\s*(\S+)',text,re.M);scoring=match[1] if match else 'vina'
        if not minimize and policy=='fast' and scoring in ('vina','vinardo'):effective='gasteiger'
    engine='antechamber' if minimize or engines!='vina' or effective=='am1bcc' else 'rdkit'
    return dict(requested=requested,effective=effective,engine=engine,scoring=scoring)


def validate_charges(molecule,values,tolerance=1e-3):
    values=np.asarray(values,float)
    if values.shape!=(molecule.GetNumAtoms(),) or not np.isfinite(values).all():
        raise ValueError('Partial charges are missing or non-finite.')
    formal=Chem.GetFormalCharge(molecule)
    if abs(float(values.sum())-formal)>tolerance+1e-10:
        raise ValueError(f'Partial charge sum {values.sum():.8f} does not match formal charge {formal}.')
    return values


@dataclass
class PreparedMolecule:
    molecule: object
    partial_charges: np.ndarray
    metadata: dict=field(default_factory=dict)

    def __post_init__(self):
        self.partial_charges=validate_charges(self.molecule,self.partial_charges)
        if self.molecule.GetNumConformers()!=1:raise ValueError('A single molecular conformer is required.')
        if not np.isfinite(self.coordinates).all():raise ValueError('Non-finite molecular coordinates.')
        for index,atom in enumerate(self.molecule.GetAtoms()):
            atom.SetDoubleProp('BioMolPartialCharge',float(self.partial_charges[index]))
            if not atom.HasProp('_TriposAtomName'):atom.SetProp('_TriposAtomName',atom_name(atom,index))
        residue_ids=[]
        for atom in self.molecule.GetAtoms():
            info=atom.GetPDBResidueInfo()
            residue_ids.append(dict(chain=info.GetChainId(),residue_name=info.GetResidueName().strip(),
                residue_number=info.GetResidueNumber(),insertion_code=info.GetInsertionCode().strip()) if info else None)
        self.metadata.update(formal_charge=Chem.GetFormalCharge(self.molecule),
            partial_charge_sum=float(self.partial_charges.sum()),atom_count=self.molecule.GetNumAtoms(),
            heavy_atom_count=self.molecule.GetNumHeavyAtoms(),
            canonical_structure=Chem.MolToSmiles(self.molecule,isomericSmiles=True),
            partial_charges=self.partial_charges.tolist(),
            formal_charges=[a.GetFormalCharge() for a in self.molecule.GetAtoms()],
            elements=[a.GetSymbol() for a in self.molecule.GetAtoms()],residue_ids=residue_ids,
            bond_orders=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in self.molecule.GetBonds()],
            atom_names=self.atom_names,rdkit_version=rdBase.rdkitVersion,
            representation='RDKit graph + conformer + named atom/residue metadata; PDBQT/MOL2 are exports only')
        self.metadata.setdefault('atom_mapping',list(range(self.molecule.GetNumAtoms())))

    @property
    def coordinates(self):return self.molecule.GetConformer().GetPositions()

    @property
    def atom_names(self):return [a.GetProp('_TriposAtomName') for a in self.molecule.GetAtoms()]


class GasteigerEngine:
    name='rdkit'
    def identity(self):return dict(engine=self.name,version=rdBase.rdkitVersion)
    def calculate(self,molecule,work):
        rdPartialCharges.ComputeGasteigerCharges(molecule,throwOnParamFailure=True)
        return validate_charges(molecule,[a.GetDoubleProp('_GasteigerCharge') for a in molecule.GetAtoms()])


class MolecularPreparationCore:
    """One normalized/protonated graph and charge array for every backend.

    No tautomer enumeration, fragment stripping or hidden minimization occurs.
    Coordinate-bearing inputs retain heavy-atom coordinates. pH handling uses
    the existing Open Babel dependency in process, not its CLI.
    """
    def __init__(self,cache_directory,ambertools_path='',ambertools_env=''):
        self.cache=CalculationCache('molecular-charges',cache_directory)
        self.ambertools_path=ambertools_path
        self.ambertools_env=ambertools_env

    def load(self,input_data,*,index=1):
        if type(index) is not int or index<1:raise ValueError('conformer_index deve ser um inteiro positivo.')
        if isinstance(input_data,Chem.Mol):return Chem.Mol(input_data),'rdkit'
        path=Path(input_data)
        if not path.is_file():
            raise ValueError('Molecular source file does not exist: '+str(path))
        suffix=path.suffix.lower()
        if suffix=='.sdf':
            supplier=Chem.SDMolSupplier(str(path),removeHs=False)
            molecule=supplier[index-1] if 1<=index<=len(supplier) else None
            if molecule is not None:
                for atom in molecule.GetAtoms():
                    if atom.HasProp('molFileAlias'):atom.SetProp('_TriposAtomName',atom.GetProp('molFileAlias'))
        elif suffix=='.pdbqt':
            # Recover only the explicit SMILES/topological index mapping that
            # our Meeko exporter embeds. Never infer bond orders from a pose.
            text=path.read_text()
            if 'REMARK SMILES ' not in text or 'REMARK SMILES IDX ' not in text:
                raise ValueError('PDBQT has no trusted chemical graph/atom mapping; provide the corresponding SDF/MOL2. PDBQT-to-MOL2 inference is disabled.')
            from meeko import PDBQTMolecule,RDKitMolCreate
            poses=PDBQTMolecule.from_file(str(path),poses_to_read=index,skip_typing=True)
            molecules=RDKitMolCreate.from_pdbqt_mol(poses)
            if len(molecules)!=1 or molecules[0] is None:raise ValueError('Ambiguous chemical graph in pose metadata.')
            molecule=Chem.Mol(molecules[0])
            if molecule.GetNumConformers()<index:raise ValueError('Selected docking pose does not exist.')
            conformer=Chem.Conformer(molecule.GetConformer(index-1));molecule.RemoveAllConformers();molecule.AddConformer(conformer)
            return molecule,'pdbqt/explicit-Meeko-graph'
        elif suffix=='.mol2':
            # SYBYL charge/aromatic extensions are not all supported by RDKit's
            # MOL2 reader (notably N.4 and O.co2). Parse with Open Babel, then
            # transfer the graph in memory through its richer SDF representation.
            from openbabel import openbabel as ob
            entries=path.read_text().split('@<TRIPOS>MOLECULE')
            selected='@<TRIPOS>MOLECULE'+entries[index] if index<len(entries) else ''
            converter=ob.OBConversion();converter.SetInAndOutFormats('mol2','sdf');native=ob.OBMol()
            molecule=Chem.MolFromMolBlock(converter.WriteString(native),removeHs=False) if selected and converter.ReadString(native,selected) else None
        else:
            raise ValueError('Native preparation requires SDF, MOL2 or RDKit input; PDBQT is an export, not a molecular graph.')
        if molecule is None:raise ValueError('Cannot read a complete molecular graph: '+path.name)
        if suffix=='.mol2':
            block=selected.split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
            rows=[line.split() for line in block.splitlines() if line.strip()]
            if len(rows)!=molecule.GetNumAtoms():raise ValueError('MOL2 export changed atom count.')
            for atom,row in zip(molecule.GetAtoms(),rows):
                info=Chem.AtomPDBResidueInfo();info.SetName(row[1]);info.SetResidueName(row[7]);info.SetResidueNumber(int(row[6]));info.SetIsHeteroAtom(True)
                atom.SetMonomerInfo(info);atom.SetProp('_TriposAtomName',row[1])
        return molecule,suffix[1:]

    def prepare(self,input_data,*,method='gasteiger',ph=7.4,add_hydrogens=True,remove_hydrogens=True,index=1,minimize=False,charge_engine='rdkit'):
        method=charge_method(method)
        if not math.isfinite(float(ph)) or not 0<=float(ph)<=14:raise ValueError('Preparation pH must be between 0 and 14.')
        molecule,source_format=self.load(input_data,index=index)
        Chem.SanitizeMol(molecule)
        if any(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()):raise ValueError('Radical ligands require explicit spin and dedicated parameterization.')
        if abs(Chem.GetFormalCharge(molecule))>3:raise ValueError('Highly charged ligands require explicit validation outside the automatic preparation envelope.')
        residues={(a.GetPDBResidueInfo().GetChainId(),a.GetPDBResidueInfo().GetResidueNumber()) for a in molecule.GetAtoms()
            if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip().upper() in protein_letters_3to1}
        if len(residues)>1:raise ValueError('Protein residues require receptor force-field templates, not whole-protein ligand charges.')
        if molecule.GetNumConformers()>1:raise ValueError('A single molecular conformer is required.')
        if molecule.GetNumConformers() and not np.isfinite(molecule.GetConformer().GetPositions()).all():raise ValueError('Non-finite molecular coordinates.')
        if molecule.GetNumConformers() and not molecule.GetConformer().Is3D():
            raise ValueError('Supply a 3D conformer or SMILES for conformer generation; 2D coordinates are not a bound pose.')
        if len(Chem.GetMolFrags(molecule))!=1:raise ValueError('Prepare molecular fragments independently; automatic salt stripping is disabled.')
        if any(a.GetAtomicNum() not in {1,6,7,8,9,15,16,17,35,53} for a in molecule.GetAtoms()):
            raise ValueError('Metal or unsupported element requires dedicated parameterization, not generic ligand charges.')
        before=Chem.MolToSmiles(molecule,isomericSmiles=True)
        generated=molecule.GetNumConformers()==0
        if molecule.GetNumConformers()==0:
            # Reuse the established seed42/UFF initial-conformer protocol.
            # This is not a replacement for the separately selected AMBER
            # minimization, which still follows its validated legacy path.
            from .visualizations import molecule_sdf
            molecule=Chem.MolFromMolBlock(molecule_sdf(Chem.MolToSmiles(molecule,canonical=False,isomericSmiles=True)),removeHs=False)
            if molecule is None:raise ValueError('Cannot generate a molecular conformer.')
        original_atoms=[a for a in molecule.GetAtoms() if a.GetAtomicNum()!=1]
        original_heavy=np.array([molecule.GetConformer().GetAtomPosition(a.GetIdx()) for a in original_atoms])
        original_info=[]
        for atom in original_atoms:
            info=atom.GetPDBResidueInfo()
            original_info.append(dict(input_index=atom.GetIdx(),element=atom.GetAtomicNum(),
                name=info.GetName().strip() if info else atom.GetProp('_TriposAtomName') if atom.HasProp('_TriposAtomName') else None,
                residue_name=info.GetResidueName().strip() if info else None,
                residue_number=info.GetResidueNumber() if info else None,chain=info.GetChainId() if info else None))
        if remove_hydrogens:molecule=Chem.RemoveHs(molecule)
        if add_hydrogens:
            from openbabel import openbabel as ob
            converter=ob.OBConversion();converter.SetInAndOutFormats('sdf','sdf')
            native=ob.OBMol()
            if not converter.ReadString(native,Chem.MolToMolBlock(molecule)):raise ValueError('Cannot normalize ligand protonation.')
            native.CorrectForPH(float(ph));native.AddHydrogens()
            molecule=Chem.MolFromMolBlock(converter.WriteString(native),removeHs=False)
            if molecule is None:raise ValueError('Protonation produced an invalid molecular graph.')
        if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in molecule.GetAtoms()):
            # RDKit bracket-H counts are not graph atoms; charges must include
            # every hydrogen rather than silently losing implicit H charge.
            molecule=Chem.AddHs(molecule,addCoords=True) if add_hydrogens else molecule
            if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in molecule.GetAtoms()):
                raise ValueError('Explicit hydrogen atoms are required before assigning partial charges.')
        heavy=np.array([molecule.GetConformer().GetAtomPosition(a.GetIdx()) for a in molecule.GetAtoms() if a.GetAtomicNum()!=1])
        if heavy.shape!=original_heavy.shape or not np.allclose(heavy,original_heavy,atol=2e-4,rtol=0):
            raise ValueError('Protonation changed heavy-atom mapping or coordinates.')
        molecule.SetProp('_Name','LIG')
        used=set();heavy_index=0;mapping=[]
        for i,a in enumerate(molecule.GetAtoms()):
            original=original_info[heavy_index] if a.GetAtomicNum()!=1 else None
            if original:heavy_index+=1
            name=original['name'] if original and original['name'] and len(original['name'])<=4 and original['name'] not in used else atom_name(a,i)
            if name in used:name=atom_name(a,i)
            if name in used:raise ValueError('Antechamber changed or duplicated the ligand atom mapping.')
            a.SetProp('_TriposAtomName',name);used.add(name);mapping.append(original['input_index'] if original else None)
        if charge_engine not in ('rdkit','antechamber'):raise ValueError('Unsupported ligand charge engine: '+str(charge_engine))
        if method=='gasteiger' and charge_engine=='rdkit':engine=GasteigerEngine()
        else:
            from .molecular_adapters import AM1BCCEngine,AntechamberGasteigerEngine
            engine=(AM1BCCEngine if method=='am1bcc' else AntechamberGasteigerEngine)(self.ambertools_path,self.ambertools_env)
        recipe=dict(method=method,engine=engine.identity(),rdkit=rdBase.rdkitVersion,
            adapter={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('molecular_core.py','molecular_adapters.py')},
            structure=Chem.MolToSmiles(molecule,isomericSmiles=True),
            atoms=[(a.GetAtomicNum(),a.GetFormalCharge(),int(a.GetChiralTag())) for a in molecule.GetAtoms()],
            bonds=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType()),str(b.GetStereo())) for b in molecule.GetBonds()],
            coordinates=np.round(molecule.GetConformer().GetPositions(),6).tolist() if method=='am1bcc' else None)
        key=identity(recipe)
        with preparation_workspace(self.cache.root) as work:
            artifact=work/'charges.json'
            with self.cache.locked(key):
                cached=self.cache.restore(key,{'charges':artifact})
                reused=cached is not None
                if reused:
                    try:
                        charge_record=json.loads(artifact.read_text())
                        values=validate_charges(molecule,charge_record['charges'])
                    except (ValueError,KeyError,TypeError):reused=False
                if not reused:
                    values=engine.calculate(molecule,work)
                    charge_record=dict(charges=values.tolist(),execution=getattr(engine,'execution',engine.name))
                    write_json(charge_record,artifact)
                    self.cache.remember(key,{'charges':artifact},dict(recipe=recipe))
        metadata=dict(ph=float(ph),protonation_engine='Open Babel' if add_hydrogens else 'input state preserved',
            conformer_source='existing seed42/UFF initial-conformer protocol' if generated else 'input 3D coordinates',
            original_heavy_atoms=original_info,
            atom_mapping=mapping,
            protonation_before=before,protonation_after=Chem.MolToSmiles(molecule,isomericSmiles=True),
            charge_method=method,charge_engine=engine.identity(),charge_execution=charge_record.get('execution',engine.name),charge_cache_hit=reused,
            source_format=source_format,minimization=False,backend='python',charge_cache_key=key)
        from .diagnostics import event,get_logger
        event(get_logger('scientific_tools'),'molecular.prepared','Molecular preparation complete',
            tool=engine.name,configuration=key,reused=reused)
        prepared=PreparedMolecule(molecule,values,metadata)
        if minimize:
            from .chimera_minimization import minimize_ligand
            with preparation_workspace(self.cache.root) as work:
                prepared=minimize_ligand(prepared,work,cache_directory=self.cache.root.parent)
        return prepared
