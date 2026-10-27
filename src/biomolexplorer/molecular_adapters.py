"""Backend exports and Amber parameterization of an already charged graph."""
from pathlib import Path

import numpy as np
from rdkit import Chem,rdBase

from .molecular_core import PreparedMolecule, validate_charges,preparation_workspace
from .processes import run_command, ScientificToolError
from .stage_cache import file_digest
from .storage import write_text
from .preparation_tools import amber_tool,run_amber_tool


class AutoDockAdapter:
    def export_string(self,prepared):
        try:from meeko import MoleculePreparation,PDBQTWriterLegacy
        except ImportError as error:
            raise ScientificToolError('Meeko/Gemmi are required for native AutoDock export.','TOOL_NOT_FOUND','Update the scientific environment using requirements.yml.') from error
        # Meeko merges nonpolar H charges into their parent atoms; it must not
        # silently replace AM1-BCC/Gasteiger with a new charge calculation.
        molecule=Chem.Mol(prepared.molecule)
        for index,atom in enumerate(molecule.GetAtoms()):
            info=atom.GetPDBResidueInfo()
            if info is None:
                info=Chem.AtomPDBResidueInfo();info.SetResidueName('LIG');info.SetResidueNumber(1);info.SetIsHeteroAtom(True)
                atom.SetMonomerInfo(info)
                info=atom.GetPDBResidueInfo()
            info.SetName(prepared.atom_names[index].ljust(4))
        setups=MoleculePreparation(charge_model='read',charge_atom_prop='BioMolPartialCharge',add_index_map=True).prepare(molecule)
        if len(setups)!=1:raise ValueError('AutoDock export requires exactly one prepared setup.')
        text,ok,error=PDBQTWriterLegacy.write_string(setups[0],add_index_map=True)
        if not ok:raise ValueError('Meeko PDBQT export failed: '+str(error))
        values=[float(line.split()[-2]) for line in text.splitlines() if line.startswith(('ATOM  ','HETATM'))]
        if not values or not np.isfinite(values).all() or abs(sum(values)-Chem.GetFormalCharge(prepared.molecule))>max(.002,.00051*len(values)):
            raise ValueError('PDBQT export lost or changed molecular charges.')
        return text

    def export(self,prepared,path):
        path=Path(path);write_text(self.export_string(prepared),path)
        from .ligand_cache import validate_prepared_ligand
        return validate_prepared_ligand(path,'pdbqt')

    def set_ligand(self,vina,prepared):
        """Direct Python API path; exporting an audit file remains optional."""
        text=self.export_string(prepared);vina.set_ligand_from_string(text)
        return text


class Dock6Adapter:
    def export_string(self,prepared):
        from openbabel import openbabel as ob
        converter=ob.OBConversion();converter.SetInAndOutFormats('sdf','mol2')
        molecule=ob.OBMol()
        if not converter.ReadString(molecule,Chem.MolToMolBlock(prepared.molecule)):
            raise ValueError('Cannot export the prepared molecular graph to MOL2.')
        if molecule.NumAtoms()!=prepared.molecule.GetNumAtoms():raise ValueError('MOL2 export changed atom count.')
        for index,atom in enumerate(ob.OBMolAtomIter(molecule)):
            if atom.GetAtomicNum()!=prepared.molecule.GetAtomWithIdx(index).GetAtomicNum():raise ValueError('MOL2 export changed atom mapping.')
            atom.SetPartialCharge(float(prepared.partial_charges[index]))
        molecule.SetPartialChargesPerceived()
        text=converter.WriteString(molecule)
        lines=text.splitlines();start=lines.index('@<TRIPOS>ATOM')+1
        # Open Babel otherwise labels even supplied AM1-BCC values GASTEIGER.
        header=lines.index('@<TRIPOS>MOLECULE');lines[header+4]='USER_CHARGES'
        index=0
        for i in range(start,len(lines)):
            if lines[i].startswith('@<TRIPOS>'):break
            if not lines[i].strip():continue
            cells=lines[i].split();cells[1]=prepared.atom_names[index]
            cells[8]=f'{prepared.partial_charges[index]:.8f}'
            lines[i]=' '.join(cells);index+=1
        if index!=prepared.molecule.GetNumAtoms():raise ValueError('Incomplete MOL2 atom export.')
        # Preserve the graph's bond orders. Open Babel may otherwise label
        # non-ring carboxylate bonds "ar", changing subsequent BCC typing.
        bonds={frozenset((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1)):b for b in prepared.molecule.GetBonds()}
        start=lines.index('@<TRIPOS>BOND')+1
        for i in range(start,len(lines)):
            if lines[i].startswith('@<TRIPOS>'):break
            if not lines[i].strip():continue
            cells=lines[i].split();bond=bonds[frozenset((int(cells[1]),int(cells[2])))]
            cells[3]='ar' if bond.GetIsAromatic() else 'am' if cells[3]=='am' and bond.GetBondTypeAsDouble()==1 else str(int(bond.GetBondTypeAsDouble()))
            lines[i]=' '.join(cells)
        return '\n'.join(lines)+'\n'

    def export(self,prepared,path):
        path=Path(path);write_text(self.export_string(prepared),path)
        from .ligand_cache import validate_prepared_ligand
        return validate_prepared_ligand(path,'mol2')


class AM1BCCEngine:
    name='antechamber'
    method='am1bcc'
    def __init__(self,prefix='',environment=''):
        self.prefix=prefix;self.environment=environment
        self.tool=amber_tool('antechamber',prefix,environment=environment)
    def identity(self):
        if hasattr(self,'_identity'):return self._identity
        root=self.tool.parent.parent
        files={name:file_digest(root/'bin'/name) for name in ('antechamber','sqm','am1bcc','atomtype','bondtype') if (root/'bin'/name).is_file()}
        data=root/'dat/antechamber'
        files.update({str(p.relative_to(root)):file_digest(p) for p in sorted(data.rglob('*')) if p.is_file()})
        self._identity=dict(engine=self.name,method=self.method,files=files,protocol=dict(j=5,acdoctor=False,qm_theory='AM1' if self.method=='am1bcc' else None))
        return self._identity

    def calculate(self,molecule,work):
        neutral=PreparedMolecule(Chem.Mol(molecule),np.full(molecule.GetNumAtoms(),Chem.GetFormalCharge(molecule)/molecule.GetNumAtoms()))
        Dock6Adapter().export(neutral,work/'input.mol2')
        # Amber18's reader mistakes SYBYL UNITY attributes for atom records.
        # The formal molecular charge is already supplied explicitly via -nc;
        # keep this compatibility normalization local to Antechamber input.
        lines=[];skip=False
        for line in (work/'input.mol2').read_text().splitlines():
            if line.startswith('@<TRIPOS>'):skip=line.startswith('@<TRIPOS>UNITY_')
            if not skip:lines.append(line)
        write_text('\n'.join(lines)+'\n',work/'input.mol2')
        command=['-i','input.mol2','-fi','mol2','-o','charged.mol2','-fo','mol2',
            '-c','bcc' if self.method=='am1bcc' else 'gas','-nc',str(Chem.GetFormalCharge(molecule)),
            '-at','gaff','-an','n','-du','n','-j','5','-dr','n','-s','2','-pf','y']
        if self.method=='am1bcc':command+=['-ek',"qm_theory='AM1',"]
        run_amber_tool('antechamber',command,cwd=work,prefix=self.prefix,environment=self.environment)
        return read_mapped_charges(work/'charged.mol2',molecule)


class AntechamberGasteigerEngine(AM1BCCEngine):
    """The legacy Chimera charge family, without starting Chimera."""
    method='gasteiger'

    def identity(self):
        identity=dict(super().identity())
        identity['formal_charge_fallback']='RDKit Gasteiger when Antechamber gas returns a different net charge'
        identity['rdkit']=rdBase.rdkitVersion
        return identity

    def calculate(self,molecule,work):
        self.execution='Antechamber Gasteiger'
        try:return super().calculate(molecule,work)
        except ValueError as error:
            # Antechamber's ``gas`` typer can report total formal charge zero
            # for charged imine/amidine patterns even when ``-nc`` is supplied.
            # Never renormalize a whole-electron error.  The in-process RDKit
            # Gasteiger implementation preserves the requested graph charge and
            # is the explicit, provenance-recorded fallback for this method.
            if Chem.GetFormalCharge(molecule)==0 or 'does not match formal charge' not in str(error):raise
            from .molecular_core import GasteigerEngine
            self.execution='RDKit Gasteiger (formal-charge fallback)'
            return GasteigerEngine().calculate(molecule,work)


def charge_existing_mol2(path, method='am1bcc', prefix=''):
    """Assign charges only: preserve the input geometry, H, types and records.

    Unlike molecular preparation, this operation must not reprotonate or
    minimize an already converted pose. Publish only after mapped validation.
    """
    import re
    from .molecular_core import MolecularPreparationCore, charge_method, atom_name
    from .calculation_cache import CalculationCache, identity
    from .ligand_cache import validate_prepared_ligand

    path=Path(path)
    validate_prepared_ligand(path,'mol2')
    text=path.read_text()
    if text.count('@<TRIPOS>MOLECULE')!=1:
        raise ValueError('Charge-only preparation requires one MOL2 molecule.')
    molecule,_=MolecularPreparationCore(path.parent).load(path)
    names=[a.GetProp('_TriposAtomName') for a in molecule.GetAtoms()]
    if len(set(names))!=len(names):
        # Duplicate source labels are legal in some imported structures.
        # Use unambiguous names only inside Antechamber, then transfer the
        # validated vector back by original atom index without renaming output.
        for index,atom in enumerate(molecule.GetAtoms()):
            atom.SetProp('_TriposAtomName',atom_name(atom,index))
    lines=text.splitlines(keepends=True)
    marker=next(i for i,line in enumerate(lines) if line.strip()=='@<TRIPOS>MOLECULE')
    normalized=list(lines);normalized[marker+4]='USER_CHARGES\n'
    in_atoms=False;records=[]
    for number,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            in_atoms=line.strip()=='@<TRIPOS>ATOM'
        elif in_atoms and line.strip():
            fields=list(re.finditer(r'\S+',line))
            if len(fields)<9:raise ValueError('MOL2 atom record lacks partial charges.')
            field=fields[8];records.append((number,field.start(),field.end()))
            normalized[number]=line[:field.start()]+'0'+line[field.end():]
    if len(records)!=molecule.GetNumAtoms():raise ValueError('MOL2 charge mapping is incomplete.')
    method=charge_method(method)
    engine=(AM1BCCEngine if method=='am1bcc' else AntechamberGasteigerEngine)(prefix)
    cache=CalculationCache('mol2-charge-only',path.parent)
    # Existing partial charges are the output being replaced, not an input.
    # Ignoring only those columns makes this operation idempotent as well as
    # reusable across input copies; geometry and chemical records remain keyed.
    recipe=dict(source=identity(''.join(normalized)),method=method,engine=engine.identity(),
                adapter=file_digest(Path(__file__)),
                core=file_digest(Path(__file__).with_name('molecular_core.py')))
    key=identity(recipe)
    with cache.locked(key):
        if cache.restore(key,{'ligand':path},lambda role,p:validate_prepared_ligand(p,'mol2')) is not None:
            return path
        with preparation_workspace(cache.root) as work:
            charges=validate_charges(molecule,engine.calculate(molecule,work))
        for value,(number,start,end) in zip(charges,records):
            line=lines[number]
            lines[number]=line[:start]+f'{value:.8f}'+line[end:]
        # MOL2 charge model occupies the fifth line including the marker.
        lines[marker+4]='USER_CHARGES\n'
        write_text(''.join(lines),path)
        cache.remember(key,{'ligand':path},dict(recipe=recipe,method=method,atom_names=names,
                                              partial_charge_sum=float(charges.sum())))
    return path


def read_mapped_charges(path,molecule):
    """Match Antechamber output by stable atom IDs, not unverified row order."""
    text=Path(path).read_text();block=text.split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
    rows=[line.split() for line in block.splitlines() if line.strip()]
    names=[a.GetProp('_TriposAtomName') for a in molecule.GetAtoms()]
    if len(rows)!=len(names) or len({r[1] for r in rows})!=len(rows) or {r[1] for r in rows}!=set(names):
        raise ValueError('Antechamber changed or duplicated the ligand atom mapping.')
    table={r[1]:r for r in rows};values=[]
    ids={r[0]:r[1] for r in rows}
    bonds=text.split('@<TRIPOS>BOND',1)[1].split('@<TRIPOS>',1)[0]
    edges={frozenset((ids[r[1]],ids[r[2]])) for r in [line.split() for line in bonds.splitlines() if line.strip()]}
    expected={frozenset((names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()])) for b in molecule.GetBonds()}
    if edges!=expected:raise ValueError('Antechamber changed ligand connectivity.')
    for i,name in enumerate(names):
        row=table[name]
        if not np.allclose([float(v) for v in row[2:5]],molecule.GetConformer().GetAtomPosition(i),atol=.002,rtol=0):
            raise ValueError('Antechamber changed input coordinates unexpectedly.')
        values.append(float(row[8]))
    values=np.asarray(values,float)
    # Some AmberTools releases serialize AM1-BCC charges with only four
    # decimals.  The accumulated rounding residual can exceed the strict
    # topology tolerance for larger ligands even though SQM reports the exact
    # requested integer charge.  Correct rounding only; a chemically different
    # charge state (for example 0 instead of +1) remains a hard failure.
    residual=Chem.GetFormalCharge(molecule)-float(values.sum())
    rounding_limit=max(.0011,.00011*len(values))
    if abs(residual)<=rounding_limit:
        values=values+residual/len(values)
    return validate_charges(molecule,values)


class AmberAdapter:
    """GAFF artifacts from supplied charges; do not recalculate AM1-BCC."""
    def __init__(self,prefix='',environment=''):self.prefix=prefix;self.environment=environment
    def parameterize(self,prepared,directory,forcefield='gaff'):
        from .calculation_cache import CalculationCache,identity,copy_artifact
        from .storage import write_json
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        cache=CalculationCache('amber-ligand-adapter',directory)
        tools={name:amber_tool(name,self.prefix,environment=self.environment) for name in ('antechamber','parmchk2','tleap')}
        root=tools['tleap'].parent.parent
        library=root/'dat/leap'
        recipe=dict(molecule=Chem.MolToMolBlock(prepared.molecule),charges=prepared.partial_charges.tolist(),
            atom_names=prepared.atom_names,
            forcefield=forcefield,tools={name:file_digest(path) for name,path in tools.items()},
            parameters={str(p.relative_to(library)):file_digest(p) for p in sorted(library.rglob('*')) if p.is_file() and (forcefield in p.name)},
            adapter=file_digest(Path(__file__)))
        key=identity(recipe)
        names=('ligand.gaff.mol2','ligand.frcmod','ligand.prmtop','ligand.inpcrd','amber_preparation.json','leap.in','leap.log')
        targets={name:directory/name for name in names}
        with cache.locked(key):
            if cache.restore(key,targets):return targets
            # A failed tleap can exit successfully without writing a topology.
            # An isolated workspace prevents stale outputs from an earlier
            # parameterization being mistaken for this scientific request.
            with preparation_workspace(cache.root) as work:
                self._compute(prepared,work,forcefield)
                for name,target in targets.items():
                    if name!='amber_preparation.json':copy_artifact(work/name,target)
            write_json(dict(recipe=recipe,preparation=prepared.metadata,cache_key=key),directory/'amber_preparation.json')
            cache.remember(key,targets,dict(recipe=recipe))
        return targets

    def _compute(self,prepared,directory,forcefield):
        if forcefield not in ('gaff','gaff2'):raise ValueError('Select GAFF or GAFF2 explicitly.')
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        tool=amber_tool('antechamber',self.prefix,environment=self.environment)
        Dock6Adapter().export(prepared,directory/'ligand.mol2')
        write_text('\n'.join(str(float(v)) for v in prepared.partial_charges)+'\n',directory/'charges.dat')
        run_amber_tool('antechamber',['-i','ligand.mol2','-fi','mol2','-o','ligand.gaff.mol2','-fo','mol2',
            '-at',forcefield,'-c','rc','-cf','charges.dat','-nc',str(Chem.GetFormalCharge(prepared.molecule)),'-an','n','-du','n','-s','2','-pf','y'],
            cwd=directory,prefix=self.prefix,environment=self.environment)
        validate_charges(prepared.molecule,read_mapped_charges(directory/'ligand.gaff.mol2',prepared.molecule))
        run_amber_tool('parmchk2',['-i','ligand.gaff.mol2','-f','mol2','-o','ligand.frcmod','-s','2' if forcefield=='gaff2' else '1'],cwd=directory,prefix=self.prefix,environment=self.environment)
        write_text(f'source leaprc.{forcefield}\nloadamberparams ligand.frcmod\nx=loadmol2 ligand.gaff.mol2\ncheck x\nsaveamberparm x ligand.prmtop ligand.inpcrd\nquit\n',directory/'leap.in')
        run_amber_tool('tleap',['-f','leap.in'],cwd=directory,prefix=self.prefix,environment=self.environment)
        from .md_analysis import topology
        info=topology(directory/'ligand.prmtop')
        if info['natom']!=prepared.molecule.GetNumAtoms() or info['names']!=prepared.atom_names or not np.allclose(info['charges'],prepared.partial_charges,atol=1e-4,rtol=0):
            raise ValueError('Amber topology changed the supplied atom charges or mapping.')
        for name in ('ligand.frcmod','ligand.inpcrd'):
            if not (directory/name).is_file() or not (directory/name).stat().st_size:raise ValueError('Incomplete Amber parameterization: '+name)
        return {name:directory/name for name in ('ligand.gaff.mol2','ligand.frcmod','ligand.prmtop','ligand.inpcrd')}
