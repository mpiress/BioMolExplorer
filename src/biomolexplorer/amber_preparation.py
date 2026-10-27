"""Shared AmberTools topology/minimization for docking and molecular preparation.

This is an explicit ff14SB/GAFF vacuum protocol, not a numerical MMTK emulation.
The historical DOCK6 MD bundle's ff94 topology is intentionally not reused for
ff14SB receptors. Ligand GAFF artifacts use the existing AmberAdapter cache.
"""
from pathlib import Path
import re

import numpy as np
from rdkit import Chem

from .calculation_cache import CalculationCache,identity
from .md_analysis import topology,PROTEIN
from .molecular_adapters import AmberAdapter
from .molecular_core import PreparedMolecule,preparation_workspace
from .preparation_tools import amber_tool,run_amber_tool
from .stage_cache import file_digest
from .storage import write_text,write_json

MINIMIZATION=dict(imin=1,maxcyc=110,ncyc=100,ntmin=1,ntb=0,igb=0,cut=999.,ntc=1,ntf=1,drms=1e-4)


def parameterization_identity(prefix='',minimize=True,environment=''):
    tools={name:amber_tool(name,prefix,environment=environment) for name in ('tleap','antechamber','parmchk2')+ (('sander',) if minimize else ())}
    root=tools['tleap'].parent.parent
    return dict(tools={name:file_digest(path) for name,path in tools.items()},
                parameters={str(p.relative_to(root)):file_digest(p) for p in sorted((root/'dat/leap').rglob('*')) if p.is_file()},
                protocol=MINIMIZATION if minimize else None,implementation=file_digest(Path(__file__)))


def coordinates(path,count=None):
    lines=Path(path).read_text().splitlines();n=int(lines[1].split()[0])
    values=[float(line[i:i+12]) for line in lines[2:] for i in range(0,len(line),12) if line[i:i+12].strip()]
    if count is not None and n!=count:raise ValueError('Amber coordinate/topology atom count mismatch.')
    xyz=np.asarray(values[:3*n]).reshape(n,3)
    if not np.isfinite(xyz).all():raise ValueError('Amber produced non-finite coordinates.')
    return xyz


def peptide_breaks(text):
    """Keep missing peptide segments explicit in Amber and downstream PDBs."""
    from .structure_selection import CAPS
    new=[];previous=None;breaks=[]
    for line in text.splitlines():
        if line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() in PROTEIN|CAPS:
            if line[12:16].strip()=='N' and previous is not None:
                position=np.array([float(line[k:k+8]) for k in (30,38,46)])
                if np.linalg.norm(position-previous)>2.2:new.append('TER');breaks.append(line[21:27])
            if line[12:16].strip()=='C':previous=np.array([float(line[k:k+8]) for k in (30,38,46)])
        if line.startswith('TER'):previous=None
        new.append(line)
    return '\n'.join(new)+'\n',breaks


def minimize_topology(prmtop,inpcrd,directory,prefix='',environment=''):
    """Same 100 SD + 10 CG cycle budget as the legacy default, vacuum, no freeze.

    Sander's line search and convergence criteria differ from MMTK. No solvent,
    restraints, equilibration or force-field switch is silently introduced.
    """
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    tool=amber_tool('sander',prefix,environment=environment);info=topology(prmtop);initial=coordinates(inpcrd,info['natom'])
    cache=CalculationCache('amber-minimization',directory)
    recipe=dict(topology=file_digest(prmtop),coordinates=file_digest(inpcrd),protocol=MINIMIZATION,
                tool=file_digest(tool),implementation=file_digest(Path(__file__)))
    key=identity(recipe)
    outputs={'min.rst7':directory/'min.rst7','min.out':directory/'min.out','minimization.json':directory/'minimization.json'}
    with cache.locked(key):
        reused=cache.restore(key,outputs)
        if reused is None:
            with preparation_workspace(cache.root) as work:
                mdin='Native docking minimization: ff14SB/GAFF vacuum\n&cntrl\n '+',\n '.join(f'{k}={v}' for k,v in MINIMIZATION.items())+',\n ntpr=1, ntxo=1,\n/\n'
                write_text(mdin,work/'min.in')
                run_amber_tool('sander',['-O','-i','min.in','-p',str(Path(prmtop).resolve()),
                             '-c',str(Path(inpcrd).resolve()),'-o','min.out','-r','min.rst7'],cwd=work,
                             prefix=prefix,environment=environment)
                final=coordinates(work/'min.rst7',info['natom']);text=(work/'min.out').read_text()
                energies=[float(m[1].replace('D','E')) for m in re.finditer(r'NSTEP\s+ENERGY\s+RMS\s+GMAX[^\n]*\n\s*\d+\s+([-+\d.EeDd]+)',text)]
                if not energies or not np.isfinite(energies).all() or energies[-1]>energies[0]+1e-3:
                    raise ValueError('Amber minimization did not provide a finite, non-increasing potential energy.')
                heavy=info['masses']>2.5;displacement=np.linalg.norm(final-initial,axis=1)
                bonds=info['bonds']
                final_bonds=[float(np.linalg.norm(final[a]-final[b])) for a,b in bonds]
                if any(length<.5 or length>2.8 for length in final_bonds):
                    raise ValueError('Amber minimization introduced an implausible covalent bond length.')
                from scipy.spatial import cKDTree
                bonded={tuple(sorted((a,b))) for a,b in bonds}
                def clashes(xyz):
                    return sum(tuple(sorted(pair)) not in bonded for pair in cKDTree(xyz).query_pairs(.8))
                if float(displacement[heavy].max())>3.:
                    raise ValueError('Unexpected heavy-atom displacement (>3 A) during short vacuum minimization; inspect this system explicitly.')
                report=dict(protocol=recipe,methodological_change='MMTK -> AmberTools/sander; same 100 SD + 10 CG budget, different optimizer implementation',
                            initial_energy=energies[0],final_energy=energies[-1],
                            heavy_atom_rmsd=float(np.sqrt(np.mean(displacement[heavy]**2))),
                            max_heavy_displacement=float(displacement[heavy].max()),
                            initial_severe_clashes=clashes(initial),final_severe_clashes=clashes(final),
                            minimum_bond_length=min(final_bonds) if final_bonds else None,
                            maximum_bond_length=max(final_bonds) if final_bonds else None,
                            restrained_atoms=[],solvent='vacuum',converged='CONVERGENCE CRITERION' in text.upper())
                write_json(report,work/'minimization.json')
                cache.remember(key,{name:work/name for name in outputs},dict(recipe=recipe))
                from .calculation_cache import copy_artifact
                for name,target in outputs.items():copy_artifact(work/name,target)
    return coordinates(outputs['min.rst7'],info['natom']),outputs['minimization.json']


def minimize_ligand(prepared,directory,prefix='',environment=''):
    artifacts=AmberAdapter(prefix,environment).parameterize(prepared,directory,forcefield='gaff')
    xyz,report=minimize_topology(artifacts['ligand.prmtop'],artifacts['ligand.inpcrd'],directory,prefix,environment)
    molecule=Chem.Mol(prepared.molecule)
    for index,position in enumerate(xyz):molecule.GetConformer().SetAtomPosition(index,position)
    import json
    metadata=dict(prepared.metadata,minimization=json.loads(report.read_text()))
    return PreparedMolecule(molecule,prepared.partial_charges.copy(),metadata)


def protonation(structure,ph,remove_hydrogens=True):
    """Preserve named variants/H evidence; otherwise record a model-compound rule.

    These reference-pKa choices are not a site-specific pKa prediction. Histidine
    defaults to HIE above pH 6, HIP below; disulfides are detected separately.
    """
    result=structure.clone();states=[];disulfides=[];cysteines=[]
    for chain in result[0]:
        for residue in chain:
            names={a.name for a in residue};old=residue.name;chosen=old;method='explicit residue variant'
            if old=='HIS':
                chosen='HIP' if {'HD1','HE2'}<=names else 'HID' if 'HD1' in names else 'HIE' if 'HE2' in names else 'HIP' if ph<6. else 'HIE'
                method='input H evidence' if {'HD1','HE2'} & names else 'reference pKa 6.0; neutral tautomer HIE'
            elif old in ('ASP','GLU','CYS','LYS'):
                choices={'ASP':('ASH',3.9,'HD2'),'GLU':('GLH',4.2,'HE2'),'CYS':('CYM',8.3,'HG'),'LYS':('LYN',10.5,'HZ3')}
                variant,pka,hydrogen=choices[old]
                if old in ('ASP','GLU'):chosen=variant if hydrogen in names or ph<pka else old
                else:chosen=old if hydrogen in names else variant if ph>pka else old
                method='input H evidence' if hydrogen in names else f'reference model-compound pKa {pka}'
            residue.name=chosen
            states.append(dict(chain=chain.name,number=residue.seqid.num,insertion_code=residue.seqid.icode.strip(),
                               input=old,state=chosen,method=method,ph=float(ph)))
            if old in ('CYS','CYX'):
                sg=next((a for a in residue if a.name=='SG'),None)
                if sg:cysteines.append((len(states)-1,residue,sg))
    for i,(index,residue,sg) in enumerate(cysteines):
        partners=[other for other in cysteines if other[1] is not residue and sg.pos.dist(other[2].pos)<2.3]
        if len(partners)>1:raise ValueError('Ambiguous cysteine coordination; explicit disulfide/metal parameters required.')
        if residue.name=='CYX' and not partners:
            raise ValueError('Explicit CYX has no disulfide partner; provide a dedicated bonded topology.')
        if partners and any(atom.name=='HG' for atom in residue):
            raise ValueError('Explicit thiol hydrogen conflicts with a disulfide-distance contact; resolve the input state.')
        if partners:
            residue.name='CYX';states[index].update(state='CYX',method='SG-SG distance <2.3 A; explicit disulfide bond')
            j=partners[0][0]
            if index<j:
                identity=lambda state:dict(chain=state['chain'],number=state['number'],
                                           insertion_code=state['insertion_code'])
                # Keep both the crystallographic identity (for auditability)
                # and the dense LEaP address used after writing a deliberately
                # renumbered input PDB.  Neither PDB sequence numbers nor raw
                # Gemmi positions alone are reliable across missing segments.
                disulfides.append(dict(residues=[identity(states[index]),identity(states[j])],
                                         amber_residues=[index+1,j+1]))
    if remove_hydrogens:result.remove_hydrogens()
    return result,states,disulfides


def receptor_parameters(structure,directory,ph=7.4,options=None,prefix='',chemistry_folder=None,config=None,environment=''):
    """ff14SB protein/water; separately parameterized GAFF organic cofactors."""
    from .structure_selection import WATER,METALS,CAPS,write_structure
    options=options or {};directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    cofactors={};config=config or {}
    for chain in structure[0]:
        for residue in chain:
            if any(a.element.name in METALS for a in residue):raise ValueError('Selected receptor contains metal/ion coordination; supply a dedicated parameterized system or explicitly choose legacy_chimera.')
            if residue.name not in PROTEIN|CAPS|WATER:
                if not chemistry_folder:raise ValueError('Nonstandard receptor component requires a trusted chemical template: '+residue.name)
                from .structure_selection import crystal_molecule
                from .molecular_core import MolecularPreparationCore
                record=['cofactor',residue.name,residue.seqid.num,chain.name]
                local=dict(config,ligand_chain=chain.name,ligand_insertion_code=residue.seqid.icode.strip())
                molecule,mapping,template=crystal_molecule(structure,record,local,chemistry_folder)
                prepared=MolecularPreparationCore(directory,prefix,environment).prepare(molecule,method=options.get('charge_type','gas'),
                    ph=ph,charge_engine='antechamber',add_hydrogens=options.get('add_hydrogens',True),
                    remove_hydrogens=options.get('remove_hydrogens',True))
                component=directory/('cofactor-'+residue.name+'-'+str(residue.seqid.num)+chain.name)
                artifacts=AmberAdapter(prefix,environment).parameterize(prepared,component)
                if residue.name in cofactors and file_digest(cofactors[residue.name]['ligand.gaff.mol2'])!=file_digest(artifacts['ligand.gaff.mol2']):
                    raise ValueError('Different states of the same cofactor require distinct residue identifiers.')
                cofactors[residue.name]=artifacts
    from .structure_selection import atom_identity
    input_identities={(chain.name,residue.seqid.num,residue.seqid.icode.strip(),atom.name):atom_identity(chain,residue,atom,config.get('structure',{}).get('model',1))
                      for chain in structure[0] for residue in chain for atom in residue}
    normalized,states,disulfides=protonation(structure,float(ph),options.get('remove_hydrogens',True))
    residues=[r for c in normalized[0] for r in c]
    # Amber reads one chosen alternate only; identity is retained in the map.
    for residue in residues:
        if residue.name in WATER:residue.name='WAT'
        for atom in residue:atom.altloc='\x00'
    leap=amber_tool('tleap',prefix,environment=environment);root=leap.parent.parent
    recipe=dict(structure=normalized.make_pdb_string(),original_identities=list(input_identities.values()),ph=ph,states=states,disulfides=disulfides,
                forcefield='ff14SB',options=options,tool=file_digest(leap),implementation=file_digest(Path(__file__)),
                cofactors={name:{key:file_digest(value) for key,value in artifacts.items()} for name,artifacts in cofactors.items()},
                parameters={str(p.relative_to(root)):file_digest(p) for p in sorted((root/'dat/leap').rglob('*')) if p.is_file()})
    cache=CalculationCache('receptor-ff14SB',directory);key=identity(recipe)
    targets={name:directory/name for name in ('receptor.prmtop','receptor.inpcrd','receptor.amber.pdb','receptor.parameters.json','leap.log')}
    with cache.locked(key):
        if cache.restore(key,targets) is None:
            with preparation_workspace(cache.root) as work:
                # Report crystallographic break identifiers even though the
                # private LEaP input below uses dense residue numbering.
                _,breaks=peptide_breaks(normalized.make_pdb_string())
                amber_input=normalized.clone()
                from .structure_selection import gemmi_module
                for index,residue in enumerate(r for chain in amber_input[0] for r in chain):
                    residue.seqid=gemmi_module().SeqId(index+1,' ')
                write_structure(amber_input,work/'input.pdb')
                # Missing residues are not fabricated. Distant peptide termini
                # become explicit chain breaks instead of accidental bonds.
                pdb,_=peptide_breaks((work/'input.pdb').read_text())
                write_text(pdb,work/'input.pdb')
                commands=['source leaprc.protein.ff14SB','source leaprc.water.tip3p']
                if cofactors:commands.append('source leaprc.gaff')
                for name,artifacts in cofactors.items():
                    text=artifacts['ligand.gaff.mol2'].read_text()
                    lines=text.splitlines();inside=False
                    for number,line in enumerate(lines):
                        if line.startswith('@<TRIPOS>'):inside=line=='@<TRIPOS>ATOM'
                        elif inside and line.strip():
                            fields=line.split();fields[7]=name;lines[number]=' '.join(fields)
                    write_text('\n'.join(lines)+'\n',work/(name+'.mol2'))
                    from .calculation_cache import copy_artifact
                    copy_artifact(artifacts['ligand.frcmod'],work/(name+'.frcmod'))
                    commands += [f'loadamberparams {name}.frcmod',f'{name}=loadmol2 {name}.mol2']
                commands.append('x=loadpdb input.pdb')
                commands += [f'bond x.{bond["amber_residues"][0]}.SG x.{bond["amber_residues"][1]}.SG'
                             for bond in disulfides]
                commands += ['check x','saveamberparm x receptor.prmtop receptor.inpcrd','savepdb x receptor.amber.pdb','quit']
                write_text('\n'.join(commands)+'\n',work/'leap.in')
                result=run_amber_tool('tleap',['-f','leap.in'],cwd=work,prefix=prefix,environment=environment)
                info=topology(work/'receptor.prmtop');xyz=coordinates(work/'receptor.inpcrd',info['natom'])
                if info['nres']!=len(residues):raise ValueError('Amber changed the receptor residue mapping.')
                original=[];added=[]
                for i,residue in enumerate(residues):
                    source={a.name:a for a in residue}
                    indices=np.flatnonzero(info['residue']==i);names={info['names'][j] for j in indices}
                    if source.keys()-names:raise ValueError('Amber lost input atoms: '+str(source.keys()-names))
                    for j in indices:
                        name=info['names'][j]
                        if name in source:
                            atom=source[name];pos=np.array([atom.pos.x,atom.pos.y,atom.pos.z])
                            if np.linalg.norm(xyz[j]-pos)>.002:raise ValueError('Amber unexpectedly moved input receptor atoms before minimization.')
                            original.append(dict(index=int(j),original=input_identities[(states[i]['chain'],states[i]['number'],states[i]['insertion_code'],name)]))
                        else:added.append(dict(index=int(j),residue=i+1,name=name,heavy=bool(info['masses'][j]>2.5)))
                if not options.get('add_hydrogens',True) and any(not a['heavy'] for a in added):
                    raise ValueError('Amber templates require absent hydrogen atoms; enable hydrogen addition or provide a complete explicitly protonated receptor.')
                write_json(dict(recipe=recipe,states=states,original_atom_mapping=original,added_atoms=added,
                                chain_breaks=breaks,missing_residues='not modeled; chain breaks retained',
                                net_charge=float(info['charges'].sum()),forcefield='ff14SB'),work/'receptor.parameters.json')
                cache.remember(key,{name:work/name for name in targets},dict(recipe=recipe))
                from .calculation_cache import copy_artifact
                for name,target in targets.items():copy_artifact(work/name,target)
    return targets,normalized


def export_receptor(parameters,original,output_mol2,output_pdb,*,minimize=False,prefix='',directory=None,environment=''):
    """Use Amber topology charges, never ligand charge models on protein."""
    from openbabel import openbabel as ob
    from .structure_selection import gemmi_module,write_structure
    info=topology(parameters['receptor.prmtop'])
    xyz=coordinates(parameters['receptor.inpcrd'],info['natom'])
    if minimize:xyz,_=minimize_topology(parameters['receptor.prmtop'],parameters['receptor.inpcrd'],directory,prefix,environment)
    import json
    states=json.loads(parameters['receptor.parameters.json'].read_text())['states']
    amber=gemmi_module().read_structure(str(parameters['receptor.amber.pdb']))
    amber_residues=[r for c in amber[0] for r in c];result=original.clone();index=0;mapping=[]
    for chain in result[0]:
        for ri in range(len(chain)):
            source=chain[ri];prepared=amber_residues[index].clone();prepared.seqid=source.seqid;index+=1
            prepared.name=states[index-1]['state']
            # Keep user atom labels and altloc identity for original atoms.
            lookup={a.name:a for a in source}
            for atom in prepared:
                if atom.name in lookup:atom.altloc=lookup[atom.name].altloc
            del chain[ri]
            chain.add_residue(prepared,ri)
    flat=[a for c in result[0] for r in c for a in r]
    if len(flat)!=info['natom'] or [a.name for a in flat]!=info['names']:raise ValueError('Amber PDB/topology mapping mismatch.')
    for atom,pos in zip(flat,xyz):atom.pos=gemmi_module().Position(*pos)
    write_structure(result,output_pdb)
    pdb,_=peptide_breaks(Path(output_pdb).read_text());write_text(pdb,output_pdb)
    # Open Babel supplies SYBYL typing for the standard residue geometry; its
    # charge model is disabled by setting the exact topology vector afterward.
    conversion=ob.OBConversion();conversion.SetInAndOutFormats('pdb','mol2');molecule=ob.OBMol()
    if not conversion.ReadString(molecule,pdb):raise ValueError('Cannot export parameterized receptor.')
    if molecule.NumAtoms()!=info['natom']:raise ValueError('Receptor export changed atom count.')
    for atom,charge,pos in zip(ob.OBMolAtomIter(molecule),info['charges'],xyz):
        if np.linalg.norm(np.array([atom.GetX(),atom.GetY(),atom.GetZ()])-pos)>.002:raise ValueError('Receptor export atom-order mismatch.')
        atom.SetPartialCharge(float(charge))
    molecule.SetPartialChargesPerceived(True)
    text=conversion.WriteString(molecule);lines=text.splitlines();marker=lines.index('@<TRIPOS>MOLECULE')
    lines[marker+1]=Path(output_mol2).stem;lines[marker+3]='PROTEIN';lines[marker+4]='USER_CHARGES'
    residues=[(chain,residue) for chain in result[0] for residue in chain]
    counts=lines[marker+2].split();counts[2]=str(len(residues));lines[marker+2]=' '.join(counts)
    atom_section=lines.index('@<TRIPOS>ATOM');in_atoms=False;j=0
    for number,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):in_atoms=line=='@<TRIPOS>ATOM'
        elif in_atoms and line.strip():
            fields=line.split();r=int(info['residue'][j]);fields[6]=str(r+1);fields[7]=residues[r][1].name
            lines[number]=' '.join(fields);j+=1
    # DOCK6's historical grid reader consumes both optional MOL2 header
    # records. A single blank line makes it eat ATOM as the comment record.
    while atom_section>0 and lines[atom_section-1]=='':
        del lines[atom_section-1];atom_section-=1
    lines[atom_section:atom_section]=['','']
    lines.append('@<TRIPOS>SUBSTRUCTURE')
    for index,(chain,residue) in enumerate(residues):
        first=int(np.flatnonzero(info['residue']==index)[0])+1
        label=f'{residue.name}{residue.seqid.num}{residue.seqid.icode.strip()}'
        lines.append(f'{index+1} {label} {first} RESIDUE 1 {chain.name or "****"} {residue.name} 0 ROOT')
    write_text('\n'.join(lines)+'\n',output_mol2)
    return result
