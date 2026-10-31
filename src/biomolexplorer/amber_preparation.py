"""Legacy structure utilities; execution uses Chimera exclusively."""
from pathlib import Path
import numpy as np
from .md_analysis import PROTEIN
from .chimera_minimization import minimize_ligand

def parameterization_identity(*args,**kwargs):
    raise ValueError('The standalone ff14SB protocol is retired; configure Chimera preparation.')

def minimize_topology(*args,**kwargs):
    raise ValueError('Standalone minimization is retired; use Chimera without silently changing force fields.')

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


def receptor_parameters(*args,**kwargs):
    raise ValueError('The experimental ff14SB backend is retired; native DOCK6 ff94 is not substituted.')

def export_receptor(*args,**kwargs):
    raise ValueError('Use Chimera receptor preparation.')
