"""Conservative, geometry-based interactions in the docked coordinate frame.

Thresholds follow PLIP's published configuration: hydrophobic 4 Å; H-bond
donor/acceptor 4.1 Å with D-H-A >100°; aromatic centres 5.5 Å, 30° angular
tolerance, 2 Å offset. This is a limited classifier, not the PLIP engine.
No hydrogen positions or protonation states are invented for classification.
"""
import math
from collections import defaultdict

from .docking_scene import atoms,PROTEIN

RINGS={
    'PHE':[('CG','CD1','CE1','CZ','CE2','CD2')],
    'TYR':[('CG','CD1','CE1','CZ','CE2','CD2')],
    'TRP':[('CG','CD1','NE1','CE2','CD2'),('CD2','CE2','CZ2','CH2','CZ3','CE3')],
    'HIS':[('CG','ND1','CE1','NE2','CD2')],
}
for _name in ('HID','HIE','HIP'):RINGS[_name]=RINGS['HIS']
ACCEPTORS={'ASP':{'OD1','OD2'},'GLU':{'OE1','OE2'},'ASN':{'OD1'},'GLN':{'OE1'},
           'SER':{'OG'},'THR':{'OG1'},'TYR':{'OH'},'HID':{'NE2'},'HIE':{'ND1'}}
DONORS={'LYS':{'NZ'},'ARG':{'NE','NH1','NH2'},'ASN':{'ND2'},'GLN':{'NE2'},
        'SER':{'OG'},'THR':{'OG1'},'TYR':{'OH'},'TRP':{'NE1'},
        'HIS':{'ND1','NE2'},'HID':{'ND1'},'HIE':{'NE2'},'HIP':{'ND1','NE2'}}
COLORS={'pi_parallel':'#7c3aed','pi_t':'#db2777','hydrogen_bond':'#0284c7','hydrophobic':'#d97706'}


def ligand_molecule(layer,smiles=None):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    text,fmt=layer['data'],layer['format']
    if fmt=='mol2':mol=Chem.MolFromMol2Block(text,sanitize=True,removeHs=True)
    elif fmt=='sdf':mol=Chem.MolFromMolBlock(text.split('$$$$')[0],removeHs=True)
    else:mol=None
    if smiles:
        template=Chem.MolFromSmiles(smiles)
        if template is None:return None
        if mol is None:
            records=atoms(text,fmt)
            block='\n'.join(f'HETATM{i:5d} {a["atom"][:4]:<4} LIG A   1    '+
                ''.join(f'{v:8.3f}' for v in a['xyz'])+'  1.00  0.00          '+f'{a["elem"]:>2}'
                for i,a in enumerate(records,1))+'\nEND\n'
            mol=Chem.MolFromPDBBlock(block,sanitize=False,removeHs=False)
        if mol is None or mol.GetNumAtoms()!=template.GetNumHeavyAtoms():return None
        mol=AllChem.AssignBondOrdersFromTemplate(template,mol)
        # Stereo labels may be absent in docking files; connectivity must agree.
        if Chem.MolToSmiles(mol,isomericSmiles=False)!=Chem.MolToSmiles(template,isomericSmiles=False):return None
    if mol is None or not mol.GetNumConformers():return None
    Chem.SanitizeMol(mol)
    if len(Chem.GetMolFrags(mol))!=1:return None
    return mol


def ring_geometry(points):
    import numpy as np
    points=np.asarray(points,dtype=float);centre=points.mean(axis=0)
    _,singular,vectors=np.linalg.svd(points-centre)
    if singular[1]<.1 or singular[-1]/math.sqrt(len(points))>.15:return None
    return centre,vectors[-1]


def interaction_diagram(receptor,pose,smiles=None):
    """Return disabled capability for incomplete or ambiguous chemistry."""
    try:
        mol=ligand_molecule(pose,smiles)
        if mol is None:return dict(available=False,reason='chemical_topology')
        return _diagram(receptor,pose,mol)
    except (ValueError,RuntimeError,IndexError,KeyError):
        return dict(available=False,reason='chemical_topology')


def _diagram(receptor,pose,mol):
    import numpy as np
    from rdkit import Chem,RDConfig
    from rdkit.Chem import ChemicalFeatures,rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from pathlib import Path
    protein=atoms(receptor['data'],receptor['format'],include_hydrogens=True)
    residues=defaultdict(list)
    for atom in protein:
        if atom['resn'] in PROTEIN:
            residues[tuple(atom[k] for k in ('chain','resi','icode','resn'))].append(atom)
    if not residues:return dict(available=False,reason='residue_identity')
    conformer=mol.GetConformer()
    xyz=np.asarray(conformer.GetPositions())
    if not np.isfinite(xyz).all():return dict(available=False,reason='coordinates')
    ligand_rings=[]
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            geometry=ring_geometry(xyz[list(ring)])
            if geometry:ligand_rings.append((list(ring),geometry))
    features=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef')).GetFeaturesForMol(mol)
    acceptors={i for f in features if f.GetFamily()=='Acceptor' for i in f.GetAtomIds()}
    donors={i for f in features if f.GetFamily()=='Donor' for i in f.GetAtomIds()}
    hydrophobic={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6 and
                 all(n.GetAtomicNum() in (1,6) for n in a.GetNeighbors())}
    ligand_h=[a for a in atoms(pose['data'],pose['format'],include_hydrogens=True) if a['elem'] in ('H','D')]
    interactions=[];seen=set()
    def add(kind,key,indices,distance,protein_atom='',angle=None):
        identity=(kind,key,tuple(indices))
        if identity in seen:return
        seen.add(identity)
        row=dict(kind=kind,color=COLORS[kind],chain=key[0],resi=key[1],icode=key[2],resn=key[3],
            ligand_atoms=list(indices),receptor_atom=protein_atom,distance=round(float(distance),3))
        if angle is not None:row['angle']=round(float(angle),1)
        interactions.append(row)
    def hydrogen_bond(donor,hydrogens,acceptor):
        distance=math.dist(donor,acceptor)
        if not .5<distance<=4.1:return None
        for hydrogen in hydrogens:
            if not .6<math.dist(donor,hydrogen)<=1.25:continue
            # Match each H only to its unique nearest donor; callers supply
            # attached H positions rather than generating a hydrogen geometry.
            first=np.asarray(donor)-hydrogen;second=np.asarray(acceptor)-hydrogen
            norm=np.linalg.norm(first)*np.linalg.norm(second)
            if not norm:continue
            angle=math.degrees(math.acos(float(np.clip(np.dot(first,second)/norm,-1,1))))
            if angle>100 and math.dist(hydrogen,acceptor)<distance:return distance,angle
        return None
    for key,records in residues.items():
        heavy=[a for a in records if a['elem'] not in ('H','D')]
        if not heavy or min(math.dist(a['xyz'],p) for a in heavy for p in xyz)>7.5:continue
        names={a['atom']:a for a in heavy}
        for ring in RINGS.get(key[3],[]):
            if not all(name in names for name in ring):continue
            geometry=ring_geometry([names[name]['xyz'] for name in ring])
            if geometry is None:continue
            centre,normal=geometry
            for indices,(lc,ln) in ligand_rings:
                delta=lc-centre;distance=np.linalg.norm(delta)
                angle=math.degrees(math.acos(float(np.clip(abs(np.dot(normal,ln)),0,1))))
                offset=min(np.linalg.norm(delta-np.dot(delta,normal)*normal),np.linalg.norm(delta-np.dot(delta,ln)*ln))
                if .5<distance<=5.5 and offset<=2 and (angle<=30 or angle>=60):
                    add('pi_parallel' if angle<=30 else 'pi_t',key,indices,distance,'/'.join(ring),angle)
        # Keep only the closest hydrophobic atom pair per residue.
        pairs=[(math.dist(a['xyz'],xyz[i]),a,i) for a in heavy if a['elem']=='C' and a['atom'] not in ('C','CA')
               and not any(b['elem'] in ('N','O','S') and math.dist(a['xyz'],b['xyz'])<1.9 for b in heavy)
               for i in hydrophobic]
        if pairs:
            distance,a,i=min(pairs,key=lambda item:item[0])
            if .5<distance<=4:add('hydrophobic',key,[i],distance,a['atom'])
        for a in heavy:
            is_donor=(a['atom']=='N' and key[3]!='PRO') or a['atom'] in DONORS.get(key[3],set())
            is_acceptor=a['atom'] in ('O','OXT') or a['atom'] in ACCEPTORS.get(key[3],set())
            attached=[np.asarray(h['xyz']) for h in records if h['elem'] in ('H','D') and
                min(heavy,key=lambda other:math.dist(other['xyz'],h['xyz'])) is a]
            if is_donor:
                for i in acceptors:
                    bond=hydrogen_bond(a['xyz'],attached,xyz[i])
                    if bond:add('hydrogen_bond',key,[i],bond[0],a['atom'],bond[1])
            if is_acceptor:
                for i in donors:
                    attached=[np.asarray(h['xyz']) for h in ligand_h if int(np.argmin(np.linalg.norm(xyz-h['xyz'],axis=1)))==i]
                    bond=hydrogen_bond(xyz[i],attached,a['xyz'])
                    if bond:add('hydrogen_bond',key,[i],bond[0],a['atom'],bond[1])
    drawing=Chem.Mol(mol);drawing.RemoveAllConformers();rdDepictor.Compute2DCoords(drawing)
    drawer=rdMolDraw2D.MolDraw2DSVG(600,440)
    drawer.DrawMolecule(drawing);drawer.FinishDrawing()
    positions=[[drawer.GetDrawCoords(i).x,drawer.GetDrawCoords(i).y] for i in range(mol.GetNumAtoms())]
    return dict(available=True,molecule_svg=drawer.GetDrawingText(),atom_positions=positions,
        interactions=interactions,types=list(COLORS),method='geometry',
        hydrogen_bonds_require_explicit_hydrogens=True)
