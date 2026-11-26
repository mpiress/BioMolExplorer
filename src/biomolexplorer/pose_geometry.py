"""Physical plausibility of docked poses, inspired by the PoseBusters criteria.

Buttenschoen, M., Morris, G. M., Deane, C. M. (2024). PoseBusters: AI-based
docking methods fail to generate physically valid poses or generalise to novel
sequences. Chemical Science, 15, 3130-3139.

Every scientific dependency (RDKit, NumPy, SciPy) is imported inside the
functions so that importing this module, the UI or the catalog never loads the
scientific stack.  Each criterion reports its own measurement and flag; a
measurement that cannot be obtained is reported as unavailable instead of
silently failing a pose.
"""
import math
import re
from pathlib import Path

# Bondi van der Waals radii (Angstrom); unlisted elements use the carbon radius.
VDW_RADII = {'H':1.10,'HE':1.40,'LI':1.81,'BE':1.53,'B':1.92,'C':1.70,'N':1.55,'O':1.52,'F':1.47,
             'NE':1.54,'NA':2.27,'MG':1.73,'AL':1.84,'SI':2.10,'P':1.80,'S':1.80,'CL':1.75,
             'AR':1.88,'K':2.75,'CA':2.31,'MN':1.97,'FE':1.94,'CO':1.92,'NI':1.84,'CU':1.40,
             'ZN':1.39,'SE':1.90,'BR':1.85,'I':1.98,'CD':1.58,'HG':1.55}
DEFAULT_RADIUS = 1.70
SOLVENT_RESIDUES = {'HOH','WAT','DOD','SOL','TIP3','TIP','H2O'}
# Deviations are compared against a relaxed conformer of the same molecule.
# A contact closer than ``vdw_clash_scale`` times the sum of the van der Waals
# radii is a clash; 0.75 keeps hydrogen bonds and tight polar contacts valid,
# as in the PoseBusters intermolecular distance criterion.
TOLERANCES = {'vdw_clash_scale':0.75,'delta_energy_limit':100.0,'bond_tolerance':0.25,
              'angle_tolerance':25.0,'planarity_tolerance':0.25,'minimization_steps':400}
REASONS = {
    'unreadable':'Arquivo da pose ilegível ou sem coordenadas.',
    'sanitization':'A pose não passa na sanitização química (valência ou aromaticidade inválida).',
    'fragments':'A pose contém fragmentos desconexos.',
    'coordinates':'Coordenadas ausentes, repetidas ou não finitas.',
    'bond_lengths':'Comprimentos de ligação fora da tolerância física.',
    'bond_angles':'Ângulos de ligação fora da tolerância física.',
    'planarity':'Anéis aromáticos ou ligações duplas não planares.',
    'internal_energy':'Energia interna muito acima do confôrmero relaxado.',
    'steric_clash':'Choque estérico com o receptor ou cofatores.',
}


def element_of(name, hint=''):
    """Element symbol from a PDB/PDBQT atom name, preferring the explicit column."""
    hint = (hint or '').strip().upper()
    if hint and hint.isalpha():
        return hint if hint in VDW_RADII else hint[:1]
    token = re.sub(r'[^A-Za-z]','',(name or '').strip()).upper()
    if not token:
        return 'C'
    if token[:2] in VDW_RADII and not token[:1] in ('C','N','O','S','P','H'):
        return token[:2]
    return token[:1]


def radius_of(element):
    return VDW_RADII.get((element or '').upper(),DEFAULT_RADIUS)


def structure_atoms(path, *, heavy_only=True, exclude_residues=()):
    """Heavy-atom coordinates and radii of a PDB/PDBQT structure, without RDKit."""
    path = Path(path)
    excluded = {str(name).strip().upper() for name in exclude_residues} | SOLVENT_RESIDUES
    coordinates = []; radii = []; elements = []
    with path.open(encoding='utf-8',errors='replace') as stream:
        for line in stream:
            if not line.startswith(('ATOM  ','HETATM')) or len(line) < 54:
                continue
            if line[17:20].strip().upper() in excluded:
                continue
            element = element_of(line[12:16],line[76:78] if len(line) >= 78 else '')
            if heavy_only and element == 'H':
                continue
            try:
                point = (float(line[30:38]),float(line[38:46]),float(line[46:54]))
            except ValueError:
                continue
            if not all(map(math.isfinite,point)):
                continue
            coordinates.append(point); radii.append(radius_of(element)); elements.append(element)
    return coordinates,radii,elements


# AutoDock atom types carry the element plus an interaction role.
AUTODOCK_ELEMENTS = {'A':'C','NA':'N','NS':'N','OA':'O','OS':'O','SA':'S','HD':'H','HS':'H','G0':'C','CG0':'C'}


def pdbqt_to_pdb(text, model=1):
    """One PDB block from a PDBQT file, with the element column restored.

    The AutoDock type in columns 78-79 is the only reliable element source for
    a docked ligand, so it is translated instead of discarded with the charges.
    """
    lines = []; current = 0; selected = model or 1
    for line in text.splitlines():
        if line.startswith('MODEL'):
            current += 1; continue
        if line.startswith('ENDMDL'):
            if current >= selected:break
            continue
        if line.startswith(('ATOM  ','HETATM')) and (current in (0,selected)):
            autodock = line[77:79].strip().upper() if len(line) > 77 else ''
            element = AUTODOCK_ELEMENTS.get(autodock,autodock if autodock in VDW_RADII else '')
            element = element or element_of(line[12:16])
            lines.append(line[:66].ljust(76)+element.rjust(2))
    return '\n'.join(lines)+'\nEND\n'


def read_pose(path, *, model=1, template_smiles=None):
    """Molecule with 3D coordinates and, when possible, explicit bond orders."""
    from rdkit import Chem
    path = Path(path); suffix = path.suffix.lower(); molecule = None; bonds = 'explicit'
    if suffix == '.sdf':
        supplier = Chem.SDMolSupplier(str(path),sanitize=False,removeHs=False)
        molecule = next((m for m in supplier if m is not None),None)
    elif suffix == '.mol2':
        molecule = Chem.MolFromMol2File(str(path),sanitize=False,removeHs=False)
    elif suffix == '.pdbqt':
        molecule = Chem.MolFromPDBBlock(pdbqt_to_pdb(path.read_text(errors='replace'),model),
                                        sanitize=False,removeHs=False,proximityBonding=True)
        bonds = 'proximity'
    elif suffix == '.pdb':
        molecule = Chem.MolFromPDBFile(str(path),sanitize=False,removeHs=False,proximityBonding=True)
        bonds = 'proximity'
    if molecule is None or not molecule.GetNumAtoms() or not molecule.GetNumConformers():
        return None,bonds
    if bonds == 'proximity' and template_smiles:
        from rdkit.Chem import AllChem
        from .molecular_parsing import mol_from_smiles
        template = mol_from_smiles(template_smiles)
        if template is not None:
            # A docking pose carries only the polar hydrogens its engine
            # added; the heavy-atom skeleton is what matches the template.
            for candidate in (molecule,_without_hydrogens(molecule)):
                if candidate is None:
                    continue
                try:
                    molecule = AllChem.AssignBondOrdersFromTemplate(template,candidate)
                    bonds = 'template'
                    break
                except (ValueError,RuntimeError):
                    continue
    return molecule,bonds


def _without_hydrogens(molecule):
    from rdkit import Chem
    try:
        editable = Chem.RWMol(molecule)
        for index in sorted((atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1),reverse=True):
            editable.RemoveAtom(index)
        probe = editable.GetMol()
        return probe if probe.GetNumAtoms() and probe.GetNumConformers() else None
    except (ValueError,RuntimeError):
        return None


def _sanitized(molecule):
    from rdkit import Chem
    probe = Chem.Mol(molecule)
    try:
        Chem.SanitizeMol(probe)
    except (ValueError,RuntimeError) as error:
        return None,str(error)
    return probe,None


def _coordinates(molecule):
    conformer = molecule.GetConformer()
    return [(conformer.GetAtomPosition(index).x,conformer.GetAtomPosition(index).y,conformer.GetAtomPosition(index).z)
            for index in range(molecule.GetNumAtoms())]


def _relaxed(molecule, steps):
    """Relaxed copy of the same connectivity: MMFF94 with a UFF fallback."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    probe = Chem.Mol(molecule)
    try:
        if AllChem.MMFFHasAllMoleculeParams(probe):
            properties = AllChem.MMFFGetMoleculeProperties(probe)
            field = AllChem.MMFFGetMoleculeForceField(probe,properties)
            method = 'MMFF94'
        elif AllChem.UFFHasAllMoleculeParams(probe):
            field = AllChem.UFFGetMoleculeForceField(probe)
            method = 'UFF'
        else:
            return None,None,None,None
    except (ValueError,RuntimeError):
        return None,None,None,None
    try:
        energy = float(field.CalcEnergy())
        field.Minimize(maxIts=max(1,int(steps)))
        relaxed_energy = float(field.CalcEnergy())
    except (ValueError,RuntimeError):
        return None,None,None,None
    if not (math.isfinite(energy) and math.isfinite(relaxed_energy)):
        return None,None,None,None
    return probe,method,energy,relaxed_energy


def _internal_deviations(molecule, relaxed):
    """Largest relative bond-length and absolute bond-angle deviations."""
    from rdkit.Chem import rdMolTransforms as transforms
    pose = molecule.GetConformer(); reference = relaxed.GetConformer()
    bond = angle = 0.0
    for link in molecule.GetBonds():
        first,second = link.GetBeginAtomIdx(),link.GetEndAtomIdx()
        expected = transforms.GetBondLength(reference,first,second)
        measured = transforms.GetBondLength(pose,first,second)
        if expected > 0:
            bond = max(bond,abs(measured-expected)/expected)
    for atom in molecule.GetAtoms():
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        for index,first in enumerate(neighbours):
            for second in neighbours[index+1:]:
                try:
                    expected = transforms.GetAngleDeg(reference,first,atom.GetIdx(),second)
                    measured = transforms.GetAngleDeg(pose,first,atom.GetIdx(),second)
                except ValueError:
                    continue
                angle = max(angle,abs(measured-expected))
    return bond,angle


def _planarity(molecule):
    """Largest out-of-plane deviation of aromatic rings and double bonds."""
    import numpy
    points = numpy.asarray(_coordinates(molecule),dtype=float)
    groups = []
    for ring in molecule.GetRingInfo().AtomRings():
        if all(molecule.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            groups.append(list(ring))
    from rdkit import Chem
    for link in molecule.GetBonds():
        if link.GetBondType() != Chem.BondType.DOUBLE or link.IsInRing():
            continue
        first,second = link.GetBeginAtom(),link.GetEndAtom()
        if first.GetAtomicNum() == 1 or second.GetAtomicNum() == 1:
            continue
        members = {first.GetIdx(),second.GetIdx()}
        for atom in (first,second):
            members.update(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
        if len(members) >= 4:
            groups.append(sorted(members))
    worst = 0.0
    for group in groups:
        selected = points[group]
        centred = selected-selected.mean(axis=0)
        # The smallest singular direction is the best-fit plane normal.
        normal = numpy.linalg.svd(centred)[2][-1]
        worst = max(worst,float(numpy.abs(centred @ normal).max()))
    return (worst if groups else None),len(groups)


def _clash(molecule, receptor_atoms, ligand_elements=None):
    """Vectorised worst van der Waals contact between heavy atoms.

    Returns the smallest distance/radii-sum ratio (the decision metric), the
    absolute overlap of that pair, the closest heavy-atom distance and how many
    pairs are closer than the sum of their radii.
    """
    import numpy
    coordinates,radii,_ = receptor_atoms
    if not coordinates:
        return None,None,None,0
    points = numpy.asarray(_coordinates(molecule),dtype=float)
    elements = ligand_elements or [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    heavy = numpy.asarray([index for index,symbol in enumerate(elements) if symbol != 'H'],dtype=int)
    if not heavy.size:
        return None,None,None,0
    points = points[heavy]
    ligand_radii = numpy.asarray([radius_of(elements[index]) for index in heavy],dtype=float)
    target = numpy.asarray(coordinates,dtype=float)
    target_radii = numpy.asarray(radii,dtype=float)
    limit = float(ligand_radii.max()+target_radii.max())
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        cKDTree = None
    gap = math.inf; ratio = math.inf; distance = math.inf; contacts = 0
    if cKDTree is not None:
        tree = cKDTree(target)
        distance = float(numpy.asarray(tree.query(points)[0]).min())
        for index,neighbours in enumerate(tree.query_ball_point(points,limit)):
            if not neighbours:
                continue
            neighbours = numpy.asarray(neighbours,dtype=int)
            separation = numpy.linalg.norm(target[neighbours]-points[index],axis=1)
            total = target_radii[neighbours]+ligand_radii[index]
            overlap = separation-total
            gap = min(gap,float(overlap.min()))
            ratio = min(ratio,float((separation/total).min()))
            contacts += int((overlap < 0).sum())
        if not math.isfinite(gap):
            # No atom pair is even within the largest possible contact radius,
            # so these lower bounds are already permissive and sufficient.
            gap = distance-limit; ratio = distance/limit
    else:
        # Broadcasting keeps the computation matricial without SciPy installed.
        separation = numpy.linalg.norm(points[:,None,:]-target[None,:,:],axis=2)
        total = target_radii[None,:]+ligand_radii[:,None]
        overlap = separation-total
        gap = float(overlap.min()); ratio = float((separation/total).min())
        distance = float(separation.min()); contacts = int((overlap < 0).sum())
    if not math.isfinite(gap):
        return None,None,None,0
    return ratio,gap,distance,contacts


def evaluate_pose(record, **options):
    """Full criteria battery for one pose; never raises for scientific content."""
    settings = dict(TOLERANCES,**{k:v for k,v in options.items() if v is not None})
    result = {'pose_valid':False,'valid_chemistry':None,'valid_connectivity':None,'valid_geometry':None,
              'valid_planarity':None,'valid_energy':None,'valid_clash':None,'min_vdw_ratio':None,
              'min_vdw_overlap':None,'min_heavy_distance':None,'clashing_contacts':None,
              'delta_energy':None,'energy_method':None,
              'max_bond_deviation':None,'max_angle_deviation':None,'max_planarity_deviation':None,
              'planar_groups':None,'bond_orders':None,'atoms':None,'rejection_reasons':[],'error':''}
    pose_path = record.get('pose')
    if not pose_path or not Path(pose_path).is_file():
        result['rejection_reasons'] = ['unreadable']; result['error'] = 'Pose indisponível.'
        return result
    try:
        molecule,bonds = read_pose(pose_path,model=int(record.get('model') or 1),
                                   template_smiles=record.get('canonical_smiles'))
    except (OSError,ValueError,RuntimeError) as error:
        result['rejection_reasons'] = ['unreadable']; result['error'] = str(error)[:400]
        return result
    if molecule is None:
        result['rejection_reasons'] = ['unreadable']; result['error'] = 'A pose não contém coordenadas legíveis.'
        return result
    result['bond_orders'] = bonds; result['atoms'] = molecule.GetNumAtoms()
    reasons = []
    points = _coordinates(molecule)
    finite = all(all(map(math.isfinite,point)) for point in points)
    distinct = len(set(points)) == len(points)
    if not finite or not distinct or not points:
        reasons.append('coordinates')
    sanitized,error = _sanitized(molecule)
    result['valid_chemistry'] = sanitized is not None
    if sanitized is None:
        reasons.append('sanitization'); result['error'] = (error or '')[:400]
    from rdkit import Chem
    probe = sanitized if sanitized is not None else molecule
    result['valid_connectivity'] = len(Chem.GetMolFrags(probe)) == 1
    if not result['valid_connectivity']:
        reasons.append('fragments')
    if sanitized is not None and finite and bonds != 'proximity':
        # Without bond orders every bond would be read as single: the force
        # field and the planarity test would measure a different molecule.
        relaxed,method,energy,relaxed_energy = _relaxed(sanitized,settings['minimization_steps'])
        if relaxed is not None:
            result['energy_method'] = method
            result['delta_energy'] = round(energy-relaxed_energy,4)
            result['valid_energy'] = result['delta_energy'] <= settings['delta_energy_limit']
            if not result['valid_energy']:
                reasons.append('internal_energy')
            bond,angle = _internal_deviations(sanitized,relaxed)
            result['max_bond_deviation'] = round(bond,4); result['max_angle_deviation'] = round(angle,3)
            result['valid_geometry'] = bond <= settings['bond_tolerance'] and angle <= settings['angle_tolerance']
            if bond > settings['bond_tolerance']:
                reasons.append('bond_lengths')
            if angle > settings['angle_tolerance']:
                reasons.append('bond_angles')
        deviation,groups = _planarity(sanitized)
        result['planar_groups'] = groups
        if deviation is not None:
            result['max_planarity_deviation'] = round(deviation,4)
            result['valid_planarity'] = deviation <= settings['planarity_tolerance']
            if not result['valid_planarity']:
                reasons.append('planarity')
    receptor = record.get('receptor')
    if receptor and Path(receptor).is_file():
        try:
            excluded = record.get('exclude_residues') or [record.get('ligand_id') or '']
            atoms = structure_atoms(receptor,exclude_residues=excluded)
            ratio,gap,distance,contacts = _clash(molecule,atoms)
        except (OSError,ValueError,MemoryError) as error:
            ratio = gap = distance = contacts = None
            result['error'] = result['error'] or str(error)[:400]
        if ratio is not None:
            result['min_vdw_ratio'] = round(ratio,4); result['min_vdw_overlap'] = round(gap,4)
            result['min_heavy_distance'] = round(distance,4); result['clashing_contacts'] = contacts
            result['valid_clash'] = ratio >= settings['vdw_clash_scale']
            if not result['valid_clash']:
                reasons.append('steric_clash')
    result['rejection_reasons'] = reasons
    result['pose_valid'] = not reasons
    return result


def describe(reasons):
    return '; '.join(REASONS.get(reason,reason) for reason in reasons)
