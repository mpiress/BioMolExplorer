"""Pose selection and geometric residue contacts in the original coordinate frame."""
import math
import re
from collections import defaultdict

PROTEIN=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP CYX SEC PYL'.split())


def best_pose(text, fmt):
    if fmt=='pdbqt' and re.search(r'^MODEL\s',text,re.M):
        blocks=re.findall(r'^MODEL\s.*?(?:^ENDMDL.*?$|\Z)',text,re.M|re.S)
    elif fmt=='mol2':
        blocks=re.split(r'(?=^##########\s+Name:)',text,flags=re.M)
        blocks=[b for b in blocks if '@<TRIPOS>MOLECULE' in b]
        if len(blocks)<=1:
            blocks=[b for b in re.split(r'(?=^@<TRIPOS>MOLECULE)',text,flags=re.M) if '@<TRIPOS>MOLECULE' in b]
            if len(blocks)==1:blocks=[text]
    else:blocks=[text]
    if not blocks:raise ValueError('A conformação selecionada não contém átomos.')
    label='REMARK VINA RESULT' if fmt=='pdbqt' else 'Grid_Score'
    scored=[]
    for index,block in enumerate(blocks):
        match=re.search(re.escape(label)+r':\s*([-+\d.eE]+)',block)
        value=float(match[1]) if match else None
        if value is not None and math.isfinite(value):scored.append((value,index,block))
    if scored:
        value,index,block=min(scored,key=lambda item:item[0])
        return block,dict(score=value,model=index+1,selection='lowest_score')
    return blocks[0],dict(score=None,model=1,selection='first_unscored')


def atoms(text,fmt,include_hydrogens=False):
    result=[];active=False
    for line in text.splitlines():
        if fmt=='mol2':
            if line.startswith('@<TRIPOS>'):active=line=='@<TRIPOS>ATOM';continue
            if not active or not line.strip():continue
            fields=line.split()
            if len(fields)<8:continue
            name=fields[7];match=re.fullmatch(r'([A-Za-z]+)(-?\d+)?',name)
            resn=match[1].upper() if match else name
            resi=int(match[2]) if match and match[2] else int(fields[6])
            xyz=[float(v) for v in fields[2:5]];element=fields[5].split('.')[0].upper()
            record=dict(atom=fields[1],resn=resn,resi=resi,chain='',icode='',elem=element,xyz=xyz)
        else:
            if line.startswith('ENDMDL'):break
            if not line.startswith(('ATOM  ','HETATM')):continue
            if line[16:17] not in (' ','A',''):continue
            xyz=[float(line[start:start+8]) for start in (30,38,46)]
            element=line[76:78].strip().upper()
            if fmt=='pdbqt':element={'A':'C','OA':'O','NA':'N','SA':'S','HD':'H','HS':'H'}.get(line.split()[-1],line.split()[-1]).upper()
            if not element:element=re.sub('[^A-Za-z]','',line[12:16]).strip()[:1].upper()
            record=dict(atom=line[12:16].strip(),resn=line[17:20].strip(),resi=int(line[22:26]),
                chain=line[21:22].strip(),icode=line[26:27].strip(),elem=element,xyz=xyz)
        if not all(math.isfinite(v) for v in xyz):raise ValueError('A estrutura contém coordenadas inválidas.')
        if include_hydrogens or record['elem'] not in ('H','D'):result.append(record)
    if not result:raise ValueError('A conformação selecionada não contém átomos.')
    return result


def contacts(receptor,ligand,cutoff=4.):
    """Minimum heavy-atom distance per protein residue, using spatial bins."""
    cutoff=float(cutoff)
    if not math.isfinite(cutoff) or not 2<=cutoff<=8:raise ValueError('Informe uma distância de contato entre 2 e 8 Å.')
    bins=defaultdict(list)
    for atom in ligand:bins[tuple(math.floor(v/cutoff) for v in atom['xyz'])].append(atom)
    nearby={}
    for atom in receptor:
        if atom['resn'] not in PROTEIN:continue
        cell=tuple(math.floor(v/cutoff) for v in atom['xyz'])
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for other in bins[(cell[0]+dx,cell[1]+dy,cell[2]+dz)]:
                        distance=math.dist(atom['xyz'],other['xyz'])
                        if distance>cutoff:continue
                        key=(atom['chain'],atom['resi'],atom['icode'],atom['resn'])
                        if key not in nearby or distance<nearby[key]['distance']:
                            nearby[key]={k:atom[k] for k in ('chain','resi','icode','resn')}
                            nearby[key].update(distance=distance,receptor_atom=atom['atom'],ligand_atom=other['atom'])
    result=sorted(nearby.values(),key=lambda r:(r['distance'],r['chain'],r['resi'],r['icode']))
    for row in result:row['distance']=round(row['distance'],3)
    return result


def scene_payload(store,token,pid,spec,cutoff=4.,include_interactions=False):
    layers=[];total=0
    from .visualizations import MAX_VIEW_BYTES
    for entry in spec['layers']:
        data=store.read_file(token,pid,entry['path'],MAX_VIEW_BYTES+1)
        total+=len(data)
        if total>MAX_VIEW_BYTES:raise ValueError('Visualização maior que o limite de 32 MB. Baixe o arquivo para consultá-lo localmente.')
        text=data.decode('utf-8');fmt=entry['format'];info={}
        if entry['role']=='pose':text,info=best_pose(text,fmt)
        if entry.get('exclude_residue'):
            name,number,chain=entry['exclude_residue']
            text='\n'.join(line for line in text.splitlines() if not (line.startswith(('ATOM  ','HETATM')) and
                line[17:20].strip()==name and int(line[22:26])==int(number) and line[21:22].strip()==chain))+'\n'
        if entry.get('residue'):
            name,number,chain=entry['residue']
            text='\n'.join(line for line in text.splitlines() if line.startswith(('ATOM  ','HETATM')) and
                line[17:20].strip()==name and int(line[22:26])==int(number) and line[21:22].strip()==chain)+'\n'
        layers.append(dict(role=entry['role'],format=fmt,data=text,**info))
    receptor=next((layer for layer in layers if layer['role']=='receptor'),None)
    pose=next(layer for layer in layers if layer['role']=='pose')
    if receptor is None:raise ValueError('O receptor associado a esta pose não está disponível.')
    protein=atoms(receptor['data'],receptor['format'])
    result=contacts(protein,atoms(pose['data'],pose['format']),cutoff)
    reference=next((layer for layer in layers if layer['role']=='reference'),None)
    complex_layer=next((layer for layer in layers if layer['role']=='complex'),None)
    native_protein=atoms(complex_layer['data'],complex_layer['format']) if complex_layer else protein
    native=contacts(native_protein,atoms(reference['data'],reference['format']),cutoff) if reference else []
    payload=dict(name=spec['name'],layers=layers,contacts=result,reference_contacts=native,cutoff=float(cutoff),
        reference_kind=spec.get('reference_kind'),coordinate_frame='unchanged')
    if include_interactions:
        from .docking_interactions import interaction_diagram
        payload['interactions']=interaction_diagram(receptor,pose,spec.get('ligand_smiles'))
    return payload
