"""Target-independent, integrity-checked parameter profiles for native AMBER."""
import json
import math
import os
import re
from pathlib import Path
import shutil

from .stage_cache import file_digest
from .storage import write_json,write_text

DATA=Path(__file__).with_name('resources')/'amber_support/ffptm_ff03'
POLICIES=('auto_curated','native_ff94','ff03_ffptm')


def profile_data():
    """Load bundled and optional project-owned profiles with provenance."""
    roots=[DATA]
    configured=os.environ.get('BIOMOL_AMBER_PTM_PROFILES','')
    if configured:roots.append(Path(configured).expanduser().resolve())
    manifests=[]
    for root in roots:
        manifest=json.loads((root/'manifest.json').read_text())
        if manifest.get('protein_force_field')!='ff03' or not manifest.get('residues'):
            raise ValueError('Protein profiles must declare a coherent ff03 force field')
        source=manifest.get('source',{})
        if not source.get('license') or not source.get('url'):raise ValueError('Parameter origin and license are mandatory')
        for name,sha in manifest['files'].items():
            if Path(name).name!=name or file_digest(root/name)!=sha:
                raise ValueError('PTM parameter integrity check failed: '+name)
        for code,model in manifest['residues'].items():
            if not re.fullmatch('[A-Z0-9]{1,3}',code):raise ValueError('Invalid parameter residue name')
            if not all(code+suffix in manifest['files'] for suffix in ('.prepin','.frcmod')):raise ValueError('Incomplete parameter files: '+code)
            if set(model['heavy_atoms'])!={name for name,atom in model['atoms'].items() if atom['element'] not in ('H','D')}:raise ValueError('Inconsistent parameter atom graph: '+code)
            charges=[atom['charge'] for atom in model['atoms'].values()]
            if not all(map(math.isfinite,charges)) or abs(sum(charges)-model['net_charge'])>.001:raise ValueError('Invalid parameter charge: '+code)
        manifests.append((root,manifest))
    return manifests


def bundle():
    manifests=profile_data();manifest=json.loads(json.dumps(manifests[0][1]))
    manifest['origins']=[item['source'] for _,item in manifests]
    manifest['file_sources']={name:str(root/name) for root,item in manifests for name in item['files']}
    for _,extra in manifests[1:]:
        if set(manifest['residues'])&set(extra['residues']):raise ValueError('Custom profiles cannot override validated bundled residues')
        manifest['residues'].update(extra['residues']);manifest['files'].update(extra['files'])
        manifest['id']+='+'+extra['id']
    return manifest


def profile_identity():
    return [{key:manifest[key] for key in ('id','files','source','residues')} for _,manifest in profile_data()]


def pdb_residues(path):
    residues={};serials={};bonds=set()
    for line in Path(path).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            key=line[21:27];name=line[12:16].strip()
            item=residues.setdefault(key,dict(name=line[17:20].strip(),atoms={}))
            if item['name']!=line[17:20].strip() or name in item['atoms']:
                raise ValueError('Duplicate PDB atom or unresolved alternate conformation: '+key.strip()+' '+name)
            xyz=tuple(float(line[i:i+8]) for i in (30,38,46))
            if not all(map(math.isfinite,xyz)):raise ValueError('Nonfinite receptor coordinates')
            element=line[76:78].strip().upper() or name.lstrip('0123456789')[0].upper()
            charge_text=line[78:80].strip()
            formal_charge=int(charge_text[0])*(1 if charge_text[-1]=='+' else -1) if len(charge_text)==2 and charge_text[0].isdigit() and charge_text[-1] in '+-' else None
            item['atoms'][name]=dict(xyz=xyz,element=element,formal_charge=formal_charge)
            serial=int(line[6:11])
            if serial in serials:raise ValueError('Duplicate PDB atom serial')
            serials[serial]=(key,name)
        elif line.startswith('CONECT'):
            ids=[int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()]
            if ids:bonds.update(tuple(sorted((ids[0],other))) for other in ids[1:] if other!=ids[0])
    return residues,serials,bonds


def plan(path,policy='auto_curated'):
    if policy not in POLICIES:raise ValueError('Unknown residue parameter policy')
    residues,serials,bonds=pdb_residues(path);manifest=bundle();registry=manifest['residues']
    modified={key:value for key,value in residues.items() if value['name'] in registry}
    selected='ff03_ffptm' if policy=='ff03_ffptm' or policy=='auto_curated' and modified else 'native_ff94'
    report=dict(schema_version=2,policy=policy,profile=selected,protein_force_field='ff03' if selected=='ff03_ffptm' else 'ff94',
        source_sha256=file_digest(path),modified_residues={key:value['name'] for key,value in modified.items()},
        rebuilt_heavy_atoms={},chemical_state_policy='Published template states only; no implicit pH, tautomer or formal-charge inference')
    if selected=='native_ff94':return report
    ordered=list(residues)
    for key,value in modified.items():
        code=value['name'];model=registry[code];atoms=value['atoms'];heavy={name:atom for name,atom in atoms.items() if atom['element'] not in ('H','D')}
        expected=model['heavy_atoms']
        if set(heavy)-set(expected) or not {'N','CA','C','O'}<=set(heavy) or any(atom['element']!=expected[name] for name,atom in heavy.items()):
            raise ValueError(code+' heavy-atom graph differs from the published parameter model: '+key.strip())
        declared=[atom['formal_charge'] for atom in atoms.values() if atom.get('formal_charge') is not None]
        if declared and sum(declared)!=model['net_charge']:raise ValueError('Explicit PDB charge conflicts with the published PTM state: '+code)
        report['rebuilt_heavy_atoms'][key]=sorted(set(expected)-set(heavy))
        if set(atoms)-set(model['atoms']):raise ValueError('Unsupported atom naming or protonation state: '+code+' '+key.strip())
        for first,second in model['bonds']:
            if first in heavy and second in heavy and not .8<math.dist(atoms[first]['xyz'],atoms[second]['xyz'])<2.5:
                raise ValueError('Invalid covalent geometry: '+code+' '+key.strip()+' '+first+'-'+second)
        index=ordered.index(key)
        before=residues[ordered[index-1]]['atoms'] if index else {}
        after=residues[ordered[index+1]]['atoms'] if index+1<len(ordered) else {}
        if not index or index+1==len(ordered) or key[0]!=ordered[index-1][0] or key[0]!=ordered[index+1][0] or \
                'C' not in before or 'N' not in after or math.dist(before['C']['xyz'],atoms['N']['xyz'])>2.0 or math.dist(atoms['C']['xyz'],after['N']['xyz'])>2.0:
            raise ValueError('Modified terminal or disconnected residue requires terminal parameters: '+code+' '+key.strip())
    for first,second in bonds:
        a=serials.get(first);b=serials.get(second)
        if not a or not b:continue
        for local,other in ((a,b),(b,a)):
            if local[0] in modified and other[0]!=local[0] and (local[1],other[1]) not in (('N','C'),('C','N')):
                raise ValueError('Crosslinked modified residue requires a bonded model: '+modified[local[0]]+' '+local[0].strip())
    codes=sorted({value['name'] for value in modified.values()})
    files={name:sha for name,sha in manifest['files'].items() if name.split('.')[0] in codes}
    report.update(bundle_id=manifest['id'],bundle_files=files,citation=manifest['source'],parameter_origins=manifest['origins'],parameter_name_adaptations=manifest.get('adaptations',[]))
    return report


def install(directory,report,root):
    directory=Path(directory);root=Path(root)
    if report['profile']=='ff03_ffptm':
        for name in ('cmd/leaprc.ff03','lib/all_amino03.lib','parm/frcmod.ff03'):
            if not (root/'parameters/leap'/name).is_file():raise ValueError('DOCK6 ff03 library is unavailable: '+name)
        manifest=bundle()
        for name in report['bundle_files']:shutil.copy2(manifest['file_sources'][name],directory/name)
        commands=[]
        for code in sorted(set(report['modified_residues'].values())):
            commands.extend(['loadAmberPrep '+code+'.prepin','loadAmberParams '+code+'.frcmod'])
    else:commands=[]
    from .md_cofactor_support import parameterize
    residues,_,_=pdb_residues(directory/'rec.pdb')
    commands.extend(parameterize(directory,root,residues,report))
    write_text('\n'.join(commands)+'\n',directory/'biomol_residue_support.leap')
    write_json(report,directory/'residue_support.json')


def audit_ligand(path):
    """Validate MOL2 connectivity before invoking any charge/parameter tool."""
    from .md_inputs import validate_structure
    path=validate_structure(path,'ligand');text=path.read_text().splitlines()
    counts=text[text.index('@<TRIPOS>MOLECULE')+2].split();atoms={};bonds={};section=None
    for line in text:
        if line.startswith('@<TRIPOS>'):section=line;continue
        if not line.strip() or line.startswith('#'):continue
        row=line.split()
        if section=='@<TRIPOS>ATOM':
            index=int(row[0])
            if index in atoms:raise ValueError('Duplicate MOL2 atom identifier')
            if row[5].lower() in ('du','lp','any','x') or not re.fullmatch(r'[A-Za-z][A-Za-z0-9]*(?:\.[A-Za-z0-9]+)?',row[5]):
                raise ValueError('Unresolved MOL2 atom type: '+row[5])
            atoms[index]=row
        elif section=='@<TRIPOS>BOND':
            index,a,b=int(row[0]),int(row[1]),int(row[2])
            if index in bonds or a==b or row[3] not in ('1','2','3','ar','am'):raise ValueError('Invalid or unresolved MOL2 bond')
            bonds[index]=(a,b)
    if len(bonds)!=int(counts[1]) or any(a not in atoms or b not in atoms for a,b in bonds.values()):
        raise ValueError('MOL2 bond count or atom references are inconsistent')
    if len({tuple(sorted(edge)) for edge in bonds.values()})!=len(bonds):raise ValueError('Duplicate MOL2 bond')
    return dict(atoms=len(atoms),bonds=len(bonds),coordinate_and_charge_finite=True)


def validate_prepared(directory,report):
    """Verify all retained PTM graphs, published charges and original coordinates."""
    from .md_analysis import topology
    directory=Path(directory);validated=[]
    modified=report.get('modified_residues',{})
    if report['profile']!='ff03_ffptm' or not modified:return validated
    registry=bundle()['residues'];codes=set(modified.values())
    before,_,_=pdb_residues(directory/'rec.pdb');after,_,_=pdb_residues(directory/'rec.amber.pdb')
    old=[value for value in before.values() if value['name'] in codes];new=[value for value in after.values() if value['name'] in codes]
    if len(old)!=len(new):raise ValueError('Modified residue was lost during AMBER preparation')
    for a,b in zip(old,new):
        if a['name']!=b['name'] or set(a['atoms'])-set(b['atoms']):raise ValueError('Modified input atoms were lost during AMBER preparation')
        if any(math.dist(a['atoms'][name]['xyz'],b['atoms'][name]['xyz'])>.002 for name in registry[a['name']]['heavy_atoms'] if name in a['atoms']):
            raise ValueError('Modified heavy-atom coordinates changed during parameterization')
    for prefix in ('rec','rec.lig.1'):
        info=topology(directory/(prefix+'.prmtop'));indices=[i for i,label in enumerate(info['labels']) if label in codes]
        if [info['labels'][i] for i in indices]!=[value['name'] for value in old]:raise ValueError('Modified residue sequence changed in the AMBER topology')
        for index in indices:
            code=info['labels'][index];model=registry[code];selected={i for i,residue in enumerate(info['residue']) if residue==index}
            names={info['names'][i] for i in selected};charge=float(info['charges'][list(selected)].sum())
            if set(model['atoms'])!=names or not math.isfinite(charge) or abs(charge-model['net_charge'])>.001:
                raise ValueError('Modified topology atom graph or charge is invalid: '+code)
            if len(info['atom_types'])!=info['natom'] or any(info['atom_types'][i]!=model['atoms'][info['names'][i]]['type'] for i in selected):
                raise ValueError('Modified residue AMBER atom types differ from the validated model: '+code)
            edges={frozenset((info['names'][a],info['names'][b])) for a,b in info['bonds'] if a in selected and b in selected}
            if not {frozenset(edge) for edge in model['bonds']}<=edges:raise ValueError('Modified residue bond was lost in the AMBER topology: '+code)
            validated.append(dict(topology=prefix,residue=code,residue_index=index+1,charge=charge,bonds_validated=True))
    return validated
