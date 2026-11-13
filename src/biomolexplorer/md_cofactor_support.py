"""CCD-backed GAFF parameterization of isolated organic target components.

Covalent PTMs are handled by complete protein profiles, never by assigning
small-molecule charges to a detached protein fragment. Metal coordination and
unknown chemistry fail explicitly instead of receiving invented parameters.
"""
import math
from pathlib import Path
import re
import shutil

from .calculation_cache import CalculationCache,identity
from .stage_cache import file_digest
from .storage import write_json


def native_residues(root):
    names=set()
    for name in ('leaprc.dock.cofactors','leaprc.dock.ions'):
        path=Path(root)/'parameters/leap/cmd'/name
        if path.is_file():names.update(re.findall(r'^\s*([A-Za-z0-9]+)\s*=\s*(?:loadMol2|loadOff)',path.read_text(),re.M|re.I))
    names.update('LI NA NA+ K RB CS CL CL- MG CA ZN SR BA'.split())
    return names


def parameterize(directory,root,residues,report):
    """Return local LEaP commands; a validated checkpoint avoids all tools."""
    from rdkit import Chem
    from .md_residue_support import pdb_residues,audit_ligand
    from .structure_sanitation import STANDARD
    from .structure_selection import component_path,chemical_template
    from .preparation_tools import resolve_amber_tool,run_amber_tool
    directory=Path(directory)
    previous=report.get('dynamic_cofactors')
    if previous is not None and all(file_digest(directory/name)==sha for item in previous for name,sha in item['files'].items()):
        return [command for item in previous for command in (item['component']+' = loadMol2 '+item['component']+'.mol2','loadAmberParams '+item['component']+'.frcmod')]
    supported=native_residues(root)|STANDARD|set(report.get('modified_residues',{}).values())
    # Native nucleic-acid and terminal library names stay on the original path.
    supported.update('A C G U T DA DC DG DT RA RC RG RU A3 A5 C3 C5 G3 G5 U3 U5 DA3 DA5 DC3 DC5 DG3 DG5 DT3 DT5 Na+ Cl-'.split())
    unknown={key:value for key,value in residues.items() if value['name'] not in supported}
    commands=[];metadata=[];roles={row['residue']:row['role'] for row in report.get('structural_sanitation',{}).get('residue_roles',[])}
    for key,value in unknown.items():
        code=value['name'];atoms=value['atoms'];heavy={name:atom for name,atom in atoms.items() if atom['element'] not in ('H','D')}
        if {'N','CA','C'}<=set(heavy) or roles.get(key)=='covalently_connected':
            raise ValueError('Unparameterized covalent/polymer component '+code+' '+key.strip()+': provide a validated bonded protein profile; detached GAFF charges would be invalid.')
        if any(atom['element'] not in ('C','N','O','S','P','F','CL','BR','I','H') for atom in atoms.values()):
            raise ValueError('Component '+code+' requires a dedicated metal/element parameter model')
        definition=component_path(code,directory,{})
        template,stereo=chemical_template(definition,code);template=Chem.RemoveHs(template)
        names=[atom.GetProp('_TriposAtomName') if atom.HasProp('_TriposAtomName') else '' for atom in template.GetAtoms()]
        if len(set(names))!=len(names) or set(names)!=set(heavy):raise ValueError('CCD/target atom inventory mismatch for '+code)
        conformer=Chem.Conformer(len(names));template.RemoveAllConformers()
        for index,atom in enumerate(template.GetAtoms()):
            observed=heavy[names[index]]
            if atom.GetSymbol().upper()!=observed['element']:raise ValueError('CCD/target element mismatch for '+code)
            conformer.SetAtomPosition(index,observed['xyz'])
        template.AddConformer(conformer);Chem.AssignStereochemistryFrom3D(template);Chem.AssignStereochemistry(template,cleanIt=True,force=True)
        for atom,name in zip(template.GetAtoms(),names):
            if name in stereo and (not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=stereo[name]):raise ValueError('CCD/target stereochemistry mismatch for '+code+' '+name)
        if any(a.GetNumRadicalElectrons() for a in template.GetAtoms()):raise ValueError('Open-shell cofactor requires dedicated parameters: '+code)
        molecule=Chem.AddHs(template,addCoords=True);charge=Chem.GetFormalCharge(molecule)
        observed_hydrogens=sum(atom['element'] in ('H','D') for atom in atoms.values())
        if observed_hydrogens and observed_hydrogens!=molecule.GetNumAtoms()-len(names):raise ValueError('Observed cofactor protonation differs from the CCD state: '+code)
        declared=[atom['formal_charge'] for atom in atoms.values() if atom.get('formal_charge') is not None]
        if declared and sum(declared)!=charge:raise ValueError('Observed cofactor formal charge differs from the CCD state: '+code)
        chemical_identity=Chem.MolToSmiles(molecule)
        previous_component=next((item for item in metadata if item['component']==code),None)
        if previous_component:
            if previous_component['chemical_identity']!=chemical_identity:raise ValueError('Conflicting chemical states for the same component: '+code)
            previous_component['instances'].append(key)
            continue
        tools={name:resolve_amber_tool(name) for name in ('antechamber','parmchk2')}
        recipe=dict(schema=1,definition=file_digest(definition),charge=charge,force_field='gaff',
            smiles=Chem.MolToSmiles(molecule),coordinates=[heavy[name]['xyz'] for name in names],
            tools={name:file_digest(tool.executable) for name,tool in tools.items()})
        cache=CalculationCache('amber-cofactor',directory);cache_key=identity(recipe)
        folder=directory/'cofactor_parameters'/cache_key;folder.mkdir(parents=True,exist_ok=True)
        mol2=folder/(code+'.mol2');frcmod=folder/(code+'.frcmod');files={'mol2':mol2,'frcmod':frcmod}
        with cache.locked(cache_key):
            if not cache.restore(cache_key,files):
                sdf=folder/'component.sdf'
                with Chem.SDWriter(str(sdf)) as writer:writer.write(molecule)
                run_amber_tool('antechamber',['-i',sdf.name,'-fi','sdf','-o',mol2.name,'-fo','mol2','-at','gaff','-c','bcc','-nc',str(charge),'-rn',code,'-s','2'],cwd=folder,timeout=600)
                lines=mol2.read_text().splitlines();section=None;output=[];index=0
                for line in lines:
                    if line.startswith('@<TRIPOS>'):section=line
                    elif section=='@<TRIPOS>ATOM' and line.strip():
                        fields=line.split();fields[1]=names[index] if index<len(names) else 'H'+str(index-len(names)+1);fields[7]=code;line=' '.join(fields);index+=1
                    output.append(line)
                mol2.write_text('\n'.join(output)+'\n');audit_ligand(mol2)
                atom_lines=mol2.read_text().split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].splitlines()
                rows=[line.split() for line in atom_lines if line.strip()]
                if len(rows)!=molecule.GetNumAtoms() or abs(sum(float(row[8]) for row in rows)-charge)>.02:
                    raise ValueError('Generated cofactor charges or atom count failed validation: '+code)
                if any(math.dist(tuple(float(v) for v in rows[i][2:5]),heavy[names[i]]['xyz'])>.002 for i in range(len(names))):
                    raise ValueError('Generated cofactor atom order or coordinates changed: '+code)
                if any(rows[i][5].split('.')[0].lower()!=molecule.GetAtomWithIdx(i).GetSymbol().lower() and not rows[i][5].lower().startswith(molecule.GetAtomWithIdx(i).GetSymbol().lower()) for i in range(len(names))):
                    raise ValueError('Generated cofactor atom order differs from CCD: '+code)
                run_amber_tool('parmchk2',['-i',mol2.name,'-f','mol2','-o',frcmod.name,'-s','gaff'],cwd=folder,timeout=120)
                if not frcmod.is_file() or re.search('ATTN|need revision',frcmod.read_text(),re.I):raise ValueError('Unresolved cofactor force-field terms: '+code)
                cache.remember(cache_key,files,dict(recipe=recipe))
        for path in (mol2,frcmod):shutil.copy2(path,directory/path.name)
        commands.extend([code+' = loadMol2 '+mol2.name,'loadAmberParams '+frcmod.name])
        metadata.append(dict(residue=key,instances=[key],chemical_identity=chemical_identity,component=code,charge=charge,force_field='GAFF',charge_method='AM1-BCC',
            definition_sha256=file_digest(definition),definition_source='wwPDB CCD or explicit local chemical definition',
            definition_license='Local chemical input; license supplied by project owner' if definition.parent==directory or definition.parent==directory/'chemistry' else 'wwPDB CCD: CC0; see https://www.wwpdb.org/about/usage-policies',cache_key=cache_key,
            files={path.name:file_digest(path) for path in (mol2,frcmod)},tools=recipe['tools']))
    report['dynamic_cofactors']=metadata
    return commands
