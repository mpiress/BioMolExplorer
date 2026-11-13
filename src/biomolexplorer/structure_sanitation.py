"""Conservative structural sanitation; coordinate files are never edited in place.

A record tag is not a biological role. Covalent polymer residues and named
cofactors survive cleanup; only recognized, unbound additives are removed.
Unknown chemistry is retained for the parameter validator to diagnose.
"""
from collections import defaultdict
import json
import math
from pathlib import Path
import re

from .stage_cache import file_digest
from .storage import write_json

STANDARD=set('ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY HIS HID HIE HIP ILE LEU LYS LYN MET PHE PRO SER THR TRP TYR VAL ACE NME'.split())
COFACTORS=set('ADP ATP AMP FAD FMN NAD NAP NDP NADP HEM HEC HEA HEB H4B PLP PMP TPP COA SAM SAH FES SF4 F3S B12 CBL PQQ MGD MG MN ZN FE CU CA CO NI MO W SE'.split())
ADDITIVES=set('HOH WAT DOD GOL EDO PEG PGE PG4 PG5 PG6 MPD DMS EOH IPA ACT ACE SO4 PO4 CL BR IOD NH4 TRS MES HEP BME'.split())


def sanitize(source,destination,*,keep=(),remove_additives=True):
    source=Path(source);destination=Path(destination)
    if source.resolve()==destination.resolve():raise ValueError('Structural sanitation requires an execution-local copy')
    lines=source.read_text().splitlines();groups=defaultdict(list);models=0
    for index,line in enumerate(lines):
        if line.startswith('MODEL'):models+=1
        if line.startswith(('ATOM  ','HETATM')):
            if len(line)<54:raise ValueError('Truncated PDB atom record')
            xyz=[float(line[i:i+8]) for i in (30,38,46)]
            if not all(map(math.isfinite,xyz)):raise ValueError('Nonfinite structural coordinates')
            groups[line[21:27]].append((index,line))
    if models>1:raise ValueError('Select one PDB model before receptor preparation')
    selected={};alternates=[];renamed=[]
    for key,records in groups.items():
        occupancies=defaultdict(float)
        for _,line in records:
            if line[16].strip():occupancies[line[16]]+=float(line[54:60].strip() or 0)
        chosen=min(occupancies,key=lambda alt:(-occupancies[alt],alt)) if occupancies else ''
        if chosen:alternates.append(dict(residue=key,selected=chosen,occupancy_sums=dict(occupancies)))
        names=set()
        for index,line in records:
            if line[16].strip() and line[16]!=chosen:continue
            name=line[12:16].strip()
            residue=line[17:20].strip()
            alias={'OT1':'O','OT2':'OXT','HN':'H'}.get(name) if residue in STANDARD else None
            if residue=='ILE' and name=='CD':alias='CD1'
            if alias:
                renamed.append(dict(residue=key,input=name,output=alias,evidence='Unambiguous standard-protein atom alias'))
                line=line[:12]+f'{alias:>4s}'+line[16:];name=alias
            if name in names:raise ValueError('Duplicate atom in selected alternate: '+key+' '+name)
            names.add(name);selected[index]=line[:16]+' '+line[17:]
    residue_atoms=defaultdict(dict);serials={}
    for index,line in selected.items():
        serial=int(line[6:11]);key=line[21:27];name=line[12:16].strip()
        if serial in serials:raise ValueError('Duplicate PDB serial')
        serials[serial]=key;residue_atoms[key][name]=tuple(float(line[i:i+8]) for i in (30,38,46))
    covalent=set();ordered=list(residue_atoms)
    for a,b in zip(ordered,ordered[1:]):
        if a[0]==b[0] and 'C' in residue_atoms[a] and 'N' in residue_atoms[b] and math.dist(residue_atoms[a]['C'],residue_atoms[b]['N'])<1.9:
            covalent.update((a,b))
    for line in lines:
        if line.startswith('CONECT'):
            ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            if ids and ids[0] in serials:
                for other in ids[1:]:
                    if other in serials and serials[other]!=serials[ids[0]]:covalent.update((serials[ids[0]],serials[other]))
        elif line.startswith('LINK'):
            a=line[21:27];b=line[51:57]
            if a in residue_atoms and b in residue_atoms:covalent.update((a,b))
    keep=set(keep);removed=set();roles=[]
    for key,records in groups.items():
        name=records[0][1][17:20].strip();hetero=all(line.startswith('HETATM') for _,line in records)
        if key in covalent:role='covalently_connected'
        elif name in keep:role='explicit_cofactor'
        elif name in COFACTORS:role='cofactor_or_metal'
        elif hetero and name in ADDITIVES:role='unbound_additive'
        elif hetero and name in STANDARD:role='isolated_amino_acid'
        elif name in STANDARD:role='polymer'
        else:role='retained_unclassified'
        if remove_additives and role in ('unbound_additive','isolated_amino_acid'):removed.add(key)
        roles.append(dict(residue=key,name=name,role=role,action='removed' if key in removed else 'retained'))
    retained_serials={serial for serial,key in serials.items() if key not in removed};output=[]
    for index,line in enumerate(lines):
        if line.startswith(('ATOM  ','HETATM')):
            if index not in selected or line[21:27] in removed:continue
            line=selected[index]
        elif line.startswith(('ANISOU','SIGATM','SIGUIJ')):continue
        elif line.startswith('CONECT'):
            ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            if not ids or ids[0] not in retained_serials:continue
            ids=[ids[0]]+[v for v in ids[1:] if v in retained_serials]
            if len(ids)<2:continue
            line='CONECT'+''.join(f'{v:5d}' for v in ids)
        output.append(line)
    destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text('\n'.join(output)+'\n')
    report=dict(schema_version=1,source_sha256=file_digest(source),sanitized_sha256=file_digest(destination),
        alternate_conformations=alternates,atom_names_normalized=renamed,residue_roles=roles,coordinates_preserved=True,
        policy='Only recognized unbound additives and isolated amino-acid HETATMs removed; unknowns retained')
    write_json(report,destination.with_suffix('.sanitation.json'))
    return report


def receptor_script(text,script_path,record=None):
    """Apply the same role policy before Chimera's broad protein selection."""
    match=re.search(r'^open (.+\.pdb)$',text,re.M)
    if not match or 'select protein' not in text:return text
    if not Path(match[1]).is_file():return text
    selection=next(line for line in text.splitlines() if line.startswith('select protein'))
    keep=re.findall(r':([A-Z0-9]{1,5})',selection)
    path=Path(script_path).with_suffix('.sanitized.pdb')
    if Path(match[1]).resolve()==path.resolve():return text
    report=sanitize(match[1],path,keep=keep,remove_additives='| solvent' not in selection)
    if record:
        excluded={(str(record[3]),int(record[2]),str(record[1]))}
        raw=path.read_text().splitlines();drop={int(line[6:11]) for line in raw if line.startswith(('ATOM  ','HETATM')) and (line[21],int(line[22:26]),line[17:20].strip()) in excluded}
        if drop:
            filtered=[]
            for line in raw:
                if line.startswith(('ATOM  ','HETATM')) and int(line[6:11]) in drop:continue
                if line.startswith('CONECT'):
                    ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
                    if ids and any(value in drop for value in ids):continue
                filtered.append(line)
            path.write_text('\n'.join(filtered)+'\n')
            for row in report['residue_roles']:
                if (row['residue'][0],int(row['residue'][1:5]),row['name']) in excluded:
                    if row['role']=='covalently_connected':raise ValueError('Covalent reference ligand requires an explicitly bonded redocking protocol')
                    row.update(action='removed',role='selected_reference_ligand')
            report['sanitized_sha256']=file_digest(path)
            write_json(report,path.with_suffix('.sanitation.json'))
    retained={row['name'] for row in report['residue_roles'] if row['action']=='retained' and row['role'] in ('cofactor_or_metal','covalently_connected','explicit_cofactor','retained_unclassified')}
    text=text[:match.start()]+('open '+str(path))+text[match.end():]
    text=text.replace(selection,selection+''.join(' | :'+code for code in sorted(retained) if code not in keep),1)
    return text
