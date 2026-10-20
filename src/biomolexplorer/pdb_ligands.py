"""Explainable ligand triage from deposited coordinates, without extra requests.

Contacts support a binding-pocket hypothesis; they do not establish catalytic
activity. SITE/REMARK 800 annotations supply separate, explicitly named evidence.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from io import StringIO
import math
import re

from Bio.PDB import NeighborSearch, PDBParser, Polypeptide

VERSION = 1
CUTOFF = 4.0
MIN_HEAVY_ATOMS = 6
MIN_CONTACT_RESIDUES = 3
MIN_CONTACT_FRACTION = .5
WATER = {'HOH', 'WAT', 'DOD'}
# Conservative list of common crystallization solvents/buffers. It is not a
# chemical definition: excluded components remain available for manual curation.
ADDITIVES = {'GOL', 'EDO', 'PEG', 'PGE', 'PG4', 'PE4', 'PE5', 'PE8', 'DMS',
             'MPD', 'MES', 'TRS', 'BME', 'DTT', 'ACT', 'FMT', 'EOH', 'IPA'}
COFACTORS = {'FAD', 'FMN', 'NAD', 'NAP', 'NDP', 'HEM', 'HEC', 'HEA', 'PLP', 'SAM', 'SAH', 'COA'}
REASONS = {
    'recommended': 'Candidato em contato com a proteína',
    'active_site': 'Contato com sítio ativo/catalítico anotado',
    'equivalent_copy': 'Cópia equivalente de um candidato priorizado',
    'additive': 'Possível solvente ou aditivo de cristalização',
    'nonorganic': 'Sem carbono observado na estrutura',
    'small': 'Poucos átomos pesados observados',
    'contacts': 'Contato insuficiente com a proteína',
    'cofactor': 'Cofator; há outro candidato orgânico no mesmo PDB',
    'other_site': 'Outro candidato possui contato com sítio ativo anotado',
    'protein_linked': 'Ligação anotada à proteína; revisar preparo do ligante',
    'unsupported_identity': 'Cadeia vazia ou código de inserção; requer seleção manual',
    'different_chain': 'Contato principal com outra cadeia; requer revisão do receptor',
}
COLUMNS = ['PDB_CODE', 'LIGAND', 'LIGAND_NAME', 'RESNUM', 'CHAIN', 'RESOLUTION', 'ICODE',
           'HEAVY_ATOMS', 'CARBON_ATOMS', 'CONTACT_RESIDUES', 'CONTACT_FRACTION',
           'MIN_CONTACT_DISTANCE', 'RECEPTOR_CHAINS', 'CONTACT_RESIDUE_IDS', 'SITE_IDS', 'ACTIVE_SITE_CONTACTS',
           'PROTEIN_LINKED', 'COFACTOR', 'RECOMMENDED', 'REASON', 'REPRESENTATIVE']


def policy():
    return dict(version=VERSION, contact_distance_angstrom=CUTOFF,
        minimum_heavy_atoms=MIN_HEAVY_ATOMS, minimum_contact_residues=MIN_CONTACT_RESIDUES,
        minimum_contact_fraction=MIN_CONTACT_FRACTION, models='first',
        equivalent_copies='same component, protein sequence/numbering and contacting residues',
        active_site='deposited SITE/REMARK 800 annotation; proximity alone is not proof')


def identity(row):
    return f"{row['LIGAND']}:{row['CHAIN']}:{row['RESNUM']}{row.get('ICODE', '')}"


def _key(residue):
    return residue.resname.strip(), residue.parent.id.strip(), residue.id[1], residue.id[2].strip()


def _heavy(residue):
    return [a for a in residue.get_atoms() if a.element not in ('H', 'D')
            and (a.get_occupancy() is None or a.get_occupancy()>0)
            and all(math.isfinite(float(v)) for v in a.coord)]


def _annotations(text):
    sites=defaultdict(set); descriptions=defaultdict(str); current=[]
    links=[]; connections=[]
    for line in text.splitlines():
        if line.startswith('SITE  '):
            site=line[11:14].strip()
            for offset in (18,29,40,51):
                try:
                    key=(line[offset:offset+3].strip(),line[offset+4:offset+5].strip(),
                         int(line[offset+5:offset+9]),line[offset+9:offset+10].strip())
                except ValueError:continue
                sites[site].add(key)
        elif line.startswith('REMARK 800'):
            value=line[10:].strip()
            if value.startswith('SITE_IDENTIFIER:'):
                current=re.findall(r'[A-Za-z0-9]+',value.partition(':')[2])
            elif value.startswith('SITE_DESCRIPTION:'):
                for site in current:descriptions[site]+=value.partition(':')[2].strip()+' '
            elif value and not value.startswith(('EVIDENCE_CODE:', 'SITE')):
                for site in current:descriptions[site]+=value+' '
        elif line.startswith('LINK  '):
            try:
                links.append(tuple((line[o:o+3].strip(),line[o+4:o+5].strip(),
                    int(line[o+5:o+9]),line[o+9:o+10].strip()) for o in (17,47)))
            except ValueError:continue
        elif line.startswith('CONECT'):
            try:connections.append([int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()])
            except ValueError:continue
    return sites,descriptions,links,connections


def analyze_ligands(text, pdb_id=''):
    names=defaultdict(list)
    for line in text.splitlines():
        if line.startswith('HETNAM'):names[line[11:14].strip()].append(line[15:70].strip())
    structure=PDBParser(QUIET=True).get_structure(pdb_id or 'structure',StringIO(text))
    model=next(iter(structure),None)
    if model is None:raise ValueError('Estrutura PDB sem modelo com coordenadas.')
    resolution=structure.header.get('resolution')
    protein=[r for r in model.get_residues() if Polypeptide.is_aa(r,standard=False)]
    protein_keys={_key(r) for r in protein}
    atoms=[a for residue in protein for a in _heavy(residue)]
    neighbors=NeighborSearch(atoms) if atoms else None
    sites,descriptions,links,connections=_annotations(text)
    serials={a.serial_number:r for r in model.get_residues() for a in r.get_unpacked_list()}
    for connection in connections:
        if not connection:continue
        source=serials.get(connection[0])
        if source is None:continue
        for serial in connection[1:]:
            target=serials.get(serial)
            if target is None:continue
            if source is not target:links.append((_key(source),_key(target)))
    # Follow glycosylation branches as well as direct protein-ligand bonds.
    # Do not propagate through ions: LINK can also describe metal coordination.
    organic={_key(r) for r in model.get_residues() if any(a.element=='C' for a in r.get_atoms())}
    graph=defaultdict(set)
    for left,right in links:
        if left in organic and right in organic:
            graph[left].add(right);graph[right].add(left)
    reached=set(protein_keys);pending=list(protein_keys & graph.keys())
    while pending:
        for linked in graph[pending.pop()]-reached:
            reached.add(linked);pending.append(linked)
    covalent=reached-protein_keys
    chain_signatures={}
    for chain in model:
        sequence=[(_key(r)[0],r.id[1],r.id[2]) for r in chain if _key(r) in protein_keys]
        chain_signatures[chain.id]=sha256(repr(sequence).encode()).hexdigest()
    candidates=[];pockets={}
    for residue in model.get_residues():
        name,chain,number,icode=_key(residue)
        if residue.id[0]==' ' or _key(residue) in protein_keys or name in WATER:continue
        heavy=_heavy(residue);contacts={};contact_atoms=0;minimum=None
        for atom in heavy:
            near=neighbors.search(atom.coord,CUTOFF) if neighbors else []
            if near:contact_atoms+=1
            for other in near:
                parent=other.parent;contacts[_key(parent)]=parent
                distance=float(atom-other)
                minimum=distance if minimum is None else min(minimum,distance)
        carbon=sum(a.element=='C' for a in heavy)
        fraction=contact_atoms/len(heavy) if heavy else 0
        chain_counts=Counter(key[1] for key in contacts)
        site_ids=[];active_contacts=0
        for site,keys in sites.items():
            overlap=set(contacts)&keys
            description=descriptions[site]
            ligand=re.search(r'FOR RESIDUE\s+(\w+)\s+(\w)\s+(-?\d+)([A-Za-z]?)',description)
            if ligand and (ligand[1],ligand[2],int(ligand[3]),ligand[4])!=_key(residue):continue
            if len(overlap)>=2:
                site_ids.append(site)
                if re.search(r'\b(?:ACTIVE|CATALYTIC) SITE\b',description,re.I):
                    active_contacts=max(active_contacts,len(overlap))
        reason='recommended'
        if not carbon:reason='nonorganic'
        elif name in ADDITIVES:reason='additive'
        elif len(heavy)<MIN_HEAVY_ATOMS:reason='small'
        elif _key(residue) in covalent:reason='protein_linked'
        elif len(contacts)<MIN_CONTACT_RESIDUES or fraction<MIN_CONTACT_FRACTION:reason='contacts'
        elif icode or not chain:reason='unsupported_identity'
        elif not any(key[1]==chain for key in contacts):reason='different_chain'
        row=dict(PDB_CODE=pdb_id.upper(),LIGAND=name,LIGAND_NAME=' '.join(names[name]),RESNUM=number,CHAIN=chain,
            RESOLUTION=resolution,ICODE=icode,HEAVY_ATOMS=len(heavy),CARBON_ATOMS=carbon,
            CONTACT_RESIDUES=len(contacts),CONTACT_FRACTION=round(fraction,4),
            MIN_CONTACT_DISTANCE=round(minimum,3) if minimum is not None else None,
            RECEPTOR_CHAINS=','.join(sorted(chain_counts)),SITE_IDS=','.join(sorted(site_ids)),
            CONTACT_RESIDUE_IDS=';'.join(f'{n}:{c}:{i}{code}' for n,c,i,code in sorted(contacts)),
            ACTIVE_SITE_CONTACTS=active_contacts,PROTEIN_LINKED=_key(residue) in covalent,
            COFACTOR=name in COFACTORS,RECOMMENDED=reason=='recommended',REASON=reason,REPRESENTATIVE='')
        pockets[identity(row)]=(name,tuple(sorted((chain_signatures[r.parent.id],r.resname,r.id[1],r.id[2]) for r in contacts.values())))
        candidates.append(row)
    # Retain alternatives; the shortlist favors annotated active sites, then
    # non-cofactor binders. No affinity or biological activity is inferred.
    eligible=[r for r in candidates if r['RECOMMENDED']]
    active=any(r['ACTIVE_SITE_CONTACTS'] for r in eligible)
    preferred=any(not r['COFACTOR'] and (not active or r['ACTIVE_SITE_CONTACTS']) for r in eligible)
    for row in eligible:
        if active and not row['ACTIVE_SITE_CONTACTS']:row.update(RECOMMENDED=False,REASON='other_site')
        elif preferred and row['COFACTOR']:row.update(RECOMMENDED=False,REASON='cofactor')
        elif row['ACTIVE_SITE_CONTACTS']:row['REASON']='active_site'
    candidates.sort(key=lambda r:(not r['RECOMMENDED'],-r['ACTIVE_SITE_CONTACTS'],not bool(r['SITE_IDS']),
        -r['CONTACT_FRACTION'],-r['CONTACT_RESIDUES'],r['LIGAND'],r['CHAIN'],r['RESNUM'],r['ICODE']))
    representatives={}
    for row in candidates:
        if not row['RECOMMENDED']:continue
        pocket=pockets[identity(row)]
        if pocket in representatives:
            row.update(RECOMMENDED=False,REASON='equivalent_copy',REPRESENTATIVE=representatives[pocket])
        else:representatives[pocket]=identity(row)
    return dict(resolution=resolution,candidates=candidates,policy=policy())
