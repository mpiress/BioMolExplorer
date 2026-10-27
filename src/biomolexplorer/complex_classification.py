"""Evidence-based receptor/ligand interaction classification.

Explicit crystallographic records outrank geometry.  Geometry can make a case
ambiguous, but never invents a covalent bond.  Missing CONECT partners (common
after selecting one alternate location) are ignored instead of being mistaken
for an inter-residue linkage.
"""
from rdkit import Chem

NON_COVALENT='NON_COVALENT'
COVALENT='COVALENT'
AMBIGUOUS='AMBIGUOUS'


class CovalentComplexError(ValueError):
    def __init__(self,report):
        self.report=report
        attachment=report['attachments'][0] if report.get('attachments') else {}
        label=(f"{attachment.get('receptor_residue','?')} {attachment.get('receptor_chain','?')} "
               f"{attachment.get('receptor_number','?')}:{attachment.get('receptor_atom','?')} ↔ "
               f"{attachment.get('ligand_residue','?')}:{attachment.get('ligand_atom','?')}")
        evidence=attachment.get('evidence','explicit crystallographic bond')
        super().__init__('Covalent crystal complex requires a validated bonded Amber topology '
                         f'({evidence}): '+label)


def _residue_key(chain,residue):
    return chain.name,residue.name,residue.seqid.num,residue.seqid.icode.strip()


def _partner_key(partner):
    return partner.chain_name,partner.res_id.name,partner.res_id.seqid.num,partner.res_id.seqid.icode.strip()


def _attachment(ligand_key,ligand_atom,chain,residue,atom,source,distance=None,bond_type='single'):
    return dict(receptor_chain=chain.name,receptor_residue=residue.name,
        receptor_number=residue.seqid.num,receptor_insertion_code=residue.seqid.icode.strip(),
        receptor_atom=atom.name,ligand_chain=ligand_key[0],ligand_residue=ligand_key[1],
        ligand_number=ligand_key[2],ligand_insertion_code=ligand_key[3],ligand_atom=ligand_atom,
        bond_type=bond_type,evidence=source,distance_angstrom=float(distance) if distance is not None else None)


def classify_complex(structure,record,config=None):
    """Classify one selected ligand and retain auditable attachment atoms."""
    from .structure_selection import protein_residue,METALS
    config=config or {};chain_id=config.get('ligand_chain',record[3]);icode=config.get('ligand_insertion_code','')
    matches=[(chain,residue) for chain in structure[0] for residue in chain
             if chain.name==chain_id and residue.name==record[1] and residue.seqid.num==int(record[2])
             and residue.seqid.icode.strip()==icode]
    if len(matches)!=1:raise ValueError('Crystallographic ligand selection is missing or ambiguous; specify chain and insertion code.')
    ligand_chain,ligand=matches[0];ligand_key=_residue_key(ligand_chain,ligand)
    residue_by_key={_residue_key(chain,residue):(chain,residue) for chain in structure[0] for residue in chain}
    atom_by_serial={atom.serial:(chain,residue,atom) for chain in structure[0] for residue in chain for atom in residue}
    ligand_serials={atom.serial for atom in ligand};attachments=[];ambiguous=[]

    for connection in structure.connections:
        partners=(connection.partner1,connection.partner2);keys=tuple(_partner_key(p) for p in partners)
        selected=[i for i,key in enumerate(keys) if key==ligand_key]
        if len(selected)!=1:continue
        li=selected[0];other_key=keys[1-li];other=residue_by_key.get(other_key)
        if other is None:continue
        chain,residue=other;distance=getattr(connection,'reported_distance',None)
        if protein_residue(residue):
            attachments.append(_attachment(ligand_key,partners[li].atom_name,chain,residue,
                next((a for a in residue if a.name==partners[1-li].atom_name),type('Atom',(),{'name':partners[1-li].atom_name})()),
                'mmCIF struct_conn/PDB LINK',distance,'single'))
        else:
            ambiguous.append(dict(evidence='mmCIF struct_conn/PDB LINK',partner=dict(chain=chain.name,
                residue=residue.name,number=residue.seqid.num,atom=partners[1-li].atom_name)))

    seen=set()
    for serial,neighbors in structure.conect_map.items():
        if serial not in ligand_serials:continue
        ligand_atom=atom_by_serial.get(serial)
        if ligand_atom is None:continue
        for neighbor in neighbors:
            other=atom_by_serial.get(neighbor)
            # An absent neighbor is normally an unselected alternate conformer.
            # It is not evidence of an inter-residue bond.
            if other is None:continue
            chain,residue,atom=other
            if _residue_key(chain,residue)==ligand_key:continue
            if protein_residue(residue):
                key=tuple(sorted((serial,neighbor)))
                if key in seen:continue
                seen.add(key);distance=ligand_atom[2].pos.dist(atom.pos)
                attachments.append(_attachment(ligand_key,ligand_atom[2].name,chain,residue,atom,'PDB CONECT',distance,'single'))
            elif atom.element.name in METALS:
                ambiguous.append(dict(evidence='PDB CONECT metal/coordination',partner=atom.name))

    # Distance is supporting evidence only.  A sub-covalent contact without an
    # explicit record is ambiguous and must be curated, never called covalent.
    ligand_atoms=[a for a in ligand if not a.element.is_hydrogen]
    periodic=Chem.GetPeriodicTable()
    for chain,residue in residue_by_key.values():
        if not protein_residue(residue):continue
        for atom in residue:
            if atom.element.is_hydrogen:continue
            for ligand_atom in ligand_atoms:
                cutoff=1.15*(periodic.GetRcovalent(atom.element.atomic_number)+
                             periodic.GetRcovalent(ligand_atom.element.atomic_number))
                distance=atom.pos.dist(ligand_atom.pos)
                if distance<cutoff and not attachments:
                    ambiguous.append(dict(evidence='geometric contact only',distance_angstrom=float(distance),
                        receptor=dict(chain=chain.name,residue=residue.name,number=residue.seqid.num,atom=atom.name),
                        ligand=dict(chain=ligand_key[0],residue=ligand_key[1],number=ligand_key[2],atom=ligand_atom.name)))

    # LINK and CONECT can describe the same bond; collapse them but retain all evidence.
    merged={}
    for bond in attachments:
        key=tuple(bond[k] for k in ('receptor_chain','receptor_number','receptor_atom','ligand_chain','ligand_number','ligand_atom'))
        if key in merged:
            evidence=set(merged[key]['evidence'].split(' + '));evidence.add(bond['evidence'])
            merged[key]['evidence']=' + '.join(sorted(evidence))
        else:merged[key]=bond
    attachments=list(merged.values())
    interaction=COVALENT if attachments else AMBIGUOUS if ambiguous else NON_COVALENT
    return dict(version=1,interaction_type=interaction,attachments=attachments,ambiguities=ambiguous,
        states=dict(crystal_complex='COVALENT_CRYSTAL_COMPLEX' if interaction==COVALENT else interaction,
            vina='COVALENT_DOCKING_UNAVAILABLE' if interaction==COVALENT else 'STANDARD_DOCKING_AVAILABLE',
            dock6='COVALENT_DOCKING_REQUIRES_DOCK_6_12' if interaction==COVALENT else 'STANDARD_DOCKING_AVAILABLE',
            amber='COVALENT_MD_REQUIRES_VALIDATED_PARAMETERS' if interaction==COVALENT else 'NON_COVALENT_MD_SUPPORTED'))
