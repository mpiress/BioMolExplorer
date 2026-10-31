"""Owned classic-Chimera preparation helpers and narrowly scoped recovery."""
import json
from pathlib import Path
import re


def ligand_hydrogenation(script):
    """Separate native AddH from charging so the observed macrostate is validated."""
    script=Path(script);text=script.read_text()
    charge=re.search(r'^addcharge all spec #0 chargeModel (\S+) method (am1|gas)$',text,re.M)
    output=re.search(r'^write format mol2 #0 (.+\.lig\.mol2)$',text,re.M)
    if not charge or not output:return None
    destination=Path(output[1]);intermediate=destination.parent/'cache'/'hydrogenated'/destination.name
    intermediate.parent.mkdir(parents=True,exist_ok=True)
    first=script.with_suffix('.hydrogenation.com')
    first.write_text(text[:charge.start()]+'write format mol2 #0 '+str(intermediate)+'\nclose session\nclose all\n')
    return dict(script=first,intermediate=intermediate,output=destination,model=charge[1],method=charge[2],tail=text[charge.end():])


def ligand_charging(script,phase):
    """Use the formal charge of the actual, hydrogenated graph with native AM1/gas."""
    from .ligand_mol2 import ensure_ligand_mol2
    from .dock6_validation import validate_mol2
    ensure_ligand_mol2(phase['intermediate'])
    molecule=validate_mol2(phase['intermediate'])
    if any(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
        raise ValueError('O ligante hidrogenado contém um grafo molecular de camada aberta.')
    charge=sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    script=Path(script);helper=script.with_suffix('.charge.py')
    # addNonstandardResCharges propagates ChargeError, unlike the command wrapper.
    helper.write_text('import chimera\nimport sys\nfrom chimera import replyobj\n'
        'from AddCharge import addNonstandardResCharges, _cmdChargeModel\n'
        'def biomol_info(message, *args, **kwargs):\n'
        '    sys.__stdout__.write(str(message))\n    sys.__stdout__.flush()\n'
        '    return original_info(message, *args, **kwargs)\n'
        'original_info = replyobj.info\nreplyobj.info = biomol_info\ntry:\n'
        '    for model in chimera.openModels.list(modelTypes=[chimera.Molecule]):\n'
        '        model.chargeModel = _cmdChargeModel('+repr(phase['model'])+')\n'
        '        addNonstandardResCharges(model.residues, '+str(charge)+', method='+repr(phase['method'])+', gaffType=True)\n'
        'finally:\n    replyobj.info = original_info\n')
    # Native minimize defaults to running DockPrep AGAIN, which can add H and
    # replace the explicitly selected charges. Minimize the already prepared graph.
    def prepared_minimize(match):
        options=re.sub(r'\s+prep\s+(?:true|false)\b','',match[2],flags=re.I)
        return match[1]+options+' prep false'
    tail=re.sub(r'^(minimize spec #0)(.*)$',prepared_minimize,phase['tail'],flags=re.M)
    script.write_text('open '+str(phase['intermediate'])+'\nopen '+str(helper)+'\n'+tail.lstrip('\n'))
    report=phase['output'].with_suffix('.preparation.json')
    hydrogen_report=phase['intermediate'].with_suffix('.preparation.json')
    metadata=json.loads(hydrogen_report.read_text()) if hydrogen_report.is_file() else {}
    metadata.update(net_charge=charge,charge_source='validated hydrogenated observed molecular graph',charge_method=phase['method'])
    report.write_text(json.dumps(metadata,indent=2)+'\n')


def validate_ligand_preparation(phase):
    """Charging/minimization must preserve the validated graph and net charge."""
    from rdkit import Chem
    from .dock6_validation import validate_mol2
    reference=validate_mol2(phase['intermediate']);prepared=validate_mol2(phase['output'])
    if Chem.MolToSmiles(reference,isomericSmiles=False)!=Chem.MolToSmiles(prepared,isomericSmiles=False):
        raise ValueError('O preparo do ligante alterou o grafo molecular validado.')
    expected=sum(atom.GetFormalCharge() for atom in reference.GetAtoms())
    total=sum(atom.GetDoubleProp('_TriposPartialCharge') for atom in prepared.GetAtoms())
    if abs(total-expected)>max(.02,prepared.GetNumAtoms()*.001):
        raise ValueError('As cargas parciais do ligante preparado não correspondem à carga líquida validada.')
    report=phase['output'].with_suffix('.preparation.json');metadata=json.loads(report.read_text())
    metadata.update(validated_graph_preserved=True,partial_charge_sum=total,minimization_preparation='already prepared; no repeated DockPrep')
    report.write_text(json.dumps(metadata,indent=2)+'\n')


def receptor_cleanup(script):
    """Normalize alternate locations with DockPrep, without completing missing atoms."""
    script=Path(script);text=script.read_text()
    helper=script.with_suffix('.cleanup.py')
    source=re.search(r'^open (.+\.pdb)$',text,re.M)
    keep=set(re.findall(r':([A-Z0-9]{1,5})',next((line for line in text.splitlines() if line.startswith('select protein')),'')))
    amino=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
    hetero=set()
    if source:
        for line in Path(source[1]).read_text().splitlines():
            if line.startswith('HETATM') and line[17:20].strip() in amino-keep:
                hetero.add((line[17:20].strip(),int(line[22:26]),line[21:22].strip()))
    cleanup='''import chimera
import DockPrep
excluded = []
for model in chimera.openModels.list(modelTypes=[chimera.Molecule]):
    for residue in list(model.residues):
        identity = (residue.type, residue.id.position, residue.id.chainId.strip())
        if identity not in isolated_candidates:
            continue
        peptide = any(atom.name in ('N', 'C') and any(other.residue != residue and other.name in ('N', 'C') and other.element.number == 6 + (atom.name == 'C') for other in atom.neighbors) for atom in residue.atoms)
        if not peptide:
            excluded.append(identity)
            model.deleteResidue(residue)
DockPrep.prep(chimera.openModels.list(modelTypes=[chimera.Molecule]), addHFunc=None, addCharges=False, delSolvent=False, delIons=False, delLigands=False, delAltLocs=True, incompleteSideChains=False, mutateMSE=False, mutate5BU=False, mutateUMS=False, mutateCSL=False, nogui=True)
import json
open(exclusion_report, 'w').write(json.dumps({'isolated_hetero_amino_acids_not_selected_as_cofactors': excluded}))
'''
    output=re.search(r'^write format mol2 #0 (.+\.mol2)$',text,re.M)
    report=Path(output[1]).with_suffix('.preparation.json') if output else script.with_suffix('.cleanup.json')
    helper.write_text('isolated_candidates = '+repr(hetero)+'\nexclusion_report = '+repr(str(report))+'\n'+cleanup)
    lines=text.splitlines();index=next((i for i,line in enumerate(lines) if line=='addh' or line.startswith('addcharge')),None)
    if index is not None:lines.insert(index,'open '+str(helper))
    script.write_text('\n'.join(lines)+'\n')


def recover_hydrogen_placement(script,error):
    """Recover a known Chimera H-bond code bug using its native geometric AddH."""
    if 'accPhiPsi() takes exactly' not in str(error):return False
    script=Path(script);text=script.read_text()
    replacement,count=re.subn(r'^addh$', 'addh hbond false',text,flags=re.M)
    if not count:return False
    script.write_text(replacement)
    return True


def record_recovery(script):
    script=Path(script)
    match=re.search(r'^write format mol2 #0 (.+\.mol2)$',script.read_text(),re.M)
    if match:
        path=Path(match[1]).with_suffix('.preparation.json')
        metadata=json.loads(path.read_text()) if path.is_file() else {}
        metadata.update(hydrogen_placement='Chimera geometric AddH',
            reason='Native FindHBond accPhiPsi signature mismatch; H-bond-guided placement failed',
            charge_and_minimization_settings='unchanged')
        path.write_text(json.dumps(metadata,indent=2)+'\n')
