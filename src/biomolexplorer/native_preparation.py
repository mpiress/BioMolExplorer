"""Main structure-preparation orchestration, without Chimera subprocesses."""
import json
from pathlib import Path

import numpy as np

from .redocking_config import DEFAULTS,pair_key
from .storage import write_json,write_dataframe
from .structure_selection import (read_selected,subset,source_structure,write_structure,
                                  crystal_molecule,protein_residue,PROTEIN,WATER,METALS)


def receptor_pdbqt(mol2,path):
    """Rigid AutoDock export from assigned receptor charges, no reprotonation."""
    from openbabel import openbabel as ob
    from .storage import write_text
    conversion=ob.OBConversion();conversion.SetInAndOutFormats('mol2','pdbqt')
    conversion.AddOption('r',ob.OBConversion.OUTOPTIONS);molecule=ob.OBMol()
    if not conversion.ReadFile(molecule,str(mol2)):raise ValueError('Cannot read charged receptor MOL2.')
    text=conversion.WriteString(molecule)
    if not any(line.startswith(('ATOM  ','HETATM')) for line in text.splitlines()):raise ValueError('Rigid receptor PDBQT contains no atoms.')
    write_text(text,path)


def prepare_structures(dock,records,charge_type,ph,redefine_center=True,settings=None,ready=None):
    from .amber_preparation import receptor_parameters,export_receptor
    from .molecular_core import MolecularPreparationCore,charge_method
    from .molecular_adapters import Dock6Adapter,AutoDockAdapter
    from .diagnostics import event
    from .stage_cache import file_digest
    from .paths import resolve_path
    RESOURCE_ROOT=Path(__file__).parent/'resources'
    import pandas as pd

    settings=settings or {};ready=ready if ready is not None else set();result=[]
    output=Path(dock.outputpath);folder=Path(dock.complexpath)
    centers_path=output/'centers.csv'
    centers=pd.read_csv(centers_path).to_dict('list') if centers_path.is_file() else {}
    for record in records:
        config=settings.get(pair_key(record),{})
        if config.get('preparation_backend','auto')=='legacy_chimera':
            result.extend(dock.prepare_for_docking([record],charge_type,ph,redefine_center,
                preparation_pairs=settings,_prepared_receptors=ready));continue
        # Arbitrary Chimera command edits cannot be reinterpreted as Python.
        for name in ('prepare_complex.template','prepare_receptor.template','prepare_ligand.template'):
            from .redocking_config import configure_template
            actual=resolve_path('src/scripts/chimera/'+name).read_text()
            canonical=(RESOURCE_ROOT/'chimera'/name).read_text()
            if name=='prepare_ligand.template':
                actual=actual.replace('close session\nopen {output_ligand}.lig.pdb\n','')
            if configure_template(actual,name,config,record).strip()!=configure_template(canonical,name,config,record).strip():
                raise ValueError('Custom Chimera template detected: select legacy_chimera explicitly, or migrate its choices to native settings.')
        structure,identities=read_selected(source_structure(folder,record[0]),config.get('structure',{}))
        cofactors=set(config.get('cofactors',[]));chain_id=record[3]
        selected=subset(structure,lambda chain,residue:chain.name==chain_id or residue.name in cofactors)
        prefix=f'{record[0]}_{chain_id}';ligand_prefix=f'{record[0]}_{record[1]}_{record[2]}{chain_id}'
        write_structure(selected,output/(prefix+'.complex.pdb'))
        from .complex_classification import classify_complex
        interaction=classify_complex(structure,record,config)
        ligand,mapping,template=crystal_molecule(structure,record,config,folder)
        if (record[0],chain_id) not in ready:
            options=dict(DEFAULTS,**config.get('receptor',{}))
            receptor=subset(selected,lambda chain,residue:protein_residue(residue) or residue.name in cofactors
                            or (residue.name in WATER and not options['remove_solvent'])
                            or any(atom.element.name in METALS for atom in residue))
            work=output/'amber'/prefix
            parameters,normalized=receptor_parameters(receptor,work,ph,options,config.get('ambertools_path',''),folder,config,
                config.get('ambertools_env',''))
            prepared_receptor=export_receptor(parameters,receptor,output/(prefix+'.dockprep.mol2'),
                output/(prefix+'.prepared.pdb'),minimize=options['minimize'],prefix=config.get('ambertools_path',''),directory=work,
                environment=config.get('ambertools_env',''))
            write_structure(prepared_receptor,output/(prefix+'.noH.pdb'),remove_hydrogens=True)
            from .amber_preparation import peptide_breaks
            from .storage import write_text
            no_h=output/(prefix+'.noH.pdb');text,_=peptide_breaks(no_h.read_text());write_text(text,no_h)
            receptor_pdbqt(output/(prefix+'.dockprep.mol2'),output/(prefix+'.dockprep.pdbqt'))
            write_json(dict(backend='python/AmberTools',source=file_digest(source_structure(folder,record[0])),
                            selected_atoms=identities,parameterization=json.loads(parameters['receptor.parameters.json'].read_text()),
                            minimization=json.loads((work/'minimization.json').read_text()) if options['minimize'] else False),output/(prefix+'.preparation.json'))
            ready.add((record[0],chain_id))
            checkpoint=getattr(dock,'receptor_checkpoint',None)
            if checkpoint:checkpoint(record)
        options=dict(DEFAULTS,**{'charge_type':charge_type,**config.get('ligand',{})})
        initial=ligand.GetConformer().GetPositions().copy()
        prepared=MolecularPreparationCore(output,config.get('ambertools_path',''),config.get('ambertools_env','')).prepare(ligand,
            method=charge_method(options['charge_type']),ph=ph,charge_engine='antechamber',
            add_hydrogens=options['add_hydrogens'],remove_hydrogens=options['remove_hydrogens'],minimize=options['minimize'])
        Dock6Adapter().export(prepared,output/(ligand_prefix+'.lig.mol2'))
        AutoDockAdapter().export(prepared,output/(ligand_prefix+'.lig.pdbqt'))
        from rdkit import Chem
        from .storage import write_text
        write_text(Chem.MolToPDBBlock(prepared.molecule),output/(ligand_prefix+'.lig.pdb'))
        before_min=np.asarray([a['coordinates'] for a in [m['original'] for m in mapping]])
        if not np.array_equal(initial,before_min):raise ValueError('Crystal extraction moved heavy atoms.')
        write_json(dict(backend='python/AmberTools',interaction=interaction,template=str(template),template_sha256=file_digest(template),
                        extraction_heavy_rmsd=0.,original_atom_mapping=mapping,preparation=prepared.metadata),output/(ligand_prefix+'.preparation.json'))
        if redefine_center:
            # Preserve existing mass-weighted center semantics, not centroid.
            masses=np.array([a.GetMass() for a in prepared.molecule.GetAtoms()])
            centers[ligand_prefix]=np.average(prepared.coordinates,axis=0,weights=masses).tolist()
            write_dataframe(pd.DataFrame(centers),centers_path)
        result.append(tuple(record))
        event(dock.logger,'preparation.native.completed','Prepared without Chimera',pair=pair_key(record))
    return result
