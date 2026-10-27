"""DOCK6 redocking of the selected reference ligands, using shared preparation."""
import math
import shutil
from pathlib import Path


def heavy_atom_rmsd(reference,pose):
    """Symmetry-aware heavy-atom RMSD in the receptor coordinate frame."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from .dock6_validation import validate_mol2
    reference=Chem.RemoveHs(validate_mol2(reference))
    pose=Chem.RemoveHs(validate_mol2(pose))
    if reference.GetNumAtoms()!=pose.GetNumAtoms():raise ValueError('A pose DOCK6 não preservou os átomos pesados do ligante.')
    return float(rdMolAlign.CalcRMS(pose,reference))


def redock(records, prepared, output, dock6_path, search_type, checkpoint, settings=None):
    from caad.docking import Dock6
    from .docking_data import binding_center, convert_structure
    from .storage import write_dataframe
    import pandas as pd
    prepared,output=Path(prepared),Path(output)
    rows=[]
    for record in records:
        code,ligand,number,chain=record[:4]
        prefix=f'{code}_{ligand}_{number}{chain}'
        root=output/'dock6'/prefix
        inputs=root/'ligands';inputs.mkdir(parents=True,exist_ok=True)
        source=prepared/(prefix+'.lig.mol2')
        reference=inputs/(prefix+'.mol2')
        shutil.copy2(source,reference)
        from .chimera_mol2 import normalize_carboxylates
        from .dock6_validation import validate_mol2
        normalize_carboxylates(reference)
        try:validate_mol2(reference)
        except ValueError:
            from .chimera_mol2 import canonicalize_component
            from .redocking_config import pair_key
            canonicalize_component(reference,prepared.parent,record,(settings or {}).get(pair_key(record),{}))
            validate_mol2(reference)
        dock=Dock6(dock6_path=dock6_path,ligand_input_path=str(inputs),
            receptor_input_path=str(prepared),base_output_path=str(root/'results'),
            pdb_code=f'{code}_{chain}',conformer_search_type=search_type,
            binding_site_center=binding_center(prepared,record))
        try:
            dock.validate_setup();dock.prepare_site()
            dock.prepare_minimization();dock.perform_dock6_evaluation()
        finally:
            dock.export_results()
        pose=root/'results'/search_type/(prefix+'_scored.mol2')
        if not pose.is_file():raise ValueError('Saída de redocking DOCK6 ausente: '+prefix)
        target=output/(prefix+'.lig.pdbqt')
        convert_structure(pose,target)
        rmsd=heavy_atom_rmsd(reference,pose)
        if not math.isfinite(rmsd):raise ValueError('RMSD inválido no redocking DOCK6: '+prefix)
        checkpoint(list(record),rmsd)
        row=dict(zip(['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION'],record))
        row['RMSD']=rmsd;rows.append(row)
    write_dataframe(pd.DataFrame(rows),prepared.parent/'pdb_codes.csv')
