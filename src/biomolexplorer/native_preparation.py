"""Chimera preparation compatibility entry points."""
from pathlib import Path

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
    """Compatibility entry point; use the same Chimera route as current blocks."""
    from .redocking_config import chimera_pairs
    return dock.prepare_for_docking(records,charge_type,ph,redefine_center,
        preparation_pairs=chimera_pairs(records,settings or {}),_prepared_receptors=ready)
