"""SMILES parsing with explicit handling of intentionally retained hydrogens."""
import re

_EXPLICIT_H = re.compile(r'\[(?:\d+)?H(?![a-z])')


def mol_from_smiles(smiles, *, retained_hydrogens=None):
    """Keep RDKit's molecular identity, including isolated H, ions and isotopes.

    RDKit's default SMILES parser tries to remove H atoms and warns for every
    unbonded H. Parse explicit H first, then remove only the usual removable H.
    This is local to this molecule; parsing errors and other RDKit logs remain
    enabled. Callers may record retained atoms in their quality report.
    """
    from rdkit import Chem
    if not isinstance(smiles, str) or not _EXPLICIT_H.search(smiles):
        return Chem.MolFromSmiles(smiles)
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    molecule = Chem.MolFromSmiles(smiles, parser)
    if molecule is None:return None
    isolated = sum(atom.GetAtomicNum() == 1 and atom.GetDegree() == 0 for atom in molecule.GetAtoms())
    removal = Chem.RemoveHsParameters()
    removal.removeDegreeZero = False
    # Retention is intentional and reported by the input quality boundary.
    # Keep other removal warnings enabled when there are no isolated H atoms.
    removal.showWarnings = not bool(isolated)
    molecule = Chem.RemoveHs(molecule, removal)
    if isolated and retained_hydrogens is not None:retained_hydrogens(isolated)
    return molecule
