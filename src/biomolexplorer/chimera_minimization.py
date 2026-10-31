"""Minimize an already charged ligand with classic Chimera, without DockPrep."""
from pathlib import Path
import shutil

import numpy as np
from rdkit import Chem

from .calculation_cache import CalculationCache,identity,copy_artifact
from .molecular_core import PreparedMolecule,preparation_workspace
from .stage_cache import file_digest
from .storage import write_text
from .processes import run_command,ScientificToolError


def minimize_ligand(prepared,directory,*legacy,cache_directory=None):
    from .molecular_adapters import Dock6Adapter,read_mapped_charges
    from .preparation_tools import run_amber_tool,amber_tool
    executable=shutil.which('chimera')
    if not executable:raise ScientificToolError('Chimera is required for minimization.','TOOL_NOT_FOUND','Make classic UCSF Chimera available on PATH.')
    executable=Path(executable).resolve()
    library=executable.parent.parent/'share'
    recipe=dict(engine='Chimera MMTK',steps=100,cgsteps=10,prep=False,
        molecule=Chem.MolToMolBlock(prepared.molecule),coordinates=prepared.coordinates.tolist(),charges=prepared.partial_charges.tolist(),
        implementation=file_digest(Path(__file__)),executable=file_digest(executable),
        typing=file_digest(amber_tool('antechamber')),
        modules={str(p.relative_to(library)):file_digest(p) for name in ('MMMD','MMTK','AddCharge','AmberInfo')
                 for p in sorted((library/name).rglob('*')) if p.is_file() and p.suffix in ('.py','.so','.dat')})
    cache=CalculationCache('chimera-ligand-minimization',cache_directory or directory);key=identity(recipe)
    target=Path(directory)/'minimized.mol2'
    with cache.locked(key):
        if cache.restore(key,{'ligand':target}) is None:
            with preparation_workspace(cache.root) as work:
                Dock6Adapter().export(prepared,work/'input.mol2')
                write_text('\n'.join(str(float(v)) for v in prepared.partial_charges)+'\n',work/'charges.dat')
                run_amber_tool('antechamber',['-i','input.mol2','-fi','mol2','-o','typed.mol2','-fo','mol2',
                    '-at','gaff','-c','rc','-cf','charges.dat','-nc',str(Chem.GetFormalCharge(prepared.molecule)),
                    '-an','n','-du','n','-j','5','-dr','n','-s','2','-pf','y'],cwd=work)
                read_mapped_charges(work/'typed.mol2',prepared.molecule)
                block=(work/'typed.mol2').read_text().split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
                types={row[1]:row[5] for row in [line.split() for line in block.splitlines() if line.strip()]}
                # Python repr is safe for filenames containing spaces/quotes.
                script=('import chimera\nfrom chimera import runCommand\n'
                    'models=chimera.openModels.open('+repr(str(work/'input.mol2'))+')\n'
                    'from AddCharge import _cmdChargeModel\n'
                    'types='+repr(types)+'\n'
                    'for model in models:\n    model.chargeModel=_cmdChargeModel(\'14sb\')\n    for atom in model.atoms:\n        atom.gaffType=types[atom.name]\n'
                    "runCommand('minimize spec #0 nsteps 100 cgsteps 10 prep false')\n"
                    'from WriteMol2 import writeMol2\nwriteMol2(chimera.openModels.list(modelTypes=[chimera.Molecule]), '+repr(str(work/'minimized.mol2'))+')\n')
                write_text(script,work/'minimize.py')
                run_command([str(executable),'--nogui','--silent','minimize.py'],cwd=work)
                _read(work/'minimized.mol2',prepared)
                cache.remember(key,{'ligand':work/'minimized.mol2'},dict(recipe=recipe))
                copy_artifact(work/'minimized.mol2',target)
    molecule=_read(target,prepared)
    return PreparedMolecule(molecule,prepared.partial_charges.copy(),dict(prepared.metadata,
        minimization=dict(engine='Chimera MMTK',ligand_force_field='GAFF',solvent='vacuum',steps=100,cgsteps=10,prep=False,cache_key=key)))


def _read(path,prepared):
    text=Path(path).read_text();rows=[r.split() for r in text.split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0].splitlines() if r.strip()]
    table={r[1]:r for r in rows}
    if len(table)!=len(rows) or set(table)!=set(prepared.atom_names):raise ValueError('Chimera minimization changed atom mapping.')
    names=prepared.atom_names;ids={r[0]:r[1] for r in rows}
    edges={frozenset((ids[r[1]],ids[r[2]])) for r in [line.split() for line in text.split('@<TRIPOS>BOND',1)[1].split('@<TRIPOS>',1)[0].splitlines() if line.strip()]}
    expected={frozenset((names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()])) for b in prepared.molecule.GetBonds()}
    if edges!=expected:raise ValueError('Chimera minimization changed connectivity.')
    xyz=np.array([[float(v) for v in table[name][2:5]] for name in names])
    charges=np.array([float(table[name][8]) for name in names])
    if not np.isfinite(xyz).all() or not np.allclose(charges,prepared.partial_charges,atol=.00011,rtol=0):raise ValueError('Chimera minimization changed charges or produced invalid coordinates.')
    molecule=Chem.Mol(prepared.molecule)
    for i,position in enumerate(xyz):molecule.GetConformer().SetAtomPosition(i,position)
    return molecule
