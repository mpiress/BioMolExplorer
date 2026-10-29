"""DOCK6 6.11 redocking in isolated runtimes, with fixed-frame heavy-atom RMSD."""
import json
import math
from pathlib import Path
import re
import shutil
import tempfile

from .redocking_config import pair_key
from .redocking_preparation import recoverable


def crystal_reference(prepared, pdb, record, destination):
    """Keep the prepared graph/charges, restoring original crystallographic heavy atoms."""
    text=Path(prepared).read_text();coordinates={}
    if Path(pdb).is_file():
        for line in Path(pdb).read_text().splitlines():
            if line.startswith('ENDMDL'):break
            if not line.startswith(('ATOM  ','HETATM')):continue
            if (line[17:20].strip(),line[21:22].strip(),line[22:26].strip()) != (str(record[1]),str(record[3]),str(record[2])):continue
            if line[16:17] not in (' ','A'):continue
            name=line[12:16].strip()
            if name in coordinates:raise ValueError('Referência cristalográfica contém nomes atômicos ambíguos.')
            coordinates[name]=tuple(float(line[i:i+8]) for i in (30,38,46))
        if not coordinates:raise ValueError('Ligante cristalográfico ausente na referência de RMSD.')
    lines=[];active=False
    for line in text.splitlines():
        if line.startswith('@<TRIPOS>'):active=line=='@<TRIPOS>ATOM'
        elif active and line.strip():
            fields=line.split()
            if fields[5].split('.')[0]!='H' and coordinates:
                if fields[1] not in coordinates:raise ValueError('Mapeamento dos átomos pesados da referência incompleto.')
                xyz=coordinates[fields[1]]
                if not all(math.isfinite(v) for v in xyz):raise ValueError('Coordenadas da referência de RMSD inválidas.')
                fields[2:5]=[f'{v:.6f}' for v in xyz];line=' '.join(fields)
        lines.append(line)
    Path(destination).write_text('\n'.join(lines)+'\n')
    return 'crystal' if coordinates else 'prepared'


def heavy_atom_rmsd(reference, pose):
    """Symmetry-aware RMSD without rotation, translation, or hydrogen permutations."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from .dock6_validation import validate_mol2
    # RemoveHs deliberately retains some nonstandard-valence H (e.g. phosphoryl
    # fragments). RMSD is defined over all and only heavy atoms in both engines.
    ref=Chem.RemoveAllHs(validate_mol2(reference))
    if Path(pose).suffix.lower()=='.pdbqt':
        from .docking_scene import best_pose
        text,_=best_pose(Path(pose).read_text(),'pdbqt')
        names={a.GetProp('_TriposAtomName') for a in ref.GetAtoms()};positions={}
        for line in text.splitlines():
            if not line.startswith(('ATOM  ','HETATM')):continue
            name=line[12:16].strip()
            if name not in names:continue
            if name in positions:raise ValueError('A pose contém nomes atômicos ambíguos.')
            positions[name]=tuple(float(line[i:i+8]) for i in (30,38,46))
        if set(positions)!=names:raise ValueError('Mapeamento dos átomos pesados da pose incompleto.')
        probe=Chem.Mol(ref)
        for atom in probe.GetAtoms():probe.GetConformer().SetAtomPosition(atom.GetIdx(),positions[atom.GetProp('_TriposAtomName')])
    else:probe=Chem.RemoveAllHs(validate_mol2(pose,scored=True))
    if Chem.MolToSmiles(ref)!=Chem.MolToSmiles(probe):raise ValueError('A pose DOCK6 alterou a identidade química do ligante.')
    value=rdMolAlign.CalcRMS(probe,ref,maxMatches=100000)
    if not math.isfinite(value) or value<0:raise ValueError('RMSD inválido no redocking DOCK6.')
    return float(value)


def run_pairs(source, output, records, *, dock6_app_path, density=.5, radius=1.4,
              distance=10., conformer_search_type='flex', max_orientations=1000,
              num_scored_conformers=1, dock6_seed=73):
    from caad.docking import Dock6
    from .docking_scene import best_pose
    from .dock6_validation import validate_mol2
    source,output=Path(source),Path(output);prepared=source/'Prepared';output.mkdir(parents=True,exist_ok=True)
    rows={};failures=[]
    with tempfile.TemporaryDirectory(prefix='bme-redock-') as temp:
        surface_cache=Path(temp)/'surfaces';surface_cache.mkdir()
        for record in records:
            prefix=f'{record[0]}_{record[1]}_{record[2]}{record[3]}';receptor=f'{record[0]}_{record[3]}'
            stage=Path(temp)/prefix;ligands=stage/'ligands';ligands.mkdir(parents=True)
            shutil.copy2(prepared/(prefix+'.lig.mol2'),ligands/(prefix+'.lig.mol2'))
            reference=output/(prefix+'.reference.mol2');dock=None;run=None
            try:
                kind=crystal_reference(ligands/(prefix+'.lig.mol2'),source/(record[0]+'.pdb'),record,reference)
                validate_mol2(reference)
                receptor_files=stage/'receptor';receptor_files.mkdir()
                for suffix in ('.dockprep.mol2','.noH.pdb'):
                    shutil.copy2(prepared/(receptor+suffix),receptor_files/(receptor+suffix))
                dock=Dock6(dock6_path=dock6_app_path,ligand_input_path=str(ligands),
                    receptor_input_path=str(receptor_files),base_output_path=str(stage/'protocol'),pdb_code=receptor,
                    density=density,radius=radius,distance=distance,conformer_search_type=conformer_search_type)
                # Native full-receptor surface and spheres are independent of the selected ligand.
                dock.precomputed_surface=surface_cache/receptor
                dock.validate_setup();dock.prepare_surface()
                surface=Path(dock.runtime_path)/'out/surface'
                cached=surface_cache/receptor;cached.mkdir(exist_ok=True)
                for extension in ('.dms','.sph'):shutil.copy2(surface/(receptor+extension),cached/(receptor+extension))
                dock.prepare_showbox();dock.prepare_gridbox()
                run=Path(dock.runtime_path)/'out/search';dock.set_outputpath(str(run))
                local_reference=Path(dock.runtime_path)/'reference.mol2';shutil.copy2(reference,local_reference)
                ligand=Path(dock.ligandpath)/(prefix+'.lig.mol2')
                script=run/'redocking.in';outprefix=run/prefix
                dock.generate_docking_script('src/scripts/dock6/redocking.template',str(script),
                    conformer_search_type=conformer_search_type,ligand_atom_file=str(ligand),
                    rmsd_reference_filename=str(local_reference),
                    receptor_site_file=str(surface/'Molecules'/(prefix+'.sph')),
                    ligand_sphere_file=str(surface/'Molecules'/(prefix+'.sph')),
                    grid_score_grid_prefix=dock._grid_prefix(prefix),dock6_path=str(Path(dock6_app_path).resolve())+'/',
                    ligand_outfile_prefix=str(outprefix),max_orientations=max_orientations,
                    num_scored_conformers=num_scored_conformers,dock6_seed=dock6_seed)
                dock.perform_subprocess(['dock6','-i',str(script),'-o',str(run/'redocking.out')],str(run))
                log=(run/'redocking.out').read_text()
                version=re.search(r'DOCK\s+v?(6\.\d+)',log)
                if not version or version[1]!='6.11':
                    from .processes import ScientificToolError
                    raise ScientificToolError('O redocking requer uma execução validada com DOCK6 6.11.','TOOL_VERSION_UNSUPPORTED','Configure DOCK6 6.11.')
                raw=outprefix.with_name(prefix+'_scored.mol2')
                text,selection=best_pose(raw.read_text(),'mol2')
                pose=output/(prefix+'.lig.mol2');pose.write_text(text)
                value=heavy_atom_rmsd(reference,pose)
                if selection['score'] is None:raise ValueError('Score ausente no resultado: Grid_Score')
                if raw.read_text().count('@<TRIPOS>MOLECULE')>1:
                    shutil.copy2(raw,output/(prefix+'.poses.mol2'))
                shutil.copy2(run/'redocking.out',output/(prefix+'.dock6.out'))
                shutil.copy2(script,output/(prefix+'.dock6.in'))
                info=dict(ENGINE='dock6',ENGINE_VERSION=version[1],RMSD=value,SCORE=selection['score'],
                    RMSD_REFERENCE=kind,RMSD_METHOD='heavy_atom_symmetry_no_fit',POSE_FILE=pose.name,
                    max_symmetry_matches=100000,pose_selection=selection)
                (output/(prefix+'.dock6.json')).write_text(json.dumps(info,indent=2)+'\n')
                rows[pair_key(record)]=info
            except Exception as error:
                if not recoverable(error):raise
                failures.append(dict(zip(('PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION'),list(record[:4])+[record[4] if len(record)>4 else None]),
                    STATUS='preparation_failed',FAILED_STEP='Redocking DOCK6',ERROR=f'{type(error).__name__}: {error}'))
                for path in output.glob(prefix+'.*'):
                    if path.suffix not in ('.out','.in'):path.unlink()
            finally:
                # Preserve native diagnostics even when the pair fails before
                # scoring or RMSD validation; the temporary runtime is removed.
                if run is not None:
                    for name,suffix in (('redocking.in','.dock6.in'),('redocking.out','.dock6.out')):
                        if (run/name).is_file():shutil.copy2(run/name,output/(prefix+suffix))
                if dock is not None:shutil.rmtree(dock.runtime_path,ignore_errors=True)
    return rows,failures
