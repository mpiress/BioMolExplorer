"""Physical screening of docked poses and identifier-based input recovery."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from biomolexplorer import pose_sources
from biomolexplorer.docking_data import write_csv
from biomolexplorer.operations import OPERATIONS, execute_operation, validate_operation
from biomolexplorer.pose_geometry import evaluate_pose, structure_atoms
from wrappers.pose_validation import REPORT, SUMMARY, validate_poses

SMILES = 'c1ccc(cc1)C(=O)NCCO'


def relaxed_molecule(smiles=SMILES, seed=11):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule,randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(molecule,maxIters=2000)
    return molecule


def write_pose(path, molecule):
    from rdkit import Chem
    writer = Chem.SDWriter(str(path)); writer.write(molecule); writer.close()
    return path


def write_mol2(path, molecule, name='LIG'):
    """Minimal Sybyl MOL2, the format DOCK6 publishes its scored poses in."""
    from rdkit import Chem
    kinds = {Chem.BondType.SINGLE:'1',Chem.BondType.DOUBLE:'2',Chem.BondType.TRIPLE:'3',Chem.BondType.AROMATIC:'ar'}
    conformer = molecule.GetConformer(); atoms = []; bonds = []
    for index,atom in enumerate(molecule.GetAtoms(),1):
        point = conformer.GetAtomPosition(index-1); symbol = atom.GetSymbol()
        atoms.append(f'{index:7d} {symbol}{index:<6} {point.x:9.4f} {point.y:9.4f} {point.z:9.4f} '
                     f'{symbol+(".ar" if atom.GetIsAromatic() else ""):<6} 1 {name:<8} 0.0000')
    for index,bond in enumerate(molecule.GetBonds(),1):
        bonds.append(f'{index:6d} {bond.GetBeginAtomIdx()+1:5d} {bond.GetEndAtomIdx()+1:5d} {kinds[bond.GetBondType()]}')
    Path(path).write_text('@<TRIPOS>MOLECULE\n'+name+f'\n{len(atoms):5d}{len(bonds):6d}     1     0     0\n'
        'SMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n'+'\n'.join(atoms)+'\n@<TRIPOS>BOND\n'+'\n'.join(bonds)+'\n')
    return Path(path)


def write_pdbqt(path, molecule, score=-7.5, name='LIG'):
    """Vina-style PDBQT: one model, PDB coordinate columns and the score remark."""
    conformer = molecule.GetConformer(); lines = ['MODEL 1',f'REMARK VINA RESULT:   {score:.1f}      0.000      0.000']
    for index,atom in enumerate(molecule.GetAtoms(),1):
        point = conformer.GetAtomPosition(index-1); symbol = atom.GetSymbol()
        lines.append(f'ATOM  {index:5d} {symbol+str(index):<4} {name} A   1    '
                     f'{point.x:8.3f}{point.y:8.3f}{point.z:8.3f}  1.00  0.00    +0.000 {symbol:<2}')
    lines += ['ENDMDL','']
    Path(path).write_text('\n'.join(lines))
    return Path(path)


def translated(molecule, offset):
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    import numpy
    probe = Chem.Mol(molecule)
    transform = numpy.eye(4)
    transform[:3,3] = offset
    rdMolTransforms.TransformConformer(probe.GetConformer(),transform)
    return probe


def receptor_pdb(path, points, residue='ALA'):
    lines = []
    for index,point in enumerate(points,1):
        lines.append(f'ATOM  {index:5d}  CA  {residue} A{index:4d}    '
                     f'{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00           C')
    lines.append('END')
    Path(path).write_text('\n'.join(lines)+'\n')
    return Path(path)


def heavy_points(molecule, count=3):
    """Coordinates of the first heavy atoms, to place receptor atoms on them."""
    conformer = molecule.GetConformer()
    points = [conformer.GetAtomPosition(atom.GetIdx()) for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    return [[point.x,point.y,point.z] for point in points[:count]]


def centroid(molecule):
    conformer = molecule.GetConformer()
    points = [conformer.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())]
    return [sum(p[axis] for p in points)/len(points) for axis in range(3)]


class PoseGeometryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.molecule = relaxed_molecule()
        self.pose = write_pose(self.root/'1ABC_LIG_301A_MOL1.sdf',self.molecule)
        self.receptor = receptor_pdb(self.root/'1ABC_A.dockprep.pdb',
                                     [[x+40.0,0.0,0.0] for x in (0.0,3.5,7.0)])

    def test_relaxed_pose_far_from_the_receptor_is_valid(self):
        result = evaluate_pose({'pose':str(self.pose),'receptor':str(self.receptor),'ligand_id':'LIG'})
        self.assertEqual(result['rejection_reasons'],[])
        self.assertTrue(result['pose_valid'])
        self.assertTrue(result['valid_chemistry'])
        self.assertTrue(result['valid_connectivity'])
        self.assertTrue(result['valid_geometry'])
        self.assertTrue(result['valid_planarity'])
        self.assertTrue(result['valid_energy'])
        self.assertTrue(result['valid_clash'])
        self.assertGreater(result['min_vdw_overlap'],0)
        self.assertEqual(result['energy_method'],'MMFF94')

    def test_pose_displaced_onto_receptor_residues_is_rejected_for_steric_clash(self):
        # The receptor residues sit exactly where the displaced pose's heavy
        # atoms are, so the overlap is the full sum of the van der Waals radii.
        occupied = receptor_pdb(self.root/'1ABC_A.occupied.pdb',heavy_points(self.molecule))
        displaced = write_pose(self.root/'clash.sdf',self.molecule)
        result = evaluate_pose({'pose':str(displaced),'receptor':str(occupied),'ligand_id':'LIG'})
        self.assertIn('steric_clash',result['rejection_reasons'])
        self.assertFalse(result['pose_valid'])
        self.assertFalse(result['valid_clash'])
        self.assertLess(result['min_vdw_overlap'],-0.5)
        self.assertLess(result['min_heavy_distance'],2.0)
        self.assertGreater(result['clashing_contacts'],0)
        # The clash is intermolecular: the molecule itself stays valid.
        self.assertTrue(result['valid_chemistry'])
        self.assertTrue(result['valid_geometry'])

    def test_distorted_internal_geometry_and_energy_are_reported(self):
        from rdkit import Chem
        probe = Chem.Mol(self.molecule)
        conformer = probe.GetConformer()
        position = conformer.GetAtomPosition(0)
        conformer.SetAtomPosition(0,(position.x+1.6,position.y,position.z))
        distorted = write_pose(self.root/'distorted.sdf',probe)
        result = evaluate_pose({'pose':str(distorted),'receptor':str(self.receptor)},delta_energy_limit=5.0)
        self.assertFalse(result['pose_valid'])
        self.assertTrue({'bond_lengths','internal_energy'} & set(result['rejection_reasons']))
        self.assertGreater(result['max_bond_deviation'],0.25)
        self.assertGreater(result['delta_energy'],5.0)

    def test_disconnected_fragments_are_rejected(self):
        molecule = relaxed_molecule('CCO.CCO')
        pose = write_pose(self.root/'fragments.sdf',molecule)
        result = evaluate_pose({'pose':str(pose)})
        self.assertIn('fragments',result['rejection_reasons'])
        self.assertFalse(result['valid_connectivity'])

    def test_receptor_parsing_skips_hydrogens_and_solvent(self):
        path = self.root/'solvated.pdb'
        path.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
                        'ATOM      2  H   ALA A   1       1.000   0.000   0.000  1.00  0.00           H\n'
                        'HETATM    3  O   HOH A   2       2.000   0.000   0.000  1.00  0.00           O\n'
                        'HETATM    4 ZN    ZN A   3       3.000   0.000   0.000  1.00  0.00          ZN\nEND\n')
        coordinates,radii,elements = structure_atoms(path)
        self.assertEqual(elements,['C','ZN'])
        self.assertEqual(len(coordinates),2)
        self.assertAlmostEqual(radii[1],1.39)

    def test_vina_pdbqt_recovers_bond_orders_from_the_table_smiles(self):
        pose = write_pdbqt(self.root/'MOL1.pdbqt',self.molecule)
        result = evaluate_pose({'pose':str(pose),'receptor':str(self.receptor),'canonical_smiles':SMILES})
        self.assertEqual(result['bond_orders'],'template')
        self.assertTrue(result['pose_valid'],result['rejection_reasons'])
        self.assertIsNotNone(result['max_planarity_deviation'])
        # Without a reference SMILES only connectivity and clash are decided.
        proximity = evaluate_pose({'pose':str(pose),'receptor':str(self.receptor)})
        self.assertEqual(proximity['bond_orders'],'proximity')
        self.assertIsNone(proximity['max_planarity_deviation'])
        self.assertTrue(proximity['valid_clash'])

    def test_unreadable_pose_is_reported_without_raising(self):
        broken = self.root/'broken.sdf'; broken.write_text('not a molecule\n')
        result = evaluate_pose({'pose':str(broken)})
        self.assertEqual(result['rejection_reasons'],['unreadable'])
        self.assertFalse(result['pose_valid'])


class PoseSourceRecoveryTests(unittest.TestCase):
    """A consensus payload may arrive with identifiers only, never a pose path."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        # Outputs live outside the input tree: a validation run must never
        # rediscover the poses it has just published.
        self.base = Path(self.temp.name); self.root = self.base/'inputs'; self.root.mkdir()
        self.molecule = relaxed_molecule()
        self.poses = self.root/'Dock6'; self.poses.mkdir()
        self.pose = write_mol2(self.poses/'1ABC_LIG_301A_MOL1_scored.mol2',self.molecule)
        write_mol2(self.poses/'1ABC_LIG_301A_MOL2_scored.mol2',relaxed_molecule('CCOCC'))
        self.prepared = self.root/'Prepared'; self.prepared.mkdir()
        self.receptor = receptor_pdb(self.prepared/'1ABC_A.dockprep.pdb',
                                     [[x+40.0,0.0,0.0] for x in (0.0,3.5,7.0)])
        # A consensus ranking with references only: no pose or receptor column.
        write_csv(self.root/'ranking.csv',[
            dict(molecule_chembl_id='MOL1',receptor_id='1ABC_A',canonical_smiles=SMILES,
                 vina=-7.5,dock6=-42.0,normalized_score=.9),
            dict(molecule_chembl_id='MOL2',receptor_id='1ABC_A',canonical_smiles='CCOCC',
                 vina=-5.5,dock6=-30.0,normalized_score=.4)])

    def files(self):
        return [p for p in sorted(self.root.rglob('*')) if p.is_file()]

    def test_poses_and_receptors_are_recovered_by_reference_identifier(self):
        records = pose_sources.collect(self.files())
        by_compound = {record['compound_id']:record for record in records}
        self.assertEqual(set(by_compound),{'MOL1','MOL2'})
        self.assertEqual(by_compound['MOL1']['pose'],str(self.pose))
        self.assertEqual(by_compound['MOL1']['receptor'],str(self.receptor))
        self.assertEqual(by_compound['MOL1']['resolved_by'],'identifier')
        self.assertEqual(by_compound['MOL1']['receptor_id'],'1ABC_A')
        self.assertEqual(by_compound['MOL1']['canonical_smiles'],SMILES)

    def test_identifiers_are_read_from_the_pose_filename_when_no_table_exists(self):
        records = pose_sources.collect([self.pose,self.receptor])
        self.assertEqual(len(records),1)
        self.assertEqual(records[0]['receptor_id'],'1ABC_A')
        self.assertEqual(records[0]['ligand_id'],'LIG')
        self.assertEqual(records[0]['receptor'],str(self.receptor))

    def test_direct_pose_column_is_preferred_over_recovery(self):
        write_csv(self.root/'docking_results.csv',[
            dict(molecule_chembl_id='MOL1',receptor_id='1ABC_A',engine='dock6',score=-42.0,
                 canonical_smiles=SMILES,conformer_file='Dock6/1ABC_LIG_301A_MOL1_scored.mol2',
                 receptor_file='Prepared/1ABC_A.dockprep.pdb')])
        records = [r for r in pose_sources.collect(self.files()) if r['table'].endswith('docking_results.csv')]
        self.assertEqual(len(records),1)
        self.assertEqual(records[0]['pose'],str(self.pose))
        self.assertEqual(records[0]['receptor'],str(self.receptor))
        self.assertNotEqual(records[0].get('resolved_by'),'identifier')

    def test_validation_runs_end_to_end_from_a_consensus_payload(self):
        output = self.base/'validation'
        summary = validate_poses(str(self.root),str(output),target='Alvo')
        self.assertEqual(summary['poses'],2)
        self.assertEqual(summary['identifier_resolved_inputs'],2)
        self.assertEqual(summary['unresolved_inputs'],0)
        with (output/REPORT).open(encoding='utf-8-sig',newline='') as stream:
            rows = {row['compound_id']:row for row in csv.DictReader(stream)}
        self.assertEqual(set(rows),{'MOL1','MOL2'})
        self.assertEqual(rows['MOL1']['receptor_id'],'1ABC_A')
        self.assertEqual(rows['MOL1']['pose_valid'],'True')
        self.assertTrue(float(rows['MOL1']['min_vdw_overlap'])>0)
        self.assertEqual(json.loads((output/SUMMARY).read_text())['valid_poses'],2)
        self.assertEqual({p.name for p in (output/'poses').iterdir()},
                         {self.pose.name,'1ABC_LIG_301A_MOL2_scored.mol2'})
        self.assertEqual({p.name for p in (output/'receptors').iterdir()},{self.receptor.name})

    def test_strict_filter_decides_which_poses_reach_the_next_block(self):
        write_mol2(self.poses/'1ABC_LIG_301A_MOL1_scored.mol2',relaxed_molecule('CCO.CCO'))
        strict = validate_poses(str(self.root),str(self.base/'strict'),target='Alvo')
        self.assertEqual(strict['valid_poses'],1)
        self.assertEqual(strict['published_poses'],1)
        self.assertEqual(strict['rejection_counts'].get('fragments'),1)
        self.assertNotIn('1ABC_LIG_301A_MOL1_scored.mol2',{p.name for p in (self.base/'strict/poses').iterdir()})
        annotated = validate_poses(str(self.root),str(self.base/'annotated'),target='Alvo',strict_filter=False)
        self.assertEqual(annotated['published_poses'],2)
        self.assertEqual(annotated['valid_poses'],1)

    def test_empty_input_explains_the_missing_connection(self):
        with self.assertRaises(ValueError) as error:
            validate_poses(str(self.base/'empty'),str(self.base/'out'),target='Alvo')
        self.assertIn('Vina',str(error.exception))


class PoseValidationContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_operation_is_registered_for_the_cli_and_the_runner(self):
        spec = OPERATIONS['pose_validation']
        self.assertEqual((spec.module,spec.function),('wrappers.pose_validation','validate_poses'))
        validate_operation('pose_validation',{'base_input_path':'poses','target':'Alvo','strict_filter':True})
        for invalid in ({'vdw_clash_scale':-1.0},{'vdw_clash_scale':2.0},{'minimization_steps':0},
                        {'strict_filter':'yes'},{'selected_receptors':['../escape']},
                        {'pose_inputs':[{'receptor':'r'}]}):
            with self.assertRaises(ValueError):
                validate_operation('pose_validation',{'base_input_path':'poses','target':'Alvo',**invalid})

    def test_cli_execution_publishes_the_report_and_the_accepted_poses(self):
        poses = self.root/'input'; poses.mkdir()
        write_pose(poses/'1ABC_LIG_301A_MOL1.sdf',relaxed_molecule())
        receptor_pdb(poses/'1ABC_A.dockprep.pdb',[[40.0,0.0,0.0]])
        result = execute_operation('pose_validation',
            {'base_input_path':str(poses),'target':'Alvo','strict_filter':True},self.root/'output')
        names = {Path(path).name for path in result.artifacts}
        self.assertIn(REPORT,names); self.assertIn(SUMMARY,names)
        self.assertIn('1ABC_LIG_301A_MOL1.sdf',names)
        self.assertEqual(result.details['valid_poses'],1)


if __name__ == '__main__':
    unittest.main()


class PublishedPoseIdentityTests(unittest.TestCase):
    """One compound docked to several receptors keeps one file per pose."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)

    def test_poses_of_different_receptors_are_never_overwritten(self):
        from wrappers.pose_validation import _publish
        molecule = relaxed_molecule('CCO')
        source = self.base/'dock6';source.mkdir()
        write_mol2(source/'MOL1_scored.mol2',molecule)
        rows = [dict(pose_valid=True,compound_id='MOL1',pose_id=code+'_MOL1',receptor_id=code+'_A',
                     pose_file=str(source/'MOL1_scored.mol2'),receptor_file='')
                for code in ('1ABC','2DEF','3GHI')]
        output = self.base/'validation'
        published = _publish(rows,output,True)
        self.assertEqual(len(published),3)
        self.assertEqual({path.name for path in published},
                         {'MOL1_scored.mol2','2DEF_A_MOL1_scored.mol2','3GHI_A_MOL1_scored.mol2'})
        self.assertEqual(len({path.name for path in (output/'poses').iterdir()}),3)
        # The same pose listed twice is still a single file.
        again = _publish(rows[:1]+rows[:1],output,True)
        self.assertEqual({path.name for path in again},{'MOL1_scored.mol2'})


class DeclaredParameterTests(unittest.TestCase):
    """Every parameter the block declares must reach its wrapper."""

    def test_the_wrapper_accepts_the_declared_parameters(self):
        import inspect
        from biomolexplorer.operations import OPERATIONS
        from wrappers.pose_validation import validate_poses
        spec = OPERATIONS['pose_validation']
        accepted = set(inspect.signature(validate_poses).parameters)
        self.assertTrue(set(spec.required+spec.optional).issubset(accepted),
                        sorted(set(spec.required+spec.optional)-accepted))

    def test_a_removed_complex_is_not_validated(self):
        from wrappers.pose_validation import validate_poses, REPORT
        temporary = tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        base = Path(temporary.name)
        molecule = relaxed_molecule('CCO')
        poses = base/'poses';poses.mkdir()
        receptor = base/'1ABC_A.noH.pdb'
        receptor_pdb(receptor,[[80.0,0.0,0.0],[83.5,0.0,0.0]])
        records = []
        for code in ('1ABC_A','2DEF_A'):
            path = poses/(code+'_MOL1.lig.mol2')
            write_mol2(path,molecule)
            records.append(dict(pose=str(path),receptor=str(receptor),pose_id=code+'_MOL1',
                                receptor_id=code,compound_id='MOL1',engine='vina'))
        summary = validate_poses(str(base),str(base/'out'),target='Alvo',pose_inputs=records,
                                 excluded_pairs=['MOL1:2DEF_A'])
        self.assertEqual(summary['poses'],1)
        with (base/'out'/REPORT).open(encoding='utf-8-sig',newline='') as stream:
            self.assertEqual([row['receptor_id'] for row in csv.DictReader(stream)],['1ABC_A'])
