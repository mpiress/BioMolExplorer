"""Analytical trajectories: alignment, physical masks and occupancies."""
import unittest
import numpy as np
from biomolexplorer.md_analysis import analyze_frames,fit_frame,radius


class TrajectoryAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.xyz=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[2.,0.,0.]])
        self.info=dict(natom=5,masses=np.ones(5),residue=np.array([0,0,0,0,1]),labels=['ALA','LIG'])
        self.masks=dict(backbone=np.arange(4),protein=np.arange(4),ca=np.array([1]),ligand=np.array([4]))

    def test_rigid_rotation_and_translation_have_zero_fitted_rmsd(self):
        rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        moved=self.xyz@rotation+np.array([5.,-3.,2.])
        np.testing.assert_allclose(fit_frame(moved,self.xyz,np.arange(4)),self.xyz,atol=1e-12)
        data=analyze_frames([self.xyz,moved],self.info,self.masks,.01)
        self.assertLess(data['summary']['mean_ligand_rmsd'],1e-12)
        self.assertLess(data['rmsf'][0]['rmsf'],1e-12)

    def test_ligand_displacement_is_not_removed_by_ligand_fitting(self):
        moved=self.xyz.copy();moved[4,0]+=10
        data=analyze_frames([self.xyz,moved],self.info,self.masks,.02)
        self.assertAlmostEqual(data['series'][1]['ligand_rmsd'],10.)
        self.assertAlmostEqual(data['contacts'][0]['occupancy'],.5)
        self.assertAlmostEqual(data['summary']['last_frame_ps'],.04)
        self.assertAlmostEqual(data['summary']['sampled_duration_ps'],.02)
        self.assertLess(data['series'][1]['protein_backbone_rmsd'],1e-12)

    def test_mass_weighted_radius_and_invalid_frames(self):
        self.assertAlmostEqual(radius(np.array([[0.,0.,0.],[2.,0.,0.]]),np.ones(2)),1.)
        for frames in ([self.xyz],[self.xyz,self.xyz*np.nan]):
            with self.assertRaises(ValueError):analyze_frames(frames,self.info,self.masks,.01)



class FinalStateCompatibilityTests(unittest.TestCase):
    def test_chimera_histidine_variants_follow_explicit_hydrogens_without_moving_atoms(self):
        from pathlib import Path
        import tempfile
        from biomolexplorer.md_inputs import normalize_receptor_variants
        def atom(index,name,residue):
            return f'ATOM  {index:5d} {name:>4s} HIS A{residue:4d}    {index:8.3f}{1.:8.3f}{2.:8.3f}  1.00 20.00           H\n'
        text=atom(1,'HD1',1)+atom(2,'HE2',2)+atom(3,'HD1',3)+atom(4,'HE2',3)+atom(5,'N',4)+'END\n'
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'source.pdb';target=Path(folder)/'rec.pdb';source.write_text(text)
            changes=normalize_receptor_variants(source,target)
            output=target.read_text().splitlines()
            self.assertEqual([line[17:20] for line in output[:-1]],['HID','HIE','HIP','HIP','HIS'])
            self.assertEqual([line[30:54] for line in output[:-1]],[line[30:54] for line in text.splitlines()[:-1]])
            self.assertEqual(source.read_text(),text)
            self.assertEqual(len(changes),3)

    def test_static_analysis_never_invents_a_trajectory_or_rmsf(self):
        from pathlib import Path
        import tempfile
        from unittest.mock import patch
        from biomolexplorer.md_analysis import analyze_final_state
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            info=dict(natom=1,masses=np.array([1.]),names=['H'])
            with patch('biomolexplorer.md_analysis.selection',return_value=(info,{})),patch('biomolexplorer.md_analysis.coordinates',return_value=np.zeros((1,3))):
                data=analyze_final_state(root)
            self.assertEqual(data['summary']['hydrogen_atoms'],1)
            self.assertNotIn('final_ligand_rmsd',data['summary'])
            self.assertNotIn('rmsf',data)
            self.assertNotIn('series',data)
            self.assertIn('No trajectory or RMSF',data['unavailable'])


if __name__=='__main__':unittest.main()
