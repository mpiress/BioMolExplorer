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


if __name__=='__main__':unittest.main()
