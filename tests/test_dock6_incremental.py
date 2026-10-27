"""DOCK6 external calls are replaced; cache, dependencies and outputs are real."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from caad.docking import Dock6

MOL2='@<TRIPOS>MOLECULE\n{code}\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C 0.0 0.0 0.0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n1 LIG 1\n'


class Dock6IncrementalTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.number=0;self.calls=[];self.fail=set()
        self.env=patch.dict(os.environ,{'BIOMOL_PREPARATION_CACHE':str(self.root/'shared/ligands'),
                                      'BIOMOL_REUSE_PREPARATIONS':'1'})
        self.env.start();self.addCleanup(self.env.stop)
        self.install=self.root/'dock6'
        for folder,names in (('bin',('dock6','sphgen','sphere_selector','showbox','grid')),
                             ('parameters',('vdw_AMBER_parm99.defn','flex.defn','flex_drive.tbl'))):
            for name in names:
                path=self.install/folder/name;path.parent.mkdir(exist_ok=True,parents=True);path.write_text('tool '+name)
        self.rec=self.root/'rec';self.rec.mkdir()
        (self.rec/'1ABC_A.dockprep.mol2').write_text(MOL2.format(code='REC'))
        (self.rec/'1ABC_A.noH.pdb').write_text('receptor')

    def run_engine(self,codes,mode='flex',center=None,expected=None):
        self.number+=1;inputs=self.root/('ligands'+str(self.number));inputs.mkdir()
        for code in codes:(inputs/(code+'.lig.mol2')).write_text(MOL2.format(code=code))
        engine=Dock6(str(self.install),str(inputs),str(self.rec),str(self.root/('run'+str(self.number))),
                     '1ABC_A',conformer_search_type=mode,binding_site_center=center or [0,0,0])
        output=Path(engine._Dock6__base_output_path)
        def site():
            self.calls.append('site')
            for relative in ('surface/Molecules/binding_site.sph','showbox/binding_site.box.pdb',
                             'gridbox/binding_site.nrg','gridbox/binding_site.bmp'):
                path=output/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('valid site '+relative)
        def minimize(ligand,tid):
            code=ligand.removesuffix('.lig.mol2');self.calls.append('min:'+code)
            (Path(engine.outputpath)/(code+'.lig.min_scored.mol2')).write_text(MOL2.format(code=code)+'########## Grid_Score: -5\n')
        def dock(tid):
            text=(Path(engine.outputpath)/(str(tid)+'_docking.in')).read_text()
            prefix=next(line.split()[-1] for line in text.splitlines() if line.strip().startswith('ligand_outfile_prefix'))
            code=Path(prefix).name;self.calls.append('dock:'+code)
            if code in self.fail:raise RuntimeError('controlled ligand failure')
            Path(prefix+'_scored.mol2').write_text(MOL2.format(code=code)+'########## Grid_Score: -7\n')
        try:
            with patch.object(engine,'prepare_surface',side_effect=site),patch.object(engine,'prepare_showbox'), \
                 patch.object(engine,'prepare_gridbox'),patch.object(engine,'perform_parallel_minimization',side_effect=minimize), \
                 patch.object(engine,'perform_parallel_docking',side_effect=dock):
                engine.prepare_site();engine.prepare_minimization();engine.perform_dock6_evaluation()
            self.assertEqual(len(list((output/mode).glob('*_scored.mol2'))),len(codes) if expected is None else expected)
        finally:engine.export_results()

    def test_same_request_new_ligand_and_mode_change_reuse_only_compatible_work(self):
        self.run_engine(['A','B']);self.assertEqual(self.calls.count('site'),1)
        before=list(self.calls);self.run_engine(['A','B']);self.assertEqual(self.calls,before)
        self.run_engine(['A','B','C'])
        self.assertEqual(self.calls.count('site'),1)
        self.assertEqual(self.calls.count('min:A'),1);self.assertEqual(self.calls.count('dock:A'),1)
        self.assertEqual(self.calls.count('min:C'),1);self.assertEqual(self.calls.count('dock:C'),1)
        self.run_engine(['A','B','C'],mode='rigid')
        self.assertEqual(self.calls.count('site'),1);self.assertEqual(self.calls.count('min:A'),1)
        self.assertEqual(self.calls.count('dock:A'),2)

    def test_changed_site_and_force_rerun_invalidate_affected_calculations(self):
        self.run_engine(['A']);self.run_engine(['A'],center=[1,0,0])
        self.assertEqual(self.calls.count('site'),2);self.assertEqual(self.calls.count('dock:A'),2)
        with patch.dict(os.environ,{'BIOMOL_REUSE_PREPARATIONS':'0'}):self.run_engine(['A'],center=[1,0,0])
        self.assertEqual(self.calls.count('site'),3);self.assertEqual(self.calls.count('dock:A'),3)

    def test_failure_resumes_only_missing_compound(self):
        self.fail={'B'};self.run_engine(['A','B'],expected=1)
        self.fail.clear();self.run_engine(['A','B'])
        self.assertEqual(self.calls.count('site'),1)
        self.assertEqual(self.calls.count('min:A'),1);self.assertEqual(self.calls.count('min:B'),1)
        self.assertEqual(self.calls.count('dock:A'),1);self.assertEqual(self.calls.count('dock:B'),2)


if __name__=='__main__':unittest.main()
