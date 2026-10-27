"""Preparation reuse must preserve science, partial work, and failure isolation."""
from contextlib import ExitStack
import os
from pathlib import Path
import tempfile
import threading
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from biomolexplorer.docking_data import read_compounds, write_csv
from biomolexplorer.docking_preparation import prepare_candidates
from biomolexplorer.molecule_quality import CURRENT_REPORT as QUALITY, QualityReport
from biomolexplorer.preparation_failures import CURRENT_REPORT, PreparationFailures


class PreparationPerformanceTests(unittest.TestCase):
    PDBQT='ROOT\nATOM      1  C   UNL     1       0.000   1.000   2.000  0.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n'
    MOL2='@<TRIPOS>MOLECULE\nA\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C 0.0 1.0 2.0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n'
    def setUp(self):
        self.stack=ExitStack();self.addCleanup(self.stack.close)
        self.root=Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        self.dataset=self.root/'input.csv'
        self.cache=self.root/'cache'
        self.stack.enter_context(patch.dict(os.environ,{
            'BIOMOL_PREPARATION_CACHE':str(self.cache),'BIOMOL_REUSE_PREPARATIONS':'1',
            'BIOMOL_CPU_WORKERS':'2'}))
        self.stack.enter_context(patch('biomolexplorer.ligand_cache.preparation_identity',return_value={'version':'test'}))
        self.stack.enter_context(patch('biomolexplorer.docking_preparation.molecule_sdf',side_effect=lambda smiles:smiles))
        self.stack.enter_context(patch('biomolexplorer.docking_preparation.convert_structure',side_effect=self.convert))
        self.calls=[]
        self.runner=self.stack.enter_context(patch('biomolexplorer.docking_preparation.subprocess.run',side_effect=self.prepare))
        self.write([('A','CCO'),('B','CCC')])

    @staticmethod
    def convert(source,destination,**kwargs):
        destination=Path(destination)
        destination.write_text(PreparationPerformanceTests.PDBQT if destination.suffix=='.pdbqt' else Path(source).read_text())

    def prepare(self,command,**kwargs):
        work=Path(kwargs['cwd'])
        self.calls.append((work/'prepare.com').read_text())
        (work/'prepared.mol2').write_text(self.MOL2)
        return SimpleNamespace(returncode=0,stdout='',stderr='')

    def write(self,entries):
        write_csv(self.dataset,[dict(molecule_chembl_id=code,canonical_smiles=smiles) for code,smiles in entries])

    def run_batch(self,name,**kwargs):
        kwargs['config']=dict(preparation_backend='legacy_chimera',**kwargs.get('config',{}))
        return read_compounds(prepare_candidates(self.dataset,self.root/name,engines=kwargs.pop('engines','vina'),**kwargs))

    def test_receptors_and_overlapping_subsets_reuse_only_selected_compounds(self):
        first=self.run_batch('receptor1')
        self.write([('B','CCC'),('C','CCN')])
        second=self.run_batch('receptor2')
        self.assertEqual(self.runner.call_count,3)
        self.assertEqual([r['molecule_chembl_id'] for r in second],['B','C'])
        self.assertEqual(Path(first[1]['prepared_pdbqt']).read_bytes(),Path(second[0]['prepared_pdbqt']).read_bytes())
        self.assertNotEqual(first[1]['prepared_pdbqt'],second[0]['prepared_pdbqt'])

    def test_scientific_settings_and_smiles_invalidate_cache(self):
        self.write([('A','CCO')])
        self.run_batch('original')
        for name,kwargs in [('ph',{'ph':6.5}),('charge',{'config':{'ligand':{'charge_type':'am1'}}}),
                            ('minimize',{'config':{'ligand':{'minimize':False}}}),('format',{'engines':'dock6'})]:
            with self.subTest(name=name):self.run_batch(name,**kwargs)
        self.write([('A','CCN')]);self.run_batch('changed-smiles')
        self.assertEqual(self.runner.call_count,6)
        self.assertTrue(any('method am1' in command for command in self.calls))
        self.assertTrue(any('minimize spec #0' not in command for command in self.calls))

    def test_pose_content_and_conformer_index_invalidate_cache(self):
        pose=self.root/'pose.mol2';pose.write_text('original pose')
        def table(index=1):
            write_csv(self.dataset,[dict(molecule_chembl_id='A',canonical_smiles='CCO',pose_file=str(pose),conformer_index=index)])
        table();self.run_batch('a');self.run_batch('b')
        pose.write_text('different coordinates');self.run_batch('c')
        table(2);self.run_batch('d')
        self.assertEqual(self.runner.call_count,3)

    def test_changed_tool_identity_invalidates_cache(self):
        self.run_batch('a')
        with patch('biomolexplorer.ligand_cache.preparation_identity',return_value={'version':'changed'}):
            self.run_batch('b')
        self.assertEqual(self.runner.call_count,4)

    def test_corrupt_and_missing_cached_outputs_are_recomputed(self):
        self.run_batch('a')
        entries=sorted(self.cache.glob('*/ligand.pdbqt'))
        entries[0].write_text('corruption');entries[1].unlink()
        self.run_batch('b')
        self.assertEqual(self.runner.call_count,4)

    def test_explicit_rerun_bypasses_reads(self):
        self.run_batch('a')
        with patch.dict(os.environ,{'BIOMOL_REUSE_PREPARATIONS':'0'}):self.run_batch('b')
        self.assertEqual(self.runner.call_count,4)

    def test_failed_compounds_are_not_cached_or_confused_with_other_rows(self):
        failure=PreparationFailures();quality=QualityReport()
        prep_token=CURRENT_REPORT.set(failure);quality_token=QUALITY.set(quality)
        self.addCleanup(CURRENT_REPORT.reset,prep_token);self.addCleanup(QUALITY.reset,quality_token)
        original=self.prepare
        def fail(command,**kwargs):
            if (Path(kwargs['cwd'])/'input.mol2').read_text()=='CCC':raise ValueError('bad ligand')
            return original(command,**kwargs)
        self.runner.side_effect=fail
        rows=self.run_batch('a')
        self.assertEqual([r['molecule_chembl_id'] for r in rows],['A'])
        self.assertEqual(failure.records[0]['compound_identifier'],'B')
        self.assertEqual(failure.records[0]['input_row'],3)
        self.assertEqual(quality.records[0]['identifiers'],['B'])
        self.assertEqual(len(list(self.cache.glob('*/manifest.json'))),1)
        self.runner.side_effect=original
        self.assertEqual(len(self.run_batch('b')),2)

    def test_partial_work_survives_fatal_interruption(self):
        original=self.prepare
        def interrupt(command,**kwargs):
            if (Path(kwargs['cwd'])/'input.mol2').read_text()=='CCC':
                raise FileNotFoundError(2,'No such file','chimera')
            return original(command,**kwargs)
        self.runner.side_effect=interrupt
        with patch.dict(os.environ,{'BIOMOL_CPU_WORKERS':'1'}):
            from biomolexplorer.processes import ScientificToolError
            with self.assertRaises(ScientificToolError) as failure:self.run_batch('interrupted')
        self.assertEqual(failure.exception.error_code,'TOOL_NOT_FOUND')
        self.assertEqual([r['molecule_chembl_id'] for r in read_compounds(self.root/'interrupted/compounds.csv')],['A'])
        self.runner.reset_mock();self.runner.side_effect=original
        self.assertEqual(len(self.run_batch('resumed')),2)
        self.runner.assert_called_once()

    def test_parallelism_obeys_budget_and_preserves_input_order(self):
        self.write([(str(i),'C'*(i+2)) for i in range(6)])
        lock=threading.Lock();active=peak=0
        original=self.prepare
        def slow(command,**kwargs):
            nonlocal active,peak
            with lock:active+=1;peak=max(peak,active)
            try:
                time.sleep(.04)
                return original(command,**kwargs)
            finally:
                with lock:active-=1
        self.runner.side_effect=slow
        rows=self.run_batch('parallel')
        self.assertEqual(peak,2)
        self.assertEqual([r['molecule_chembl_id'] for r in rows],list(map(str,range(6))))
        self.assertTrue(all('method gas' in c and 'minimize spec #0' in c for c in self.calls))

    def test_already_prepared_inputs_skip_chimera(self):
        source=self.root/'existing.pdbqt';source.write_text(self.PDBQT)
        write_csv(self.dataset,[dict(molecule_chembl_id='A',canonical_smiles='CCO',prepared_pdbqt=str(source))])
        row=self.run_batch('exported')[0]
        self.runner.assert_not_called()
        self.assertEqual(Path(row['prepared_pdbqt']),source)

    def test_nonfinite_pdbqt_is_rejected_before_cache_or_vina(self):
        from biomolexplorer.ligand_cache import validate_prepared_ligand
        invalid=self.root/'nan.pdbqt';invalid.write_text(self.PDBQT.replace('0.000   1.000   2.000','   nan     nan     nan'))
        with self.assertRaisesRegex(ValueError,'coordenadas'):
            validate_prepared_ligand(invalid,'pdbqt')


if __name__=='__main__':unittest.main()
