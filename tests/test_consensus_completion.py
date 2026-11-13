"""Consensus must complete without duplicate input transport or GUI plotting."""
import csv
import json
import tempfile
import unittest
import threading
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.consensus_io import read_scores, materialize_scores, write_table
from biomolexplorer.docking_data import write_csv


class ConsensusCompletionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def inputs(self):
        for engine in ('vina','dock6'):
            folder=self.root/engine
            folder.mkdir()
            pose=folder/'pose.mol2'
            pose.write_text('original pose')
            write_csv(folder/'docking_results.csv',[
                dict(molecule_chembl_id='M1',score=-7,receptor_id='T',engine=engine,
                     canonical_smiles='CCO',conformer_file='pose.mol2')])

    def test_authoritative_summary_never_walks_pose_tree(self):
        self.inputs()
        with patch.object(Path,'rglob',side_effect=AssertionError('Recursive scan')):
            rows=list(read_scores(self.root/'vina','vina'))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['conformer_file'],str(self.root/'vina/pose.mol2'))

    def test_materialization_preserves_authorized_paths_without_copying_poses(self):
        self.inputs()
        rows=list(read_scores(self.root/'vina','vina'))
        destination=self.root/'transport'
        materialize_scores(rows,[self.root/'vina/docking_results.csv'],destination)
        with (destination/'docking_results.csv').open() as stream:
            result=next(csv.DictReader(stream))
        self.assertEqual(result['conformer_file'],str(self.root/'vina/pose.mol2'))
        self.assertFalse((destination/'poses').exists())
        previous=(destination/'consensus_input_files.json').read_text()
        (self.root/'vina/pose.mol2').write_text('modified pose')
        materialize_scores(rows,[self.root/'vina/docking_results.csv'],destination)
        self.assertNotEqual((destination/'consensus_input_files.json').read_text(),previous)

    def test_table_is_complete_before_plotting_and_plot_failure_is_nonfatal(self):
        from wrappers.docking import generate_consensus
        self.inputs()
        output=self.root/'output'
        def fail_plot(frame,directory):
            self.assertTrue((output/'ranking.csv').is_file())
            raise RuntimeError('No graphical backend')
        with patch('wrappers.docking.plot_scatter_comparison',side_effect=fail_plot):
            frame=generate_consensus(str(self.root),str(output),'ranking',
                base_vina_path=str(self.root/'vina'),base_dock6_path=str(self.root/'dock6'))
        self.assertEqual(len(frame),1)
        self.assertEqual(json.loads((output/'consensus_report.json').read_text())['rows'],1)
        self.assertEqual((self.root/'vina/pose.mol2').read_text(),'original pose')

    def test_atomic_csv_failure_preserves_previous_table(self):
        import pandas as pd
        path=self.root/'ranking.csv'
        path.write_text('previous complete table')
        with patch.object(pd.DataFrame,'to_csv',side_effect=OSError('write failure')):
            with self.assertRaises(OSError):write_table(pd.DataFrame([{'x':1}]),path)
        self.assertEqual(path.read_text(),'previous complete table')
        self.assertEqual(list(self.root.iterdir()),[path])

    def test_plot_uses_agg_without_pyplot_gui_factory(self):
        import pandas as pd
        from wrappers.docking import plot_scatter_comparison
        frame=pd.DataFrame([dict(vina=1.,dock6=1.,normalized_score=1.,consensus_model='test')])
        with patch('wrappers.docking.plt.subplots',side_effect=AssertionError('GUI figure')):
            plot_scatter_comparison(frame,str(self.root))
        self.assertGreater((self.root/'correlation.png').stat().st_size,0)

    def test_terminal_job_does_not_rewrite_running_history_for_final_progress(self):
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.catalog import new_stage
        service=object.__new__(PipelineService)
        service.store=SimpleNamespace(project_dir=lambda pid:self.root,_require_user=lambda *args:None)
        service._lock=threading.RLock()
        service.cpu_workers=1
        service.worker_python=None
        service.job_timeout=None
        service._resolve=lambda *args,**kwargs:{}
        stage=new_stage('consensus')
        item={'configuration':stage}
        job=dict(id='job',status='succeeded',log_path='log',
                 progress={'message':'worker done'},result={'artifacts':[],'details':{}})
        manager=unittest.mock.Mock()
        manager.submit.return_value=job
        manager.get.return_value=job
        with patch('biomolexplorer.pipeline.JobManager',return_value=manager), \
             patch('biomolexplorer.pipeline.materialize_templates',return_value=self.root), \
             patch('biomolexplorer.stage_cache.runtime_signature',return_value={}), \
             patch('biomolexplorer.pipeline.input_manifests',return_value={}), \
             patch('biomolexplorer.artifact_naming.publish',return_value=[]), \
             patch.object(service,'_update',create=True) as update:
            service._execute_single('run','project','user',item,[item],self.root,{},
                                    {'event':threading.Event()},stage=stage)
        self.assertEqual(update.call_count,1)
        self.assertEqual(item['progress'],job['progress'])
        manager.close.assert_called_once()

    def test_history_is_replayable_without_retaining_all_decoded_runs(self):
        from biomolexplorer.published_results import ReviewHistory
        history=ReviewHistory(object(),list(range(12)))
        with patch('biomolexplorer.result_catalog.read_proofs',side_effect=lambda store,rid:[{'id':str(rid)}]):
            self.assertEqual([rid for rid,stages in history],list(range(12)))
            self.assertLessEqual(len(history.loaded),2)
            self.assertEqual([rid for rid,stages in history],list(range(12)))
            self.assertLessEqual(len(history.loaded),2)

    def test_reuse_stops_history_scan_at_first_valid_match(self):
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.catalog import new_stage
        from biomolexplorer.published_results import ReviewHistory
        stage=new_stage('consensus')
        history=ReviewHistory(object(),['new','must-not-open'])
        service=object.__new__(PipelineService)
        service.store=object()
        previous={'id':stage['id'],'configuration':stage,'status':'succeeded','artifacts':['table.csv'],
                  'artifact_manifest':{'table.csv':'sha'},'artifact_records':[{'path':'table.csv'}]}
        with patch('biomolexplorer.result_catalog.read_proofs',return_value=[previous]) as read, \
             patch('biomolexplorer.published_results.retained_configuration',return_value=stage), \
             patch('biomolexplorer.published_results._hydrate',side_effect=lambda store,rid,entry:entry):
            reused=service._reuse_stage('project','user',stage,{},'software',[],history)
        self.assertIsNotNone(reused)
        self.assertEqual(read.call_count,1)
        self.assertEqual(history.loaded[0][0],'new')


if __name__=='__main__':
    unittest.main()
