"""Report overview avoids result I/O and chart limits bound actual work."""
import unittest
from unittest.mock import patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.scientific_report import analyze_run
from biomolexplorer.experimental_report import build_report_preview
import test_semantic_artifacts as fixtures


class ReportPerformanceTests(unittest.TestCase):
    setUp=fixtures.SemanticArtifactTests.setUp
    tearDown=fixtures.SemanticArtifactTests.tearDown

    def test_preview_queries_counts_without_loading_history_or_reading_artifacts(self):
        with patch.object(self.store,'_run_row',side_effect=AssertionError('No execution deserialization')), patch.object(self.store,'scoped_path',side_effect=AssertionError('No artifact reads')):
            preview=build_report_preview(self.store,self.token,self.pid)
        self.assertEqual(preview['project']['id'],self.pid)
        self.assertEqual(preview['document_scope'],{'execution_count':0,'artifact_count':0})

    def test_chart_budget_prevents_rendering_discarded_charts_and_keeps_metrics(self):
        stage=new_stage('admet');stage.update(status='succeeded',artifacts=[])
        for index in range(20):
            path=self.store.project_dir(self.pid)/f'metrics-{index}.csv'
            path.write_text('compound,ADMET_selected,MW\nA,true,100\nB,false,200\nC,true,300\n')
            stage['artifacts'].append(str(path))
        with patch('biomolexplorer.scientific_report._histogram',return_value=b'chart') as draw:
            analyses=analyze_run(self.store,self.pid,{'stages':[stage]},[stage])
        self.assertEqual(draw.call_count,8)
        self.assertEqual(len(analyses[0]['charts']),8)
        self.assertEqual(len({metric['dataset'] for metric in analyses[0]['metrics']}),20)


if __name__=='__main__':unittest.main()
