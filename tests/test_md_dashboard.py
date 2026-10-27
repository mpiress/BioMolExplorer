"""Result-panel construction and equivalent-protocol comparison contracts."""
from types import SimpleNamespace
import unittest
from unittest.mock import Mock

from biomolexplorer.md_results import MDResults
from biomolexplorer.ui.md_results import MDResultsTable


class ComparisonTests(unittest.TestCase):
    def test_comparison_groups_by_receptor_content_and_protocol_not_labels(self):
        service=MDResults.__new__(MDResults)
        service.pairs=Mock(return_value=[dict(id=k,receptor='same label',ligand=k,status='succeeded') for k in ('A','B','C')])
        def dashboard(*args):
            key=args[-1]
            return dict(protocol=dict(inputs={'receptor':{'sha256':'a' if key!='C' else 'b'}},settings={'temperature':300},software={}),
                        analysis={'summary':{'mean_ligand_rmsd':1.}},energies={'amber_score':-1.})
        service.dashboard=dashboard
        rows=service.comparison('token','pid','rid','sid')
        self.assertEqual(rows[0]['group'],rows[1]['group'])
        self.assertNotEqual(rows[0]['group'],rows[2]['group'])
        self.assertEqual([row['ligand'] for row in rows],['A','B','C'])


class DashboardTests(unittest.IsolatedAsyncioTestCase):
    async def test_mobile_dashboard_and_comparison_build_with_only_calculated_metrics(self):
        dialogs=[]
        async def call(fn,*args):return fn(*args)
        table=MDResultsTable.__new__(MDResultsTable)
        table.active=True;table.children=[]
        table.ui=SimpleNamespace(token='token',page=SimpleNamespace(width=360,show_dialog=dialogs.append),call=call)
        table.results=SimpleNamespace(token='token',project_id='pid',run_id='rid',stage={'id':'sid'},valid=lambda:True)
        data=dict(analysis=dict(summary=dict(frames=4,mean_ligand_rmsd=.2,
            most_flexible_residues=[dict(name='ALA',residue=2,rmsf=.1,binding_site=True)],
            most_persistent_contacts=[dict(name='ALA',residue=2,occupancy=.5)]),definitions={'ligand_rmsd':'After protein fitting'}),
            protocol={'settings':{'md_movable_mask':':1'}},energies={'amber_score':-2.})
        table.service=SimpleNamespace(files=lambda *a:[],dashboard=lambda *a:data,comparison=lambda *a:[])
        await table.open('a'*64)
        self.assertEqual(dialogs[-1].content.width,280)
        await table.compare()
        self.assertEqual(len(dialogs),2)


if __name__=='__main__':unittest.main()
