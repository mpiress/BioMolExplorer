"""Authorization, pagination, previews and atomic dataset curation."""
import asyncio
import csv
import json
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.compound_tables import CompoundTables,compound_links
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.visualizations import molecule_conformer,molecule_sdf
from biomolexplorer.workspace import AccessDenied,WorkspaceStore

try:
    from biomolexplorer.ui.compound_table import CompoundTableViewer
except ImportError:CompoundTableViewer=None


class CompoundTablesTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.temp.name)
        self.owner=self.store.register('Owner','owner@example.org','owner-password')
        self.guest=self.store.register('Guest','guest@example.org','guest-password')
        self.project=self.store.create_project(self.owner,'Curation');self.pid=self.project['id']
        self.rid=uuid4().hex;self.stage=new_stage('retrieve_compounds');self.sid=self.stage['id']
        root=self.store.project_dir(self.pid)/'runs'/self.rid/self.sid/'artifacts'
        self.path=root/'compounds'/'CHEMBL220'/'compounds.csv';self.path.parent.mkdir(parents=True)
        self.original='molecule_chembl_id,canonical_smiles,source,extra\nCHEMBL1,CCO,ChEMBL,"original, metadata"\nPUBCHEM3,O,PubChem,kept\n'
        self.path.write_text(self.original)
        artifacts=[str(self.path)]
        for suffix in ('FULL','MOLS','SIMS'):
            path=root/'ChEMBL'/'DrugBank'/('CHEMBL220_'+suffix+'.csv');path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(self.original);artifacts.append(str(path))
        pubchem=root/'PubChem'/'similars'/'CHEMBL220'/'compounds.csv';pubchem.parent.mkdir(parents=True);pubchem.write_text(self.original);artifacts.append(str(pubchem))
        stage=dict(id=self.sid,name=self.stage['name'],operation=self.stage['operation'],configuration=self.stage,status='succeeded',
                   artifacts=artifacts,artifact_manifest=artifact_manifest(artifacts),finished_at=time.time())
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(self.rid,self.pid,self.store.user(self.owner)['id'],'succeeded',json.dumps([stage]),time.time(),time.time(),None))
        self.tables=CompoundTables(self.store)

    def tearDown(self):self.temp.cleanup()

    def page(self,token=None,**kwargs):
        return self.tables.page(token or self.owner,self.pid,self.rid,self.sid,str(self.path),**kwargs)

    def remove(self,token=None,index=0,version=None):
        return self.tables.remove(token or self.owner,self.pid,self.rid,self.sid,str(self.path),index,version or self.page()['version'])

    def share(self,role):
        self.store.invite(self.owner,self.pid,'guest@example.org',role);self.store.accept_invitation(self.guest,self.pid)

    def test_only_summary_and_integrated_tables_are_available(self):
        choices=self.tables.tables(self.owner,self.pid,self.rid,self.sid)
        self.assertEqual([c['name'] for c in choices],['compounds.csv','CHEMBL220_FULL.csv','CHEMBL220_MOLS.csv','CHEMBL220_SIMS.csv'])
        self.assertTrue(choices[0]['integrated'])

    def test_pagination_search_and_original_row_identity(self):
        page=self.page(limit=1,offset=1)
        self.assertEqual(page['total'],2);self.assertEqual(page['rows'][0]['id'],'PUBCHEM3')
        filtered=self.page(query='pubchem')
        self.assertEqual(filtered['matched'],1);self.assertEqual(filtered['rows'][0]['index'],1)
        self.assertEqual(self.page(query='cco')['rows'][0]['id'],'CHEMBL1')
        self.assertEqual(self.page(offset=100)['rows'],[])

    def test_provider_links_use_compound_ids_and_optional_pubchem_metadata(self):
        rows=self.page()['rows']
        self.assertEqual(rows[0]['links'],[{'provider':'ChEMBL','url':'https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL1'}])
        self.assertEqual(rows[1]['links'],[{'provider':'PubChem','url':'https://pubchem.ncbi.nlm.nih.gov/compound/3'}])
        self.assertEqual(compound_links({'molecule_chembl_id':'Imported1','PubChem_CID':'2244.0'}),
            [{'provider':'PubChem','url':'https://pubchem.ncbi.nlm.nih.gov/compound/2244'}])
        self.assertEqual(len(compound_links({'molecule_chembl_id':'CHEMBL25','PubChem_CID':'2244'})),2)
        for identifier in ('Imported1','CHEMBL1/other','PUBCHEM0','PUBCHEM3?cid=4','2244'):
            with self.subTest(identifier=identifier):
                self.assertEqual(compound_links({'molecule_chembl_id':identifier,'PubChem_CID':'invalid'}),[])

    def test_removal_preserves_metadata_other_tables_and_audit_backup(self):
        original_manifest=self.store.get_run(self.owner,self.rid)['stages'][0]['artifact_manifest']
        self.remove(index=0)
        with self.path.open() as stream:rows=list(csv.DictReader(stream))
        self.assertEqual(rows,[{'molecule_chembl_id':'PUBCHEM3','canonical_smiles':'O','source':'PubChem','extra':'kept'}])
        self.assertEqual((self.path.parents[2]/'ChEMBL'/'DrugBank'/'CHEMBL220_FULL.csv').read_text(),self.original)
        with self.store.connect() as db:edit=dict(db.execute('SELECT * FROM compound_edits').fetchone())
        self.assertEqual(Path(edit['backup']).read_text(),self.original)
        self.assertEqual(edit['compound_id'],'CHEMBL1')
        current=self.store.get_run(self.owner,self.rid)['stages'][0]
        self.assertNotEqual(current['artifact_manifest'],original_manifest)
        self.assertEqual(current['artifact_manifest'],artifact_manifest(current['artifacts']))

    def test_project_rollback_restores_curated_csv_and_its_run_manifest(self):
        self.remove(index=1)
        event=self.store.history(self.owner,self.pid)[0]
        self.assertEqual(event['action'],'curation')
        self.store.rollback(self.owner,self.pid,event['id'])
        self.assertEqual(self.path.read_text(),self.original)
        item=self.store.get_run(self.owner,self.rid)['stages'][0]
        self.assertEqual(item['artifact_manifest'],artifact_manifest(item['artifacts']))
        with self.store.connect() as db:self.assertEqual(db.execute('SELECT COUNT(*) FROM compound_edits WHERE project_id=?',(self.pid,)).fetchone()[0],0)

    def test_viewer_can_read_and_preview_but_cannot_remove(self):
        self.share('viewer')
        self.assertFalse(self.page(self.guest)['can_edit'])
        with self.assertRaises(AccessDenied):self.remove(self.guest)
        self.assertEqual(self.path.read_text(),self.original)

    def test_editor_can_remove_and_uninvited_user_cannot_read(self):
        with self.assertRaises(AccessDenied):self.page(self.guest)
        self.share('editor');self.assertTrue(self.page(self.guest)['can_edit'])
        self.remove(self.guest,index=1);self.assertEqual(self.page()['total'],1)

    def test_active_pipeline_blocks_edit_and_stale_version_cannot_delete(self):
        version=self.page()['version'];self.remove(index=1)
        with self.assertRaisesRegex(ValueError,'Atualize'):self.remove(index=0,version=version)
        with self.store.connect() as db:
            db.execute('UPDATE runs SET status=? WHERE id=?',('running',self.rid))
        self.assertFalse(self.page()['can_edit'])
        with self.assertRaisesRegex(ValueError,'Aguarde'):self.remove()

    def test_selected_file_must_be_registered_and_within_this_project(self):
        other=self.store.create_project(self.owner,'Other')
        with self.assertRaises(AccessDenied):self.tables.page(self.owner,other['id'],self.rid,self.sid,str(self.path))
        with self.assertRaises(AccessDenied):self.tables.page(self.owner,self.pid,self.rid,self.sid,str(self.path.parent/'unknown.csv'))
        with self.assertRaises(ValueError):self.page(limit=1000)

    def test_write_failure_restores_original_dataset_and_rolls_back_audit(self):
        from unittest.mock import patch
        with patch('biomolexplorer.compound_tables.artifact_manifest',side_effect=OSError('disk error')):
            with self.assertRaises(OSError):self.remove()
        self.assertEqual(self.path.read_text(),self.original)
        with self.store.connect() as db:
            self.assertIsNone(db.execute("SELECT name FROM sqlite_master WHERE name='compound_edits'").fetchone())

    def test_3d_conformer_has_reproducible_coordinates_and_valid_bonds(self):
        model=molecule_conformer('CCO')
        self.assertEqual(model,molecule_conformer('CCO'))
        self.assertEqual(len(model['atoms']),9)
        self.assertTrue(all(0<=b['a']<9 and 0<=b['b']<9 for b in model['bonds']))
        for smiles in ('',None,'invalid'):
            with self.subTest(smiles=smiles),self.assertRaises(ValueError):molecule_conformer(smiles)

    def test_sdf_preserves_three_dimensions_and_bond_orders(self):
        from rdkit import Chem
        smiles='CC(=O)O'
        sdf=molecule_sdf(smiles)
        self.assertEqual(sdf,molecule_sdf(smiles))
        molecule=Chem.MolFromMolBlock(sdf,removeHs=False)
        self.assertTrue(molecule.GetConformer().Is3D())
        self.assertEqual(molecule.GetNumAtoms(),len(molecule_conformer(smiles)['atoms']))
        self.assertIn(2.,[b.GetBondTypeAsDouble() for b in molecule.GetBonds()])
        self.assertTrue(sdf.endswith('$$$$\n'))

    @unittest.skipIf(CompoundTableViewer is None,'Install ui extra')
    def test_native_table_buttons_show_generated_views_and_enforce_read_only(self):
        self.share('viewer')
        dialogs=[]
        page=SimpleNamespace(width=1440,update=lambda:None,show_dialog=dialogs.append,pop_dialog=lambda:None)
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):await action()
        ui=SimpleNamespace(page=page,token=self.guest,current=self.store.project(self.guest,self.pid),store=self.store,call=call,guard=guard)
        viewer=CompoundTableViewer(ui,self.pid,self.rid,self.sid,self.tables.tables(self.guest,self.pid,self.rid,self.sid));ui.compound_viewer=viewer
        viewer.build();asyncio.run(viewer.load())
        self.assertEqual(len(viewer.table.rows),2)
        self.assertTrue(viewer.table.rows[0].cells[3].content.disabled)
        for row,provider,url in zip(viewer.table.rows,('ChEMBL','PubChem'),
                ('https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL1','https://pubchem.ncbi.nlm.nih.gov/compound/3')):
            link=row.cells[0].content.controls[1]
            self.assertEqual(link.tooltip,'Abrir no '+provider)
            self.assertEqual(link.url.url,url);self.assertEqual(link.url.target.value,'_blank')
        viewer.selector.value=next(t['path'] for t in viewer.tables if t['name']=='CHEMBL220_MOLS.csv')
        asyncio.run(viewer.change_table(None))
        self.assertIn('CHEMBL220_MOLS.csv',viewer.count.value)
        self.assertEqual(viewer.offset,0)
        asyncio.run(viewer.preview(viewer.data['rows'][0],False));self.assertTrue(dialogs[-1].title.value.endswith('2D'))
        # A success notification must not receive the molecule's close command.
        import flet as ft
        molecule_dialog=dialogs[-1];molecule_dialog.open=True
        notification=ft.SnackBar(content=ft.Text('Saved'),open=True);dialogs.append(notification)
        molecule_dialog.actions[-1].on_click(None)
        self.assertFalse(molecule_dialog.open);self.assertTrue(notification.open)
        from unittest.mock import AsyncMock,Mock,patch
        ui.compound_view_url=Mock(return_value='https://example.org/molecular-viewer/compound-ticket')
        count=len(dialogs)
        with patch('flet.UrlLauncher.launch_url',new_callable=AsyncMock) as launch:
            asyncio.run(viewer.preview(viewer.data['rows'][0],True))
            launch.assert_awaited_once()
            self.assertEqual(launch.call_args.args[0],'https://example.org/molecular-viewer/compound-ticket')
            self.assertEqual(launch.call_args.kwargs['web_only_window_name'],'_blank')
            ui.compound_view_url.assert_called_once_with(self.pid,'CCO','CHEMBL1',self.guest)
            self.assertEqual(len(dialogs),count)
            ui.token='new-session';launch.reset_mock()
            asyncio.run(viewer.preview(viewer.data['rows'][0],True));launch.assert_not_awaited()

    @unittest.skipIf(CompoundTableViewer is None,'Install ui extra')
    def test_stale_compound_preview_does_not_launch_browser(self):
        from unittest.mock import AsyncMock,patch
        async def call(function,*args):return function(*args)
        ui=SimpleNamespace(token=self.owner,current={'id':self.pid},store=self.store,call=call)
        viewer=CompoundTableViewer(ui,self.pid,self.rid,self.sid,self.tables.tables(self.owner,self.pid,self.rid,self.sid),inline=True)
        def generated(*args):
            viewer.preview_version+=1
            return 'https://example.org/molecular-viewer/stale'
        ui.compound_view_url=generated
        with patch('flet.UrlLauncher.launch_url',new_callable=AsyncMock) as launch:
            asyncio.run(viewer.preview({'id':'CHEMBL1','smiles':'CCO'},True))
            launch.assert_not_awaited()
