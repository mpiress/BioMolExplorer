"""Target selection preserves compound, pose, receptor and capability identity."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock, patch

from test_docking_scene import SceneCapabilityTests, atom
from biomolexplorer.docking_results import DockingResults
from biomolexplorer.result_state import input_digest
from biomolexplorer.workspace import AccessDenied
from biomolexplorer.pdb_view import PREFIX


class ConsensusTargetActionsTests(SceneCapabilityTests):
    def fixture(self):
        stages=[];tables={};targets=['1AAA_A','2BBB_A','3CCC_A']
        for engine in ('vina','dock6'):
            folder=self.root/'runs'/'origin-run'/engine;folder.mkdir(parents=True)
            rows=[];files=[]
            for n,target in enumerate(targets):
                receptor=folder/(target.split('_')[0]+'.pdb')
                receptor.write_text('TITLE     '+['ACETYLCHOLINESTERASE','BUTYRYLCHOLINESTERASE','BACE1'][n]+'\n'+atom(1,'CA','ALA','A',2,n*10))
                pose=folder/(target+'.pdbqt');pose.write_text('REMARK VINA RESULT: -4 0 0\n'+atom(1,'C1','LIG','A',1,n*10+3))
                pdf=folder/(target+'.pdf');pdf.write_bytes(b'%PDF-1.4\n'+target.encode())
                rows.append(dict(molecule_chembl_id='M1',canonical_smiles='CCO',receptor_id=target,score=-4-n,
                    engine=engine,conformer_file=pose.name,receptor_file=receptor.name,footprint_file=pdf.name,footprint_origin='docked_pose'))
                files.extend([pose,receptor,pdf])
            table=folder/'docking_results.csv'
            with table.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            files.append(table);tables[engine]=table
            stages.append(dict(id=engine,operation='docking_'+engine,status='succeeded',configuration={},artifacts=[str(p) for p in files],artifact_manifest={str(p):input_digest(p) for p in files}))
        folder=self.root/'consensus';folder.mkdir();table=folder/'consensus.csv'
        row=dict(molecule_chembl_id='M1',canonical_smiles='CCO',receptor_id=targets[0],vina=.8,dock6=.9,normalized_score=.85,
            target_count=3,target_scores=json.dumps([dict(receptor_id=t,vina_raw=-4-n,dock6_raw=-4-n) for n,t in enumerate(targets)]))
        with table.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader();writer.writerow(row)
        stages.append(dict(id=self.sid,operation='consensus',status='succeeded',configuration={'bindings':{
            'base_vina_path':{'stage':'vina'},'base_dock6_path':{'stage':'dock6'}}},
            input_files={'base_vina_path':[str(tables['vina'])],'base_dock6_path':[str(tables['dock6'])]},artifacts=[str(table)]))
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        service=DockingResults(self.store)
        version=service.page(self.token,self.pid,self.rid,self.sid,str(table))['version']
        return service,dict(table=str(table),index=0,version=version),tables

    def test_each_engine_and_footprint_uses_selected_target(self):
        service,selection,tables=self.fixture()
        choices=service.target_choices(self.token,self.pid,self.rid,self.sid,**selection)
        self.assertEqual([x['label'].split(' · ')[0] for x in choices],['AChE','BChE','BACE1'])
        self.assertTrue(all(x['available'] for x in choices))
        for n,choice in enumerate(choices):
            target=choice['id']
            for engine in ('vina','dock6'):
                spec=service.scene(self.token,self.pid,self.rid,self.sid,**selection,column=engine+'_pose',target=target)
                self.assertEqual(Path(spec['layers'][0]['path']).name,target.split('_')[0]+'.pdb')
                self.assertEqual(Path(spec['layers'][1]['path']).name,target+'.pdbqt')
                self.assertEqual(Path(spec['layers'][1]['path']).parent.name,engine)
            self.assertEqual(service.footprint_pdf(self.token,self.pid,self.rid,self.sid,**selection,target=target),b'%PDF-1.4\n'+target.encode())

    def test_target_survives_browser_capabilities(self):
        service,selection,tables=self.fixture();selection['target']='2BBB_A'
        key=self.viewers.issue_result(self.token,self.pid,self.rid,self.sid,'docking',dict(selection,column='vina_pose'))
        response=self.viewers.response(PREFIX+'/'+key+'/scene')
        self.assertEqual(response[0],200)
        self.assertIn('2BBB_A',response[2].decode())
        key=self.viewers.issue_document(self.token,self.pid,self.rid,self.sid,selection)
        self.assertEqual(self.viewers.response(PREFIX+'/'+key)[2],b'%PDF-1.4\n2BBB_A')
        self.store.logout(self.token)
        self.assertEqual(self.viewers.response(PREFIX+'/'+key)[0],403)

    def test_invalid_target_column_or_changed_source_is_rejected(self):
        service,selection,tables=self.fixture()
        for target,column in [('unrelated','vina_pose'),('2BBB_A','invalid')]:
            with self.assertRaises(AccessDenied):service.scene(self.token,self.pid,self.rid,self.sid,**selection,target=target,column=column)
        tables['vina'].write_text(tables['vina'].read_text().replace('M1','OTHER'))
        with self.assertRaises(ValueError):service.scene(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A',column='vina_pose')
        choices=service.target_choices(self.token,self.pid,self.rid,self.sid,**selection)
        self.assertFalse(any(x['available'] for x in choices))

    def test_modified_pose_is_unavailable_but_existing_pdf_remains_usable(self):
        service,selection,tables=self.fixture()
        pose=tables['dock6'].parent/'2BBB_A.pdbqt'
        pose.write_text(pose.read_text().replace('13.000','93.000'))
        with self.assertRaises(ValueError):
            service.scene(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A',column='dock6_pose')
        self.assertEqual(service.footprint_pdf(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A'),b'%PDF-1.4\n2BBB_A')

    def test_reused_consensus_resolves_authorized_historical_producer(self):
        service,selection,tables=self.fixture()
        with self.store.connect() as db:
            stages=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
            db.execute('INSERT INTO runs(id,project_id,user_id,status,stages,created,updated) VALUES(?,?,?,?,?,?,?)',
                ('origin-run',self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps(stages[:2]),1,1))
            for producer in stages[:2]:producer['artifacts']=[];producer['artifact_manifest']={}
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        for engine in ('vina','dock6'):
            scene=service.scene(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A',column=engine+'_pose')
            self.assertIn('origin-run',scene['layers'][1]['path'])
        self.assertEqual(service.footprint_pdf(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A'),b'%PDF-1.4\n2BBB_A')
        with self.store.connect() as db:db.execute('DELETE FROM runs WHERE id=?',('origin-run',))
        with self.assertRaises(AccessDenied):
            service.scene(self.token,self.pid,self.rid,self.sid,**selection,target='2BBB_A',column='vina_pose')

    def test_large_target_chooser_does_not_parse_poses_before_selection(self):
        service,selection,tables=self.fixture()
        with patch.object(service,'scene',side_effect=AssertionError('eager pose validation')):
            choices=service.target_choices(self.token,self.pid,self.rid,self.sid,**selection,validate=False)
        self.assertEqual(len(choices),3)
        self.assertTrue(all(choice['available'] for choice in choices))



class TargetChooserTests(__import__('unittest').IsolatedAsyncioTestCase):
    def viewer(self,multi=True):
        from biomolexplorer.ui.docking_results import DockingResultsTable
        viewer=object.__new__(DockingResultsTable)
        viewer.token='token';viewer.project_id='p';viewer.run_id='r';viewer.stage_id='s'
        viewer.sequence=1;viewer.selector=SimpleNamespace(value='table');viewer.data={'version':'sha'}
        viewer.valid=lambda:True
        viewer.service=SimpleNamespace(target_choices=lambda:None)
        async def call(*args,**kwargs):return [dict(id='A',label='AChE',available=True),dict(id='B',label='BChE',available=True)]
        page=SimpleNamespace(web=False,show_dialog=lambda dialog:setattr(viewer,'dialog',dialog),width=1000,height=800,update=lambda:None)
        async def guard(action):return await action()
        viewer.ui=SimpleNamespace(page=page,call=call,guard=guard,preview_docking=AsyncMock())
        row={'id':'M1','index':0,'properties':{'target_count':2 if multi else 1}}
        return viewer,row

    async def test_multi_target_defers_pose_until_selection(self):
        viewer,row=self.viewer()
        with patch('biomolexplorer.ui.docking_results.fit_dialog'),patch('biomolexplorer.ui.docking_results.close_dialog'):
            await viewer.pose_preview(row,'vina_pose')
            viewer.ui.preview_docking.assert_not_awaited()
            selector=viewer.dialog.content.content.controls[1]
            selector.value='B';selector.on_select(None)
            await viewer.dialog.actions[1].on_click(None)
        self.assertEqual(viewer.ui.preview_docking.await_args.args[4]['target'],'B')

    async def test_single_target_opens_directly(self):
        viewer,row=self.viewer(False)
        await viewer.pose_preview(row,'dock6_pose')
        self.assertFalse(hasattr(viewer,'dialog'))
        self.assertNotIn('target',viewer.ui.preview_docking.await_args.args[4])

    async def test_footprint_and_dock6_choose_target_and_cancel_is_inert(self):
        viewer,row=self.viewer()
        with patch('biomolexplorer.ui.docking_results.fit_dialog'),patch('biomolexplorer.ui.docking_results.close_dialog'):
            viewer.footprint=AsyncMock()
            await viewer.choose_target(row,'footprint','dock6_pose')
            selector=viewer.dialog.content.content.controls[1]
            selector.value='B';selector.on_select(None)
            await viewer.dialog.actions[1].on_click(None)
            viewer.footprint.assert_awaited_once_with(row,'B')
            await viewer.choose_target(row,'pose','dock6_pose')
            viewer.dialog.actions[0].on_click(None)
            await viewer.dialog.actions[1].on_click(None)
            viewer.ui.preview_docking.assert_not_awaited()

    async def test_browser_target_buttons_keep_direct_navigation(self):
        viewer,row=self.viewer();viewer.ui.page.web=True
        viewer.ui.docking_view_url=lambda *args,**kwargs:'https://viewer.example/'+args[4]['target']
        with patch('biomolexplorer.ui.docking_results.fit_dialog'),patch('biomolexplorer.ui.docking_results.close_dialog'):
            await viewer.pose_preview(row,'dock6_pose')
            selector=viewer.dialog.content.content.controls[1]
            self.assertTrue(selector.enable_filter and selector.enable_search)
            self.assertTrue(viewer.dialog.actions[1].disabled)
            selector.value='B';selector.on_select(None)
            self.assertEqual(viewer.dialog.actions[1].url.url,'https://viewer.example/B')
            selector=viewer.dialog.content.content.controls[1]
            selector.value='B';selector.on_select(None)
            await viewer.dialog.actions[1].on_click(None)
        viewer.ui.preview_docking.assert_not_awaited()

    async def test_hundreds_of_targets_use_one_bounded_searchable_combobox(self):
        viewer,row=self.viewer()
        viewer.ui.call=AsyncMock(return_value=[dict(id=str(n),label='Target '+str(n),available=True) for n in range(500)])
        with patch('biomolexplorer.ui.docking_results.close_dialog'):
            await viewer.pose_preview(row,'vina_pose')
            selector=viewer.dialog.content.content.controls[1]
            self.assertEqual(len(selector.options),500)
            self.assertEqual(selector.menu_height,280)
            self.assertTrue(selector.enable_filter and selector.enable_search)
            self.assertEqual(len(viewer.dialog.actions),2)
            viewer.ui.preview_docking.assert_not_awaited()
            selector.value='499';selector.on_select(None)
            await viewer.dialog.actions[1].on_click(None)
        self.assertEqual(viewer.ui.preview_docking.await_args.args[4]['target'],'499')
