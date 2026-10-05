"""Inline result exploration, file curation and automatic fingerprint selection."""
import asyncio
import csv
import io
import json
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.fingerprint_selection import generated_kind
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import AccessDenied,WorkspaceStore

try:
    import flet as ft
    from biomolexplorer.ui.guided import GuidedForm
    from biomolexplorer.ui.stage_results import StageResults
    from biomolexplorer.ui.file_table import FileTable
except ImportError:
    GuidedForm=StageResults=FileTable=None


class ResultsDisplayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@example.org','owner-password')
        self.guest=self.store.register('Guest','guest@example.org','guest-password')
        self.pid=self.store.create_project(self.token,'Results')['id']
        self.rid=uuid4().hex;self.retrieval=new_stage('retrieve_compounds');self.admet=new_stage('admet')
        root=self.store.project_dir(self.pid)/'runs'/self.rid/'artifacts';root.mkdir(parents=True)
        self.compounds=root/'compounds'/'CHEMBL220'/'compounds.csv';self.compounds.parent.mkdir(parents=True)
        self.compounds.write_text('molecule_chembl_id,canonical_smiles\n'+''.join(f'MOL{i},CCO\n' for i in range(31)))
        self.admet_csv=root/'admet.csv'
        self.admet_csv.write_text('molecule_chembl_id,canonical_smiles,TPSA,WLOGP,BBB,HIA,source\n'
            'A,CCO,20,1,BBB+,HIA+,ChEMBL\nB,CCC,10,4,BBB-,HIA+,PubChem\nC,O,150,0,BBB-,HIA-,Imported\n')
        internal=root/'internal.json';internal.write_text('{}')
        stages=[]
        for config,files in [(self.retrieval,[str(self.compounds)]),(self.admet,[str(self.admet_csv),str(internal)])]:
            stages.append(dict(id=config['id'],name=config['name'],operation=config['operation'],configuration=config,
                status='succeeded',artifacts=files,artifact_manifest=artifact_manifest(files)))
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(self.rid,self.pid,self.store.user(self.token)['id'],
                'succeeded',json.dumps(stages),time.time(),time.time(),None))
        self.service=ResultFiles(self.store)

    def tearDown(self):self.temp.cleanup()

    def ui(self):
        dialogs=[];saved=[]
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):await action()
        async def save_file(**kwargs):saved.append(kwargs)
        return SimpleNamespace(token=self.token,current=self.store.project(self.token,self.pid),store=self.store,
            page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append,width=1440),call=call,guard=guard,
            picker=SimpleNamespace(save_file=save_file),dialogs=dialogs,saved=saved,notify=lambda s:None,
            service=SimpleNamespace(dock6_path=None))

    def test_admet_subsets_and_exports_match_selected_points(self):
        self.assertEqual([f['name'] for f in self.service.files(self.token,self.pid,self.rid,self.admet['id'])],['admet.csv'])
        for subset,ids in [('all',{'A','B','C'}),('BBB+',{'A'}),('BBB-',{'B','C'}),('HIA+',{'A','B'})]:
            with self.subTest(subset=subset):
                args=(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv),subset)
                model=self.service.admet_model(*args)
                self.assertEqual({n['id'] for n in model['nodes']},ids)
                rows=list(csv.DictReader(io.StringIO(self.service.admet_csv(*args).decode())))
                self.assertEqual({r['molecule_chembl_id'] for r in rows},ids)
                self.assertTrue(all('source' in r for r in rows))
        self.assertTrue(self.service.admet_png(*args).startswith(b'\x89PNG'))

    @unittest.skipIf(StageResults is None,'Install ui extra')
    def test_exclusion_report_is_listed_and_downloadable_for_admet(self):
        from biomolexplorer.ui.app import WorkspaceUI
        path=self.admet_csv.parent/'molecule_exclusions.json'
        path.write_text(json.dumps({'version':1,'excluded_records':1,'records':[{'identifiers':['BAD'],'reason':'SMILES inválido'}]}))
        run=self.store.get_run(self.token,self.rid)
        stage=run['stages'][1];stage['artifacts'].append(str(path));stage['excluded_records']=1
        stage['artifact_manifest']=artifact_manifest(stage['artifacts'])
        with self.store.connect() as db:
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),self.rid))
        self.assertIn('molecule_exclusions.json',[f['name'] for f in self.service.files(self.token,self.pid,self.rid,self.admet['id'])])
        fixture=self.ui();ui=object.__new__(WorkspaceUI)
        ui.token=self.token;ui.current=fixture.current;ui.store=self.store;ui.call=fixture.call;ui.picker=fixture.picker
        asyncio.run(ui.download_stage_report(self.pid,str(path)))
        self.assertEqual(fixture.saved[-1]['src_bytes'],path.read_bytes())

    def test_empty_admet_subset_is_valid_and_preserves_headers(self):
        self.admet_csv.write_text('molecule_chembl_id,canonical_smiles,TPSA,WLOGP,BBB,HIA\nA,CCO,20,1,BBB+,HIA+\n')
        args=(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv),'BBB-')
        self.assertEqual(self.service.admet_model(*args)['nodes'],[])
        self.assertIn('canonical_smiles',self.service.admet_csv(*args).decode())
        self.assertTrue(self.service.admet_png(*args).startswith(b'\x89PNG'))

    def test_unclassified_uploaded_admet_uses_backend_classifier(self):
        self.admet_csv.write_text('molecule_chembl_id,canonical_smiles,TPSA,WLOGP,MW,HBD,RB\nA,CCO,20,1,,,\n')
        model=self.service.admet_model(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv),'BBB+')
        self.assertEqual(model['nodes'][0]['properties']['BBB'],'BBB+')
        self.assertEqual(model['nodes'][0]['properties']['HIA'],'HIA+')

    def test_deletion_is_audited_invalidates_cache_and_rollback_restores_bytes(self):
        original=self.admet_csv.read_bytes()
        self.service.remove(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv))
        self.assertFalse(self.admet_csv.exists())
        stage=self.store.get_run(self.token,self.rid)['stages'][1]
        self.assertTrue(stage['outputs_removed']);self.assertNotIn(str(self.admet_csv),stage['artifacts'])
        event=self.store.history(self.token,self.pid)[0];self.assertEqual(event['action'],'remove_file')
        self.store.rollback(self.token,self.pid,event['id'])
        self.assertEqual(self.admet_csv.read_bytes(),original)
        self.assertFalse(self.store.get_run(self.token,self.rid)['stages'][1].get('outputs_removed'))

    def test_failed_file_deletion_restores_files_and_metadata(self):
        with patch('biomolexplorer.project_state.record',side_effect=OSError('disk failure')):
            with self.assertRaises(OSError):self.service.remove(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv))
        self.assertTrue(self.admet_csv.exists())
        self.assertIn(str(self.admet_csv),self.store.get_run(self.token,self.rid)['stages'][1]['artifacts'])

    def test_permissions_registered_files_and_active_run_are_enforced(self):
        self.store.invite(self.token,self.pid,'guest@example.org','viewer');self.store.accept_invitation(self.guest,self.pid)
        args=(self.pid,self.rid,self.admet['id'],str(self.admet_csv))
        with self.assertRaises(AccessDenied):self.service.remove(self.guest,*args)
        with self.assertRaises(AccessDenied):self.service.remove(self.token,*args[:-1],str(self.compounds))
        with self.store.connect() as db:db.execute('UPDATE runs SET status=? WHERE id=?',('running',self.rid))
        with self.assertRaisesRegex(ValueError,'Aguarde'):self.service.remove(self.token,*args)
        self.assertTrue(self.admet_csv.exists())

    def test_input_file_deletion_requires_unlinking_and_can_be_rolled_back(self):
        ticket=self.store.prepare_upload(self.token,self.pid,'custom.csv','compounds')
        (self.store.staging/ticket).write_text('molecule_chembl_id,canonical_smiles\nOWN,CCO\n')
        asset=self.store.finish_upload(self.token,ticket)
        path=self.store.asset_path(self.token,self.pid,asset);original=path.read_bytes()
        stage=new_stage('admet');stage['bindings']['base_input_path']={'asset':asset,'selector':'auto'}
        self.store.save_pipeline(self.token,self.pid,[stage],0)
        with self.assertRaisesRegex(ValueError,'associado'):self.service.remove_asset(self.token,self.pid,asset)
        self.store.save_pipeline(self.token,self.pid,[],1)
        self.service.remove_asset(self.token,self.pid,asset)
        self.assertFalse(path.exists());self.assertEqual(self.store.assets(self.token,self.pid),[])
        event=self.store.history(self.token,self.pid)[0]
        self.store.rollback(self.token,self.pid,event['id'])
        self.assertEqual(self.store.asset_path(self.token,self.pid,asset).read_bytes(),original)

    @unittest.skipIf(StageResults is None,'Install ui extra')
    def test_retrieval_expansion_loads_inline_25_rows_and_page_size_control(self):
        ui=self.ui();stage=self.store.get_run(self.token,self.rid)['stages'][0]
        results=StageResults(ui,self.pid,self.rid,stage,True)
        self.assertFalse(results.loaded)
        asyncio.run(results.expand(SimpleNamespace(data='true')))
        self.assertTrue(results.loaded);self.assertTrue(results.child.inline)
        self.assertEqual(len(results.child.table.rows),25)
        asyncio.run(results.child.next.on_click(None));self.assertEqual(len(results.child.table.rows),6)
        results.child.page_size.value='50';asyncio.run(results.child.change_size(None))
        self.assertEqual(len(results.child.table.rows),31);self.assertEqual(results.child.offset,0)
        results.close();self.assertFalse(results.child.valid())

    @unittest.skipIf(StageResults is None,'Install ui extra')
    def test_admet_inline_filter_hover_popup_and_download_callbacks(self):
        ui=self.ui();stage=self.store.get_run(self.token,self.rid)['stages'][1]
        results=StageResults(ui,self.pid,self.rid,stage,True);asyncio.run(results.load())
        self.assertEqual(len(results.subset.options),4)
        results.subset.value='BBB+';asyncio.run(results.change_admet(None))
        viewer=results.child;self.assertEqual(set(viewer.points),{'A'})
        x,y=viewer.points['A'];event=SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))
        viewer.hover(event);self.assertEqual(viewer.tip_text.value,'A')
        asyncio.run(viewer.click(event));self.assertIn('molécula 2D',ui.dialogs[-1].title.value)
        self.assertTrue(any(isinstance(c,ft.Image) for c in ui.dialogs[-1].content.content.controls))
        asyncio.run(viewer.actions[1].on_click(None))
        self.assertEqual(ui.saved[-1]['file_name'],'admet_BBB+.csv')
        self.assertEqual(len(list(csv.DictReader(io.StringIO(ui.saved[-1]['src_bytes'].decode())))),1)

    @unittest.skipIf(FileTable is None,'Install ui extra')
    def test_file_tables_have_bounded_pages_and_actions_on_right(self):
        ui=self.ui();table=FileTable(ui,self.pid,[{'name':f'file{i}.csv','size':100} for i in range(77)])
        self.assertEqual(len(table.table.rows),25)
        table.move(1);self.assertEqual(table.offset,25)
        table.page_size.value='50';table.resize(None);self.assertEqual(len(table.table.rows),50)
        self.assertEqual(table.table.columns[-1].label.value,'Ações')

    @unittest.skipIf(StageResults is None,'Install ui extra')
    def test_runs_view_keeps_open_stages_and_refreshes_changed_manifests(self):
        from biomolexplorer.ui.app import WorkspaceUI
        fake=self.ui();ui=object.__new__(WorkspaceUI);ui.__dict__.update(fake.__dict__)
        ui.current_run=ui.run_progress=None;ui.polling=SimpleNamespace(done=lambda:False)
        asyncio.run(ui.runs_view(True))
        key=(self.rid,self.retrieval['id']);old=ui.stage_results[key]
        self.assertFalse(old.loaded)
        asyncio.run(old.expand(SimpleNamespace(data=True)))
        asyncio.run(ui.runs_view(True));self.assertIs(ui.stage_results[key],old)
        self.compounds.write_text(self.compounds.read_text()+'NEW,CCN\n')
        with self.store.connect() as db:
            stages=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
            stages[0]['artifact_manifest']=artifact_manifest(stages[0]['artifacts'])
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        asyncio.run(ui.runs_view(True))
        new=ui.stage_results[key];self.assertIsNot(new,old);self.assertFalse(old.active)
        self.assertTrue(new.expanded);asyncio.run(new.load());self.assertEqual(new.child.data['total'],32)


@unittest.skipIf(GuidedForm is None,'Install ui extra')
class FingerprintFormTests(unittest.TestCase):
    def ui(self,stages):
        return SimpleNamespace(current={'pipeline':stages},service=SimpleNamespace(dock6_path=None),page=SimpleNamespace(update=lambda:None))

    def test_single_fingerprint_and_conditional_morgan_fields(self):
        stage=new_stage('fingerprints');form=GuidedForm(self.ui([stage]),stage,[],True);form.layout()
        self.assertEqual(sum(stage['parameters'][k] for k in ('morgan','maccs','pharmacophore')),1)
        form.fingerprint_choice.value='maccs';form.fingerprint_choice.on_select(SimpleNamespace())
        self.assertFalse(form.field_controls['radius'].visible);self.assertFalse(form.conditional_cells['radius'].visible)
        form.field_controls['radius'].value='not a number'
        data=form.read()['parameters'];self.assertTrue(data['maccs']);self.assertFalse(data['morgan']);self.assertNotIn('radius',data)
        form.fingerprint_choice.value='morgan';form.fingerprint_choice.on_select(SimpleNamespace())
        self.assertTrue(form.conditional_cells['morgan_n_bits'].visible)

    def test_similarity_follows_generated_selection_and_unlocks_only_for_custom_input(self):
        producer=new_stage('fingerprints');producer['parameters'].update(morgan=False,maccs=True)
        stage=new_stage('similarity');stage['bindings']['base_input_path']={'stage':producer['id'],'selector':'auto'}
        ui=self.ui([producer,stage]);form=GuidedForm(ui,stage,[{'id':'own','name':'custom.csv','kind':'fingerprints'}],True)
        control=form.field_controls['fingerprint'];self.assertTrue(control.disabled);self.assertEqual(control.value,'maccs')
        self.assertEqual(form.read()['parameters']['fingerprint'],'maccs')
        source=form.input_editors['base_input_path'].rows[0][0]
        source.value='asset:own';source.on_select(None);self.assertFalse(control.disabled)
        control.value='pharmacophore';self.assertEqual(form.read()['parameters']['fingerprint'],'pharmacophore')

    def test_old_multi_algorithm_results_resolve_by_filename_and_reject_mixed_types(self):
        producer=new_stage('fingerprints');producer['parameters'].update(morgan=True,maccs=True,pharmacophore=True)
        stage=new_stage('similarity');stage['bindings']['base_input_path']={'stage':producer['id'],'selector':'maccs_molecules.csv'}
        self.assertEqual(generated_kind(stage,[producer,stage]),('maccs',False))
        stage['bindings']['base_input_path']={'sources':[{'stage':producer['id'],'selector':'morgan_molecules.csv'},
            {'stage':producer['id'],'selector':'maccs_molecules.csv'}]}
        with self.assertRaisesRegex(ValueError,'mesmo tipo'):generated_kind(stage,[producer,stage])

    def test_backend_overrides_stale_similarity_algorithm(self):
        with tempfile.TemporaryDirectory() as root:
            store=WorkspaceStore(root);token=store.register('Owner','owner@example.org','owner-password')
            pid=store.create_project(token,'Fingerprint')['id'];producer=new_stage('fingerprints');stage=new_stage('similarity')
            producer['parameters'].update(morgan=False,maccs=True)
            stage['bindings']['base_input_path']={'stage':producer['id'],'selector':'maccs_compounds.csv'}
            store.save_pipeline(token,pid,[producer,stage],0)
            path=store.project_dir(pid)/'maccs_compounds.csv';path.write_text('molecule_chembl_id,fingerprint\nA,"[0,1]"\n')
            service=PipelineService(store)
            try:
                params=service._resolve(pid,store.user(token)['id'],stage,{producer['id']:[str(path)]})
                self.assertEqual(params['fingerprint'],'maccs')
                self.assertEqual(params['filename'],'selected_fingerprints.csv')
                self.assertTrue((Path(params['base_input_path'])/params['filename']).is_file())
            finally:service.close()
