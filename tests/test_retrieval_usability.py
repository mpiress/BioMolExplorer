"""Curation persists scientific metadata and does not weaken access or cache checks."""
import asyncio
import csv
import io
import json
import tempfile
import time
import unittest
from types import SimpleNamespace
from uuid import uuid4
from biomolexplorer.catalog import new_stage
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import WorkspaceStore,AccessDenied

try:
    import flet as ft
except ImportError:ft=None

PDB=(b'HEADER    TEST\n'
     b'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n'
     b'ATOM      2  CA  GLY A   2       2.000   1.000   0.000  1.00 20.00           C  \n'
     b'HETATM    3  C1  LIG A 101       1.000   2.000   3.000  1.00 20.00           C  \n'
     b'HETATM    4  P   ATP B 102       3.000   2.000   1.000  1.00 20.00           P  \nEND\n')


class PDBCurationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@example.org','owner-password')
        self.guest=self.store.register('Guest','guest@example.org','guest-password')
        self.pid=self.store.create_project(self.token,'PDB results')['id'];self.rid=uuid4().hex
        self.config=new_stage('retrieve_structures');self.sid=self.config['id']
        root=self.store.project_dir(self.pid)/'runs'/self.rid/'artifacts'/'PDB';root.mkdir(parents=True)
        self.path=root/'1ABC.pdb';self.path.write_bytes(PDB)
        self.csv=root/'pdb_codes.csv';self.csv.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,LIG,101,A,2.0\n2XYZ,ATP,3,B,1.5\n')
        report=root/'retrieval_report.json';report.write_text('{}')
        files=[str(p) for p in (self.path,self.csv,report)]
        self.stage=dict(id=self.sid,name=self.config['name'],operation='retrieve_structures',configuration=self.config,
                        status='succeeded',artifacts=files,artifact_manifest=artifact_manifest(files))
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(self.rid,self.pid,self.store.user(self.token)['id'],
                'awaiting_input',json.dumps([self.stage]),time.time(),time.time(),None))
        self.service=ResultFiles(self.store);self.args=(self.token,self.pid,self.rid,self.sid,str(self.path))

    def tearDown(self):self.temp.cleanup()
    def context(self):return self.service.pdb_ligands(*self.args)
    def save(self,rows,revision=None):return self.service.set_pdb_ligands(*self.args,rows,revision or self.context()['revision'])

    def test_internal_files_hidden_but_preserved_for_pipeline(self):
        self.assertEqual([f['name'] for f in self.service.files(*self.args[:4])],['1ABC.pdb'])
        from biomolexplorer.artifact_choices import choices
        self.assertEqual(list(choices(self.stage)),['PDB/1ABC.pdb'])
        self.assertTrue(self.csv.exists())
        self.assertEqual({r['LIGAND'] for r in self.context()['available_ligands']},{'ATP','LIG'})
        self.assertEqual(len(self.store.get_run(self.token,self.rid)['stages'][0]['artifacts']),3)

    def test_replace_remove_add_preserves_other_structures_and_manifest(self):
        before=self.context();old_manifest=self.stage['artifact_manifest']
        result=self.save([{'LIGAND':'atp','RESNUM':'102','CHAIN':'B'}],before['revision'])
        self.assertEqual([(r['LIGAND'],r['RESNUM'],r['CHAIN']) for r in result['ligands']],[('ATP','102','B')])
        rows=list(csv.DictReader(io.StringIO(self.csv.read_text())))
        self.assertIn({'PDB_CODE':'2XYZ','LIGAND':'ATP','RESNUM':'3','CHAIN':'B','RESOLUTION':'1.5'},rows)
        self.assertNotEqual(old_manifest,self.store.get_run(self.token,self.rid)['stages'][0]['artifact_manifest'])
        self.assertEqual(self.save([])['ligands'],[])
        self.assertEqual(len(list(csv.DictReader(io.StringIO(self.csv.read_text())))),1)
        with self.store.connect() as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM project_history WHERE project_id=? AND action='curate_pdb_ligands'",(self.pid,)).fetchone()[0],2)

    def test_invalid_residue_stale_revision_and_active_run_do_not_mutate(self):
        original=self.csv.read_bytes();revision=self.context()['revision']
        with self.assertRaisesRegex(ValueError,'não foi encontrado'):
            self.save([{'LIGAND':'ATP','RESNUM':'999','CHAIN':'B'}])
        self.assertEqual(self.csv.read_bytes(),original)
        self.save([])
        with self.assertRaisesRegex(ValueError,'Reabra'):self.save([],revision)
        with self.store.connect() as db:db.execute("UPDATE runs SET status='running' WHERE id=?",(self.rid,))
        with self.assertRaisesRegex(ValueError,'pausar'):self.save([])

    def test_revision_describes_the_bytes_shown_even_when_read_races_with_write(self):
        from hashlib import sha256
        from unittest.mock import patch
        original=self.csv.read_bytes();read=self.store.read_file
        def raced(token,pid,path,limit=None):
            value=read(token,pid,path,limit)
            if path==str(self.csv):self.csv.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,ATP,102,B,2.0\n')
            return value
        with patch.object(self.store,'read_file',side_effect=raced):context=self.context()
        self.assertEqual(context['revision'],sha256(original).hexdigest())
        with self.assertRaisesRegex(ValueError,'Reabra'):self.save([],context['revision'])

    def test_cannot_read_other_projects_or_write_as_viewer(self):
        with self.assertRaises(AccessDenied):self.service.pdb_ligands(self.guest,*self.args[1:])
        self.store.invite(self.token,self.pid,'guest@example.org','viewer');self.store.accept_invitation(self.guest,self.pid)
        self.assertEqual(self.service.pdb_ligands(self.guest,*self.args[1:])['pdb_id'],'1ABC')
        with self.assertRaises(AccessDenied):self.service.set_pdb_ligands(self.guest,*self.args[1:],[],self.context()['revision'])
        with self.assertRaises(AccessDenied):self.service.pdb_ligands(*self.args[:4],str(self.csv))

    def test_curated_ligands_feed_selected_structure_and_change_input_cache(self):
        from biomolexplorer.pipeline import PipelineService
        service=object.__new__(PipelineService);service.store=self.store
        stage=new_stage('prepare_structures');stage['parameters']['target']='Selected'
        refs=[{'stage':self.sid,'selector':'1ABC.pdb'}];results={self.sid:self.stage['artifacts']}
        before,_,_=service._materialize_inputs(self.pid,stage,'base_input_path',refs,results)
        self.save([{'LIGAND':'ATP','RESNUM':102,'CHAIN':'B'}])
        after,_,_=service._materialize_inputs(self.pid,stage,'base_input_path',refs,results)
        self.assertNotEqual(before,after)
        rows=list(csv.DictReader(io.StringIO((after/'Selected'/'pdb_codes.csv').read_text())))
        self.assertEqual([(r['PDB_CODE'],r['LIGAND']) for r in rows],[('1ABC','ATP')])
        self.save([])
        with self.assertRaisesRegex(ValueError,'não possuem registros PDB'):
            service._materialize_inputs(self.pid,stage,'base_input_path',refs,results)

    @unittest.skipIf(ft is None,'Install ui extra')
    def test_file_selection_adds_pdb_actions_and_keeps_metadata_out_of_choices(self):
        from biomolexplorer.ui.file_selection import FileSelection
        config=new_stage('prepare_structures');config['bindings']={'base_input_path':{'stage':self.sid,'selector':'auto'}}
        pending={'configuration':config,'name':config['name']};seen=[]
        def actions(stage,path):
            import flet as ft
            seen.append(path.name);return [ft.IconButton(ft.Icons.SCIENCE_OUTLINED)]
        form=FileSelection({'stages':[self.stage]},pending,file_actions=actions)
        self.assertEqual(seen,['1ABC.pdb'])
        checks=form.rows['base_input_path'];self.assertEqual(len(checks),1)
        checks[0][0].value=True
        self.assertEqual(form.read()['bindings']['base_input_path']['selector'],'1ABC.pdb')

    @unittest.skipIf(ft is None,'Install ui extra')
    def test_ligand_dialog_removes_and_adds_existing_residue_before_saving(self):
        from biomolexplorer.ui.pdb_results import PDBActions
        dialogs=[];messages=[]
        async def call(function,*args):return function(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(call=call,guard=guard,notify=messages.append,editing_stage=True,
                           page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append))
        results=SimpleNamespace(ui=ui,token=self.token,project_id=self.pid,run_id=self.rid,
                                stage=self.stage,service=self.service,valid=lambda:True)
        actions=PDBActions(results,True)
        asyncio.run(actions.open_ligands({'path':str(self.path)}))
        dialog=dialogs[-1];controls=dialog.content.content.controls;listing=controls[1]
        listing.controls[0].controls[-1].on_click(None)
        controls[2].value='0';controls[3].on_click(None)
        self.assertEqual(listing.controls[0].controls[0].value,'ATP')
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual([r['LIGAND'] for r in self.context()['ligands']],['ATP'])
        self.assertIn('histórico',messages[0]);self.assertTrue(ui.editing_stage)

    @unittest.skipIf(ft is None,'Install ui extra')
    def test_local_preview_and_structure_specific_buttons(self):
        from biomolexplorer.pdb_view import viewer_document
        document=viewer_document('1ABC.pdb')
        self.assertIn('3Dmol-min.js',document);self.assertIn('pdb-viewer.js',document)
        from biomolexplorer.ui.pdb_results import PDBActions
        results=SimpleNamespace(ui=SimpleNamespace(),valid=lambda:True)
        actions=PDBActions(results,True).actions({'name':self.path.name,'path':str(self.path)})
        self.assertEqual(len(actions),3);self.assertEqual(actions[2].url,'https://www.rcsb.org/structure/1ABC')


@unittest.skipIf(ft is None,'Install ui extra')
class RetrievalControlsTests(unittest.TestCase):
    def test_activity_selection_survives_search_and_paging_and_custom_values(self):
        from biomolexplorer.ui.activity_measures import ActivityMeasures
        page=SimpleNamespace(update=lambda:None)
        picker=ActivityMeasures(page,['Ki','IC50'])
        self.assertGreater(len(picker.types),6000)
        picker.search.value='Inhibition';picker.filter(None)
        check=picker.grid.controls[0].content;check.value=True;check.on_change(SimpleNamespace(control=check))
        picker.move(1);picker.search.value='IC50';picker.filter(None)
        check=picker.grid.controls[0].content;check.value=False;check.on_change(SimpleNamespace(control=check))
        picker.custom.value='Future assay measure';picker.add(None)
        self.assertIn('Ki',picker.values());self.assertNotIn('IC50',picker.values());self.assertIn('Future assay measure',picker.values())
        self.assertEqual(picker.grid.controls[0].col['md'],4)

    def test_organism_and_assay_dropdown_values_serialize_scientific_codes(self):
        import flet as ft
        from biomolexplorer.ui.guided import GuidedForm
        ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None))
        form=GuidedForm(ui,new_stage('retrieve_compounds'),[],True)
        organism=next(c for c in form.filter_tiles['target'].controls if isinstance(c,ft.Dropdown) and c.label=='Organismo')
        self.assertIn('Homo sapiens',organism.hint_text);self.assertTrue(organism.editable)
        organism.text='Arabidopsis thaliana';organism.on_text_change(SimpleNamespace(control=organism))
        assay=next(c for c in form.filter_tiles['bioactivity'].controls if isinstance(c,ft.Dropdown) and c.label=='Tipo de ensaio')
        self.assertEqual({c.key for c in assay.options},{'','B','F','A','T','P','U'})
        assay.value='F';result=form.read()
        self.assertEqual(json.loads(result['templates']['crawlers/target.json'])['organism'],'Arabidopsis thaliana')
        self.assertEqual(json.loads(result['templates']['crawlers/bioactivity.json'])['assay_type'],'F')
        from biomolexplorer.ui.localization import LocalizedPage
        localized=LocalizedPage(ui.page,'en');localized.localize([organism,assay])
        self.assertEqual(next(c.text for c in assay.options if c.key=='B'),'B - Binding')
        organism.value='';organism.text='Any';organism.on_select(SimpleNamespace(control=organism))
        organism.on_text_change(SimpleNamespace(control=organism))
        self.assertNotIn('organism',json.loads(form.template_readers['crawlers/target.json']()))
