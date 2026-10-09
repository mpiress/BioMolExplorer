"""Stage handoffs show current, compatible files and require explicit choices."""
import asyncio
import copy
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock, patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
import test_pipeline_execution as execution

try:
    import flet as ft
    from biomolexplorer.ui.app import WorkspaceUI
    from biomolexplorer.ui.file_selection import FileSelection
except ImportError:
    FileSelection = None


@unittest.skipIf(FileSelection is None, 'Install the ui extra')
class FileSelectionTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_popup_resumes_same_run_and_only_passes_checked_files(self):
        source=new_stage('retrieve_compounds')
        stage=new_stage('admet')
        stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        manager,calls=self.immediate_manager()
        class MultipleFiles(manager):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    extra=manager.config.workspace/'chosen.csv'
                    extra.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL2,CCC\n')
                    report=manager.config.workspace/'report.csv'
                    report.write_text('status\nok\n')
                    job['result']['artifacts'] += [str(extra),str(report)]
                return job
        with patch('biomolexplorer.pipeline.JobManager',MultipleFiles):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                dialogs=[]
                def show(dialog):dialog.open=True;dialogs.append(dialog)
                page=SimpleNamespace(width=1440,height=1000,update=lambda:None,show_dialog=show,
                    run_task=lambda *args:SimpleNamespace(done=lambda:False))
                ui=WorkspaceUI(page,self.store,service)
                ui.token=self.token;ui.current=self.store.project(self.token,self.project_id)
                ui.draw_project=AsyncMock()
                async def local_call(function,*args,**kwargs):return function(*args,**kwargs)
                ui.call=local_call
                asyncio.run(ui.configure_pending(run['id']))
                asyncio.run(ui.configure_pending(run['id']))
                self.assertEqual(len(dialogs),1)
                dialog=dialogs[0]
                checks=[c for c in dialog.content.content.controls if isinstance(c,ft.Checkbox)]
                self.assertEqual({c.label for c in checks},{'chosen.csv','compounds.csv'})
                self.assertFalse(any(c.value for c in checks))
                # Closing the popup never advances the backend.
                dialog.actions[0].on_click(None)
                self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'awaiting_input')
                for check in checks:check.value=check.label=='chosen.csv'
                asyncio.run(dialog.actions[-1].on_click(None))
                finished=self.wait(run)
                self.assertEqual(finished['id'],run['id'])
                self.assertEqual(finished['status'],'succeeded',finished['error'])
                selected=finished['stages'][1]['input_files']['base_input_path']
                self.assertEqual([Path(p).name for p in selected],['chosen.csv'])
                stored=self.store.project(self.token,self.project_id)['pipeline'][1]
                self.assertEqual(stored['bindings']['base_input_path'],
                    {'stage':source['id'],'selector':'chosen.csv'})
                self.assertEqual(stored['input_processing'],'individual')
                self.assertEqual(len(calls),2)
                self.assertFalse(dialog.open)
            finally:service.close()

    def test_duplicate_filenames_have_distinct_selectors_and_empty_selection_is_rejected(self):
        source,stage=new_stage('retrieve_compounds'),new_stage('admet')
        stage['bindings']['base_input_path']={'stage':source['id']}
        root=self.store.project_dir(self.project_id)
        paths=[]
        for folder in ('first','second'):
            path=root/folder/'compounds.csv';path.parent.mkdir()
            path.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\n');paths.append(str(path))
        pending={'id':stage['id'],'name':stage['name'],'configuration':stage,'status':'awaiting_input'}
        run={'stages':[dict(id=source['id'],name=source['name'],status='succeeded',artifacts=paths),pending]}
        form=FileSelection(run,pending)
        rows=form.rows['base_input_path']
        self.assertEqual({ref['selector'] for _,ref in rows},{'first/compounds.csv','second/compounds.csv'})
        with self.assertRaisesRegex(ValueError,'Selecione ao menos'):form.read()
        rows[1][0].value=True
        configured=form.read()
        self.assertEqual(configured['bindings']['base_input_path']['selector'],'second/compounds.csv')
        self.assertEqual(stage['bindings']['base_input_path'],{'stage':source['id']})

    def test_existing_explicit_choices_are_checked_without_ambiguous_filename_matches(self):
        source,stage=new_stage('retrieve_compounds'),new_stage('fingerprints')
        root=self.store.project_dir(self.project_id)
        paths=[]
        for folder in ('first','second'):
            path=root/folder/'compounds.csv';path.parent.mkdir()
            path.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\n');paths.append(str(path))
        stage['bindings']['base_input_path']={'stage':source['id'],'selector':'second/compounds.csv'}
        stage['input_processing']='merge'
        pending={'id':stage['id'],'name':stage['name'],'configuration':stage,'status':'awaiting_input'}
        run={'stages':[dict(id=source['id'],name=source['name'],status='succeeded',artifacts=paths),pending]}
        form=FileSelection(run,pending)
        self.assertEqual([ref['selector'] for check,ref in form.rows['base_input_path'] if check.value],
                         ['second/compounds.csv'])
        self.assertEqual(form.read()['input_processing'],'merge')
        stage['bindings']['base_input_path']['selector']='compounds.csv'
        ambiguous=FileSelection(run,pending)
        self.assertFalse(any(check.value for check,_ in ambiguous.rows['base_input_path']))

    def test_deferring_then_configuring_keeps_files_and_processing_mode(self):
        source,stage=new_stage('retrieve_compounds'),new_stage('admet')
        stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                dialogs=[]
                def show(dialog):dialog.open=True;dialogs.append(dialog)
                page=SimpleNamespace(width=1440,height=1000,update=lambda:None,show_dialog=show,
                    run_task=lambda *args:SimpleNamespace(done=lambda:False))
                ui=WorkspaceUI(page,self.store,service)
                ui.token=self.token;ui.current=self.store.project(self.token,self.project_id)
                ui.draw_project=AsyncMock()
                async def local_call(function,*args,**kwargs):return function(*args,**kwargs)
                ui.call=local_call
                asyncio.run(ui.configure_pending(run['id']))
                first=dialogs[-1]
                checks=[c for c in first.content.content.controls if isinstance(c,ft.Checkbox)]
                checks[0].value=True
                processing=next(c for c in first.content.content.controls if isinstance(c,ft.Dropdown))
                processing.value='merge'
                first.actions[0].on_click(None)
                self.assertFalse(ui.editing_stage)
                asyncio.run(ui.configure_pending(run['id']))
                second=dialogs[-1]
                self.assertTrue(next(c for c in second.content.content.controls if isinstance(c,ft.Checkbox)).value)
                self.assertEqual(next(c for c in second.content.content.controls if isinstance(c,ft.Dropdown)).value,'merge')
                asyncio.run(second.actions[1].on_click(None))
                settings=dialogs[-1]
                self.assertTrue(ui.editing_stage)
                def walk(control):
                    yield control
                    child=getattr(control,'content',None)
                    if child is not None:yield from walk(child)
                    for child in getattr(control,'controls',[]) or []:yield from walk(child)
                dropdowns=[c for c in walk(settings.content) if isinstance(c,ft.Dropdown)]
                selected=next(c for c in dropdowns if c.label=='Resultado usado nesta entrada')
                self.assertEqual(selected.value,'compounds.csv')
                self.assertEqual(next(c for c in dropdowns if c.label=='Como processar os arquivos selecionados?').value,'merge')
                self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'awaiting_input')
                asyncio.run(settings.actions[-1].on_click(None))
                finished=self.wait(run)
                self.assertEqual(finished['status'],'succeeded',finished['error'])
                stored=self.store.project(self.token,self.project_id)['pipeline'][1]
                self.assertEqual(stored['bindings']['base_input_path']['selector'],'compounds.csv')
                self.assertEqual(stored['input_processing'],'merge')
            finally:service.close()

    def test_deferring_keeps_explicitly_unchecked_files_unchecked(self):
        source,stage=new_stage('retrieve_compounds'),new_stage('admet')
        path=self.store.project_dir(self.project_id)/'compounds.csv'
        path.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\n')
        stage['bindings']['base_input_path']={'stage':source['id'],'selector':'compounds.csv'}
        pending={'id':stage['id'],'name':stage['name'],'configuration':stage,'status':'awaiting_input'}
        run={'stages':[dict(id=source['id'],name=source['name'],status='succeeded',artifacts=[str(path)]),pending]}
        form=FileSelection(run,pending)
        form.rows['base_input_path'][0][0].value=False
        reopened=FileSelection(run,pending,state=form.state())
        self.assertFalse(reopened.rows['base_input_path'][0][0].value)
        with self.assertRaisesRegex(ValueError,'Selecione ao menos'):reopened.read()
        draft=reopened.read(require_selection=False)
        self.assertEqual(draft['bindings']['base_input_path'],{'stage':source['id'],'selector':'auto'})

    def test_poller_automatically_opens_once_for_each_pending_stage(self):
        dialogs=[]
        page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append)
        ui=WorkspaceUI(page,None,None)
        ui.token='token';ui.current={'id':'project','role':'editor'}
        ui.current_run='run';ui.tab='Pipeline'
        first,second=new_stage('admet'),new_stage('fingerprints')
        def waiting(stage):return {'id':'run','project_id':'project','status':'awaiting_input',
                                  'stages':[{'id':stage['id'],'status':'awaiting_input'}]}
        responses=iter([waiting(first),waiting(first),waiting(second),
                       {'id':'run','project_id':'project','status':'succeeded','stages':[]}])
        ui.store=SimpleNamespace(get_run=lambda *args:next(responses))
        async def local_call(function,*args):return function(*args)
        ui.call=local_call;ui.configure_pending=AsyncMock()
        with patch('biomolexplorer.ui.app.asyncio.sleep',new=AsyncMock()):
            asyncio.run(ui.poll_runs())
        self.assertEqual(ui.configure_pending.await_count,2)
        self.assertIsNone(ui.current_run)

    def test_prepared_selection_keeps_only_one_chain_and_its_metadata(self):
        source,stage=new_stage('prepare_structures'),new_stage('docking_vina')
        root=self.store.project_dir(self.project_id)/'prepared'/'MeuAlvo'
        prepared=root/'Prepared';prepared.mkdir(parents=True)
        atom='ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        for chain in ('A','B'):
            (prepared/f'1ABC_{chain}.dockprep.pdbqt').write_text(atom)
            (prepared/f'1ABC_{chain}.dockprep.mol2').write_text('@<TRIPOS>MOLECULE\nreceptor\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n')
            (prepared/f'1ABC_LIG_1_{chain}.lig.pdbqt').write_text(atom)
        (root/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,1,A\n1ABC,LIG,1,B\n')
        (prepared/'centers.csv').write_text('1ABC_LIG_1_A,1ABC_LIG_1_B\n0,1\n0,1\n0,1\n')
        files=[str(p) for p in root.rglob('*') if p.is_file()]
        service=PipelineService(self.store)
        try:
            path,_,selected=service._materialize_inputs(self.project_id,stage,'base_input_path',
                [{'stage':source['id'],'selector':'1ABC_A.dockprep.pdbqt'}],{source['id']:files})
            data=path/'MeuAlvo'
            self.assertTrue((data/'Prepared'/'1ABC_A.dockprep.pdbqt').exists())
            self.assertTrue((data/'Prepared'/'1ABC_A.dockprep.mol2').exists())
            self.assertFalse((data/'Prepared'/'1ABC_LIG_1_A.lig.pdbqt').exists())
            self.assertFalse((data/'Prepared'/'1ABC_B.dockprep.pdbqt').exists())
            self.assertFalse((data/'Prepared'/'1ABC_LIG_1_B.lig.pdbqt').exists())
            self.assertNotIn('1ABC,LIG,1,B',(data/'pdb_codes.csv').read_text())
            self.assertEqual((data/'Prepared'/'centers.csv').read_text().splitlines()[0],'1ABC_LIG_1_A')
            for operation in ('docking_vina','docking_dock6'):
                path,_,selected=service._materialize_inputs(self.project_id,new_stage(operation),'base_input_path',
                    [{'stage':source['id']}],{source['id']:files})
                self.assertFalse(any('.lig.' in p.name for p in selected))
                self.assertEqual({p.name for p in selected if p.suffix=='.pdbqt'},
                                 {'1ABC_A.dockprep.pdbqt','1ABC_B.dockprep.pdbqt'})
                self.assertTrue((path/'MeuAlvo'/'Prepared'/'1ABC_A.dockprep.mol2').exists())
            redocking=new_stage('redocking');redocking['parameters']['prepare_complex']=False
            path,_,_=service._materialize_inputs(self.project_id,redocking,'base_input_path',
                [{'stage':source['id'],'selector':'1ABC_A.dockprep.pdbqt'}],{source['id']:files})
            self.assertTrue((path/'MeuAlvo'/'Prepared'/'1ABC_LIG_1_A.lig.pdbqt').exists())
        finally:service.close()

    def test_docking_target_labels_and_configuration_only_offer_receptors(self):
        from biomolexplorer.catalog import operation_fields
        from biomolexplorer.flow import input_ports
        from biomolexplorer.ui.input_editor import InputEditor
        from biomolexplorer.ui.localization import Translator
        source=new_stage('prepare_structures')
        names={'Prepared/1ABC_A.dockprep.pdbqt','Prepared/1ABC_A.dockprep.mol2',
               'Prepared/1ABC_LIG_1_A.lig.pdbqt','Prepared/1ABC_LIG_1_A.lig.mol2','1ABC.pdb','pdb_codes.csv'}
        assets=[dict(id=str(i),name=name,kind='prepared_structures') for i,name in enumerate(sorted(names))]
        ui=SimpleNamespace(current={'pipeline':[source]},artifact_choices={source['id']:names},
                           page=SimpleNamespace(update=lambda:None))
        for operation in ('docking_vina','docking_dock6'):
            with self.subTest(operation=operation):
                stage=new_stage(operation)
                field=next(f for f in operation_fields(operation) if f['name']=='base_input_path')
                self.assertEqual(field['label'],'Alvo');self.assertEqual(Translator('en')(field['label']),'Target')
                self.assertEqual(input_ports(stage)[0]['label'],'Alvo')
                stage['bindings']['base_input_path']={'stage':source['id'],'selector':'Prepared/1ABC_LIG_1_A.lig.pdbqt'}
                editor=InputEditor(ui,stage,field,assets,True)
                self.assertEqual([o.key for o in editor.selectors('stage:'+source['id'])],
                                 ['auto','Prepared/1ABC_A.dockprep.pdbqt'])
                self.assertEqual([o.text for o in editor.source_options() if o.key.startswith('asset:')],
                                 ['Prepared/1ABC_A.dockprep.pdbqt'])
                # Old ligand selections must not be added back to the menu.
                self.assertEqual(editor.rows[0][1].value,'auto')
                editor.rows[0][1].value='Prepared/1ABC_LIG_1_A.lig.pdbqt'
                with self.assertRaisesRegex(ValueError,'receptor preparado'):editor.read()

    def test_docking_handoff_only_offers_prepared_receptors_and_labels_target(self):
        source=new_stage('prepare_structures')
        root=self.store.project_dir(self.project_id)/'files'/'Prepared'
        files=[str(root/name) for name in ('1ABC_A.dockprep.pdbqt','1ABC_A.dockprep.mol2',
               '1ABC_LIG_1_A.lig.pdbqt','1ABC_LIG_1_A.lig.mol2','1ABC.pdb','centers.csv')]
        for operation in ('docking_vina','docking_dock6'):
            with self.subTest(operation=operation):
                stage=new_stage(operation);stage['bindings']['base_input_path']={'stage':source['id']}
                pending={'configuration':stage,'name':stage['name']}
                run={'stages':[dict(id=source['id'],name=source['name'],operation=source['operation'],status='succeeded',artifacts=files)]}
                form=FileSelection(run,pending)
                self.assertEqual([r['selector'] for _,r in form.rows['base_input_path']],['1ABC_A.dockprep.pdbqt'])
                self.assertTrue(any(isinstance(c,ft.Text) and c.value=='Alvo' for c in form.control.controls))
                receptor={'id':'receptor','name':'1ABC_A.dockprep.pdbqt'}
                ligand={'id':'ligand','name':'1ABC_LIG_1_A.lig.pdbqt'}
                stage['bindings']['base_input_path']={'sources':[{'asset':'receptor'},{'asset':'ligand'}]}
                form=FileSelection(run,pending,assets=[receptor,ligand])
                self.assertEqual([r for _,r in form.rows['base_input_path']],[{'asset':'receptor'}])

    def test_docking_backend_rejects_ligands_as_targets(self):
        source=new_stage('prepare_structures')
        root=self.store.project_dir(self.project_id)/'prepared';root.mkdir()
        ligand=root/'1ABC_LIG_1_A.lig.pdbqt'
        ligand.write_text('ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n')
        asset=self.store.import_local_file(self.token,self.project_id,str(ligand),'prepared_structures')
        service=PipelineService(self.store)
        try:
            for operation in ('docking_vina','docking_dock6'):
                with self.subTest(operation=operation):
                    with self.assertRaisesRegex(ValueError,'receptor preparado'):
                        service._materialize_inputs(self.project_id,new_stage(operation),'base_input_path',
                            [{'stage':source['id'],'selector':ligand.name}],{source['id']:[str(ligand)]})
                    with self.assertRaisesRegex(ValueError,'receptor preparado'):
                        service._materialize_inputs(self.project_id,new_stage(operation),'base_input_path',[{'asset':asset}],{})
        finally:service.close()
