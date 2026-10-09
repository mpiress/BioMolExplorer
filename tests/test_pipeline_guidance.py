"""User-facing connection rejection, block guidance and compound selection."""
import asyncio
import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

import test_pipeline_execution as execution
import test_docking_handoffs as handoffs
from biomolexplorer.catalog import TITLES,new_stage,operation_fields
from biomolexplorer.pipeline import PipelineService,validate_pipeline
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.block_help import block_information
from biomolexplorer.ui.input_editor import InputEditor
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.docking_data import read_compounds,write_csv


class ConnectionGuidanceTests(unittest.TestCase):
    def test_rejected_click_and_drag_are_safe_and_next_valid_link_works(self):
        import flet as ft
        for operation in ('retrieve_compounds','retrieve_pubchem'):
            source,target,fingerprints=new_stage(operation),new_stage('similarity'),new_stage('fingerprints')
            stages=[source,target,fingerprints];notifications=[]
            ui=SimpleNamespace(current={'pipeline':stages},selected=None,dirty=False,
                page=SimpleNamespace(update=Mock()),notify=notifications.append)
            editor=FlowCanvas(ui,True);editor.pending=source['id']
            before=copy.deepcopy(stages)
            # Invoke the actual DragTarget callback rendered for the similarity input.
            node=editor.node(target,1)
            def walk(control):
                yield control
                for child in getattr(control,'controls',[]) or []:yield from walk(child)
                if getattr(control,'content',None) is not None:yield from walk(control.content)
            port=next(c for c in walk(node) if isinstance(c,ft.DragTarget))
            asyncio.run(port.on_accept(SimpleNamespace(src=SimpleNamespace(data=source['id']))))
            self.assertEqual(stages,before)
            self.assertEqual(notifications,['Dados incompatíveis. Escolha uma saída do tipo esperado por esta entrada.'])
            self.assertEqual(editor.history,[]);self.assertIsNone(editor.pending)
            asyncio.run(editor.link(fingerprints['id'],target['id'],'base_input_path'))
            self.assertEqual(target['bindings']['base_input_path']['stage'],fingerprints['id'])
            self.assertEqual(len(editor.history),1)

    def test_every_block_has_information_and_popup_is_available_to_viewers(self):
        from biomolexplorer.ui.localization import Translator
        from biomolexplorer.ui.block_help import TYPE_DESCRIPTIONS,NOTES
        translate=Translator('en')
        for text in list(TYPE_DESCRIPTIONS.values())+list(NOTES.values()):
            self.assertNotEqual(translate(text),text)
        page=SimpleNamespace(update=Mock(),show_dialog=Mock())
        ui=SimpleNamespace(current={'pipeline':[]},selected=None,page=page)
        editor=FlowCanvas(ui,False)
        for operation in TITLES:
            stage=new_stage(operation)
            description,inputs,outputs,_=block_information(stage)
            self.assertTrue(description);self.assertTrue(inputs);self.assertTrue(outputs)
            editor.show_information(stage)
            dialog=page.show_dialog.call_args.args[0]
            self.assertEqual(dialog.title.value,TITLES[operation][0])
        self.assertIn('Gerar fingerprints',block_information(new_stage('similarity'))[3])


class CompoundSelectionUITests(unittest.TestCase):
    def editor(self,operation='docking_vina'):
        origin,stage=new_stage('retrieve_compounds'),new_stage(operation)
        stage['bindings']['base_selected_mols']={'stage':origin['id'],'selector':'first.csv','compound_id':'B'}
        ui=SimpleNamespace(current={'pipeline':[origin,stage]},page=SimpleNamespace(update=Mock()),
            artifact_choices={origin['id']:{'first.csv','second.csv'}},
            docking_compounds={'stage:'+origin['id']:{'/results/first.csv':{'A':'CCO','B':'CCN'},
                '/results/second.csv':{'C':'CCC'}}})
        field=next(f for f in operation_fields(operation) if f['name']=='base_selected_mols')
        return InputEditor(ui,stage,field,[],True)

    def test_compound_choice_tracks_file_and_empty_means_all_in_both_engines(self):
        for operation in ('docking_vina','docking_dock6'):
            editor=self.editor(operation);source,selector,_=editor.rows[0]
            choice=editor.compound_fields[id(source)]
            self.assertTrue(choice.visible)
            self.assertEqual({o.key for o in choice.options},{'','A','B'})
            self.assertEqual(editor.read()['compound_id'],'B')
            choice.value=''
            self.assertNotIn('compound_id',editor.read())
            choice.value='B';selector.value='second.csv';selector.on_select(None)
            self.assertEqual(choice.value,'')
            self.assertEqual({o.key for o in choice.options},{'','C'})
            choice.value='C';self.assertEqual(editor.read()['compound_id'],'C')

    def test_file_popup_preserves_compound_through_state_and_settings(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'compounds.csv'
            write_csv(path,[dict(molecule_chembl_id='A',canonical_smiles='CCO'),
                            dict(molecule_chembl_id='B',canonical_smiles='CCN')])
            source,stage=new_stage('retrieve_compounds'),new_stage('docking_vina')
            stage['bindings']['base_selected_mols']={'stage':source['id'],'selector':'compounds.csv','compound_id':'B'}
            run={'stages':[dict(source,status='succeeded',artifacts=[str(path)])]}
            pending=dict(stage,configuration=stage)
            form=FileSelection(run,pending)
            check,_=form.rows['base_selected_mols'][0]
            self.assertEqual(form.compound_fields[id(check)].value,'B')
            state=form.state();form=FileSelection(run,pending,state=state)
            self.assertEqual(form.read()['bindings']['base_selected_mols']['compound_id'],'B')
            check,_=form.rows['base_selected_mols'][0]
            form.compound_fields[id(check)].value=''
            self.assertNotIn('compound_id',form.read()['bindings']['base_selected_mols'])


class CompoundSelectionPipelineTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager
    prepared=handoffs.DockingHandoffTests.prepared

    def test_select_one_compound_or_all_and_isolate_materialized_inputs(self):
        asset=self.upload('compounds.csv','name,smiles\nA,CCO\nB,CCN\n')
        service=PipelineService(self.store)
        try:
            for operation in ('docking_vina','docking_dock6'):
                stage=new_stage(operation)
                chosen=[]
                for code in ('A','B',None):
                    ref={'asset':asset}
                    if code:ref['compound_id']=code
                    stage['bindings']['base_selected_mols']=ref
                    folder,stem,_=service._materialize_inputs(self.project_id,stage,'base_selected_mols',[ref],{})
                    chosen.append(folder)
                    self.assertEqual({r['molecule_chembl_id'] for r in read_compounds(folder/(stem+'.csv'))},
                        {code} if code else {'A','B'})
                self.assertEqual(len(set(chosen)),3)
                with self.assertRaisesRegex(ValueError,'não está no arquivo'):
                    service._materialize_inputs(self.project_id,stage,'base_selected_mols',[{'asset':asset,'compound_id':'missing'}],{})
        finally:service.close()

    def test_selection_is_per_file_and_preserves_prepared_compound_files(self):
        root=self.store.project_dir(self.project_id);folder=root/'candidates';folder.mkdir()
        pose=folder/'B.lig.mol2';pose.write_text(handoffs.MOL2)
        table=folder/'compounds.csv'
        write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='CCO',prepared_mol2='',prepared_origin=''),
            dict(molecule_chembl_id='B',canonical_smiles='CCN',prepared_mol2=str(pose),prepared_origin='library')])
        asset=self.upload('other.csv','name,smiles\nA,CCO\nC,CCC\n')
        origin,stage=new_stage('prepare_structures'),new_stage('docking_dock6')
        refs=[{'stage':origin['id'],'selector':'compounds.csv','compound_id':'B'},
              {'asset':asset,'compound_id':'C'}]
        stage['bindings']['base_selected_mols']={'sources':refs}
        service=PipelineService(self.store)
        try:
            output,stem,_=service._materialize_inputs(self.project_id,stage,'base_selected_mols',refs,
                {origin['id']:[str(table),str(pose)]})
            rows=read_compounds(output/(stem+'.csv'))
            self.assertEqual({r['molecule_chembl_id'] for r in rows},{'B','C'})
            chosen=next(r for r in rows if r['molecule_chembl_id']=='B')
            self.assertEqual(Path(chosen['prepared_mol2']).read_text(),pose.read_text())
            self.assertEqual(chosen['prepared_origin'],'library')
        finally:service.close()

    def test_preconfigured_docking_executes_without_popup_and_auto_still_asks(self):
        self.prepared(None)
        root=self.store.project_dir(self.project_id)
        dock6=root/'dock6-bin';dock6.mkdir()
        asset=self.upload('compounds.csv','name,smiles\nA,CCO\nB,CCN\n')
        for operation in ('docking_vina','docking_dock6'):
            for explicit in (True,False):
                with self.subTest(operation=operation,explicit=explicit):
                    source=new_stage('import_results');source['parameters']['asset_ids']=[asset]
                    stage=new_stage(operation)
                    stage['parameters']['base_input_path']=str(root/'1ABC_A')
                    if operation=='docking_dock6':stage['parameters'].update(dock6_app_path=str(dock6),charge_type='gas')
                    ref={'stage':source['id'],'selector':'compounds.csv' if explicit else 'auto'}
                    if explicit:ref['compound_id']='B'
                    stage['bindings']['base_selected_mols']=ref
                    self.save([source,stage]);manager,calls=self.immediate_manager()
                    with patch('biomolexplorer.pipeline.JobManager',manager):
                        service=PipelineService(self.store,dock6_path=dock6)
                        try:
                            run=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                            if run['status']=='awaiting_input':service.cancel(self.token,run['id'])
                        finally:service.close()
                    self.assertEqual(run['status'],'succeeded' if explicit else 'awaiting_input',run['error'])
                    self.assertEqual(run['stages'][1]['requires_curation'],not explicit)
                    if explicit:
                        self.assertEqual([op for op,_ in calls],[operation])
                        params=calls[0][1]
                        table=Path(params['base_selected_mols'])/(params['mol_filename']+'.csv')
                        self.assertEqual([r['molecule_chembl_id'] for r in read_compounds(table)],['B'])
                    else:self.assertEqual(calls,[])

    def test_binding_validation_rejects_compound_selection_on_other_inputs(self):
        stage=new_stage('admet')
        stage['bindings']['base_input_path']={'asset':'f'*32,'compound_id':'A'}
        with self.assertRaisesRegex(ValueError,'identificador de composto'):validate_pipeline([stage])
