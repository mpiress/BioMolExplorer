"""Large analysis inputs must not be prepared by the confirmation request."""
import asyncio
import json
import time
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.ui.guided import GuidedForm
import test_pipeline_execution as execution


class AnalysisHandoffTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown

    def pending_run(self,operation,content):
        producer=new_stage({'fingerprints':'retrieve_compounds','similarity':'fingerprints','graphs':'similarity'}[operation])
        stage=new_stage(operation)
        field='similarity_path' if operation=='graphs' else 'base_input_path'
        paths=[]
        for name in ('alpha','beta'):
            path=self.store.project_dir(self.project_id)/(name+'.csv')
            path.write_text(content);paths.append(str(path))
        stage['bindings'][field]={'sources':[{'stage':producer['id'],'selector':name+'.csv'} for name in ('alpha','beta')]}
        stages=[dict(id=producer['id'],name=producer['name'],operation=producer['operation'],
                     configuration=producer,status='succeeded',artifacts=paths),
                dict(id=stage['id'],name=stage['name'],operation=operation,
                     configuration=stage,status='awaiting_input',artifacts=[])]
        rid=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                       (rid,self.project_id,self.store.user(self.token)['id'],'awaiting_input',
                        json.dumps(stages),time.time(),time.time(),None))
        return rid,stage

    def test_confirmation_queues_multiple_files_without_cleaning_or_materializing(self):
        cases={'fingerprints':'molecule_chembl_id,canonical_smiles\nA,CCO\n',
               'similarity':'molecule_chembl_id,fingerprint\nA,"[0, 1]"\n',
               'graphs':'source,target,value\nA,B,0.9\n'}
        service=PipelineService(self.store)
        try:
            for operation,content in cases.items():
                for mode in ('individual','merge'):
                    with self.subTest(operation=operation,mode=mode):
                        rid,stage=self.pending_run(operation,content)
                        stage['input_processing']=mode
                        with patch.object(service,'_resolve',side_effect=AssertionError('preparing inputs on confirmation')), \
                             patch.object(service._executor,'submit') as submit:
                            resumed=service.resume(self.token,rid,stage)
                        self.assertEqual(resumed['status'],'queued')
                        self.assertTrue(resumed['stages'][1]['curation_confirmed'])
                        submit.assert_called_once()
                        self.assertFalse((self.store.project_dir(self.project_id)/'.input-cache').exists())
        finally:service.close()

    def test_invalid_file_remains_pending(self):
        rid,stage=self.pending_run('similarity','molecule_chembl_id,canonical_smiles\nA,CCO\n')
        service=PipelineService(self.store)
        try:
            with patch.object(service._executor,'submit') as submit:
                with self.assertRaisesRegex(ValueError,'fingerprint'):service.resume(self.token,rid,stage)
                submit.assert_not_called()
            self.assertEqual(self.store.get_run(self.token,rid)['status'],'awaiting_input')
        finally:service.close()

    def test_graph_complement_is_optional_when_unselected(self):
        rid,stage=self.pending_run('graphs','source,target,value\nA,B,0.9\n')
        stage['bindings']['base_input_path']={'asset':uuid4().hex}
        pending=dict(name=stage['name'],configuration=stage)
        form=FileSelection(self.store.get_run(self.token,rid),pending)
        form.rows['base_input_path'][0][0].value=False
        self.assertNotIn('base_input_path',form.read()['bindings'])

    def test_picker_changes_refresh_inferred_fingerprint_and_preserve_multiple_files(self):
        producer=new_stage('fingerprints');producer['parameters'].update(morgan=True,maccs=True)
        stage=new_stage('similarity');stage['bindings']['base_input_path']={'stage':producer['id']}
        paths=[]
        for name in ('morgan_alpha','maccs_alpha','maccs_beta'):
            path=self.store.project_dir(self.project_id)/(name+'.csv')
            path.write_text('molecule_chembl_id,fingerprint\nA,"[0, 1]"\n');paths.append(str(path))
        item=dict(id=producer['id'],operation='fingerprints',name=producer['name'],
                  status='succeeded',artifacts=paths,configuration=producer)
        ui=SimpleNamespace(current={'pipeline':[producer,stage]},artifact_items={producer['id']:item},
                           page=SimpleNamespace(update=lambda:None),service=SimpleNamespace(dock6_path=None))
        form=GuidedForm(ui,stage,[],True)
        picker=form.input_editors['base_input_path'].selection.pickers['base_input_path']
        picker.set_selected([i for i,(_,ref) in enumerate(picker.rows) if ref['selector'].startswith('maccs_')],True)
        self.assertEqual(form.field_controls['fingerprint'].value,'maccs')
        result=form.read()
        self.assertEqual(len(result['bindings']['base_input_path']['sources']),2)
        self.assertEqual(result['parameters']['fingerprint'],'maccs')

    def test_repeated_click_is_ignored_and_failure_reenables_button(self):
        async def scenario():
            ui=WorkspaceUI(SimpleNamespace(update=lambda:None),None,None)
            ui.notify=lambda message:None
            entered=asyncio.Event();release=asyncio.Event();calls=[]
            async def handler(e):
                calls.append(True);entered.set();await release.wait()
                raise ValueError('Seleção inválida')
            control=SimpleNamespace(disabled=False)
            event=SimpleNamespace(control=control)
            action=ui.event_action(handler)
            first=asyncio.create_task(action(event));await entered.wait()
            self.assertTrue(control.disabled)
            await action(event)
            self.assertEqual(len(calls),1)
            release.set();await first
            self.assertFalse(control.disabled)
            await action(event)
            self.assertEqual(len(calls),2)
        asyncio.run(scenario())
