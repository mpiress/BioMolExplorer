"""Saving PDB ligands must change the executable plan, not only its metadata."""
import copy
import asyncio
import json
from pathlib import Path
import shutil
import unittest
from types import SimpleNamespace
from unittest.mock import patch, Mock, AsyncMock

from biomolexplorer.bindings import sources
from biomolexplorer.catalog import new_stage
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.redocking_results import RedockingResults
from biomolexplorer.operations import execute_operation
from biomolexplorer.pipeline import PipelineService
import test_pipeline_execution as execution
import test_redocking_incremental as incremental
from test_redocking_partial_failures import ATOM
import test_retrieval_usability as retrieval


class CuratedPairPlanTests(unittest.TestCase):
    setUp=retrieval.PDBCurationTests.setUp
    tearDown=retrieval.PDBCurationTests.tearDown
    context=retrieval.PDBCurationTests.context
    save=retrieval.PDBCurationTests.save

    def consumer(self,selector='auto'):
        stage=new_stage('redocking')
        stage['bindings']['base_input_path']={'stage':self.sid,'selector':selector}
        stage['parameters']['pdb_codes']=[['2XYZ','ATP',3,'B',1.5]]
        revision=self.store.project(self.token,self.pid)['revision']
        self.store.save_pipeline(self.token,self.pid,[self.config,stage],revision)
        return stage

    def test_new_curated_pair_is_added_to_connected_execution_plan(self):
        consumer=self.consumer('2XYZ.pdb')
        self.save([])
        self.save([{'LIGAND':'LIG','RESNUM':101,'CHAIN':'A'}])
        stage=self.store.project(self.token,self.pid)['pipeline'][1]
        self.assertEqual([r[:4] for r in stage['parameters']['pdb_codes']],
                         [['2XYZ','ATP',3,'B'],['1ABC','LIG',101,'A']])
        refs=sources(stage['bindings']['base_input_path'])
        self.assertEqual({Path(r['selector']).name for r in refs},{'2XYZ.pdb','1ABC.pdb'})
        self.assertEqual(stage['id'],consumer['id'])

    def test_already_saved_ligand_repairs_missing_plan_and_repeated_save_is_idempotent(self):
        self.consumer()
        records=[{'LIGAND':'LIG','RESNUM':101,'CHAIN':'A'}]
        first=self.save(records);second=self.save(records)
        stage=self.store.project(self.token,self.pid)['pipeline'][1]
        self.assertEqual(len(stage['parameters']['pdb_codes']),2)
        self.assertEqual(first['added_redocking_pairs'],1)
        self.assertEqual(second['added_redocking_pairs'],0)
        self.assertEqual(second['updated_redocking_stages'],[])

    def test_placeholder_for_structure_without_ligand_does_not_block_first_pair(self):
        self.consumer()
        self.csv.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,,,,\n')
        self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A')])
        stage=self.store.project(self.token,self.pid)['pipeline'][1]
        self.assertEqual(stage['parameters']['pdb_codes'][-1][:4],['1ABC','LIG',101,'A'])

    def test_new_ligand_keeps_shared_receptor_options_and_other_pairs(self):
        self.path.write_bytes(self.path.read_bytes().replace(b'ATP B 102',b'ATP A 102'))
        consumer=self.consumer()
        consumer['parameters']['pdb_codes'].append(['1ABC','LIG',101,'A'])
        options={'receptor':{'minimize':False},'ligand':{'minimize':False},'cofactors':[]}
        consumer['parameters']['preparation_pairs']={'1ABC|LIG|101|A':options}
        project=self.store.project(self.token,self.pid)
        self.store.save_pipeline(self.token,self.pid,[self.config,consumer],project['revision'])
        self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A'),dict(LIGAND='ATP',RESNUM=102,CHAIN='A')])
        stage=self.store.project(self.token,self.pid)['pipeline'][1]
        settings=stage['parameters']['preparation_pairs']
        self.assertEqual(settings['1ABC|LIG|101|A'],options)
        self.assertEqual(settings['1ABC|ATP|102|A'],options)
        self.assertEqual(len(stage['parameters']['pdb_codes']),3)

    def test_new_receptor_inherits_general_preparation_without_other_receptor_cofactors(self):
        consumer=self.consumer()
        general={'receptor':{'minimize':False},'ligand':{'charge_type':'am1'},
                 'preparation_backend':'legacy_chimera','ambertools_env':'BioMolExplorer-AmberTools'}
        consumer['parameters']['preparation_pairs']={'2XYZ|ATP|3|B':dict(general,cofactors=['FAD'])}
        project=self.store.project(self.token,self.pid)
        self.store.save_pipeline(self.token,self.pid,[self.config,consumer],project['revision'])
        self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A')])
        settings=self.store.project(self.token,self.pid)['pipeline'][1]['parameters']['preparation_pairs']
        self.assertEqual(settings['1ABC|LIG|101|A'],general)
        self.assertEqual(settings['2XYZ|ATP|3|B']['cofactors'],['FAD'])

    def test_removal_prunes_only_the_edited_metadata_pair(self):
        self.consumer()
        self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A')]);self.save([])
        stage=self.store.project(self.token,self.pid)['pipeline'][1]
        self.assertEqual([row[:4] for row in stage['parameters']['pdb_codes']],[['2XYZ','ATP',3,'B']])

    def test_failed_commit_restores_metadata_and_execution_plan(self):
        self.consumer();project=self.store.project(self.token,self.pid);content=self.csv.read_bytes()
        with patch('biomolexplorer.stage_cache.artifact_manifest',side_effect=RuntimeError('test rollback')):
            with self.assertRaisesRegex(RuntimeError,'test rollback'):
                self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A')])
        self.assertEqual(self.csv.read_bytes(),content)
        current=self.store.project(self.token,self.pid)
        self.assertEqual(current['pipeline'],project['pipeline'])
        self.assertEqual(current['revision'],project['revision'])

    def test_unconnected_redocking_is_untouched(self):
        consumer=self.consumer();consumer['bindings']['base_input_path']={'stage':'other_source'}
        # Use the pure synchronizer: this deliberately unrelated fixture source
        # need not exist in a valid persisted pipeline.
        from biomolexplorer.redocking_selection import sync_curated_pairs
        before=copy.deepcopy(consumer)
        changes=sync_curated_pairs([consumer],self.sid,self.path,self.stage['artifacts'],[],
            [dict(PDB_CODE='1ABC',LIGAND='LIG',RESNUM=101,CHAIN='A')])
        self.assertEqual(changes['updated_redocking_stages'],[])
        self.assertEqual(consumer,before)

    def test_frontend_refresh_keeps_canvas_edits_and_adopts_new_pairs(self):
        self.consumer();current=self.store.project(self.token,self.pid)
        base=copy.deepcopy(current['pipeline'])
        current['pipeline'][0]['position']={'x':600,'y':700}
        self.save([dict(LIGAND='LIG',RESNUM=101,CHAIN='A')])
        async def call(fn,*args):return fn(*args)
        ui=SimpleNamespace(current=current,base_pipeline=base,dirty=True,call=call,
            schedule_save=Mock(),tab='Pipeline',flow_editor=SimpleNamespace(refresh=Mock()))
        results=SimpleNamespace(ui=ui,project_id=self.pid,token=self.token,service=self.service,valid=lambda:True)
        from biomolexplorer.ui.pdb_results import PDBActions
        asyncio.run(PDBActions(results,True).refresh_pipeline())
        self.assertIs(ui.current,current)
        self.assertEqual(current['pipeline'][0]['position'],{'x':600,'y':700})
        self.assertIn(['1ABC','LIG',101,'A'],[r[:4] for r in current['pipeline'][1]['parameters']['pdb_codes']])
        self.assertTrue(ui.dirty);ui.schedule_save.assert_called_once()
        ui.flow_editor.refresh.assert_called_once()

    def test_saving_dialog_refreshes_the_plan_before_run_pipeline(self):
        self.consumer();project=self.store.project(self.token,self.pid)
        dialogs=[];messages=[]
        async def call(fn,*args):return fn(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(current=project,base_pipeline=copy.deepcopy(project['pipeline']),dirty=False,
            call=call,guard=guard,save_draft=AsyncMock(),notify=messages.append,editing_stage=False,
            page=SimpleNamespace(update=Mock(),show_dialog=dialogs.append))
        result=SimpleNamespace(ui=ui,project_id=self.pid,token=self.token,run_id=self.rid,
            stage=self.stage,service=self.service,valid=lambda:True)
        from biomolexplorer.ui.pdb_results import PDBActions
        asyncio.run(PDBActions(result,True).open_ligands({'path':str(self.path)}))
        asyncio.run(dialogs[-1].actions[-1].on_click(None))
        self.assertEqual([row[:4] for row in ui.current['pipeline'][1]['parameters']['pdb_codes']],
            [['2XYZ','ATP',3,'B'],['1ABC','LIG',101,'A']])
        self.assertFalse(ui.dirty);self.assertFalse(ui.editing_stage)
        self.assertEqual(ui.current['revision'],self.store.project(self.token,self.pid)['revision'])
        ui.save_draft.assert_awaited_once_with(silent=True)
        self.assertIn('redocking',messages[0])


class CuratedRedockingExecutionTests(unittest.TestCase):
    def test_1p0i_bua_saved_after_first_run_is_docked_and_accumulated(self):
        fixture=execution.PipelineExecutionTests();fixture.setUp();self.addCleanup(fixture.tearDown)
        science=incremental.IncrementalRedockingTests();science.setUp();self.addCleanup(science.doCleanups)
        bua='HETATM    3  C   BUA A   3       2.000   0.000   0.000  1.00 20.00           C\n'
        (science.source/'1P0I.pdb').write_text(ATOM+bua)
        producer=new_stage('retrieve_structures');producer['name']='BChE PDB'
        producer['parameters'].update(target='Target',pdb_ids=['1ABC','1P0I'])
        consumer=new_stage('redocking');consumer['input_processing']='individual'
        consumer['parameters'].update(target='Target',pdb_codes=[['1ABC','LIG',1,'A']],
            preparation_pairs={'1ABC|LIG|1|A':{'preparation_backend':'legacy_chimera'}})
        consumer['bindings']['base_input_path']={'stage':producer['id'],'selector':'1ABC.pdb'}
        fixture.save([producer,consumer]);calls=[]
        class Manager:
            def __init__(self,config):self.config=config
            def submit(self,operation,parameters):
                calls.append((operation,copy.deepcopy(parameters)))
                output=self.config.workspace/'artifacts';output.mkdir(parents=True,exist_ok=True)
                if operation=='retrieve_structures':
                    target=output/'Target';target.mkdir()
                    for code in ('1ABC','1P0I'):shutil.copy2(science.source/(code+'.pdb'),target/(code+'.pdb'))
                    (target/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,LIG,1,A,2.0\n')
                    result=dict(artifacts=[str(p) for p in target.iterdir()],details={})
                else:result=execute_operation(operation,parameters,output).to_dict()
                self.job=dict(id='worker',status='succeeded',log_path=str(output/'execution.log'),progress=None,error=None,result=result)
                return self.job
            def get(self,job):return self.job
            def cancel(self,job):self.job['status']='cancelled'
            def close(self):pass
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            service=PipelineService(fixture.store);self.addCleanup(service.close)
            first=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
            self.assertEqual(first['status'],'succeeded',first.get('error'))
            raw=next(p for p in first['stages'][0]['artifacts'] if Path(p).name=='1P0I.pdb')
            results=ResultFiles(fixture.store)
            args=(fixture.token,fixture.project_id,first['id'],producer['id'],raw)
            context=results.pdb_ligands(*args);self.assertEqual(context['ligands'],[])
            results.set_pdb_ligands(*args,[dict(LIGAND='BUA',RESNUM=3,CHAIN='A')],context['revision'])
            second=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
        self.assertEqual(second['status'],'succeeded',second.get('error'))
        self.assertEqual(science.docked,[('1ABC','LIG'),('1P0I','BUA')])
        rows=RedockingResults(fixture.store).simulations(fixture.token,fixture.project_id,second['id'],consumer['id'])
        self.assertEqual({(r['pdb'],r['ligand']) for r in rows},{('1ABC','LIG'),('1P0I','BUA')})
        self.assertEqual(len(rows),2)
        self.assertEqual(sum(op=='retrieve_structures' for op,_ in calls),1)


if __name__=='__main__':unittest.main()
