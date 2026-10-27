"""Preparation failures stay visible; deleting them removes exact input origins."""
import asyncio
import csv
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import pandas as pd
import test_results_display as display
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_results import RedockingResults
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import AccessDenied
from biomolexplorer.ui.stage_results import StageResults
from wrappers.redocking import perform_redocking

ATOM='ATOM      1  CA  ALA A   2       0.000   0.000   0.000  1.00 20.00           C\n'
LIG='HETATM    2  C   LIG A   1       1.000   0.000   0.000  1.00 20.00           C\n'
BAD='HETATM    3  C   BAD A   3       2.000   0.000   0.000  1.00 20.00           C\n'


class PartialPreparationTests(unittest.TestCase):
    def run_redocking(self,records,tool_error=False):
        from caad.docking import Docking
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        root=Path(self.temp.name);source=root/'input'/'Target';source.mkdir(parents=True)
        for code in {r[0] for r in records}:(source/(code+'.pdb')).write_text(ATOM+LIG+BAD)
        self.scripts=[]
        def chimera(dock,filename):
            self.scripts.append(filename)
            if tool_error:raise FileNotFoundError(2,'tool missing','chimera')
            if filename.startswith('prepare_complex_2ABC'):raise RuntimeError('invalid receptor')
            if '_BAD_' in filename:
                (Path(dock.outputpath)/'1ABC_BAD_3A.lig.mol2').write_text('partial')
                raise RuntimeError('invalid ligand')
            directory=Path(dock.outputpath)
            if filename.startswith('prepare_complex_'):
                name=filename.removeprefix('prepare_complex_').removesuffix('.com')+'.complex.pdb'
                (directory/name).write_text(ATOM+LIG+BAD)
            elif filename.startswith('prepare_receptor_'):
                name=filename.removeprefix('prepare_receptor_').removesuffix('.com')+'.dockprep.mol2'
                (directory/name).write_text('@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n')
            else:
                name=filename.removeprefix('prepare_ligand_').removesuffix('.com')
                (directory/(name+'.lig.mol2')).write_text('@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n')
                (directory/(name+'.lig.pdb')).write_text(LIG)
        def obabel(dock,inputfile,outputfile,*args,**kwargs):(Path(dock.outputpath)/outputfile).write_text(ATOM)
        def vina_run(vina):
            frame=vina.call_args.kwargs['pdb_codes'].copy();frame['RMSD']=.2
            frame.to_csv(source/'pdb_codes.csv',index=False)
        with patch.object(Docking,'prepare_on_chimera',chimera),patch.object(Docking,'prepare_on_obabel',obabel), \
             patch.object(Docking,'calculate_ligand_centerofmass',return_value=[0,0,0]), \
             patch('biomolexplorer.redocking_config.validate_redocking_tools'),patch('wrappers.redocking.DockVina') as vina:
            vina.return_value.redocking.side_effect=lambda **kwargs:vina_run(vina)
            result=perform_redocking(str(root/'input'),'Target',str(root/'output'),records,
                                    preparation_pairs={'|'.join(map(str,r[:4])):{'preparation_backend':'legacy_chimera'} for r in records})
            self.vina_called=vina.called
        self.source=source;self.output=root/'output'
        return result

    def test_bad_ligand_and_bad_receptor_do_not_stop_remaining_pairs(self):
        records=[['1ABC','LIG',1,'A',1.7],['1ABC','BAD',3,'A',1.7],['2ABC','LIG',1,'A',2.0],['3ABC','LIG',1,'A',2.1]]
        result=self.run_redocking(records)
        rows=list(csv.DictReader(io.StringIO((self.source/'pdb_codes.csv').read_text())))
        self.assertEqual([r['STATUS'] for r in rows],['succeeded','preparation_failed','preparation_failed','succeeded'])
        self.assertEqual([r['RMSD'] for r in rows],['0.2','','','0.2'])
        self.assertIn('invalid ligand',rows[1]['ERROR'])
        self.assertEqual(result,{'redocking_failures':2,'redocking_succeeded':2})
        self.assertEqual(self.scripts.count('prepare_receptor_1ABC_A.com'),1)
        self.assertFalse((self.source/'Prepared/1ABC_BAD_3A.lig.mol2').exists())
        failures=self.output/'Target/redocking_failures.csv'
        self.assertEqual(len(list(csv.DictReader(io.StringIO(failures.read_text())))),2)

    def test_all_preparations_can_fail_without_running_vina_or_raising(self):
        result=self.run_redocking([['1ABC','BAD',3,'A',1.7]])
        self.assertFalse(self.vina_called)
        self.assertEqual(result['redocking_succeeded'],0)
        rows=pd.read_csv(self.source/'pdb_codes.csv')
        self.assertEqual(rows.iloc[0]['STATUS'],'preparation_failed')
        self.assertTrue(pd.isna(rows.iloc[0]['RMSD']))

    def test_missing_installation_is_not_silently_classified_as_a_bad_complex(self):
        with self.assertRaises(FileNotFoundError):self.run_redocking([['1ABC','LIG',1,'A',1.7]],tool_error=True)


class FailedResultDeletionTests(unittest.TestCase):
    setUpBase=display.ResultsDisplayTests.setUp
    tearDown=display.ResultsDisplayTests.tearDown
    ui=display.ResultsDisplayTests.ui

    def setUp(self):
        self.setUpBase()
        self.producer=new_stage('retrieve_structures');self.consumer=new_stage('redocking')
        self.consumer['bindings']['base_input_path']={'stage':self.producer['id']}
        self.consumer['parameters']['pdb_codes']=[['1ABC','BAD',3,'A'],['2ABC','LIG',1,'A']]
        self.producer['parameters']['pdb_ids']=['1ABC','2ABC']
        self.store.save_pipeline(self.token,self.pid,[self.producer,self.consumer],0)
        root=self.store.project_dir(self.pid)/'runs'/self.rid/'partial'
        source=root/'retrieval'/'PDB';source.mkdir(parents=True)
        self.original=source/'1ABC.pdb';self.original.write_text(ATOM+LIG+BAD)
        self.good=source/'2ABC.pdb';self.good.write_text(ATOM+LIG)
        self.original_meta=source/'pdb_codes.csv'
        self.original_meta.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,BAD,3,A\n2ABC,LIG,1,A\n')
        self.snapshot=root/'redocking'/'artifacts'/'structures'/'Target';self.snapshot.mkdir(parents=True)
        (self.snapshot/'1ABC.pdb').write_text(self.original.read_text())
        (self.snapshot/'2ABC.pdb').write_text(self.good.read_text())
        self.metadata=self.snapshot/'pdb_codes.csv'
        self.metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD,STATUS,FAILED_STEP,ERROR\n1ABC,BAD,3,A,,preparation_failed,Chimera,invalid ligand\n2ABC,LIG,1,A,0.2,succeeded,,\n')
        self.unrelated=self.store.project_dir(self.pid)/'unrelated'/'1ABC.pdb';self.unrelated.parent.mkdir()
        self.unrelated.write_text(ATOM+LIG+BAD)
        self.cache=self.store.project_dir(self.pid)/'.input-cache/key/Target/1ABC.pdb';self.cache.parent.mkdir(parents=True)
        self.cache.write_text(self.original.read_text())
        stages=[]
        for config,files in ((self.producer,[self.original,self.good,self.original_meta]),(self.consumer,list(self.snapshot.glob('*')))):
            stages.append(dict(id=config['id'],name=config['name'],operation=config['operation'],configuration=config,
                status='succeeded',artifacts=[str(p) for p in files],artifact_manifest=artifact_manifest(files)))
        stages[1]['input_files']={'base_input_path':[str(self.original),str(self.good),str(self.original_meta)]}
        stages[1]['redocking_failures']=1
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        self.stage=stages[1];self.service=RedockingResults(self.store)
        self.args=(self.token,self.pid,self.rid,self.consumer['id'])

    def failure(self):return next(r for r in self.service.simulations(*self.args) if r['status']=='preparation_failed')

    def test_failure_is_visible_and_deletion_cleans_source_configuration_and_cache(self):
        rows=self.service.simulations(*self.args)
        self.assertEqual([r['rmsd'] for r in rows],[None,.2])
        self.service.remove_failure(*self.args,self.failure()['id'])
        self.assertFalse(self.original.exists());self.assertFalse(self.cache.exists())
        self.assertFalse((self.snapshot/'1ABC.pdb').exists())
        self.assertTrue(self.good.exists());self.assertTrue(self.unrelated.exists())
        self.assertEqual([r['pdb'] for r in self.service.simulations(*self.args)],['2ABC'])
        source_rows=list(csv.DictReader(io.StringIO(self.original_meta.read_text())))
        self.assertEqual([r['PDB_CODE'] for r in source_rows],['2ABC'])
        pipeline=self.store.project(self.token,self.pid)['pipeline']
        self.assertEqual(pipeline[0]['parameters']['excluded_pdb_ids'],['1ABC'])
        self.assertEqual(pipeline[1]['parameters']['pdb_codes'],[['2ABC','LIG',1,'A']])
        stages=self.store.get_run(self.token,self.rid)['stages']
        self.assertTrue(all(s['outputs_removed'] for s in stages))
        self.assertNotIn(str(self.original),stages[0]['artifacts'])

    def test_cannot_delete_success_forged_key_as_viewer_or_while_running(self):
        success=next(r for r in self.service.simulations(*self.args) if r['status']=='succeeded')
        with self.assertRaisesRegex(ValueError,'Somente linhas'):self.service.remove_failure(*self.args,success['id'])
        with self.assertRaises(AccessDenied):self.service.remove_failure(*self.args,'forged')
        self.store.invite(self.token,self.pid,'guest@example.org','viewer');self.store.accept_invitation(self.guest,self.pid)
        with self.assertRaises(AccessDenied):self.service.remove_failure(self.guest,*self.args[1:],self.failure()['id'])
        key=self.failure()['id']
        with self.store.connect() as db:db.execute("UPDATE runs SET status='running' WHERE id=?",(self.rid,))
        with self.assertRaisesRegex(ValueError,'Aguarde'):self.service.remove_failure(*self.args,key)
        self.assertTrue(self.original.exists())

    def test_transaction_failure_restores_files_and_metadata(self):
        before=self.metadata.read_bytes();original=self.original.read_bytes();project=self.store.project(self.token,self.pid)
        with patch('biomolexplorer.redocking_curation.artifact_manifest',side_effect=RuntimeError('test rollback')):
            with self.assertRaisesRegex(RuntimeError,'test rollback'):self.service.remove_failure(*self.args,self.failure()['id'])
        self.assertEqual(self.metadata.read_bytes(),before)
        self.assertEqual(self.original.read_bytes(),original)
        self.assertTrue(self.cache.exists())
        self.assertEqual(self.store.project(self.token,self.pid)['revision'],project['revision'])

    def test_changed_source_is_not_deleted_as_if_it_were_the_failed_input(self):
        self.original.write_text(ATOM+LIG)
        with self.assertRaisesRegex(ValueError,'origem foi alterado'):
            self.service.remove_failure(*self.args,self.failure()['id'])
        self.assertTrue(self.original.exists());self.assertTrue(self.cache.exists())

    def test_successful_sibling_from_the_same_pdb_keeps_its_prepared_results(self):
        prepared=self.snapshot/'Prepared/1ABC_A.dockprep.pdbqt';prepared.parent.mkdir()
        prepared.write_text(ATOM)
        self.metadata.write_text(self.metadata.read_text()+'1ABC,LIG,1,A,0.3,succeeded,,\n')
        run=self.store.get_run(self.token,self.rid);run['stages'][1]['artifacts'].append(str(prepared))
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),self.rid))
        self.service.remove_failure(*self.args,self.failure()['id'])
        self.assertFalse(self.original.exists());self.assertTrue(prepared.exists())
        successes=self.service.simulations(*self.args)
        self.assertEqual([(r['pdb'],r['rmsd']) for r in successes],[('2ABC',.2),('1ABC',.3)])

    def test_imported_source_asset_and_its_import_configuration_are_removed(self):
        ticket=self.store.prepare_upload(self.token,self.pid,'1ABC.pdb','structures')
        (self.store.staging/ticket).write_text(ATOM+LIG+BAD)
        asset=self.store.finish_upload(self.token,ticket);path=self.store.asset_path(self.token,self.pid,asset)
        project=self.store.project(self.token,self.pid)
        project['pipeline'][0].update(operation='import_results',parameters={'kind':'structures','asset_ids':[asset]})
        project['pipeline'][0].pop('process_all',None)
        self.store.save_pipeline(self.token,self.pid,project['pipeline'],project['revision'])
        run=self.store.get_run(self.token,self.rid)
        run['stages'][0].update(operation='import_results',artifacts=[str(path)])
        run['stages'][1]['input_files']={'base_input_path':[str(path)]}
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),self.rid))
        self.service.remove_failure(*self.args,self.failure()['id'])
        self.assertFalse(path.exists())
        self.assertNotIn(asset,[r['id'] for r in self.store.assets(self.token,self.pid)])
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'][0]['parameters']['asset_ids'],[])
        self.assertTrue(self.original.exists())  # Same name, unrelated provenance.

    def test_failed_row_is_highlighted_with_diagnostics_and_confirmed_removal(self):
        ui=self.ui();result=StageResults(ui,self.pid,self.rid,self.stage,True)
        asyncio.run(result.load());table=result.child
        failed=table.table.rows[0]
        self.assertEqual(failed.color,'#FEF2F2')
        self.assertEqual(failed.cells[4].content.value,'Falha na preparação')
        self.assertEqual(failed.cells[5].content.value,'—')
        self.assertTrue(table.warning.visible)
        self.assertTrue(failed.cells[-3].content.disabled)
        asyncio.run(failed.cells[-1].content.controls[0].on_click(None))
        self.assertTrue(any(getattr(c,'value',None)=='invalid ligand' for c in ui.dialogs[-1].content.content.controls))
        asyncio.run(failed.cells[-1].content.controls[-1].on_click(None))
        dialog=ui.dialogs[-1]
        self.assertTrue(self.original.exists())
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertFalse(self.original.exists());self.assertEqual(len(table.simulations),1)

    def test_retrieval_exclusion_prevents_a_deleted_pdb_from_being_downloaded_again(self):
        from crawlers.complex import PDBComplex
        self.service.remove_failure(*self.args,self.failure()['id'])
        excluded=self.store.project(self.token,self.pid)['pipeline'][0]['parameters']['excluded_pdb_ids']
        with tempfile.TemporaryDirectory() as directory:
            crawler=PDBComplex(directory);crawler.get_pdb_ids_with_filters=Mock(return_value=['1ABC','2ABC'])
            session=Mock();session.get.return_value=SimpleNamespace(text=ATOM+LIG,raise_for_status=lambda:None)
            with patch('requests.Session',return_value=session):crawler.get_pdb_files({'pdb_ids':['1ABC','2ABC'],'excluded_pdb_ids':excluded})
            self.assertEqual(session.get.call_count,1)
            self.assertIn('2ABC.pdb',session.get.call_args.args[0])


if __name__=='__main__':unittest.main()
