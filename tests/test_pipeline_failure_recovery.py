"""Regressions from the 5NN0 export and legacy ADMET consumer failures."""
import asyncio
import json
import shutil
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from biomolexplorer.admet_config import resolve_config
from biomolexplorer.admet_outputs import export_standard
from biomolexplorer.artifacts import validate_reference,slot_id
from biomolexplorer.catalog import new_stage
import test_pipeline_execution as execution
import test_redocking_incremental as incremental


class LegacyAdmetSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.stage='a'*32;self.root=Path(self.temp.name)/self.stage/'individual/001_B01__chembl__mol'
        self.root.mkdir(parents=True)
        self.evaluated=pd.DataFrame([
            dict(molecule_chembl_id='A',canonical_smiles='CCO',BBB='BBB+',HIA='HIA-',Lipinski_violations=0,Veber_violations=0,StructuralAlert_count=0,StructuralAlerts=''),
            dict(molecule_chembl_id='B',canonical_smiles='CCCC',BBB='BBB-',HIA='HIA+',Lipinski_violations=1,Veber_violations=0,StructuralAlert_count=1,StructuralAlerts='BRENK: test'),
        ])
        self.config=resolve_config(dict(admet_profile='cns',rule_set='oral'))
        self.stem='ACHE_MOL';self.outputs=export_standard(self.evaluated,self.config,self.root,self.stem)
        report=self.root/(self.stem+'_admet_report.json')
        report.write_text(json.dumps(dict(schema_version=4,outputs=self.outputs,output_stem=self.stem,configuration=self.config)))
        self.files=[str(p) for p in self.root.iterdir()]

    def ref(self,role):
        name='B04__admet__'+role+'.csv';selector='001_B01__chembl__mol/'+name
        return dict(stage=self.stage,selector=selector,artifact_id=slot_id(self.stage,'individual/'+selector))

    def test_legacy_endpoint_rules_and_alert_slots_resolve_to_declared_outputs(self):
        for old,role in [('bbb-pos','BBB+'),('bbb-neg','BBB-'),('hia-pos','HIA+'),('hia-neg','HIA-'),('rules-lipinski-veber','rules-lipinski-veber'),('alerts','alerts')]:
            self.assertEqual(validate_reference(self.ref(old),self.files),str(self.root/self.outputs[role]))

    def test_partitions_and_rules_do_not_inherit_the_bbb_filter(self):
        self.assertEqual(pd.read_csv(self.root/self.outputs['BBB-'])['molecule_chembl_id'].tolist(),['B'])
        self.assertEqual(pd.read_csv(self.root/self.outputs['HIA+'])['molecule_chembl_id'].tolist(),['B'])
        self.assertEqual(pd.read_csv(self.root/self.outputs['rules-lipinski-veber'])['molecule_chembl_id'].tolist(),['A'])
        alerts=pd.read_csv(self.root/self.outputs['alerts'])
        self.assertEqual(alerts['Alert_BRENK'].tolist(),[False,True])
        self.assertEqual(alerts['Alert_PAINS'].tolist(),[False,False])

    def test_forged_identity_other_batch_and_missing_rule_remain_errors(self):
        with self.assertRaises(ValueError):validate_reference(dict(self.ref('bbb-pos'),artifact_id='f'*64),self.files)
        wrong=self.ref('bbb-pos');wrong['selector']=wrong['selector'].replace('001_','002_')
        wrong['artifact_id']=slot_id(self.stage,'individual/'+wrong['selector'])
        with self.assertRaises(ValueError):validate_reference(wrong,self.files)
        with self.assertRaises(ValueError):validate_reference(self.ref('rules-veber'),self.files)
        without=[p for p in self.files if p!=str(self.root/self.outputs['BBB+'])]
        with self.assertRaises(ValueError):validate_reference(self.ref('bbb-pos'),without)


class NativeLigandExportTests(unittest.TestCase):
    def test_old_5nn0_export_with_missing_hydrogens_is_rejected_without_changing_bytes(self):
        from biomolexplorer.ligand_mol2 import ensure_ligand_mol2
        from biomolexplorer.dock6_validation import validate_mol2
        fixture=Path(__file__).parent/'fixtures/redocking'
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'5NN0_92H_627A.lig.mol2';shutil.copy2(fixture/path.name,path)
            cache=root/'cache/calculations/ccd';cache.mkdir(parents=True);shutil.copy2(fixture/'92H.cif',cache/'92H.cif')
            def atoms(text):return text.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0]
            before=atoms(path.read_text())
            with patch.dict('os.environ',{'BIOMOL_CACHE_DIR':str(root/'cache')}),patch('requests.get',side_effect=AssertionError('Cached definition must be reused')):
                with self.assertRaisesRegex(ValueError,'hidrogênios'):ensure_ligand_mol2(path)
            self.assertEqual(atoms(path.read_text()),before)
            self.assertFalse(path.with_suffix('.export.json').exists())

    def test_ccd_input_corrects_saturated_carbon_before_hydrogen_and_charge_assignment(self):
        from biomolexplorer.ligand_mol2 import prepare_chimera_ligand
        from biomolexplorer.dock6_validation import validate_mol2
        fixture=Path(__file__).parent/'fixtures/redocking'
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);cache=root/'cache/calculations/ccd';cache.mkdir(parents=True)
            shutil.copy2(fixture/'92H.cif',cache/'92H.cif')
            atoms=(fixture/'5NN0_92H_627A.lig.mol2').read_text().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0]
            lines=[]
            for line in atoms.splitlines():
                fields=line.split()
                if fields[5]=='H':continue
                i,name=int(fields[0]),fields[1];x,y,z=map(float,fields[2:5]);element=fields[5].split('.')[0]
                lines.append(f'HETATM{i:5d} {name:>4s} 92H A 627    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}')
            complex=root/'5NN0_A.complex.pdb';complex.write_text('\n'.join(lines)+'\n')
            output=root/'5NN0_92H_627A.lig.mol2';script=root/'prepare_ligand_5NN0_92H_627A.com'
            script.write_text(f'open {complex}\nselect :627.A\nselect invert\ndelete selected\ndelete element.H\naddh\naddcharge all spec #0 chargeModel 14sb method am1\nwrite format mol2 #0 {output}\n')
            with patch.dict('os.environ',{'BIOMOL_CACHE_DIR':str(root/'cache')}):
                self.assertTrue(prepare_chimera_ligand(script))
            molecule=validate_mol2(output.with_suffix('.input.mol2'))
            by_name={a.GetProp('_TriposAtomName'):a for a in molecule.GetAtoms()}
            self.assertEqual(by_name['CAM'].GetProp('_TriposAtomType'),'C.3')
            self.assertEqual(by_name['CAT'].GetProp('_TriposAtomType'),'C.3')
            self.assertIn('addcharge',script.read_text())
            self.assertNotIn('select invert',script.read_text())

    def test_pdbqt_atom_names_are_restored_before_docking_without_moving_atoms(self):
        from biomolexplorer.ligand_mol2 import restore_pdbqt_names
        fixture=Path(__file__).parent/'fixtures/redocking/5NN0_92H_627A.lig.mol2'
        atoms=fixture.read_text().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0]
        heavy=[line.split() for line in atoms.splitlines() if line.split()[5]!='H']
        with tempfile.TemporaryDirectory() as folder:
            destination=Path(folder)/'ligand.pdbqt';lines=[]
            for index,a in enumerate(reversed(heavy),1):
                x,y,z=map(float,a[2:5]);element=a[5].split('.')[0]
                lines.append(f'ATOM  {index:5d}  {element:<3s} UNL     1    {x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00     0.123 {element}')
            destination.write_text('\n'.join(lines)+'\n')
            restore_pdbqt_names(fixture,destination)
            restored=destination.read_text().splitlines()
            self.assertEqual([line[12:16].strip() for line in restored],[a[1] for a in reversed(heavy)])
            self.assertEqual([line[16:] for line in restored],[line[16:] for line in lines])
            destination.write_text('\n'.join(lines[:-1])+'\n')
            with self.assertRaisesRegex(ValueError,'Mapeamento'):restore_pdbqt_names(fixture,destination)

    def test_template_exports_charges_without_reopening_a_pdb(self):
        from biomolexplorer.templates import RESOURCE_ROOT
        text=(RESOURCE_ROOT/'chimera/prepare_ligand.template').read_text()
        self.assertLess(text.index('write format mol2'),text.index('close session'))
        self.assertNotIn('open {output_ligand}.lig.pdb',text)


class PairFailureIsolationTests(unittest.TestCase):
    setUp=incremental.IncrementalRedockingTests.setUp
    run_pairs=incremental.IncrementalRedockingTests.run_pairs

    def test_a_vina_or_rmsd_failure_does_not_abort_remaining_pairs(self):
        from wrappers import redocking
        original=redocking.DockVina.return_value.redocking.side_effect
        def run(**kwargs):
            frame=redocking.DockVina.call_args.kwargs['pdb_codes']
            if frame.iloc[0]['LIGAND']=='BAD':raise ValueError('invalid ligand topology during RMSD')
            return original(**kwargs)
        redocking.DockVina.return_value.redocking.side_effect=run
        _,_,rows,docked=self.run_pairs([['1ABC','LIG',1,'A'],['1ABC','BAD',3,'A'],['2ABC','LIG',1,'A']])
        self.assertEqual([r['STATUS'] for r in rows],['succeeded','preparation_failed','succeeded'])
        self.assertIn('RMSD',rows[1]['FAILED_STEP'])
        self.assertEqual(docked,[('1ABC','LIG'),('2ABC','LIG')])

    def test_vina_does_not_rerun_scripts_from_previously_processed_pairs(self):
        from caad.docking import DockVina
        from unittest.mock import Mock
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'previous.vina').write_text('old');(root/'current.vina').write_text('new')
            dock=object.__new__(DockVina);dock.outputpath=folder;dock.perform_subprocess=Mock()
            dock.perform_vina_evaluation(files=['current.vina'])
            self.assertEqual(dock.perform_subprocess.call_count,1)
            self.assertEqual(dock.perform_subprocess.call_args.args[0],['vina','--config','current.vina'])


class TerminalRunHistoryTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown

    def test_terminal_update_uses_its_existing_transaction_to_read_history(self):
        from biomolexplorer.pipeline import PipelineService
        stage=new_stage('retrieve_compounds');run_id='b'*32
        user=self.store.user(self.token)['id']
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(run_id,self.project_id,user,'running','[]',time.time(),time.time(),None))
        original=self.store.project_dir
        def project_dir(pid,db=None):
            self.assertIsNotNone(db,'A nested reader can deadlock after SQLite spills a large run update')
            return original(pid,db)
        service=PipelineService(self.store)
        try:
            with patch.object(self.store,'project_dir',side_effect=project_dir):
                service._update(run_id,'failed',[dict(stage,status='failed',configuration=stage,artifacts=[])],'pair error')
            result=self.store.get_run(self.token,run_id)
            self.assertEqual(result['status'],'failed')
        finally:service.close()


class IndependentBranchReuseTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_removing_redocking_outputs_preserves_pubchem_and_fingerprint_reuse(self):
        import sys
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.result_files import ResultFiles
        from biomolexplorer.stage_cache import artifact_manifest
        source=new_stage('retrieve_compounds');pubchem=new_stage('retrieve_pubchem');fingerprints=new_stage('fingerprints')
        pubchem['bindings']={'base_input_path':{'stage':source['id']}}
        fingerprints['bindings']={'base_input_path':{'stage':source['id']}}
        self.save([source,pubchem,fingerprints]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(first['status'],'succeeded',first.get('error'))
                redocking=new_stage('redocking');path=self.store.project_dir(self.project_id)/'redocking-example.pdbqt';path.write_text('pose')
                historical=dict(redocking,status='succeeded',configuration=redocking,artifacts=[str(path)],artifact_manifest=artifact_manifest([path]))
                stages=first['stages']+[historical]
                with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),first['id']))
                ResultFiles(self.store).remove_all(self.token,self.project_id,first['id'],redocking['id'])
                self.assertFalse(path.exists())
                before=len(calls)
                second=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(second['status'],'succeeded',second.get('error'))
                self.assertEqual(len(calls),before)
                for item in second['stages']:self.assertTrue(item.get('reused'),item['operation'])
            finally:service.close()
