"""Project authorization and execution contracts, independent of Flet."""
import csv
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1] / 'src'))
from biomolexplorer.catalog import new_stage, operation_fields
from biomolexplorer.operations import OPERATIONS
from biomolexplorer.pipeline import PipelineService, select_input, validate_pipeline
from biomolexplorer.templates import RESOURCE_ROOT, validate_templates, materialize_templates
from biomolexplorer.workspace import AccessDenied, WorkspaceStore


class WorkspaceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='biomol workspace ')
        self.store=WorkspaceStore(self.temp.name)
        self.owner=self.store.register('Owner','owner@example.org','owner-password')
        self.guest=self.store.register('Guest','guest@example.org','guest-password')
        self.project=self.store.create_project(self.owner,'Exploração',tags=['MAO'])
        self.pid=self.project['id']

    def tearDown(self):
        self.temp.cleanup()

    def upload(self,name='compounds.csv',content='name,smiles\nETHANOL,CCO\n',kind='compounds'):
        ticket=self.store.prepare_upload(self.owner,self.pid,name,kind)
        (self.store.staging / ticket).write_text(content)
        return self.store.finish_upload(self.owner,ticket)

    def test_passwords_sessions_and_logout(self):
        with self.store.connect() as db:
            user=db.execute('SELECT password FROM users LIMIT 1').fetchone()
            self.assertNotEqual(user['password'],'owner-password')
            self.assertIsNone(db.execute('SELECT 1 FROM sessions WHERE token=?',(self.owner,)).fetchone())
        with self.assertRaises(AccessDenied):
            self.store.login('owner@example.org','incorrect')
        self.assertEqual(self.store.login('OWNER@example.org','owner-password') is not None,True)
        self.store.logout(self.owner)
        with self.assertRaises(AccessDenied):self.store.user(self.owner)
        with self.assertRaises(AccessDenied):self.store.user(None)

    def test_isolation_invitation_roles_and_revocation(self):
        self.assertEqual(self.store.list_projects(self.guest),[])
        with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)
        self.assertEqual(len(self.store.invitations(self.guest)),1)
        self.store.accept_invitation(self.guest,self.pid)
        self.assertEqual(self.store.project(self.guest,self.pid)['role'],'viewer')
        with self.assertRaises(AccessDenied):self.store.save_pipeline(self.guest,self.pid,[],0)
        with self.assertRaises(AccessDenied):self.store.prepare_upload(self.guest,self.pid,'x.csv','compounds')
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        self.store.accept_invitation(self.guest,self.pid)
        self.store.save_pipeline(self.guest,self.pid,[],0)
        with self.assertRaises(AccessDenied):self.store.delete_project(self.guest,self.pid)
        with self.assertRaises(AccessDenied):self.store.invite(self.guest,self.pid,'owner@example.org')
        guest_id=self.store.user(self.guest)['id']
        self.store.revoke(self.owner,self.pid,guest_id)
        with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)

    def test_metadata_archive_delete_and_revision_conflict(self):
        self.store.update_project(self.owner,self.pid,'Novo','Descrição','#6366F1',['teste'],True)
        self.assertEqual(self.store.list_projects(self.owner),[])
        self.assertEqual(self.store.list_projects(self.owner,True)[0]['name'],'Novo')
        self.store.save_pipeline(self.owner,self.pid,[],0)
        with self.assertRaises(ValueError):self.store.save_pipeline(self.owner,self.pid,[],0)
        self.store.delete_project(self.owner,self.pid)
        with self.assertRaises(AccessDenied):self.store.project(self.owner,self.pid)

    def test_invitation_requires_pending_consent_and_cannot_resurrect_revoked_access(self):
        with self.assertRaises(AccessDenied):self.store.accept_invitation(self.guest,self.pid)
        for role in ('viewer','editor'):
            self.store.invite(self.owner,self.pid,' GUEST@example.org ',role)
            self.assertEqual(self.store.invitations(self.guest)[0]['role'],role)
            self.store.accept_invitation(self.guest,self.pid,False)
            self.assertEqual(self.store.invitations(self.guest),[])
            with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        self.store.accept_invitation(self.guest,self.pid)
        with self.assertRaises(ValueError):self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        self.assertEqual(self.store.project(self.guest,self.pid)['role'],'viewer')
        with self.assertRaises(AccessDenied):self.store.accept_invitation(self.guest,self.pid,False)
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)
        self.store.revoke(self.owner,self.pid,self.store.user(self.guest)['id'])
        with self.assertRaises(AccessDenied):self.store.accept_invitation(self.guest,self.pid)
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        self.store.delete_project(self.owner,self.pid)
        self.assertEqual(self.store.invitations(self.guest),[])
        with self.assertRaises(AccessDenied):self.store.accept_invitation(self.guest,self.pid)

    def test_readers_can_inspect_results_but_only_editors_can_modify_and_execute(self):
        asset=self.upload()
        first=new_stage('import_results');first['parameters']['asset_ids']=[asset]
        self.store.save_pipeline(self.owner,self.pid,[first],0)
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        self.store.accept_invitation(self.guest,self.pid)
        file=self.store.asset_path(self.guest,self.pid,asset)
        self.assertIn(b'ETHANOL',self.store.read_file(self.guest,self.pid,file))
        service=PipelineService(self.store,worker_python=Path(sys.executable))
        try:
            with self.assertRaises(AccessDenied):service.submit(self.guest,self.pid)
            with self.assertRaises(AccessDenied):self.store.update_project(self.guest,self.pid,'Bad','','#6366F1',[])
            with self.assertRaises(AccessDenied):self.store.members(self.guest,self.pid)
            with self.assertRaises(AccessDenied):self.store.revoke(self.guest,self.pid,self.store.user(self.owner)['id'])
            self.store.invite(self.owner,self.pid,'guest@example.org','editor')
            self.store.accept_invitation(self.guest,self.pid)
            run=self.wait_run(service,service.submit(self.guest,self.pid))
            self.assertEqual(run['status'],'succeeded',run['error'])
            self.assertEqual(self.store.get_run(self.guest,run['id'])['id'],run['id'])
            self.assertEqual(len(self.store.list_runs(self.guest,self.pid)),1)
            ticket=self.store.prepare_upload(self.guest,self.pid,'guest.csv','compounds')
            (self.store.staging/ticket).write_text('name,smiles\nGUEST,CCO\n')
            self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
            self.store.accept_invitation(self.guest,self.pid)
            with self.assertRaises(AccessDenied):self.store.finish_upload(self.guest,ticket)
            with self.assertRaises(AccessDenied):service.cancel(self.guest,run['id'])
            self.assertTrue(self.store.read_file(self.guest,self.pid,run['stages'][0]['artifacts'][0]))
            self.store.revoke(self.owner,self.pid,self.store.user(self.guest)['id'])
            with self.assertRaises(AccessDenied):self.store.get_run(self.guest,run['id'])
            with self.assertRaises(AccessDenied):self.store.read_file(self.guest,self.pid,file)
            with self.assertRaises(AccessDenied):self.store.assets(self.guest,self.pid)
        finally:service.close()

    def test_invite_rejects_owner_unknown_account_invalid_role_and_cross_project_access(self):
        for email,role in [('owner@example.org','viewer'),('unknown@example.org','viewer'),('guest@example.org','owner')]:
            with self.assertRaises(ValueError):self.store.invite(self.owner,self.pid,email,role)
        stranger=self.store.register('Other','other@example.org','other-password')
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        with self.assertRaises(AccessDenied):self.store.accept_invitation(stranger,self.pid)
        self.store.accept_invitation(self.guest,self.pid)
        other=self.store.create_project(self.owner,'Privado')
        with self.assertRaises(AccessDenied):self.store.project(self.guest,other['id'])

    def test_invitation_role_cannot_change_between_display_and_acceptance(self):
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        displayed=self.store.invitations(self.guest)[0]
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        for accept in (True,False):
            with self.assertRaises(ValueError):
                self.store.accept_invitation(self.guest,self.pid,accept,displayed['role'])
        with self.assertRaises(AccessDenied):self.store.project(self.guest,self.pid)
        self.assertEqual(self.store.invitations(self.guest)[0]['role'],'editor')
        self.store.accept_invitation(self.guest,self.pid,True,'editor')
        self.assertEqual(self.store.project(self.guest,self.pid)['role'],'editor')

    def test_invitation_response_and_expected_role_are_validated(self):
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        for accept,role in (('false','viewer'),(1,'viewer'),(True,'owner'),(None,None)):
            with self.assertRaises(ValueError):self.store.accept_invitation(self.guest,self.pid,accept,role)
        self.assertEqual(len(self.store.invitations(self.guest)),1)

    def test_upload_tickets_cross_project_paths_and_limits(self):
        asset=self.upload()
        file=self.store.asset_path(self.owner,self.pid,asset)
        self.assertIn(b'ETHANOL',self.store.read_file(self.owner,self.pid,file))
        with self.assertRaises(AccessDenied):self.store.asset_path(self.guest,self.pid,asset)
        with self.assertRaises(AccessDenied):self.store.read_file(self.owner,self.pid,Path(self.temp.name)/'workspace.sqlite3')
        other=self.store.create_project(self.owner,'Other')
        with self.assertRaises(AccessDenied):self.store.asset_path(self.owner,other['id'],asset)
        for name in ('../bad.csv','bad\n.csv','bad\\file.csv'):
            with self.assertRaises(ValueError):self.store.prepare_upload(self.owner,self.pid,name,'other')
        ticket=self.store.prepare_upload(self.owner,self.pid,'large.csv','other')
        (self.store.staging/ticket).write_bytes(b'12345')
        self.store.max_upload_bytes=3
        with self.assertRaises(ValueError):self.store.finish_upload(self.owner,ticket)
        ticket=self.store.prepare_upload(self.owner,self.pid,'x.csv','other')
        with self.assertRaises(AccessDenied):self.store.finish_upload(self.guest,ticket)

    def test_symlink_cannot_escape_project(self):
        alias=self.store.project_dir(self.pid)/'escape'
        alias.symlink_to(self.store.database)
        with self.assertRaises(AccessDenied):self.store.read_file(self.owner,self.pid,alias)

    def test_dag_order_cycle_and_unknown_input(self):
        a,b=new_stage('import_results'),new_stage('admet')
        b['bindings']['base_input_path']={'stage':a['id']}
        self.assertEqual(validate_pipeline([b,a]),[a['id'],b['id']])
        a['depends_on']=[b['id']]
        with self.assertRaises(ValueError):validate_pipeline([a,b])
        a['depends_on']=[]
        b['bindings']['base_input_path']={'stage':'f'*32}
        with self.assertRaises(ValueError):validate_pipeline([a,b])

    def test_catalog_matches_all_operations(self):
        for operation,spec in OPERATIONS.items():
            self.assertEqual({f['name'] for f in operation_fields(operation)},set(spec.required+spec.optional)-{'verbose'})

    def test_templates_preserve_contract_and_are_scoped(self):
        name='vina/config.template'
        original=(RESOURCE_ROOT/name).read_text()
        changed=original+'\n# ajuste deste estágio\n'
        validate_templates({name:changed})
        destination=Path(self.temp.name)/'overrides'
        materialize_templates(destination,{name:changed})
        self.assertEqual((destination/name).read_text(),changed)
        self.assertEqual((RESOURCE_ROOT/name).read_text(),original)
        for mapping in ({'../workspace.sqlite3':'bad'},{name:original+'\noutput = /tmp/private\n'},
                        {name:original+'\ncommand; rm anything\n'},
                        {'chimera/prepare_receptor.template':(RESOURCE_ROOT/'chimera/prepare_receptor.template').read_text()+'\nsystem ls\n'}):
            with self.assertRaises(ValueError):validate_templates(mapping)

    def wait_run(self,service,run):
        deadline=time.monotonic()+40
        while time.monotonic()<deadline:
            run=self.store.get_run(self.owner,run['id'])
            if run['status'] not in ('queued','running'):return run
            time.sleep(.05)
        self.fail('Pipeline did not finish within 40 seconds')

    def test_import_normalizes_own_compounds_and_runs_admet_worker(self):
        asset=self.upload('my.compounds.csv')
        first,second=new_stage('import_results'),new_stage('admet')
        first['parameters']['asset_ids']=[asset]
        second['bindings']['base_input_path']={'stage':first['id']}
        self.store.save_pipeline(self.owner,self.pid,[second,first],0)
        with patch.dict(os.environ,{'MPLCONFIGDIR':'/tmp/biomol-mpl'}):
            service=PipelineService(self.store,worker_python=Path(sys.executable))
            try:
                run=self.wait_run(service,service.submit(self.owner,self.pid))
                self.assertEqual(run['status'],'awaiting_input')
                second['bindings']['base_input_path']['selector']='my.compounds.csv'
                run=self.wait_run(service,service.resume(self.owner,run['id'],second))
                self.assertEqual(run['status'],'succeeded',run['error'])
                self.assertTrue(run['stages'][1]['artifacts'])
                imported=Path(run['stages'][0]['artifacts'][0])
                self.assertIn('canonical_smiles',imported.read_text())
                self.assertIn('molecule_chembl_id',imported.read_text())
                self.assertEqual(self.store.asset_path(self.owner,self.pid,asset).read_text(),'name,smiles\nETHANOL,CCO\n')
            finally:service.close()

    def test_prepared_structure_import_bypasses_retrieval(self):
        assets=[self.upload('1ABC_A.dockprep.pdbqt','ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n','prepared_structures'),
                self.upload('centers.csv','1ABC_LIG_1_A\n0\n0\n0\n','prepared_structures'),
                self.upload('pdb_codes.csv','PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,1,A\n','prepared_structures')]
        stage=new_stage('import_results')
        stage['parameters'].update(kind='prepared_structures',target='Target',asset_ids=assets)
        self.store.save_pipeline(self.owner,self.pid,[stage],0)
        service=PipelineService(self.store)
        try:
            run=self.wait_run(service,service.submit(self.owner,self.pid))
            self.assertEqual(run['status'],'succeeded',run['error'])
            root=select_input(run['stages'][0]['artifacts'],'base_input_path')
            self.assertTrue((root/'Target'/'Prepared'/'1ABC_A.dockprep.pdbqt').is_file())
        finally:service.close()

    def test_disabled_stage_and_cross_project_input_are_rejected(self):
        a,b=new_stage('import_results'),new_stage('admet')
        a['enabled']=False
        b['bindings']['base_input_path']={'stage':a['id']}
        self.store.save_pipeline(self.owner,self.pid,[a,b],0)
        service=PipelineService(self.store)
        try:
            with self.assertRaises(ValueError):service.submit(self.owner,self.pid)
            b['bindings']={};b['parameters']={'base_input_path':str(self.store.root)}
            current=self.store.project(self.owner,self.pid)
            self.store.save_pipeline(self.owner,self.pid,[b],current['revision'])
            with self.assertRaises(AccessDenied):service.submit(self.owner,self.pid)
            with self.assertRaises(RuntimeError):PipelineService(self.store)
        finally:service.close()

    def test_consensus_accepts_independent_vina_and_dock6_directories(self):
        from wrappers.docking import generate_consensus
        vina=Path(self.temp.name)/'vina-results'
        dock6=Path(self.temp.name)/'dock6-results'
        output=Path(self.temp.name)/'consensus'
        vina.mkdir();dock6.mkdir();output.mkdir()
        for index in range(3):
            (vina/f'MOL{index}.lig.pdbqt').write_text(f'REMARK VINA RESULT: {-5.0-index}\n')
            (dock6/f'MOL{index}_scored.mol2').write_text(f'Grid_Score: {-20.0-index*2}\nInternal_energy_repulsive: 1.0\n')
        with patch('wrappers.docking.plot_scatter_comparison'):
            generate_consensus(str(Path(self.temp.name)/'unused'),str(output),'Target',
                               base_vina_path=str(vina),base_dock6_path=str(dock6))
        with (output/'Target.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),3)
        self.assertEqual({float(row['vina']) for row in rows},{-5.,-6.,-7.})

    def test_selected_csv_and_expansion_resolve_correct_sources(self):
        root=self.store.project_dir(self.pid)/'source'
        root.mkdir()
        for name in ('molecules.csv','molecules_BBB+.csv'):
            (root/name).write_text('canonical_smiles,molecule_chembl_id\nCCO,A\n')
        upstream=new_stage('retrieve_compounds')
        downstream=new_stage('admet')
        downstream['bindings']['base_input_path']={'stage':upstream['id'],'selector':'molecules_BBB+.csv'}
        results={upstream['id']:[str(p) for p in root.glob('*.csv')]}
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.pid,self.store.user(self.owner)['id'],downstream,results)
            self.assertEqual(params['input_file'],'selected_compounds.csv')
            self.assertTrue((Path(params['base_input_path'])/params['input_file']).is_file())
            original=root/'ChEMBL'/'molecules'/'CHEMBL220'/'original.csv'
            original.parent.mkdir(parents=True)
            original.write_text('canonical_smiles,molecule_chembl_id\nCCO,A\n')
            results[upstream['id']].append(str(original))
            expansion=new_stage('expand_similar_compounds')
            expansion['bindings']['base_input_path']={'stage':upstream['id']}
            params=service._resolve(self.pid,self.store.user(self.owner)['id'],expansion,results)
            self.assertEqual(Path(params['base_input_path']),root)
        finally:service.close()


if __name__=='__main__':unittest.main()
