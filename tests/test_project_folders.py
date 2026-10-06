"""Destructive project operations are scoped, confirmed and recoverable."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.workspace import WorkspaceStore, AccessDenied


class ProjectFolderTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.store=WorkspaceStore(self.root/'workspace')
        self.token=self.store.register('Owner','owner@example.org','owner-password')
        self.other=self.store.register('Other','other@example.org','other-password')
        self.parent=self.root/'Downloads';self.parent.mkdir()
    def tearDown(self):self.temp.cleanup()
    def create(self,name='Study',**kwargs):
        return self.store.create_project_in_parent(self.token,name,parent=str(self.parent),**kwargs)

    def test_named_folder_preserves_parent_and_siblings(self):
        (self.parent/'other.txt').write_text('keep')
        project=self.create('Estudo com espaços')
        self.assertEqual(project['directory'],str(self.parent/'Estudo com espaços'))
        self.assertTrue((self.parent/'Estudo com espaços'/'project.json').exists())
        self.assertEqual((self.parent/'other.txt').read_text(),'keep')

    def test_replacement_requires_unchanged_confirmation(self):
        old=self.create();folder=Path(old['directory']);(folder/'result.csv').write_text('old')
        plan=self.store.project_destination(self.token,'Study',self.parent)
        with self.assertRaisesRegex(ValueError,'Confirme'):
            self.create()
        (folder/'result.csv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'Confirme'):
            self.create(confirmation=plan)
        plan=self.store.project_destination(self.token,'Study',self.parent)
        new=self.create(confirmation=plan)
        self.assertNotEqual(old['id'],new['id']);self.assertFalse((folder/'result.csv').exists())
        with self.assertRaises(AccessDenied):self.store.project(self.token,old['id'])
        with self.store.connect() as db:
            self.assertIsNone(db.execute('SELECT 1 FROM project_locations WHERE project_id=?',(old['id'],)).fetchone())
            self.assertEqual(db.execute('SELECT COUNT(*) FROM pending_deletions').fetchone()[0],0)

    def test_foreign_owner_active_runs_nested_projects_and_symlinks_are_protected(self):
        old=self.create()
        with self.assertRaises(AccessDenied):self.store.project_destination(self.other,'Study',self.parent)
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',('run',old['id'],self.store.user(self.token)['id'],'awaiting_input','[]',1,1,None))
        with self.assertRaisesRegex(ValueError,'execução'):self.store.project_destination(self.token,'Study',self.parent)
        with self.assertRaisesRegex(ValueError,'execução'):self.store.delete_project(self.token,old['id'])
        with self.store.connect() as db:db.execute('DELETE FROM runs')
        with self.assertRaisesRegex(ValueError,'outro projeto'):self.store.project_destination(self.token,'child',old['directory'])
        outer=self.parent/'Outer';outer.mkdir()
        self.store.create_project(self.token,'Nested',directory=outer/'Nested')
        with self.assertRaisesRegex(ValueError,'outro projeto'):self.store.project_destination(self.token,'Outer',self.parent)
        (self.parent/'Link').symlink_to(Path(old['directory']),target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'simbólicos'):self.store.project_destination(self.token,'Link',self.parent)

    def test_name_cannot_escape_selected_parent(self):
        for name in ('', '..', '../escape', '/tmp/escape', 'a/b', 'a\\b', '.hidden', 'CON', 'bad:folder'):
            with self.assertRaises(ValueError):self.create(name)
        self.assertEqual(list(self.parent.iterdir()),[])

    def test_unregistered_folder_requires_confirmation(self):
        target=self.parent/'Study';target.mkdir();(target/'old.txt').write_text('old')
        with self.assertRaises(ValueError):self.create()
        plan=self.store.project_destination(self.token,'Study',self.parent)
        self.create(confirmation=plan)
        self.assertFalse((target/'old.txt').exists())

    def test_creation_failure_restores_existing_project_and_files(self):
        old=self.create();folder=Path(old['directory']);(folder/'old.txt').write_text('old')
        ticket=self.store.prepare_upload(self.token,old['id'],'input.csv','compounds')
        (self.store.staging/ticket).write_text('pending')
        plan=self.store.project_destination(self.token,'Study',self.parent)
        with patch('biomolexplorer.project_state.record',side_effect=OSError('disk full')):
            with self.assertRaises(OSError):self.create(confirmation=plan)
        self.assertEqual((folder/'old.txt').read_text(),'old')
        self.assertEqual((self.store.staging/ticket).read_text(),'pending')
        self.assertEqual(self.store.project(self.token,old['id'])['name'],'Study')
        self.assertFalse(list(self.parent.glob('.biomol-delete-*')))

    def test_failed_creation_restores_old_folder_even_if_cleanup_fails(self):
        old=self.create();folder=Path(old['directory']);(folder/'old.txt').write_text('preserved')
        plan=self.store.project_destination(self.token,'Study',self.parent)
        with patch('biomolexplorer.project_state.record',side_effect=RuntimeError('creation failed')), \
                patch('biomolexplorer.project_folders.shutil.rmtree',side_effect=PermissionError('cleanup failed')):
            with self.assertLogs('biomolexplorer.project_folders',level='ERROR'):
                with self.assertRaisesRegex(RuntimeError,'creation failed'):self.create(confirmation=plan)
        self.assertEqual((folder/'old.txt').read_text(),'preserved')
        self.assertEqual(self.store.project(self.token,old['id'])['name'],'Study')
        WorkspaceStore(self.store.root)
        self.assertFalse(list(self.parent.glob('.biomol-delete-*')))

    def test_delete_purges_database_staging_and_folder_but_not_symlink_targets(self):
        project=self.create();folder=Path(project['directory'])
        external=self.parent/'external.txt';external.write_text('keep')
        (folder/'link').symlink_to(external)
        ticket=self.store.prepare_upload(self.token,project['id'],'input.csv','compounds')
        (self.store.staging/ticket).write_text('pending')
        with self.assertRaises(AccessDenied):self.store.delete_project(self.other,project['id'])
        self.store.delete_project(self.token,project['id'])
        self.assertFalse(folder.exists());self.assertFalse((self.store.staging/ticket).exists())
        self.assertEqual(external.read_text(),'keep')
        with self.store.connect() as db:
            for table in ('members','assets','uploads','runs','project_history','project_locations'):
                self.assertEqual(db.execute(f'SELECT COUNT(*) FROM {table} WHERE project_id=?',(project['id'],)).fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM projects').fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM pending_deletions').fetchone()[0],0)
        self.create()

    def test_database_failure_during_deletion_restores_folder(self):
        project=self.create();folder=Path(project['directory']);(folder/'result.txt').write_text('keep')
        with patch('biomolexplorer.project_folders.purge_rows',side_effect=RuntimeError('database failed')):
            with self.assertRaises(RuntimeError):self.store.delete_project(self.token,project['id'])
        self.assertEqual((folder/'result.txt').read_text(),'keep')
        self.assertEqual(self.store.project(self.token,project['id'])['directory'],str(folder))
        self.assertFalse(list(self.parent.glob('.biomol-delete-*')))

    def test_deletion_rejects_a_different_path_than_the_confirmation(self):
        project=self.create()
        with self.assertRaisesRegex(ValueError,'mudou'):
            self.store.delete_project(self.token,project['id'],expected_directory=str(self.parent/'Other'))
        self.assertTrue(Path(project['directory']).exists())

    def test_cleanup_failure_can_be_retried_and_resumed_on_restart(self):
        project=self.create()
        with patch('biomolexplorer.project_folders.shutil.rmtree',side_effect=PermissionError('denied')):
            with self.assertRaisesRegex(ValueError,'remoção física'):self.store.delete_project(self.token,project['id'])
        self.assertEqual(self.store.list_projects(self.token),[])
        with self.assertRaises(AccessDenied):self.store.delete_project(self.other,project['id'])
        self.store.delete_project(self.token,project['id'])
        self.assertFalse(list(self.parent.glob('.biomol-delete-*')))
        project=self.create()
        with patch('biomolexplorer.project_folders.shutil.rmtree',side_effect=PermissionError('denied')):
            with self.assertRaises(ValueError):self.store.delete_project(self.token,project['id'])
        WorkspaceStore(self.store.root)
        self.assertFalse(list(self.parent.glob('.biomol-delete-*')))

    def test_missing_directory_and_legacy_deleted_projects_can_be_removed(self):
        import shutil
        project=self.create();shutil.rmtree(project['directory'])
        with self.store.connect() as db:db.execute('UPDATE projects SET deleted=1 WHERE id=?',(project['id'],))
        self.store.delete_project(self.token,project['id'])
        self.create()
