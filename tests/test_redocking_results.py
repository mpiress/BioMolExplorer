"""RMSD summaries isolate simulation artifacts and preserve authorized 3D access."""
import asyncio
import io
import json
import time
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock
from uuid import uuid4
import test_results_display as display
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_results import RedockingResults
from biomolexplorer.ui.stage_results import StageResults
from biomolexplorer.ui.redocking_results import RedockingResultsTable
from biomolexplorer.ui.localization import LocalizedPage, Translator
from biomolexplorer.pdb_view import StructureViewers, PREFIX
from biomolexplorer.workspace import AccessDenied

ATOM='HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'

class RedockingResultsTests(unittest.TestCase):
    setUp=display.ResultsDisplayTests.setUp
    tearDown=display.ResultsDisplayTests.tearDown
    ui=display.ResultsDisplayTests.ui

    def stage(self):
        config=new_stage('redocking');root=self.store.project_dir(self.pid)/'runs'/self.rid/'simulations'
        files=[]
        for dataset,score in (('first',0.15),('second',0.25)):
            base=root/dataset/'artifacts';structure=base/'structures'/'Estruturas';prepared=structure/'Prepared'
            poses=base/'Estruturas';prepared.mkdir(parents=True);poses.mkdir()
            metadata=structure/'pdb_codes.csv'
            metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD\n'+
                f'1ABC,LIG,1,A,{score}\n1ABC,OTH,2,A,0.9\n1ABC,BAD,3,A,nan\n1ABC,BAD,4,A,inf\n')
            (prepared/'centers.csv').write_text('1ABC_LIG_1A\n0\n0\n0\n')
            for p in (structure/'1ABC.pdb',prepared/'1ABC_A.dockprep.pdbqt',prepared/'1ABC_LIG_1A.lig.pdb',
                      prepared/'1ABC_LIG_1A.lig.pdbqt',prepared/'1ABC_OTH_2A.lig.pdb',poses/'1ABC_LIG_1A.lig.pdbqt',
                      prepared/'1ABC_AA.dockprep.pdbqt',prepared/'1ABC_LIG_1AA.lig.pdbqt'):
                p.write_text(ATOM)
            files.extend(str(p) for p in base.rglob('*') if p.is_file())
        stage=dict(id=config['id'],name=config['name'],operation='redocking',configuration=config,status='succeeded',artifacts=files)
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        return stage

    def test_rmsd_rows_ignore_invalid_values_and_group_each_simulation(self):
        stage=self.stage();service=RedockingResults(self.store)
        rows=service.simulations(self.token,self.pid,self.rid,stage['id'])
        self.assertEqual([r['rmsd'] for r in rows],[0.15,0.9,0.25,0.9])
        simulation=service.simulation(self.token,self.pid,self.rid,stage['id'],rows[0]['id'])
        self.assertEqual(len(simulation['files']),7)
        self.assertTrue(all('/first/' in f['path'] for f in simulation['files']))
        self.assertFalse(any('OTH' in f['name'] or 'AA' in f['name'] for f in simulation['files']))
        self.assertEqual(sum(f['name']=='1ABC_LIG_1A.lig.pdbqt' for f in simulation['files']),2)
        with zipfile.ZipFile(io.BytesIO(service.archive(self.token,self.pid,self.rid,stage['id'],rows[0]['id']))) as archive:
            self.assertEqual(len(archive.namelist()),7)
            self.assertIn('structures/Estruturas/Prepared/1ABC_LIG_1A.lig.pdbqt',archive.namelist())
            self.assertIn('Estruturas/1ABC_LIG_1A.lig.pdbqt',archive.namelist())

    def test_summary_and_archive_require_current_project_and_simulation(self):
        stage=self.stage();service=RedockingResults(self.store)
        with self.assertRaises(AccessDenied):service.simulations(self.guest,self.pid,self.rid,stage['id'])
        with self.assertRaises(AccessDenied):service.simulation(self.token,self.pid,self.rid,stage['id'],'forged')
        with self.assertRaises(AccessDenied):service.simulations(self.token,'other-project',self.rid,stage['id'])
        stage['status']='failed'
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        self.assertEqual(service.simulations(self.token,self.pid,self.rid,stage['id']),[])

    def test_imported_collections_and_duplicate_artifacts_stay_separate(self):
        stage=self.stage();files=[]
        for collection in ('first','second'):
            root=self.store.project_dir(self.pid)/'imports'/collection
            prepared=root/'Prepared';prepared.mkdir(parents=True);(root/'poses').mkdir()
            metadata=root/'pdb_codes.csv'
            metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD\n1ABC,LIG,1,A,0.1\n../../,LIG,1,A,0.1\n')
            (prepared/'centers.csv').write_text('1ABC_LIG_1A\n0\n0\n0\n')
            for path in (prepared/'1ABC_A.dockprep.pdbqt',prepared/'1ABC_LIG_1A.lig.pdbqt',root/'poses'/'1ABC_LIG_1A.lig.pdbqt'):
                path.write_text(ATOM)
            files.extend(str(p) for p in root.rglob('*') if p.is_file())
        stage['artifacts']=files+files
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        service=RedockingResults(self.store);rows=service.simulations(self.token,self.pid,self.rid,stage['id'])
        self.assertEqual(len(rows),2)
        simulation=service.simulation(self.token,self.pid,self.rid,stage['id'],rows[0]['id'])
        self.assertEqual(len(simulation['files']),5)
        self.assertTrue(all('/first/' in f['path'] for f in simulation['files']))

    def test_inline_rmsd_table_popup_downloads_and_preview_callbacks(self):
        stage=self.stage();ui=self.ui();ui.preview_artifact=AsyncMock()
        result=StageResults(ui,self.pid,self.rid,stage,False);asyncio.run(result.load())
        self.assertIsInstance(result.child,RedockingResultsTable)
        table=result.child;self.assertEqual(len(table.table.rows),4)
        asyncio.run(table.table.rows[0].cells[-1].content.on_click(None))
        dialog=ui.dialogs[-1];self.assertIn('0.150',dialog.title.value)
        file_table=table.file_tables[-1]
        molecular=next(row for row in file_table.table.rows if len(row.cells[-1].content.controls)==2)
        asyncio.run(molecular.cells[-1].content.controls[-1].on_click(None))
        ui.preview_artifact.assert_awaited_once()
        asyncio.run(molecular.cells[-1].content.controls[0].on_click(None))
        self.assertTrue(ui.saved[-1]['src_bytes'])
        asyncio.run(dialog.actions[0].on_click(None))
        self.assertTrue(ui.saved[-1]['file_name'].endswith('_redocking.zip'))
        with zipfile.ZipFile(io.BytesIO(ui.saved[-1]['src_bytes'])) as archive:self.assertEqual(len(archive.namelist()),7)
        result.close();self.assertFalse(file_table.valid())
        count=len(ui.saved);asyncio.run(dialog.actions[0].on_click(None));self.assertEqual(len(ui.saved),count)

    def test_ui_localizes_summary_and_dialog_without_changing_identifiers(self):
        stage=self.stage();ui=self.ui();result=StageResults(ui,self.pid,self.rid,stage,True)
        asyncio.run(result.load());page=LocalizedPage(SimpleNamespace(), 'en');page.localize(result.root)
        table=result.child
        self.assertEqual(table.table.columns[1].label.value,'Ligand')
        self.assertEqual(table.table.rows[0].cells[0].content.value,'1ABC')
        self.assertEqual(table.table.rows[0].cells[-1].content.content,'View simulation')
        self.assertIn('simulations',table.count.value)
        asyncio.run(table.open(table.simulations[0]['id']));page.localize(ui.dialogs[-1])
        self.assertEqual(ui.dialogs[-1].actions[0].content,'Download all (ZIP)')

    def test_molecular_viewer_serves_pdbqt_and_mol2_with_correct_format(self):
        stage=self.stage();viewers=StructureViewers(self.store)
        try:
            for extension,content in (('pdbqt',ATOM),('mol2','@<TRIPOS>MOLECULE\nligand\n')):
                path=self.store.project_dir(self.pid)/('ligand.lig.'+extension);path.write_text(content)
                key=viewers.issue(self.token,self.pid,str(path),'en')
                response=viewers.response(PREFIX+'/'+key)
                self.assertEqual(response[0],200)
                self.assertIn('"format": "'+extension+'"',response[2].decode())
                self.assertIn('"representation": "sticks"',response[2].decode())
                self.assertEqual(viewers.response(PREFIX+'/'+key+'/structure')[2],content.encode())
            self.store.logout(self.token)
            self.assertEqual(viewers.response(PREFIX+'/'+key+'/structure')[0],403)
        finally:viewers.close()

    def test_row_overlay_and_residue_actions_use_the_selected_simulation(self):
        from test_docking_scene import atom
        stage=self.stage()
        raw=next(Path(p) for p in stage['artifacts'] if '/first/' in p and Path(p).name=='1ABC.pdb')
        raw.write_text(atom(1,'CA','ALA','A',10,3)+atom(2,'C1','LIG','A',1,0,record='HETATM'))
        ui=self.ui();ui.preview_docking=AsyncMock()
        result=StageResults(ui,self.pid,self.rid,stage,False);asyncio.run(result.load())
        table=result.child
        asyncio.run(table.table.rows[0].cells[-3].content.on_click(None))
        self.assertEqual(ui.preview_docking.await_args.args[3],'redocking')
        self.assertEqual(ui.preview_docking.await_args.args[4],table.simulations[0]['id'])
        asyncio.run(table.table.rows[0].cells[-2].content.on_click(None))
        dialog=ui.dialogs[-1];self.assertEqual(dialog.title.value,'Resíduos próximos ao ligante')
        asyncio.run(dialog.actions[1].on_click(None))
        self.assertEqual(ui.saved[-1]['file_name'],'residue_contacts.csv')
        self.assertIn(b'ALA',ui.saved[-1]['src_bytes'])
        result.close()
        calls=ui.preview_docking.await_count
        asyncio.run(dialog.actions[0].on_click(None))
        self.assertEqual(ui.preview_docking.await_count,calls)
