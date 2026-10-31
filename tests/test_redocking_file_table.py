"""Target grouping, durable ligand review and bounded redocking file rendering."""
import asyncio
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.ui.pdb_results import PDBActions
from biomolexplorer.ui.redocking_files import RedockingFilePicker
import test_retrieval_usability as fixtures


class RedockingTableTests(unittest.TestCase):
    def test_resizing_fills_panel_and_preserves_compact_actions_and_review(self):
        flag=SimpleNamespace(value=True,label='1ABC.pdb',physical_file='/target/1ABC.pdb',
            target_label='Target A',reviewed=True,ligand_summary='ATP · 102 · A',available=True,on_change=None)
        actions=[ft.IconButton(ft.Icons.EDIT),ft.IconButton(ft.Icons.VISIBILITY),ft.IconButton(ft.Icons.OPEN_IN_NEW)]
        page=SimpleNamespace(width=1440,update=Mock())
        picker=RedockingFilePicker([(flag,{'stage':'a','selector':'1ABC.pdb'})],{'stage:a':'PDB A'},page,{id(flag):actions})
        picker.resize_table(SimpleNamespace(width=1160))
        self.assertEqual(len(picker.table.columns),5)
        self.assertEqual(picker.table.width,1160)
        widths=sum(cell.content.width for cell in picker.table.rows[0].cells)
        self.assertAlmostEqual(widths+2*picker.table.horizontal_margin+4*picker.table.column_spacing,1160)
        picker.resize_table(SimpleNamespace(width=560))
        self.assertEqual(len(picker.table.columns),3)
        self.assertEqual(picker.table.width,560)
        self.assertEqual(picker.table.rows[0].color,'#F0FDF4')
        buttons=[button for row in picker.table.rows[0].cells[-1].content.controls for button in row.controls]
        self.assertEqual(buttons[:3],actions)
        before=page.update.call_count
        picker.resize_table(SimpleNamespace(width=560))
        self.assertEqual(page.update.call_count,before)
        self.assertTrue(flag.value)
        buttons[-1].on_click(None)
        self.assertFalse(flag.value)

    def test_selected_table_occurs_once_and_all_pages_are_accessible(self):
        rows=[(SimpleNamespace(value=True,label=f'{i:04d}.pdb',available=True,on_change=None),
               {'stage':'source','selector':f'{i:04d}.pdb'}) for i in range(61)]
        picker=RedockingFilePicker(rows,{'stage:source':'Target'},SimpleNamespace(width=800,update=Mock()))
        # Controls must have one parent: duplicate selected columns broke Flutter
        # layout and replaced the navigation row after the first 25 results.
        self.assertEqual(sum(c is picker.selected for c in picker.control.controls),1)
        navigation=next(c for c in picker.control.controls if picker.next in (getattr(c,'controls',[]) or []))
        self.assertIn(picker.previous,navigation.controls)
        seen=[]
        while True:
            seen.extend(row.cells[0].content.controls[0].value for row in picker.table.rows)
            if picker.next.disabled:break
            picker.next.on_click(None)
        self.assertEqual(len(seen),61)
        self.assertEqual(len(set(seen)),61)
        self.assertTrue(all(flag.value for flag,_ in rows))

    def test_obsolete_selections_can_be_removed_without_losing_current_inputs(self):
        rows=[(SimpleNamespace(value=True,label=f'{i:04d}.pdb',available=i>=30,on_change=None),
               {'stage':'source','selector':f'{i:04d}.pdb'}) for i in range(61)]
        picker=RedockingFilePicker(rows,{'stage:source':'Target'},SimpleNamespace(width=360,update=Mock()))
        picker.review_filter.value='unavailable';picker.filter_changed()
        self.assertEqual(len(picker.selected_rows()),30)
        picker.remove_view()
        self.assertEqual(sum(flag.value for flag,_ in rows),31)
        self.assertTrue(all(flag.value==flag.available for flag,_ in rows))
        picker.review_filter.value='all';picker.filter_changed()
        self.assertEqual(len(picker.selected_rows()),31)
        picker.resize_table(SimpleNamespace(width=248))
        self.assertGreater(picker.table.width,248)
        self.assertEqual(picker.table_host.content.scroll,ft.ScrollMode.ALWAYS)

    def test_current_and_unavailable_files_share_source_group_without_generic_target(self):
        from biomolexplorer.ui.localization import Translator
        rows=[(SimpleNamespace(value=True,label=f'{i}.pdb',available=i==0,on_change=None,
            **({'target_label':'Estruturas'} if i==0 else {})),
            {'stage':'source','selector':f'{i}.pdb'}) for i in range(2)]
        picker=RedockingFilePicker(rows,{'stage:source':'BACE1 PDB'},
            SimpleNamespace(width=800,update=Mock(),translator=Translator('en')))
        self.assertEqual(picker.group(*rows[0]),picker.group(*rows[1]))
        self.assertIn('1 targets/sources',picker.overview.value)
        self.assertIn('1 unavailable',picker.overview.value)
        self.assertIn('of 2 PDBs',picker.count.value)
        self.assertTrue(picker.unavailable_hint.visible)

    def test_large_selection_filters_targets_reviews_and_search_without_losing_files(self):
        updates=[];rows=[]
        for i in range(5000):
            target='Kinase A' if i<2500 else 'Kinase B'
            flag=SimpleNamespace(value=False,label=f'{i:04d}.pdb',physical_file=f'/{target}/{i:04d}.pdb',
                target_label=target,available=True,on_change=None)
            rows.append((flag,{'stage':'a' if i<2500 else 'b','selector':flag.label}))
        created=[]
        def action(index):
            created.append(index);return [ft.IconButton(ft.Icons.EDIT)]
        extras={id(flag):(lambda index=i:action(index)) for i,(flag,_) in enumerate(rows)}
        picker=RedockingFilePicker(rows,{'stage:a':'PDB A','stage:b':'PDB B'},SimpleNamespace(update=lambda:updates.append(1)),extras)
        self.assertEqual(created,[])
        picker.add_all_sources()
        self.assertEqual(len(created),25)
        self.assertEqual(len(picker.table.rows),25)
        self.assertEqual(sum(flag.value for flag,_ in rows),5000)
        self.assertEqual(len(picker.target.options),3)
        picker.target.value='stage:b|Kinase B';picker.filter_changed()
        self.assertEqual(len(picker.selected_rows()),2500)
        picker.search.value='4999';picker.search_changed(None)
        self.assertEqual(len(picker.table.rows),1)
        self.assertEqual(picker.selected_rows()[0][0],4999)
        picker.apply_reviews({'/Kinase B/4999.pdb':{'reviewed':True,'ligands':[{'LIGAND':'ATP','RESNUM':'102','CHAIN':'B'}]}})
        self.assertEqual(picker.table.rows[0].color,'#F0FDF4')
        picker.review_filter.value='pending';picker.filter_changed()
        self.assertEqual(picker.table.rows,[])
        picker.review_filter.value='reviewed';picker.filter_changed()
        self.assertEqual(len(picker.table.rows),1)
        self.assertEqual(sum(flag.value for flag,_ in rows),5000)
        remove=picker.table.rows[0].cells[-1].content.controls[-1]
        remove.on_click(None)
        self.assertFalse(rows[4999][0].value)
        self.assertEqual(sum(flag.value for flag,_ in rows),4999)
        picker.search.value='';picker.review_filter.value='all'
        picker.remove_view()
        self.assertEqual(sum(flag.value for flag,_ in rows),2500)
        self.assertTrue(all(flag.value for flag,_ in rows[:2500]))
        self.assertLessEqual(len(created),51)

    def test_duplicate_pdb_codes_keep_origin_actions_and_selections_when_reopened(self):
        with tempfile.TemporaryDirectory() as directory:
            sources=[new_stage('retrieve_structures') for _ in range(2)]
            run={'id':'run','stages':[]}
            for i,stage in enumerate(sources):
                stage['parameters']['target']='Target '+str(i)
                path=Path(directory)/str(i)/'1ABC.pdb';path.parent.mkdir();path.write_text('END\n')
                run['stages'].append(dict(stage,configuration=stage,status='succeeded',artifacts=[str(path)]))
            stage=new_stage('redocking');stage['bindings']={'base_input_path':{'sources':[{'stage':s['id']} for s in sources]}}
            pending={'configuration':stage,'name':stage['name']};actions=[]
            def file_actions(item,path):
                buttons=[ft.IconButton(ft.Icons.EDIT,on_click=lambda e:None),
                    ft.IconButton(ft.Icons.VISIBILITY,url='https://example.org/view'),
                    ft.IconButton(ft.Icons.OPEN_IN_NEW,url='https://www.rcsb.org/structure/1ABC')]
                actions.append(buttons);return buttons
            form=FileSelection(run,pending,file_actions=file_actions)
            picker=form.pickers['base_input_path'];self.assertIsInstance(picker,RedockingFilePicker)
            picker.add_all_sources()
            self.assertEqual(len(picker.table.rows),2)
            self.assertIs(picker.table.rows[0].cells[-1].content.controls[0],actions[0][0])
            path=run['stages'][0]['artifacts'][0]
            reviews={path:{'reviewed':True,'ligands':[]}}
            form.apply_pdb_reviews(reviews)
            picker.target.value='stage:'+sources[0]['id']+'|Target 0'
            picker.review_filter.value='reviewed';picker.filter_changed()
            state=form.state();selected=form.read()['bindings']['base_input_path']['sources']
            self.assertEqual(len(selected),2)
            self.assertEqual(len({ref['artifact_id'] for ref in selected}),2)
            reopened=FileSelection(run,pending,state=state,file_actions=file_actions)
            reopened.apply_pdb_reviews(reviews)
            restored=reopened.pickers['base_input_path']
            self.assertEqual(restored.target.value,picker.target.value)
            self.assertEqual(restored.review_filter.value,'reviewed')
            self.assertEqual(len(restored.table.rows),1)
            self.assertEqual(reopened.read()['bindings'],form.read()['bindings'])


class PDBReviewTests(unittest.TestCase):
    setUp=fixtures.PDBCurationTests.setUp
    tearDown=fixtures.PDBCurationTests.tearDown
    context=fixtures.PDBCurationTests.context
    save=fixtures.PDBCurationTests.save

    def test_review_is_saved_per_structure_and_reloaded_without_parsing_pdbs(self):
        self.assertFalse(self.service.pdb_review_statuses(self.token,self.pid,self.rid)[str(self.path)]['reviewed'])
        self.assertFalse(self.context()['reviewed'])
        self.assertTrue(self.save([{'LIGAND':'ATP','RESNUM':102,'CHAIN':'B'}])['reviewed'])
        # Editing another structure's metadata must not reset this review.
        self.csv.write_text(self.csv.read_text().replace('2XYZ,ATP,3,B','2XYZ,ATP,4,B'))
        with patch('biomolexplorer.pdb_ligands.analyze_ligands',side_effect=AssertionError('Do not parse bulk PDBs')):
            status=ResultFiles(self.store).pdb_review_statuses(self.token,self.pid,self.rid)[str(self.path)]
        self.assertTrue(status['reviewed']);self.assertEqual(status['ligands'][0]['LIGAND'],'ATP')
        self.path.write_bytes(self.path.read_bytes()+b'REMARK changed\n')
        self.assertFalse(self.service.pdb_review_statuses(self.token,self.pid,self.rid)[str(self.path)]['reviewed'])

    def test_failed_save_does_not_mark_reviewed_and_rolls_back_audit_failure(self):
        original=self.csv.read_bytes()
        with self.assertRaises(ValueError):self.save([{'LIGAND':'BAD','RESNUM':999,'CHAIN':'A'}])
        with patch('biomolexplorer.project_state.record',side_effect=OSError('disk failed')):
            with self.assertRaises(OSError):self.save([])
        self.assertEqual(self.csv.read_bytes(),original)
        self.assertFalse(self.service.pdb_review_statuses(self.token,self.pid,self.rid)[str(self.path)]['reviewed'])
        self.assertNotIn('pdb_reviews',self.store.get_run(self.token,self.rid)['stages'][0])

    def test_editor_notifies_table_only_after_successful_save(self):
        dialogs=[];callback=Mock()
        async def call(function,*args):return function(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(call=call,guard=guard,notify=Mock(),editing_stage=True,
            page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append))
        result=SimpleNamespace(ui=ui,token=self.token,project_id=self.pid,run_id=self.rid,stage=self.stage,
            service=self.service,valid=lambda:True)
        actions=PDBActions(result,True,on_saved=callback);item={'name':self.path.name,'path':str(self.path)}
        asyncio.run(actions.open_ligands(item));callback.assert_not_called()
        asyncio.run(dialogs[-1].actions[-1].on_click(None))
        callback.assert_called_once()
        self.assertTrue(callback.call_args.args[1]['reviewed'])
        self.assertTrue(ui.editing_stage)


if __name__=='__main__':unittest.main()
