"""Shared pair tables, complete receptor catalogs and editable empty inputs."""
import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock,AsyncMock,patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.artifacts import records
from biomolexplorer.md_inputs import MDInputs,receptor_input_files
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.pair_table import PairTable
from biomolexplorer.ui.localization import LocalizedPage,Translator
from test_frontend_navigation import Page


class Store:
    def __init__(self,root):self.root=root
    def project_dir(self,pid):return self.root
    def scoped_path(self,pid,path):
        result=(self.root/path).resolve() if not Path(path).is_absolute() else Path(path).resolve()
        result.relative_to(self.root)
        return result


class AmberPairSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.store=Store(self.root)
        self.producer=new_stage('redocking');self.producer['name']='Curated redocking'
        files=[];self.targets={}
        for target,count,prefix in [('BCHE',12,'B'),('ACHE',12,'A'),('BACE1',10,'C')]:
            for index in range(count):
                code=f'{prefix}{index:03d}';folder=self.root/target/code
                prepared=folder/'Prepared';prepared.mkdir(parents=True)
                for path in (folder/(code+'.pdb'),prepared/(code+'_A.noH.pdb'),
                        prepared/(code+'_A.complex.pdb'),prepared/(code+'_LIG_101A.lig.pdb')):
                    path.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
                    files.append(str(path));self.targets[str(path)]=target
        extra=self.root/'BCHE'/'B000'/'Prepared'/'B000_OTHER_102A.lig.pdb'
        extra.write_text('END\n');files.append(str(extra));self.targets[str(extra)]='BCHE'
        self.item=dict(self.producer,status='succeeded',artifacts=files)
        catalog=records(self.item)
        for record in catalog:
            record['inputs']={'base_input_path':[dict(path=str(self.root/'retrieval'/self.targets[record['path']]/record['name']),
                stage_name=self.targets[record['path']],stage_id='retrieval-'+self.targets[record['path']])]}
        self.item['artifact_records']=catalog
        self.stage=new_stage('dynamics_dock6');self.stage['bindings']={'base_input_path':{'stage':self.producer['id'],'selector':'auto'}}
        ligand=self.root/'lig.mol2';ligand.write_text('@<TRIPOS>MOLECULE\nLIG\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n@<TRIPOS>BOND\n')
        self.stage['parameters']['base_selected_mols']=str(ligand)

    def form(self,language='en',choose_receptors=True):
        self.page=Page();self.ui=SimpleNamespace(current={'id':'project','pipeline':[self.producer,self.stage]},
            store=self.store,service=SimpleNamespace(dock6_path=None),artifact_items={self.producer['id']:self.item},
            page=LocalizedPage(self.page,language),preview_artifact=AsyncMock())
        async def guard(action):await action()
        self.ui.guard=guard
        form=GuidedForm(self.ui,self.stage,[],True)
        if choose_receptors and form.input_editors['base_input_path'].selection:
            form.input_editors['base_input_path'].selection.pickers['base_input_path'].add_all_sources()
        return form

    def test_picker_and_execution_offer_all_34_receptors_without_intermediates(self):
        picker=FileSelection({'stages':[self.item]},{'configuration':self.stage,'name':self.stage['name']})
        rows=picker.rows['base_input_path']
        self.assertEqual(len(self.item['artifacts']),137)
        self.assertEqual(len(rows),34)
        self.assertEqual({self.targets[flag.physical_file] for flag,_ in rows},{'ACHE','BCHE','BACE1'})
        self.assertEqual(len(picker.pickers['base_input_path'].available.options),34)
        expected=set(receptor_input_files(self.item['artifacts']))
        actual=MDInputs(self.store,'project',{self.producer['id']:self.item['artifacts']}).files(self.stage,'base_input_path')
        self.assertEqual({Path(row['path']) for row in actual},expected)
        self.assertTrue(all(p.name.endswith('.noH.pdb') for p in expected))

    def test_preference_is_per_target_and_preserves_other_raw_receptors_and_chains(self):
        base=self.root/'mixed';base.mkdir()
        files=[base/'ACHE.pdb',base/'Prepared/ACHE_A.noH.pdb',base/'Prepared/ACHE_B.noH.pdb',
               base/'BACE1.pdb',base/'Prepared/ACHE_A.prepared.pdb',base/'Prepared/ACHE_LIG_1A.lig.pdb']
        for path in files:
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
        self.assertEqual(set(receptor_input_files(files)),{files[2],files[3],files[4]})

    def test_shared_table_filters_targets_selects_green_and_does_not_refetch(self):
        form=self.form();settings=form.md_settings
        self.assertIsInstance(settings.pair_table,PairTable)
        self.assertEqual(len(settings.pairs),34)
        self.assertEqual({p['receptor']['target'] for p in settings.pairs},{'ACHE','BCHE','BACE1'})
        settings.target.value='ACHE';settings.filter(None)
        self.assertEqual(len(settings.table.rows),12)
        with patch.object(MDInputs,'files',side_effect=AssertionError('Selection must not reread artifacts')):
            settings.all(None)
            self.assertTrue(all(row.selected and row.color=='#F0FDF4' for row in settings.table.rows))
            selected=set(form.readers['selected_pairs']())
            self.assertEqual(len(selected),12)
            settings.target.value='BACE1';settings.filter(None);settings.all(None)
            self.assertEqual(len(form.readers['selected_pairs']()),22)
            self.assertTrue(selected<=set(form.readers['selected_pairs']()))
            settings.pair_table.actions_table.rows[0].cells[0].content.controls[-1].on_click(None)
            self.assertEqual(len(form.readers['selected_pairs']()),21)
        settings.pair_table.resize(SimpleNamespace(width=1200));self.assertEqual(settings.table.width,1088)
        settings.pair_table.resize(SimpleNamespace(width=360));self.assertEqual(settings.table.width,700)

    def test_removing_last_input_is_editable_and_does_not_restore_automatic_sources(self):
        form=self.form();editor=form.input_editors['base_input_path']
        picker=editor.selection.pickers['base_input_path']
        picker.add_all_sources();form.md_settings.all(None)
        self.assertEqual(len(form.md_settings.pairs),34)
        picker.clear()
        self.assertIsNone(editor.read())
        self.assertFalse(form.md_settings.error.visible)
        self.assertEqual(form.md_settings.pairs,[])
        saved=form.read()
        self.assertNotIn('base_input_path',saved['bindings'])
        self.assertEqual(saved['parameters']['selected_pairs'],[])
        self.assertIn('base_input_path',self.stage['bindings'])
        with self.assertRaises(ValueError):MDInputs(self.store,'project',{}).resolve(saved)
        picker.set_selected([0],True)
        self.assertEqual(len(form.md_settings.pairs),1)
        self.assertEqual(len(form.readers['selected_pairs']()),1)

    def test_pagination_stays_bounded_and_search_keeps_selections(self):
        settings=self.form().md_settings
        seed=settings.pairs[0]
        settings.pairs=[dict(seed,id=f'{index:064x}',ligand=dict(seed['ligand'],name=f'Ligand {index}')) for index in range(1000)]
        settings.redraw();self.assertEqual(len(settings.table.rows),50)
        settings.move(1);self.assertEqual(len(settings.table.rows),50)
        settings.search.value='Ligand 999';settings.filter(None)
        self.assertEqual(len(settings.table.rows),1)
        settings.table.rows[0].on_select_change(SimpleNamespace(data='true'))
        self.assertEqual(len(settings.selected),1)
        settings.search.value='';settings.filter(None)
        self.assertEqual(len(settings.selected),1)

    def test_preview_actions_and_messages_follow_session_language(self):
        import asyncio
        form=self.form('en');settings=form.md_settings
        self.ui.page.add(form.layout())
        self.assertEqual(settings.table.columns[0].label.value,'Target')
        self.assertEqual(settings.target.label,'Pair target')
        self.assertEqual(settings.target.options[0].text,'All targets')
        picker=form.input_editors['base_input_path'].selection.pickers['base_input_path']
        self.assertEqual(picker.search.label,'Search files')
        self.assertTrue(picker.count.value.endswith('Page 1'))
        self.assertTrue(settings.page_count.value.startswith('Showing pairs'))
        pair=settings.pairs[0];actions=settings.pair_table.actions_table.rows[0].cells[0].content.controls
        self.assertEqual(len(actions),2)
        asyncio.run(actions[0].on_click(None))
        self.assertEqual(self.ui.preview_artifact.call_args_list[0].args,('project',pair['ligand']['path']))
        picker=form.input_editors['base_input_path'].selection.pickers['base_input_path']
        receptor_actions=picker.pair_table.actions_table.rows[0].cells[0].content.controls
        self.assertEqual(len(receptor_actions),2)
        asyncio.run(receptor_actions[0].on_click(None))
        self.assertEqual(self.ui.preview_artifact.call_args_list[1].args[0],'project')
        self.assertTrue(self.ui.preview_artifact.call_args_list[1].args[1].endswith('.noH.pdb'))
        self.assertNotIn('Residue Number',[column.label.value for column in settings.table.columns])
        self.ui.page.set_language('pt');self.ui.page.update()
        self.assertEqual(settings.table.columns[0].label.value,'Alvo')
        self.assertTrue(settings.page_count.value.startswith('Exibindo pares'))
        message='Selecione ao menos um arquivo para esta entrada.'
        self.assertEqual(Translator('en')(message),'Select at least one file for this input.')
        self.assertEqual(Translator('pt')(message),message)

    def test_removing_a_direct_path_does_not_resurrect_it_during_refresh(self):
        self.stage['parameters']['base_input_path']=str(next(p for p in map(Path,self.item['artifacts']) if p.name.endswith('.noH.pdb')))
        self.stage['bindings']={}
        form=self.form();editor=form.input_editors['base_input_path']
        self.assertEqual(len(form.md_settings.pairs),1)
        editor.rows[0][2].controls[-1].on_click(None)
        self.assertEqual(form.md_settings.pairs,[])
        self.assertNotIn('base_input_path',form.read()['parameters'])
        self.assertIn('base_input_path',self.stage['parameters'])

    def test_draft_clear_preserves_strict_execution_selection_validation(self):
        form=self.form();editor=form.input_editors['base_input_path']
        editor.selection.pickers['base_input_path'].add_all_sources()
        editor.selection.pickers['base_input_path'].clear()
        self.assertIsNone(editor.read())
        with self.assertRaises(ValueError) as error:editor.selection.read(require_selection=True)
        translated=Translator('en')(str(error.exception))
        self.assertTrue(translated.startswith('Select at least one file'),translated)
        self.assertNotIn('Selecione',translated)

    def test_redocking_uses_the_same_table_and_removing_input_updates_rows(self):
        # Raw pair curation receives retrieval structures. Redocking exports
        # only isolated receptors for downstream consumers.
        self.producer['operation']='retrieve_structures';self.item['operation']='retrieve_structures'
        self.stage=new_stage('redocking')
        self.stage['bindings']={'base_input_path':{'stage':self.producer['id'],'selector':'auto'}}
        self.stage['parameters']['pdb_codes']=[['B000','LIG',101,'A']]
        form=self.form();pairs=form.redocking_pairs
        self.assertIsInstance(pairs.pair_table,PairTable)
        editor=form.input_editors['base_input_path'];picker=editor.selection.pickers['base_input_path']
        index=next(i for i,(flag,ref) in enumerate(picker.rows) if Path(flag.physical_file).name=='B000.pdb')
        picker.set_selected([index],True)
        self.assertEqual(len(pairs.table.rows),1)
        picker.clear()
        self.assertEqual(pairs.table.rows,[])
        self.assertEqual(pairs.rows,[])
