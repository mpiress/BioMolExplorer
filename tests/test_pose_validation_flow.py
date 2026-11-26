"""Targets and compounds of Validate Pose: catalog, canvas, dialog and execution."""
import asyncio
import copy
import json
import subprocess
import sys
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch
from uuid import uuid4

import flet as ft

from biomolexplorer import flow
from biomolexplorer.artifact_registry import register_run
from biomolexplorer.artifacts import index_run
from biomolexplorer.bindings import pack, sources
from biomolexplorer.catalog import TITLES, new_stage
from biomolexplorer.docking_data import write_csv
from biomolexplorer.pipeline import PipelineService, StageAdjustment, execution_stage, requires_curation, validate_pipeline
from biomolexplorer.pose_inputs import PoseInputs
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.guided import GuidedForm
import test_pipeline_editor_regressions as editor_fixtures
from test_frontend_navigation import walk
from test_pose_validation import receptor_pdb, relaxed_molecule, write_mol2

STARTUP = ('import sys\n'
           'import biomolexplorer.ui.pose_settings\n'
           'import biomolexplorer.pose_sources, biomolexplorer.pose_geometry\n'
           'import biomolexplorer.pose_inputs, biomolexplorer.pose_catalog\n'
           'import wrappers.pose_validation\n'
           'from biomolexplorer.catalog import new_stage\n'
           'new_stage("pose_validation")\n'
           'loaded=sorted(name for name in sys.modules if name.split(".")[0] in ("rdkit","scipy","numpy","pandas","matplotlib"))\n'
           'print(";".join(loaded))\n')
SMILES = 'c1ccc(cc1)C(=O)NCCO'


class LazyImportTests(unittest.TestCase):
    def test_opening_the_block_does_not_import_the_scientific_stack(self):
        result = subprocess.run([sys.executable,'-c',STARTUP],capture_output=True,text=True,
                                env={'PYTHONPATH':'src','PYTHONDONTWRITEBYTECODE':'1','MPLBACKEND':'Agg',
                                     'PATH':'/usr/bin:/bin','HOME':str(Path.home())},timeout=300)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(result.stdout.strip(),'',
                         'Scientific modules imported at startup: '+result.stdout.strip())


class CanvasIntegrityTests(unittest.TestCase):
    setUp = editor_fixtures.PipelineEditorRegressionTests.setUp

    def canvas(self, stages):
        self.ui.current['pipeline'] = stages
        self.ui.schedule_save = Mock()
        return FlowCanvas(self.ui,True)

    def test_block_is_offered_in_the_docking_group(self):
        self.assertEqual(TITLES['pose_validation'][0],'Validar poses')
        self.assertEqual(TITLES['pose_validation'][2],'Docking')

    def test_ports_follow_the_docking_naming_and_order(self):
        ports = flow.input_ports(new_stage('pose_validation'))
        self.assertEqual([(port['field'],port['label'],port['required']) for port in ports],
                         [('base_input_path','Targets',True),('base_selected_mols','Compounds',True)])

    def test_targets_come_from_preparation_and_compounds_from_docking(self):
        validation = new_stage('pose_validation')
        for operation in ('redocking','prepare_structures','retrieve_structures'):
            self.assertTrue(flow.compatible(new_stage(operation),validation,'base_input_path'),operation)
        # Retrieval publishes raw complexes only: they are targets, not poses.
        self.assertFalse(flow.compatible(new_stage('retrieve_structures'),validation,'base_selected_mols'))
        for operation in ('docking_vina','docking_dock6','consensus'):
            self.assertTrue(flow.compatible(new_stage(operation),validation,'base_selected_mols'),operation)
        amber = new_stage('dynamics_dock6')
        self.assertTrue(flow.compatible(validation,amber,'base_selected_mols'))
        self.assertTrue(flow.compatible(validation,new_stage('consensus'),'base_vina_path'))
        redocking,vina = new_stage('redocking'),new_stage('docking_vina')
        stages = [redocking,vina,validation,amber]
        flow.connect(stages,redocking['id'],validation['id'],'base_input_path')
        flow.connect(stages,vina['id'],validation['id'],'base_selected_mols')
        flow.connect(stages,validation['id'],amber['id'],'base_selected_mols')
        self.assertEqual(len(validate_pipeline(stages)),4)
        self.assertEqual(validation['bindings']['base_input_path']['output_port'],'structures')
        self.assertEqual(validation['bindings']['base_selected_mols']['output_port'],'docking_results')
        self.assertEqual(amber['bindings']['base_selected_mols']['output_port'],'validated_poses')
        self.assertEqual(flow.stage_issues(new_stage('pose_validation')),['Conecte targets','Conecte compounds'])

    def test_repeating_a_connection_adds_no_edge_history_or_autosave(self):
        vina,validation = new_stage('docking_vina'),new_stage('pose_validation')
        editor = self.canvas([vina,validation])
        asyncio.run(editor.link(vina['id'],validation['id'],'base_selected_mols'))
        self.assertEqual(len(editor.history),1)
        self.ui.schedule_save.assert_called_once()
        asyncio.run(editor.link(vina['id'],validation['id'],'base_selected_mols'))
        self.assertEqual(len(editor.history),1)
        self.ui.schedule_save.assert_called_once()
        self.assertEqual(editor.status.value,'Esta conexão já existe.')
        self.assertEqual(len(sources(validation['bindings']['base_selected_mols'])),1)

    def test_a_curated_choice_is_not_duplicated_by_an_automatic_reconnection(self):
        vina,validation = new_stage('docking_vina'),new_stage('pose_validation')
        stages = [vina,validation]
        flow.connect(stages,vina['id'],validation['id'],'base_selected_mols')
        validation['bindings']['base_selected_mols']['selector'] = 'docking_results.csv'
        flow.connect(stages,vina['id'],validation['id'],'base_selected_mols')
        references = sources(validation['bindings']['base_selected_mols'])
        self.assertEqual(len(references),1)
        self.assertEqual(references[0]['selector'],'docking_results.csv')
        curated = {'stage':vina['id'],'selector':'docking_results.csv','output_port':'docking_results'}
        automatic = {'stage':vina['id'],'selector':'auto','output_port':'docking_results'}
        self.assertEqual(pack([curated,automatic],collapse=True),curated)
        self.assertEqual(len(sources(pack([curated,automatic]))),2)


class PublishedResultsFixture(unittest.TestCase):
    """A redocking and a Vina block published and indexed like a real run."""
    setUp = editor_fixtures.PipelineEditorRegressionTests.setUp

    def publish(self, targets=('4M0E','2WSL'), compounds=('CHEMBL92663','CHEMBL100020')):
        root = self.store.project_dir(self.pid)/'runs'/uuid4().hex
        molecule = relaxed_molecule(SMILES)
        redocking = new_stage('redocking'); vina = new_stage('docking_vina')
        prepared = root/'redocking/artifacts/structures/Estruturas/Prepared'; prepared.mkdir(parents=True)
        files = []
        for code in targets:
            receptor_pdb(prepared/f'{code}_A.dockprep.pdbqt',[[80.0,0.0,0.0],[83.5,0.0,0.0]])
            files.append(str(prepared/f'{code}_A.dockprep.pdbqt'))
        metadata = prepared.parent/'pdb_codes.csv'
        write_csv(metadata,[dict(PDB_CODE=code,LIGAND='1YL',RESNUM=604,CHAIN='A',RESOLUTION='2.0',
                                 RMSD=0.2,STATUS='succeeded') for code in targets])
        files.append(str(metadata))
        poses = root/'vina/artifacts/Estruturas/Vina'; poses.mkdir(parents=True)
        receptors = root/'vina/artifacts/receptors'; receptors.mkdir(parents=True)
        rows = []; published = []
        for code in targets:
            receptor_pdb(receptors/f'{code}_A.noH.pdb',[[80.0,0.0,0.0],[83.5,0.0,0.0]])
            published.append(str(receptors/f'{code}_A.noH.pdb'))
            for compound in compounds:
                name = f'{code}_1YL_604A_{compound}.lig.mol2'
                write_mol2(poses/name,molecule)
                published.append(str(poses/name))
                rows.append(dict(molecule_chembl_id=compound,canonical_smiles=SMILES,receptor_id=code+'_A',
                                 engine='vina',score=-10.5,conformer_file='Estruturas/Vina/'+name,
                                 receptor_file=f'receptors/{code}_A.noH.pdb'))
        table = root/'vina/artifacts/docking_results.csv'
        write_csv(table,rows)
        published.append(str(table))
        stages = [dict(id=redocking['id'],name=redocking['name'],operation='redocking',configuration=redocking,
                       status='succeeded',artifacts=files),
                  dict(id=vina['id'],name=vina['name'],operation='docking_vina',configuration=vina,
                       status='succeeded',artifacts=published)]
        run_id = uuid4().hex
        index_run(stages,run_id)
        with self.store.connect() as db:
            register_run(db,self.pid,run_id,stages)
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                       (run_id,self.pid,self.store.user(self.token)['id'],'succeeded',
                        json.dumps(stages),time.time(),time.time(),None))
        validation = new_stage('pose_validation')
        pipeline = [redocking,vina,validation]
        flow.connect(pipeline,redocking['id'],validation['id'],'base_input_path')
        flow.connect(pipeline,vina['id'],validation['id'],'base_selected_mols')
        self.ui.current['pipeline'] = pipeline
        return SimpleNamespace(redocking=redocking,vina=vina,validation=validation,run_id=run_id,
                               table=table,poses=poses,receptors=receptors,prepared=prepared)


class CatalogTests(PublishedResultsFixture):
    def test_targets_are_listed_by_pdb_code_with_their_prepared_receptor(self):
        scene = self.publish()
        targets = PoseInputs(self.store,self.pid).targets(scene.validation)
        self.assertEqual([target['pdb_code'] for target in targets],['2WSL','4M0E'])
        first = next(target for target in targets if target['pdb_code'] == '4M0E')
        self.assertEqual(first['receptor_id'],'4M0E_A')
        self.assertEqual(first['chain'],'A')
        self.assertEqual(first['ligand'],'1YL')
        self.assertEqual(first['label'],'4M0E · cadeia A')
        self.assertTrue(first['file'].endswith('4M0E_A.dockprep.pdbqt'))
        self.assertEqual(first['rmsd'],'0.2')

    def test_compounds_come_from_the_score_table_not_from_the_pose_files(self):
        scene = self.publish()
        inputs = PoseInputs(self.store,self.pid)
        compounds = inputs.compounds(scene.validation)
        self.assertEqual([entry['compound_id'] for entry in compounds],['CHEMBL100020','CHEMBL92663'])
        entry = compounds[1]
        self.assertEqual(entry['receptors'],['2WSL_A','4M0E_A'])
        self.assertEqual(entry['engines'],['vina'])
        self.assertEqual(entry['smiles'],SMILES)
        self.assertEqual(entry['poses'],2)
        self.assertEqual(entry['best_score'],-10.5)
        # Only the consolidated table is read, never the pose tree.
        self.assertEqual({Path(record['table']).name for record in inputs.records(scene.validation)},
                         {'docking_results.csv'})

    def test_compounds_are_restricted_to_the_selected_targets(self):
        scene = self.publish(compounds=('CHEMBL92663',))
        inputs = PoseInputs(self.store,self.pid)
        compounds = inputs.compounds(scene.validation,['4M0E_A'])
        self.assertEqual([entry['receptors'] for entry in compounds],[['4M0E_A']])
        self.assertEqual(inputs.compounds(scene.validation,['9ZZZ_A']),[])

    def test_the_catalog_does_not_read_execution_json(self):
        scene = self.publish()
        inputs = PoseInputs(self.store,self.pid)
        with patch.object(self.store,'list_runs',side_effect=AssertionError('run JSON loaded')), \
             patch.object(self.store,'get_run',side_effect=AssertionError('run JSON loaded')):
            self.assertTrue(inputs.targets(scene.validation))
            self.assertTrue(inputs.compounds(scene.validation))


class ConfigurationTests(PublishedResultsFixture):
    def form(self, scene, writable=True):
        self.ui.artifact_items = {}; self.ui.artifact_choices = {}
        return GuidedForm(self.ui,scene.validation,[],writable)

    def no_repaint(self):
        self.page.update = Mock(side_effect=AssertionError('The block must not repaint the page'))

    def open_dialog(self, scene):
        forms = []

        def capture(*args,**kwargs):
            form = GuidedForm(*args,**kwargs); forms.append(form); return form

        async def guard(action):
            await action()
        self.ui.guard = guard
        with patch('biomolexplorer.ui.guided.GuidedForm',side_effect=capture):
            asyncio.run(self.ui.open_stage_dialog(scene.validation['id']))
        return forms[-1],self.page.dialogs[-1]

    def test_the_dialog_opens_without_loading_every_run(self):
        scene = self.publish()
        with patch.object(self.store,'list_runs',side_effect=AssertionError('run JSON loaded')):
            form,dialog = self.open_dialog(scene)
        self.assertTrue(dialog.open)
        self.assertIsNotNone(dialog.content.height)
        self.assertTrue(dialog.content.content.expand)
        settings = form.pose_settings
        self.assertEqual([t['pdb_code'] for t in settings.targets.entries],['2WSL','4M0E'])
        self.assertEqual([c['compound_id'] for c in settings.compounds.entries],
                         ['CHEMBL100020','CHEMBL92663'])

    def test_targets_and_compounds_are_separate_sections(self):
        scene = self.publish()
        form = self.form(scene)
        titles = [card.content.controls[0].value for card in form.layout().controls
                  if getattr(getattr(card,'content',None),'controls',None)]
        self.assertEqual(titles[:5],['Identificação','Targets','Compounds','Resultados prontos','Parâmetros'])
        self.assertNotIn('Entradas',titles)
        self.assertEqual(titles.count('Targets'),1)
        self.assertEqual(titles.count('Compounds'),1)
        settings = form.pose_settings
        self.assertFalse(settings.targets.table_box.visible)
        self.assertFalse(settings.compounds.table_box.visible)
        self.assertTrue(settings.targets.empty.visible)
        self.assertTrue(settings.compounds.empty.visible)
        # The generic source rows would list producer outputs; the sections own the origin.
        self.assertFalse(form.input_editors['base_input_path'].control.visible)
        self.assertFalse(form.input_editors['base_selected_mols'].control.visible)

    def test_each_section_has_origin_combo_search_include_buttons_and_table(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        for panel,origin,keys in ((settings.targets,scene.redocking,['2WSL_A','4M0E_A']),
                                  (settings.compounds,scene.vina,['CHEMBL100020','CHEMBL92663'])):
            self.assertEqual(panel.origin.value,'stage:'+origin['id'])
            self.assertIn('stage:'+origin['id'],[option.key for option in panel.origin.options])
            self.assertEqual([option.key for option in panel.combo.options],keys)
            self.assertFalse(panel.combo.disabled)
            self.assertEqual(panel.include_one.content,'Incluir selecionado')
            self.assertEqual(panel.include_list.content,'Incluir todos (2)')
            self.assertIsNotNone(panel.search.on_change)
            self.assertIsNotNone(panel.search_add.on_click)
            self.assertFalse(panel.table_box.visible)
            self.assertTrue(panel.empty.visible)
        target = next(o for o in settings.targets.combo.options if o.key == '4M0E_A')
        self.assertEqual(target.text,'4M0E · cadeia A · ligante 1YL · RMSD 0.2 Å')
        self.assertEqual([column.label.value for column in settings.targets.table.columns],
                         ['PDB','Cadeia','Ligante de referência','Resolução','RMSD','Ações'])
        self.assertEqual([column.label.value for column in settings.compounds.table.columns],
                         ['Composto','Target','Motores','Melhor score','Poses','Ações'])

    def test_combo_includes_one_item_or_the_whole_list(self):
        scene = self.publish()
        form = self.form(scene); settings = form.pose_settings
        self.no_repaint()
        targets = settings.targets
        targets.include_selected()
        self.assertTrue(targets.feedback.visible)
        self.assertEqual(targets.selected,[])
        targets.combo.value = '4M0E_A'; targets.include_selected()
        self.assertEqual(form.readers['selected_receptors'](),['4M0E_A'])
        self.assertEqual([option.key for option in targets.combo.options],['2WSL_A'])
        self.assertEqual([cell.content.value for cell in targets.table.rows[0].cells[:5]],
                         ['4M0E','A','1YL','2 Å','0.2 Å'])
        self.assertTrue(targets.table_box.visible)
        targets.include_all()
        self.assertEqual(form.readers['selected_receptors'](),['4M0E_A','2WSL_A'])
        self.assertTrue(targets.combo.disabled)
        self.assertEqual(targets.include_list.content,'Incluir todos')
        targets.remove('2WSL_A')
        self.assertEqual(form.readers['selected_receptors'](),['4M0E_A'])
        settings.compounds.include_all()
        self.assertEqual(form.readers['selected_ligands'](),['CHEMBL100020','CHEMBL92663'])
        settings.compounds.clear()
        self.assertIsNone(form.readers['selected_ligands']())
        self.assertFalse(settings.compounds.clear_button.visible)
        self.page.update.assert_not_called()

    def test_search_autocompletes_and_its_button_includes_the_match(self):
        scene = self.publish()
        form = self.form(scene); settings = form.pose_settings
        self.no_repaint()
        targets,compounds = settings.targets,settings.compounds
        targets.search.value = '4m'; targets.filter()
        self.assertTrue(targets.suggestion_box.visible)
        self.assertEqual([row.content.controls[1].value for row in targets.suggestions.controls],['4M0E'])
        self.assertEqual([option.key for option in targets.combo.options],['4M0E_A'])
        self.assertEqual(targets.include_list.content,'Incluir todos (1)')
        targets.include_search()
        self.assertEqual(form.readers['selected_receptors'](),['4M0E_A'])
        self.assertEqual(targets.search.value,'')
        self.assertFalse(targets.suggestion_box.visible)
        # Compounds follow the included targets.
        self.assertEqual([entry['receptors'] for entry in compounds.entries],[['4M0E_A'],['4M0E_A']])
        compounds.search.value = 'CHEMBL'; compounds.filter()
        self.assertEqual(len(compounds.suggestions.controls),2)
        compounds.include_search()
        self.assertIn('2 correspondências',compounds.feedback.value)
        self.assertIsNone(form.readers['selected_ligands']())
        compounds.search.value = 'XYZ'; compounds.filter(); compounds.include_search()
        self.assertIn('Nenhum composto disponível',compounds.feedback.value)
        compounds.search.value = '92663'; compounds.filter()
        self.assertEqual([row.content.controls[1].value for row in compounds.suggestions.controls],['CHEMBL92663'])
        compounds.suggestions.controls[0].on_click(None)
        self.assertEqual(form.readers['selected_ligands'](),['CHEMBL92663'])
        self.assertFalse(compounds.feedback.visible)
        self.assertIn('1 poses serão validadas',settings.plan.content.value)
        targets.remove('4M0E_A')
        self.assertIsNone(form.readers['selected_receptors']())
        self.page.update.assert_not_called()

    def test_long_lists_stay_bounded_in_the_combo_and_the_suggestions(self):
        from biomolexplorer.ui import pose_settings
        scene = self.publish()
        panel = self.form(scene).pose_settings.compounds
        panel.entries = [dict(compound_id=f'CHEMBL{index}',receptors=['4M0E_A'],engines=['vina'],poses=1,
                              best_score=-7.0,smiles='') for index in range(1500)]
        panel.draw_choices()
        self.assertEqual(len(panel.combo.options),pose_settings.COMBO_LIMIT)
        self.assertEqual(panel.include_list.content,'Incluir todos (1500)')
        self.assertIn('pesquise para ver os demais',panel.available.value)
        panel.search.value = 'CHEMBL1'; panel.filter(patch=False)
        self.assertEqual(len(panel.suggestions.controls),pose_settings.SUGGESTION_ROWS)
        self.assertIsInstance(panel.suggestions,ft.Column)
        self.assertIn('continue digitando',panel.suggestion_count.value)
        self.assertTrue(panel.suggestions.controls[0].content.controls[1].value.startswith('CHEMBL1'))
        expected = sum(1 for entry in panel.entries if 'chembl1' in entry['compound_id'].casefold())
        panel.include_all()
        self.assertEqual(len(panel.selected),expected)
        self.assertEqual(len(panel.table.rows),pose_settings.PAGE_SIZE)
        self.assertTrue(panel.pager.visible)

    def test_search_is_debounced_and_only_the_last_term_is_rendered(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        tasks = []

        def run_task(handler,*args):
            task = asyncio.ensure_future(handler(*args)); tasks.append(task); return task

        async def scenario():
            self.page.run_task = run_task
            settings.targets.search.value = '2'; settings.targets.changed()
            settings.targets.search.value = '4M0E'; settings.targets.changed()
            self.assertTrue(tasks[0].cancelled() or tasks[0].cancelling())
            await asyncio.gather(*tasks,return_exceptions=True)
            return [row.content.controls[1].value for row in settings.targets.suggestions.controls]
        self.assertEqual(asyncio.run(scenario()),['4M0E'])

    def test_a_compound_is_listed_once_for_every_included_target(self):
        scene = self.publish()
        form = self.form(scene); settings = form.pose_settings
        self.no_repaint()
        settings.targets.include_all()
        self.assertEqual(form.readers['selected_receptors'](),['2WSL_A','4M0E_A'])
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        # One row per ligand-receptor complex, the parameter still one compound.
        self.assertEqual(form.readers['selected_ligands'](),['CHEMBL92663'])
        self.assertEqual([[cell.content.value for cell in row.cells[:5]]
                          for row in settings.compounds.table.rows],
                         [['CHEMBL92663','2WSL_A','vina','-10.5','1'],
                          ['CHEMBL92663','4M0E_A','vina','-10.5','1']])
        self.assertIn('2 complexos ligante-receptor',settings.compounds.selected_count.value)
        # Each row opens the pose of its own complex.
        poses = [row.cells[-1].content.controls[0] for row in settings.compounds.table.rows]
        self.assertEqual([button.tooltip for button in poses],['Visualizar complexo receptor-ligante em 3D']*2)
        settings.targets.remove('2WSL_A')
        self.assertEqual([[cell.content.value for cell in row.cells[:2]]
                          for row in settings.compounds.table.rows],[['CHEMBL92663','4M0E_A']])

    def test_a_complex_without_poses_is_flagged_in_the_table(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        self.no_repaint()
        # 2WSL has no result for this compound; the row stays, flagged.
        settings.pair_index = {key:value for key,value in settings.pair_index.items()
                               if key != ('CHEMBL92663','2WSL_A')}
        settings.targets.include_all()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        rows = settings.compounds.table.rows
        self.assertEqual([[cell.content.value for cell in row.cells[:5]][1:] for row in rows],
                         [['2WSL_A','—','—','0'],['4M0E_A','vina','-10.5','1']])
        self.assertEqual(rows[0].color,'#FEF3C7')
        self.assertIsNone(rows[1].color)
        self.assertIn('1 sem pose no target',settings.compounds.selected_count.value)

    def test_each_engine_of_a_complex_opens_its_own_pose(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        scenes = []

        async def preview_complex(project_id,layers,name,ligand_smiles=None,reference_kind=None):
            scenes.append(dict(layers=layers,name=name))

        async def guard(action):
            await action()
        self.ui.preview_complex = preview_complex; self.ui.guard = guard
        self.no_repaint()
        # DOCK6 refines the geometry Vina found: another file for the complex.
        vina = scene.poses/'4M0E_1YL_604A_CHEMBL92663.lig.mol2'
        dock6 = scene.poses/'4M0E_1YL_604A_CHEMBL92663.dock6.mol2'
        dock6.write_bytes(vina.read_bytes())
        summary = settings.pair_index[('CHEMBL92663','4M0E_A')]
        summary['engines'].add('dock6')
        summary['by_engine']['dock6'] = dict(summary['by_engine']['vina'],pose=str(dock6),score=-42.0)
        for engine,name in (('vina',vina.name),('dock6',dock6.name)):
            chosen = settings.complex_scene('CHEMBL92663','4M0E_A',summary,engine)
            self.assertEqual(chosen['engine'],engine)
            self.assertTrue(chosen['layers'][1]['path'].endswith(name))
            self.assertEqual(chosen['name'],'CHEMBL92663 · 4M0E_A · '+engine.upper())
        # Without an engine the Vina pose comes first; an engine without a pose
        # of its own offers no complex at all.
        self.assertEqual(settings.complex_scene('CHEMBL92663','4M0E_A',summary)['engine'],'vina')
        self.assertIsNone(settings.complex_scene('CHEMBL92663','4M0E_A',summary,'amber'))
        settings.targets.combo.value = '4M0E_A'; settings.targets.include_selected()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        card = settings.compounds.table.rows[0].cells[-1].content.controls[1]
        asyncio.run(card.on_click(None))
        dialog = self.page.dialogs[-1]
        chips = [control for control in walk(dialog.content) if isinstance(control,ft.Container)
                 and control.tooltip == 'Visualizar complexo receptor-ligante em 3D']
        self.assertEqual(len(chips),2)
        for chip in chips:
            asyncio.run(chip.on_click(None))
        self.assertEqual([item['name'] for item in scenes],
                         ['CHEMBL92663 · 4M0E_A · VINA','CHEMBL92663 · 4M0E_A · DOCK6'])
        self.assertEqual([Path(item['layers'][1]['path']).name for item in scenes],
                         [vina.name,dock6.name])

    def test_the_card_never_presents_another_target_as_its_own(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings

        async def guard(action):
            await action()
        self.ui.guard = guard
        self.no_repaint()
        # Only 2WSL keeps a result for this compound, as the Consensus publishes
        # a single target per compound.
        settings.pair_index.pop(('CHEMBL92663','4M0E_A'))
        settings.targets.include_all()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        data = settings.card_data('CHEMBL92663','4M0E_A')
        # Both engines keep their row, with no values and no target of their own.
        self.assertEqual([(item['engine'],item['available'],item['receptor'])
                          for item in data['results']],
                         [('vina',False,'4M0E_A'),('dock6',False,'4M0E_A')])
        row = next(line for line in settings.compounds.table.rows
                   if line.cells[1].content.value == '4M0E_A')
        view,card,_ = row.cells[-1].content.controls
        self.assertIsNone(view.on_click)
        asyncio.run(card.on_click(None))
        texts = [getattr(control,'value',None) for control in walk(self.page.dialogs[-1].content)]
        self.assertIn('4M0E_A',texts)
        self.assertIn('Resultados de docking em 4M0E_A',texts)
        # The complex of another target is never shown as this card's own.
        self.assertNotIn('2WSL_A',texts)
        self.assertNotIn('Outros targets deste composto',texts)

    def test_every_number_of_the_card_keeps_three_decimals(self):
        import re
        scene = self.publish()
        settings = self.form(scene).pose_settings

        async def guard(action):
            await action()
        self.ui.guard = guard
        self.no_repaint()
        settings.targets.combo.value = '4M0E_A'; settings.targets.include_selected()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        asyncio.run(settings.compounds.table.rows[0].cells[-1].content.controls[1].on_click(None))
        values = [control.value for control in walk(self.page.dialogs[-1].content)
                  if isinstance(getattr(control,'value',None),str)]
        self.assertIn('165.192',values)                      # the molecular weight, not 165.19199999999995
        long = [(value,found) for value in values
                for found in re.findall(r'-?\d+\.\d{4,}',value)
                if 'c1ccc' not in value]                     # the SMILES carries no numbers
        self.assertEqual(long,[])

    def test_removing_a_row_drops_only_that_complex(self):
        scene = self.publish()
        form = self.form(scene); settings = form.pose_settings
        self.no_repaint()
        settings.targets.include_all()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        self.assertEqual([cell.content.value for row in settings.compounds.table.rows
                          for cell in row.cells[1:2]],['2WSL_A','4M0E_A'])
        settings.compounds.remove('CHEMBL92663','4M0E_A')
        # Only that complex leaves; the compound stays on its other target.
        self.assertEqual([cell.content.value for row in settings.compounds.table.rows
                          for cell in row.cells[1:2]],['2WSL_A'])
        self.assertEqual(form.readers['selected_ligands'](),['CHEMBL92663'])
        self.assertEqual(form.readers['excluded_pairs'](),['CHEMBL92663:4M0E_A'])
        # The execution skips exactly that complex.
        resolver = PoseInputs(self.store,self.pid,pipeline=self.ui.current['pipeline'])
        stage = settings.candidate_stage()
        self.assertEqual(sorted({record['receptor_id'] for record in resolver.poses(stage)}),['2WSL_A'])
        # Removing the last complex removes the compound and its exclusions.
        settings.compounds.remove('CHEMBL92663','2WSL_A')
        self.assertEqual(settings.compounds.table.rows,[])
        self.assertIsNone(form.readers['selected_ligands']())
        self.assertIsNone(form.readers['excluded_pairs']())

    def test_the_poses_column_counts_the_poses_of_each_complex(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        self.no_repaint()
        # DOCK6 published its own pose for one of the complexes.
        dock6 = scene.poses/'4M0E_1YL_604A_CHEMBL92663.dock6.mol2'
        dock6.write_bytes((scene.poses/'4M0E_1YL_604A_CHEMBL92663.lig.mol2').read_bytes())
        settings.records = list(settings.records)+[dict(receptor_id='4M0E_A',compound_id='CHEMBL92663',
            ligand_id='CHEMBL92663',engine='dock6',score='-42.0',pose=str(dock6),
            receptor=str(scene.receptors/'4M0E_A.noH.pdb'),canonical_smiles=SMILES)]
        settings.pair_index = settings.build_pairs(settings.records)
        settings.targets.include_all()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        rows = {row.cells[1].content.value:[cell.content.value for cell in row.cells[2:5]]
                for row in settings.compounds.table.rows}
        self.assertEqual(rows['4M0E_A'],['dock6/vina','-42','2'])
        self.assertEqual(rows['2WSL_A'],['vina','-10.5','1'])

    def test_the_docking_blocks_complete_a_connected_consensus(self):
        scene = self.publish()
        resolver = PoseInputs(self.store,self.pid,pipeline=self.ui.current['pipeline'])
        # A consensus origin keeps one target per compound, republishes a copy
        # of the pose and scores it with a relevance in [0, 1].
        copied = scene.poses/'0_vina_pose.mol2'
        copied.write_bytes((scene.poses/'2WSL_1YL_604A_CHEMBL92663.lig.mol2').read_bytes())
        aggregate = [dict(receptor_id='2WSL_A',compound_id='CHEMBL92663',ligand_id='CHEMBL92663',
                          engine='vina',score='0.86',pose=str(copied),
                          receptor=str(scene.receptors/'2WSL_A.noH.pdb'),canonical_smiles=SMILES)]
        with patch.object(PoseInputs,'records',lambda self,stage: list(aggregate)), \
             patch.object(PoseInputs,'docking_stages',lambda self,stage: [scene.vina['id']]):
            records = resolver.combined(scene.validation)
        complexes = {(r['receptor_id'],r['compound_id'],r['engine']):r for r in records}
        # The block that produced the pose answers for the complex it published.
        own = complexes[('2WSL_A','CHEMBL92663','vina')]
        self.assertEqual(own['score'],'-10.5')
        self.assertTrue(own['pose'].endswith('2WSL_1YL_604A_CHEMBL92663.lig.mol2'))
        # And the targets the aggregate dropped are back in the selection.
        self.assertIn(('4M0E_A','CHEMBL92663','vina'),complexes)

    def test_3d_is_disabled_for_an_engine_without_its_own_result(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        scenes = []

        async def preview_complex(project_id,layers,name,ligand_smiles=None,reference_kind=None):
            scenes.append(name)

        async def guard(action):
            await action()
        self.ui.preview_complex = preview_complex; self.ui.guard = guard
        self.no_repaint()
        summary = settings.pair_index[('CHEMBL92663','4M0E_A')]
        summary['engines'].add('dock6')
        # The pipeline ran DOCK6, but its pose is not among the published files.
        summary['by_engine']['dock6'] = dict(summary['by_engine']['vina'],
                                             pose=str(scene.poses/'missing.mol2'),score=-42.0)
        self.assertIsNone(settings.complex_scene('CHEMBL92663','4M0E_A',summary,'dock6'))
        self.assertIsNone(settings.view_complex('CHEMBL92663','4M0E_A',summary,'dock6'))
        results = settings.card_data('CHEMBL92663')['results']
        table = settings.result_table('CHEMBL92663',results,'4M0E_A')
        chips = {}
        for line in table.content.controls[1:]:
            cells = line.content.controls
            chips[(cells[0].content.value,cells[1].content.value)] = cells[-1].content
        self.assertEqual(sorted(chips),[('dock6','4M0E_A'),('vina','2WSL_A'),('vina','4M0E_A')])
        dock6 = chips[('dock6','4M0E_A')]; vina = chips[('vina','4M0E_A')]
        self.assertIsNone(dock6.on_click)
        self.assertEqual(dock6.tooltip,
                         'Este motor não tem resultado para o complexo nos blocos de docking do pipeline.')
        self.assertIsNotNone(vina.on_click)
        self.assertEqual(vina.tooltip,'Visualizar complexo receptor-ligante em 3D')
        asyncio.run(vina.on_click(None))
        self.assertEqual(scenes,['CHEMBL92663 · 4M0E_A · VINA'])

    def test_both_tables_follow_the_dialog_width_and_label_the_viewers(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        self.no_repaint()
        settings.targets.combo.value = '4M0E_A'; settings.targets.include_selected()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        for panel in (settings.targets,settings.compounds):
            rendered = panel.pair_table
            # Nothing scrolls sideways and nothing is wider than the container.
            self.assertIsNone(getattr(panel.table_box,'scroll',None))
            self.assertIsNone(panel.table_box.width)
            self.assertEqual(len(rendered.list.controls),1+len(panel.selected))
            for line in rendered.list.controls:
                cells = line.content.controls
                self.assertEqual(len(cells),len(rendered.columns))
                for cell in cells[:-1]:
                    self.assertIsNone(cell.width); self.assertTrue(cell.expand)
                self.assertEqual(cells[-1].width,rendered.actions_width)
            # Only the actions keep a fixed width, and it fits a narrow dialog.
            self.assertLessEqual(rendered.actions_width,150)
        view,_ = settings.targets.table.rows[0].cells[-1].content.controls
        self.assertEqual(view.content.value,'3D')
        compound_view,card,_ = settings.compounds.table.rows[0].cells[-1].content.controls
        self.assertEqual(compound_view.content.value,'3D')
        self.assertEqual(card.content.value,'2D')

    def test_table_actions_open_the_target_in_3d_and_the_compound_card_in_2d(self):
        scene = self.publish()
        settings = self.form(scene).pose_settings
        previews = []

        async def preview(project_id,path):
            previews.append(path)

        scenes = []

        async def preview_complex(project_id,layers,name,ligand_smiles=None,reference_kind=None):
            scenes.append(dict(layers=layers,name=name,smiles=ligand_smiles))

        async def guard(action):
            await action()
        self.ui.preview_artifact = preview; self.ui.preview_complex = preview_complex; self.ui.guard = guard
        settings.targets.combo.value = '4M0E_A'; settings.targets.include_selected()
        settings.compounds.combo.value = 'CHEMBL92663'; settings.compounds.include_selected()
        view,remove = settings.targets.table.rows[0].cells[-1].content.controls
        self.assertEqual(view.tooltip,'Visualizar target em 3D')
        self.assertEqual(remove.tooltip,'Remover target')
        asyncio.run(view.on_click(None))
        self.assertTrue(previews[-1].endswith('4M0E_A.dockprep.pdbqt'))
        compound_view,card,remove = settings.compounds.table.rows[0].cells[-1].content.controls
        self.assertEqual(compound_view.tooltip,'Visualizar complexo receptor-ligante em 3D')
        asyncio.run(compound_view.on_click(None))
        # The complex of the row, as the docking results open it.
        self.assertEqual([layer['role'] for layer in scenes[-1]['layers']],['receptor','pose'])
        self.assertTrue(scenes[-1]['layers'][0]['path'].endswith('receptors/4M0E_A.noH.pdb'))
        self.assertTrue(scenes[-1]['layers'][1]['path'].endswith('4M0E_1YL_604A_CHEMBL92663.lig.mol2'))
        self.assertEqual(scenes[-1]['name'],'CHEMBL92663 · 4M0E_A · VINA')
        self.assertEqual(scenes[-1]['smiles'],SMILES)
        self.assertEqual(card.tooltip,'Abrir card 2D do composto')
        self.assertEqual(remove.tooltip,'Remover este complexo')
        asyncio.run(card.on_click(None))
        dialog = self.page.dialogs[-1]
        self.assertEqual(dialog.title.value,'Detalhes do composto')
        texts = [getattr(control,'value',None) for control in walk(dialog.content)]
        self.assertIn('CHEMBL92663',texts)
        self.assertIn(SMILES,texts)
        # The card is about the complex of its row: the target is explicit and
        # the other targets live in their own section.
        self.assertIn('Resultados de docking em 4M0E_A',texts)
        self.assertIn('-10.5',texts)
        self.assertTrue(any(isinstance(control,ft.Image) for control in walk(dialog.content)))
        self.assertEqual([action.content for action in dialog.actions],['Abrir complexo em 3D','Fechar'])
        for section in ('Estrutura 2D','Propriedades físico-químicas','Identidade e demais atributos'):
            self.assertIn(section,texts)
        # Two rows of two panels: drawing beside properties, core beside results.
        first,second = dialog.content.content.controls[1:3]
        self.assertEqual([[cell.col['md'] for cell in row.controls] for row in (first,second)],[[5,7],[5,7]])
        for row,expected in ((first,('Estrutura 2D','Propriedades físico-químicas')),
                             (second,('Núcleo Bemis–Murcko','Resultados de docking em 4M0E_A'))):
            titles = [getattr(control,'value',None) for control in walk(row)]
            for title in expected:
                self.assertIn(title,titles)
        # The card lists every complex of the compound, each with its own 3D.
        chips = [control for control in walk(dialog.content) if isinstance(control,ft.Container)
                 and control.tooltip == 'Visualizar complexo receptor-ligante em 3D']
        # Vina has a result for this complex; DOCK6 is not in this pipeline.
        self.assertEqual(len(chips),1)
        asyncio.run(chips[0].on_click(None))
        self.assertEqual(scenes[-1]['name'],'CHEMBL92663 · 4M0E_A · VINA')
        greyed = [line for line in walk(dialog.content) if isinstance(line,ft.Row)
                  and getattr(getattr(line.controls[0],'content',None),'value',None) == 'dock6']
        self.assertEqual([getattr(cell.content,'value',None) for cell in greyed[0].controls[:4]],
                         ['dock6','—','—','—'])

    def test_changing_the_origin_rebinds_the_input_and_reloads_the_catalog(self):
        scene = self.publish()
        form = self.form(scene); settings = form.pose_settings
        targets = settings.targets
        targets.combo.value = '4M0E_A'; targets.include_selected()
        self.assertIn('stage:'+scene.vina['id'],[option.key for option in targets.origin.options])
        targets.origin.value = 'stage:'+scene.vina['id']
        asyncio.run(targets.origin_changed(None))
        self.assertEqual(sources(form.read()['bindings']['base_input_path'])[0]['stage'],scene.vina['id'])
        self.assertEqual([entry['receptor_id'] for entry in targets.entries],['2WSL_A','4M0E_A'])
        self.assertTrue(all(entry['file'].endswith('.noH.pdb') for entry in targets.entries))
        self.assertEqual(targets.selected,['4M0E_A'])
        targets.origin.value = 'stage:'+scene.redocking['id']
        asyncio.run(targets.origin_changed(None))
        self.assertEqual(sources(form.read()['bindings']['base_input_path'])[0]['stage'],scene.redocking['id'])

    def test_apply_persists_the_selection_and_the_run_needs_no_curation(self):
        scene = self.publish()
        form,dialog = self.open_dialog(scene)
        form.pose_settings.targets.combo.value = '4M0E_A'; form.pose_settings.targets.include_selected()
        form.pose_settings.compounds.search.value = 'CHEMBL92663'; form.pose_settings.compounds.include_search()
        form.field_controls['strict_filter'].value = False
        form.field_controls['vdw_clash_scale'].value = '0.7'
        asyncio.run(dialog.actions[-1].on_click(SimpleNamespace(control=dialog.actions[-1])))
        saved = self.store.project(self.token,self.pid)['pipeline'][-1]
        self.assertFalse(dialog.open)
        self.assertEqual(saved['operation'],'pose_validation')
        self.assertEqual(saved['parameters']['selected_receptors'],['4M0E_A'])
        self.assertEqual(saved['parameters']['selected_ligands'],['CHEMBL92663'])
        self.assertFalse(saved['parameters']['strict_filter'])
        self.assertEqual(saved['parameters']['vdw_clash_scale'],.7)
        self.assertEqual(sources(saved['bindings']['base_input_path'])[0]['stage'],scene.redocking['id'])
        self.assertEqual(sources(saved['bindings']['base_selected_mols'])[0]['stage'],scene.vina['id'])
        self.assertFalse(requires_curation(saved))
        # Reopening shows the saved selection in the tables.
        form,_ = self.open_dialog(scene)
        self.assertEqual(form.pose_settings.targets.selected,['4M0E_A'])
        self.assertEqual(len(form.pose_settings.compounds.table.rows),1)
        for field in ('base_input_path','base_selected_mols'):
            self.assertEqual(len(sources(saved['bindings'][field])),1)

    def test_cancel_discards_every_change_of_the_dialog(self):
        scene = self.publish()
        before = copy.deepcopy(self.store.project(self.token,self.pid)['pipeline'])
        current = copy.deepcopy(self.ui.current['pipeline'])
        self.ui.dirty = False
        form,dialog = self.open_dialog(scene)
        settings = form.pose_settings
        settings.targets.include_all(); settings.compounds.include_all()
        settings.targets.origin.value = 'stage:'+scene.vina['id']
        asyncio.run(settings.targets.origin_changed(None))
        form.field_controls['strict_filter'].value = False
        dialog.actions[0].on_click(SimpleNamespace(control=dialog.actions[0]))
        self.assertFalse(dialog.open)
        self.assertFalse(self.ui.dirty)
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],before)
        self.assertEqual(self.ui.current['pipeline'],current)
        form,_ = self.open_dialog(scene)
        self.assertEqual(form.pose_settings.targets.selected,[])
        self.assertEqual(form.pose_settings.targets.origin.value,'stage:'+scene.redocking['id'])

    def test_a_read_only_block_cannot_change_the_selection(self):
        scene = self.publish()
        settings = self.form(scene,writable=False).pose_settings
        for panel in (settings.targets,settings.compounds):
            for control in (panel.origin,panel.combo,panel.include_one,panel.include_list,panel.search,
                            panel.search_add,panel.clear_button):
                self.assertTrue(control.disabled)
            panel.include_all()
            self.assertEqual(panel.selected,[])


class ExecutionTests(PublishedResultsFixture):
    def resolve(self, scene, **parameters):
        stage = dict(scene.validation,parameters=dict(scene.validation['parameters'],**parameters))
        service = PipelineService(self.store)
        try:
            item = {}
            params = service._resolve(self.pid,self.store.user(self.token)['id'],stage,{},item=item,
                                      pipeline=[scene.redocking,scene.vina,stage])
            return params,item
        finally:
            service.close()

    def test_execution_reuses_the_published_pose_and_prepared_receptor(self):
        scene = self.publish()
        params,item = self.resolve(scene,selected_receptors=['4M0E_A'],selected_ligands=['CHEMBL92663'])
        self.assertEqual(len(params['pose_inputs']),1)
        entry = params['pose_inputs'][0]
        self.assertEqual(entry['pose'],str(scene.poses/'4M0E_1YL_604A_CHEMBL92663.lig.mol2'))
        self.assertEqual(entry['receptor'],str(scene.receptors/'4M0E_A.noH.pdb'))
        self.assertEqual((entry['receptor_id'],entry['pdb_code'],entry['compound_id']),
                         ('4M0E_A','4M0E','CHEMBL92663'))
        self.assertEqual(entry['canonical_smiles'],SMILES)
        self.assertEqual(item['input_files']['base_input_path'],[entry['receptor']])
        self.assertEqual(item['input_files']['base_selected_mols'],[entry['pose']])
        # No input was copied or regenerated for this execution.
        self.assertFalse((self.store.project_dir(self.pid)/'.input-cache').exists())

    def test_without_a_selection_every_published_pose_is_validated(self):
        scene = self.publish()
        params,_ = self.resolve(scene)
        self.assertEqual(len(params['pose_inputs']),4)
        self.assertEqual({entry['receptor_id'] for entry in params['pose_inputs']},{'2WSL_A','4M0E_A'})

    def test_a_selection_without_results_explains_what_to_review(self):
        scene = self.publish()
        with self.assertRaises(ValueError) as error:
            self.resolve(scene,selected_receptors=['9ZZZ_A'])
        self.assertIn('targets e compostos selecionados',str(error.exception))

    def test_a_prepared_target_replaces_a_missing_table_receptor(self):
        scene = self.publish()
        for path in scene.receptors.iterdir():
            path.unlink()
        params,_ = self.resolve(scene,selected_receptors=['4M0E_A'],selected_ligands=['CHEMBL92663'])
        entry = params['pose_inputs'][0]
        self.assertEqual(entry['receptor'],str(scene.prepared/'4M0E_A.dockprep.pdbqt'))
        # A prepared receptor holds no crystal ligand, so nothing is excluded.
        self.assertNotIn('exclude_residues',entry)

    def test_a_retrieved_complex_is_used_last_and_excludes_the_crystal_ligand(self):
        """Without any prepared receptor, the raw PDB still allows the clash test."""
        scene = self.publish()
        for folder in (scene.receptors,scene.prepared):
            for path in folder.iterdir():
                path.unlink()
        retrieval = new_stage('retrieve_structures')
        root = self.store.project_dir(self.pid)/'runs'/uuid4().hex/'retrieval/artifacts/Estruturas'
        root.mkdir(parents=True)
        raw = receptor_pdb(root/'4M0E.pdb',[[80.0,0.0,0.0]])
        metadata = root/'pdb_codes.csv'
        write_csv(metadata,[dict(PDB_CODE='4M0E',LIGAND='1YL',RESNUM=604,CHAIN='A',RESOLUTION='2.0')])
        stages = [dict(id=retrieval['id'],name=retrieval['name'],operation='retrieve_structures',
                       configuration=retrieval,status='succeeded',artifacts=[str(raw),str(metadata)])]
        run_id = uuid4().hex
        index_run(stages,run_id)
        with self.store.connect() as db:
            register_run(db,self.pid,run_id,stages)
        pipeline = [retrieval,*self.ui.current['pipeline']]
        validation = pipeline[-1]
        flow.connect(pipeline,retrieval['id'],validation['id'],'base_input_path')
        scene.validation['bindings'] = validation['bindings']
        params,_ = self.resolve(scene,selected_receptors=['4M0E'],selected_ligands=['CHEMBL92663'])
        entry = params['pose_inputs'][0]
        self.assertEqual(entry['receptor'],str(raw))
        self.assertEqual(entry['exclude_residues'],['1YL'])

    def test_validation_runs_end_to_end_on_the_published_files(self):
        import csv
        from wrappers.pose_validation import REPORT, validate_poses
        scene = self.publish()
        params,_ = self.resolve(scene,selected_receptors=['4M0E_A'],selected_ligands=['CHEMBL92663'])
        output = self.store.project_dir(self.pid)/'validation'
        summary = validate_poses(base_input_path=params['base_input_path'],base_output_path=str(output),
                                 target='Estruturas',base_selected_mols=params['base_selected_mols'],
                                 pose_inputs=params['pose_inputs'],strict_filter=True)
        self.assertEqual(summary['poses'],1)
        self.assertEqual(summary['valid_poses'],1)
        self.assertEqual(summary['published_poses'],1)
        with (output/REPORT).open(encoding='utf-8-sig',newline='') as stream:
            row = next(csv.DictReader(stream))
        self.assertEqual(row['pose_valid'],'True')
        self.assertEqual(row['receptor_id'],'4M0E_A')
        self.assertEqual(row['compound_id'],'CHEMBL92663')
        self.assertEqual(row['bond_orders'],'explicit')
        self.assertGreater(float(row['min_vdw_ratio']),.75)
        self.assertEqual({p.name for p in (output/'poses').iterdir()},
                         {'4M0E_1YL_604A_CHEMBL92663.lig.mol2'})
        self.assertEqual({p.name for p in (output/'receptors').iterdir()},{'4M0E_A.noH.pdb'})


if __name__ == '__main__':
    unittest.main()


class PoseAuditTests(unittest.TestCase):
    """Saving a pose validation block never inflates the docking execution JSON."""
    from test_consensus_history import ConsensusHistoryTests
    setUp = ConsensusHistoryTests.setUp
    tearDown = ConsensusHistoryTests.tearDown
    completed = ConsensusHistoryTests.completed

    def test_pose_only_save_uses_the_compact_rollback_snapshot(self):
        from biomolexplorer.project_state import read_snapshot, translate
        rid,item = self.completed(); original = item['configuration']
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original]),self.pid))
        project = self.store.project(self.owner,self.pid)
        validation = new_stage('pose_validation')
        saved = self.store.save_pipeline(self.owner,self.pid,[original,validation],project['revision'])
        event = self.store.history(self.owner,self.pid)[0]['id']
        raw = json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertTrue(raw['_compact_runs'])
        before = read_snapshot(self.store,self.pid,event,'before')
        after = read_snapshot(self.store,self.pid,event,'after')
        self.assertEqual(before['project']['pipeline'],[original])
        self.assertEqual(after['project']['pipeline'],[original,validation])
        self.assertEqual(after['runs'][0]['stages'],translate([item],str(self.folder),'@project'))
        updated = copy.deepcopy(validation)
        updated['parameters'].update(selected_receptors=['4M0E_A'],selected_ligands=['CHEMBL92663'])
        self.store.save_pipeline(self.owner,self.pid,[original,updated],saved['revision'])
        event = self.store.history(self.owner,self.pid)[0]['id']
        raw = json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertTrue(raw['_compact_runs'])
        self.assertEqual(self.store.project(self.owner,self.pid)['pipeline'][-1]['parameters']['selected_ligands'],
                         ['CHEMBL92663'])

    def test_a_change_to_another_block_keeps_the_complete_audit_format(self):
        rid,item = self.completed(); original = item['configuration']
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original]),self.pid))
        project = self.store.project(self.owner,self.pid)
        validation = new_stage('pose_validation')
        self.store.save_pipeline(self.owner,self.pid,[dict(original,name='Renamed'),validation],project['revision'])
        event = self.store.history(self.owner,self.pid)[0]['id']
        raw = json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertNotIn('_compact_runs',raw)


class SelectionFeedbackTests(PublishedResultsFixture):
    def test_unavailable_or_repeated_items_are_explained_instead_of_ignored(self):
        scene = self.publish(compounds=('CHEMBL92663',))
        self.ui.artifact_items = {}; self.ui.artifact_choices = {}
        settings = GuidedForm(self.ui,scene.validation,[],True).pose_settings
        settings.targets.include(['4M0E_A'])
        settings.targets.include(['4M0E_A'])
        self.assertEqual(settings.targets.feedback.value,'O item escolhido já está incluído.')
        settings.compounds.include(['CHEMBL1'])
        self.assertEqual(settings.compounds.feedback.value,
                         'CHEMBL1 não tem resultado nos targets incluídos nesta origem.')
        settings.targets.clear()
        settings.compounds.include(['CHEMBL1'])
        self.assertEqual(settings.compounds.feedback.value,'CHEMBL1 não está disponível na origem selecionada.')
        self.assertEqual(settings.compounds.selected,[])


class SubmissionTests(PublishedResultsFixture):
    """A second block of the same operation must not be confused with the first."""
    def save(self, pipeline):
        revision = self.store.project(self.token,self.pid)['revision']
        self.store.save_pipeline(self.token,self.pid,pipeline,revision)

    def test_an_unconfigured_twin_block_is_identified_instead_of_guessed(self):
        folder = self.store.project_dir(self.pid)/'Estruturas'
        (folder/'Vina').mkdir(parents=True)
        configured = new_stage('pose_validation')
        configured['parameters'].update(base_input_path=str(folder),base_selected_mols=str(folder/'Vina'),
                                        selected_receptors=['4M0E_A'],selected_ligands=['CHEMBL92663'])
        spare = new_stage('pose_validation')
        # Both blocks carry the default name, so only the identifier tells them apart.
        self.assertEqual(spare['name'],configured['name'])
        self.save([configured,spare])
        service = PipelineService(self.store)
        try:
            with self.assertRaises(StageAdjustment) as refused:
                service.submit(self.token,self.pid)
            self.assertEqual(refused.exception.stage_id,spare['id'])
            self.assertEqual(refused.exception.stage_name,spare['name'])
            self.assertIn('Conecte um bloco de origem',str(refused.exception))
            self.assertEqual(self.store.list_runs(self.token,self.pid),[])
            # The configured block was never the problem: it validates on its own.
            service._validate_parameters(self.pid,execution_stage(configured,'curated'),partial=True)
        finally:
            service.close()

    def test_the_refused_block_is_opened_instead_of_only_named(self):
        scene = self.publish(compounds=('CHEMBL92663',))
        spare = new_stage('pose_validation')
        self.ui.current['pipeline'] = [scene.redocking,scene.vina,scene.validation,spare]
        opened = []
        async def open_stage_dialog(stage_id,*arguments,**options):opened.append(stage_id)
        self.ui.open_stage_dialog = open_stage_dialog
        messages = []
        self.ui.notify = messages.append
        async def refuse():
            raise StageAdjustment('O bloco “Validate poses” precisa de ajustes: Conecte um bloco de origem',spare)
        asyncio.run(self.ui.guard(refuse))
        self.assertEqual(opened,[spare['id']])
        self.assertTrue(messages and 'precisa de ajustes' in messages[0])
        # An error without a block leaves the canvas alone.
        opened.clear()
        async def fail():raise ValueError('Sem bloco associado.')
        asyncio.run(self.ui.guard(fail))
        self.assertEqual(opened,[])


class ResultsViewTests(PublishedResultsFixture):
    """The validation report is the block's result; its files sit beneath it."""

    FIELDS = ('pose_id,receptor_id,ligand_id,compound_id,engine,score,pose_valid,valid_chemistry,'
              'valid_connectivity,valid_geometry,valid_planarity,valid_energy,valid_clash,min_vdw_ratio,'
              'min_vdw_overlap,min_heavy_distance,clashing_contacts,delta_energy,energy_method,'
              'max_bond_deviation,max_angle_deviation,max_planarity_deviation,planar_groups,atoms,'
              'bond_orders,resolved_by,rejection_reasons,rejection_details,error,pose_file,receptor_file,source_table')

    def validated(self):
        """One approved and one rejected pose, published like the wrapper does."""
        from biomolexplorer.artifact_registry import register_run
        molecule = relaxed_molecule(SMILES)
        stage = new_stage('pose_validation')
        run_id = uuid4().hex
        root = self.store.project_dir(self.pid)/'runs'/run_id/'artifacts'
        (root/'poses').mkdir(parents=True);(root/'receptors').mkdir(parents=True)
        receptor_pdb(root/'receptors/4M0E_A.noH.pdb',[[80.0,0.0,0.0],[83.5,0.0,0.0]])
        write_mol2(root/'poses/4M0E_1YL_604A_CHEMBL92663.lig.mol2',molecule)
        write_mol2(root/'poses/4M0E_1YL_604A_CHEMBL100020.lig.mol2',molecule)
        approved = ('4M0E_1YL_604A_CHEMBL92663,4M0E_A,CHEMBL92663,CHEMBL92663,vina,-10.5,True,True,True,True,'
                    'True,True,True,0.87,-0.39,2.64,6,26.07,MMFF94,0.047,3.1,0.02,3,24,assigned,binding,,,,'
                    +str(root/'poses/4M0E_1YL_604A_CHEMBL92663.lig.mol2')+','
                    +str(root/'receptors/4M0E_A.noH.pdb')+',')
        rejected = ('4M0E_1YL_604A_CHEMBL100020,4M0E_A,CHEMBL100020,CHEMBL100020,dock6,-55.6,False,True,True,'
                    'False,True,True,False,0.41,-1.82,1.21,19,180.4,MMFF94,0.51,41.2,0.03,3,26,assigned,binding,'
                    'bond_lengths;steric_clash,"Comprimentos de ligação fora da tolerância física.; '
                    'Choque estérico com o receptor ou cofatores.",,'
                    +str(root/'poses/4M0E_1YL_604A_CHEMBL100020.lig.mol2')+','
                    +str(root/'receptors/4M0E_A.noH.pdb')+',')
        (root/'pose_validation.csv').write_text('\n'.join([self.FIELDS,approved,rejected])+'\n',encoding='utf-8')
        (root/'pose_validation.json').write_text(json.dumps(dict(target='Estruturas',poses=2,valid_poses=1,
            invalid_poses=1,published_poses=1,strict_filter=True,criteria='posebusters_inspired_v1',
            rejection_counts={'bond_lengths':1,'steric_clash':1},unresolved_inputs=0,
            tolerances=dict(vdw_clash_scale=0.75,delta_energy_limit=100.0))),encoding='utf-8')
        artifacts = [str(root/'pose_validation.csv'),str(root/'pose_validation.json'),
                     str(root/'poses/4M0E_1YL_604A_CHEMBL92663.lig.mol2'),
                     str(root/'poses/4M0E_1YL_604A_CHEMBL100020.lig.mol2'),
                     str(root/'receptors/4M0E_A.noH.pdb')]
        stages = [dict(id=stage['id'],name=stage['name'],operation='pose_validation',configuration=stage,
                       status='succeeded',artifacts=artifacts)]
        index_run(stages,run_id)
        with self.store.connect() as db:
            register_run(db,self.pid,run_id,stages)
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                       (run_id,self.pid,self.store.user(self.token)['id'],'succeeded',
                        json.dumps(stages),time.time(),time.time(),None))
        self.ui.current['pipeline'] = [stage]
        return run_id,stages[0],root

    def test_the_report_is_read_with_its_counts_tolerances_and_published_poses(self):
        from biomolexplorer.pose_results import PoseValidationResults
        run_id,stage,root = self.validated()
        report = PoseValidationResults(self.store).report(self.token,self.pid,run_id,stage['id'])
        self.assertEqual(report['counts'],dict(total=2,valid=1,invalid=1))
        self.assertEqual(report['summary']['published_poses'],1)
        self.assertEqual(report['summary']['tolerances']['vdw_clash_scale'],0.75)
        # The pose that needs a decision is listed first.
        self.assertEqual([row['compound_id'] for row in report['rows']],['CHEMBL100020','CHEMBL92663'])
        rejected,approved = report['rows']
        self.assertEqual(rejected['failed'],['Geometria','Choque estérico'])
        self.assertEqual(rejected['rejection_reasons'],['bond_lengths','steric_clash'])
        self.assertAlmostEqual(rejected['min_vdw_ratio'],0.41)
        self.assertEqual(rejected['clashing_contacts'],19)
        self.assertTrue(approved['pose_valid'])
        self.assertEqual(approved['failed'],[])
        for row in report['rows']:
            self.assertEqual(Path(row['pose_path']).parent,root/'poses')
            self.assertEqual(Path(row['receptor_path']),root/'receptors/4M0E_A.noH.pdb')

    def test_a_report_outside_the_project_is_never_read_as_a_pose(self):
        from biomolexplorer.pose_results import PoseValidationResults
        run_id,stage,root = self.validated()
        service = PoseValidationResults(self.store)
        rows = service._row(self.pid,{},{'pose_file':'/etc/passwd','receptor_file':'../../secret.pdb',
                                         'pose_valid':'True','pose_id':'x'})
        self.assertIsNone(rows['pose_path'])
        self.assertIsNone(rows['receptor_path'])
        with self.assertRaises(Exception):
            service.report(self.token,self.pid,run_id,'00000000000000000000000000000000')

    def test_the_table_is_the_result_and_the_files_stay_collapsed_beneath_it(self):
        from biomolexplorer.ui.pose_results import PoseValidationTable
        from biomolexplorer.ui.stage_results import StageResults
        run_id,stage,root = self.validated()
        async def guard(action):await action()
        self.ui.guard = guard
        results = StageResults(self.ui,self.pid,run_id,stage,True)
        self.addCleanup(results.close)
        asyncio.run(results.load())
        self.assertIsInstance(results.child,PoseValidationTable)
        controls = list(walk(results.root))
        tiles = [c for c in controls if isinstance(c,ft.ExpansionTile)]
        self.assertEqual([t.title.value for t in tiles],['Arquivos deste bloco'])
        self.assertFalse(any(getattr(tile,'initially_expanded',False) for tile in tiles))
        table = results.child.table
        self.assertEqual([column.label.value for column in table.columns][:6],
                         ['Pose','Alvo','Composto','Motor','Score','Situação'])
        self.assertEqual(len(table.rows),2)
        self.assertIn('2 poses avaliadas · 1 aprovadas · 1 reprovadas',results.child.count.value)
        text = ' '.join(c.value for c in controls if isinstance(c,ft.Text) and isinstance(c.value,str))
        self.assertIn('Comprimentos de ligação fora da tolerância física',text)
        self.assertIn('Somente as poses aprovadas seguiram',text)
        self.assertIn('pose_validation.csv',[item['name'] for item in results.child.files.files])
        # Reprovadas first, and the filter keeps only what was asked for.
        results.child.choice.value='valid';asyncio.run(results.child.change_filter(None))
        self.assertEqual(len(results.child.table.rows),1)
