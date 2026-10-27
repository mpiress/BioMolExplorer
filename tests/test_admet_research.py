"""Scientific screening boundaries, auditable output and fitted EGG interaction."""
import asyncio
import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pandas as pd
from rdkit import Chem
from biomolexplorer.admet_config import resolve_config,selection_reasons,PROFILES
from biomolexplorer.catalog import new_stage
from biomolexplorer.operations import validate_operation,execute_operation
from biomolexplorer.egg_plot import bounds,egg_png
from biomolexplorer.visualizations import egg_view,load_view,SUFFIX
from caad.admet import MoleculeEvaluator
from wrappers.admet import ADMETWrapper


class ResearchAdmetTests(unittest.TestCase):
    def test_descriptors_are_shared_across_outputs_but_current_filters_are_reapplied(self):
        import os
        with patch.dict(os.environ,{'BIOMOL_PREPARATION_CACHE':str(self.root/'shared/ligand-preparation')}):
            first=ADMETWrapper(str(self.root/'first'),str(self.source),'compounds.csv',admet_profile='exploratory')
            with patch.object(first,'generate_plot'),patch.object(first.evaluator,'calculate_properties',wraps=first.evaluator.calculate_properties) as calculate:
                first.run_pipeline();self.assertEqual(calculate.call_count,3)
            second=ADMETWrapper(str(self.root/'second'),str(self.source),'compounds.csv',
                admet_profile='exploratory',admet_mode='filter',property_limits={'MW':[None,1]})
            with patch.object(second,'generate_plot'),patch.object(second.evaluator,'calculate_properties',wraps=second.evaluator.calculate_properties) as calculate:
                result=second.run_pipeline()
                self.assertEqual(calculate.call_count,0)
                self.assertTrue(result.empty)

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'source';self.source.mkdir()
        pd.DataFrame([dict(molecule_chembl_id='NITRO',canonical_smiles='O=[N+]([O-])c1ccccc1',source='Imported'),
            dict(molecule_chembl_id='ETHANOL',canonical_smiles='CCO',source='Reference'),
            dict(molecule_chembl_id='CHLORO',canonical_smiles='Clc1ccccc1',source='ChEMBL')]).to_csv(self.source/'compounds.csv',index=False)

    def run_admet(self,**options):
        wrapper=ADMETWrapper(str(self.root/'output'),str(self.source),'compounds.csv',**options)
        result=wrapper.run_pipeline();return wrapper,result

    def test_new_stages_are_descriptive_and_saved_legacy_configuration_stays_legacy(self):
        stage=new_stage('admet');self.assertEqual(stage['parameters']['admet_profile'],'exploratory')
        self.assertEqual(stage['parameters']['admet_mode'],'annotate')
        self.assertEqual(resolve_config({})['admet_profile'],'legacy')
        validate_operation('admet',dict(base_input_path=str(self.source),**stage['parameters']))

    def test_exploration_preserves_alerted_molecules_metadata_and_method_provenance(self):
        with patch.object(MoleculeEvaluator,'is_toxic',side_effect=AssertionError('Legacy exclusion must not run')):
            wrapper,result=self.run_admet(admet_profile='exploratory',admet_mode='annotate')
        self.assertEqual(len(result),3);self.assertEqual(result['source'].tolist(),['Imported','Reference','ChEMBL'])
        self.assertGreater(result.loc[0,'StructuralAlert_count'],0)
        self.assertEqual(result['PGP'].unique().tolist(),['not_evaluated'])
        report=json.loads((self.root/'output/compounds_admet_report.json').read_text())
        self.assertEqual(report['output'],3);self.assertIn('not validated',report['methods']['BBB'])
        self.assertEqual(report['provider']['domain_assessment'],'not_assessed')
        self.assertTrue((self.root/'output/compounds_egg.png').read_bytes().startswith(b'\x89PNG'))
        model=load_view((self.root/'output'/('compounds_egg'+SUFFIX)).read_bytes())
        self.assertEqual(len(model['nodes']),3)

    def test_alert_exclusion_keeps_complete_audit_and_annotate_mode_keeps_failed_candidates(self):
        wrapper,filtered=self.run_admet(admet_profile='custom',admet_mode='filter',alert_catalogs=['BRENK'],alert_policy='exclude')
        self.assertNotIn('NITRO',filtered['molecule_chembl_id'].tolist())
        evaluated=pd.read_csv(self.root/'output/compounds_alerts.csv')
        rejected=evaluated[evaluated['ADMET_selected']==False]
        self.assertIn('Alertas estruturais',rejected.iloc[0]['ADMET_reasons'])
        _,annotated=self.run_admet(admet_profile='custom',admet_mode='annotate',alert_catalogs=['BRENK'],alert_policy='exclude')
        self.assertEqual(len(annotated),3);self.assertFalse(annotated.loc[0,'ADMET_selected'])

    def test_cns_and_peripheral_profiles_partition_the_same_descriptors(self):
        _,cns=self.run_admet(admet_profile='cns',admet_mode='filter')
        _,peripheral=self.run_admet(admet_profile='peripheral',admet_mode='filter')
        self.assertTrue(len(cns)>0 and len(peripheral)>0)
        self.assertEqual(set(cns['molecule_chembl_id'])|set(peripheral['molecule_chembl_id']),{'NITRO','ETHANOL','CHLORO'})
        self.assertFalse(set(cns['molecule_chembl_id'])&set(peripheral['molecule_chembl_id']))

    def test_all_rejected_and_empty_inputs_preserve_schema(self):
        _,result=self.run_admet(admet_profile='custom',admet_mode='filter',property_limits={'MW':[900,None]})
        self.assertTrue(result.empty);self.assertEqual(len(pd.read_csv(self.root/'output/compounds_alerts.csv')),3)
        self.assertIn('ADMET_selected',result.columns)
        pd.DataFrame(columns=['molecule_chembl_id','canonical_smiles']).to_csv(self.source/'compounds.csv',index=False)
        _,empty=self.run_admet(admet_profile='exploratory');self.assertTrue(empty.empty)
        self.assertTrue((self.root/'output/compounds_egg.png').exists())

    def test_limits_are_inclusive_and_missing_model_values_never_pass_a_required_limit(self):
        props=dict(BBB='BBB+',HIA='HIA+',Lipinski_violations=1,Veber_violations=0,
            StructuralAlert_count=0,MW=500,QED=.5)
        config=resolve_config(dict(admet_profile='oral',property_limits={'MW':[None,500],'QED':[.5,None]}))
        self.assertEqual(selection_reasons(props,config),[])
        props['Veber_violations']=1;self.assertIn('Veber',selection_reasons(props,config)[0])
        config=resolve_config(dict(admet_profile='custom',prediction_backend='admet_ai',property_limits={'ML_hERG':[None,.2]}))
        self.assertEqual(selection_reasons(props,config),['ML_hERG: valor indisponível'])

    def test_invalid_research_options_fail_before_execution(self):
        for options in [dict(admet_profile='invalid'),dict(property_limits={'MW':[600,500]}),
            dict(property_limits={'QED':[float('nan'),None]}),dict(property_limits={'unknown':[0,1]}),
            dict(alert_catalogs=['toxic']),dict(max_rule_violations=True),dict(property_limits={'ML_hERG':[0,1]})]:
            with self.subTest(options=options),self.assertRaises(ValueError):resolve_config(options)

    def test_descriptor_reference_and_potts_guy_units(self):
        evaluator=MoleculeEvaluator();props=evaluator.calculate_properties('CCO')
        self.assertAlmostEqual(props['MW'],46.069,places=3);self.assertAlmostEqual(props['TPSA'],20.23,places=2)
        self.assertEqual((props['HBD'],props['HBA'],props['RB']),(1,1,0))
        screening=evaluator.screening_properties(props,[])
        self.assertAlmostEqual(screening['logKp'],-6.5820149,places=5)
        self.assertTrue(0<=props['QED']<=1);self.assertEqual(props['FractionCSP3'],1)
        self.assertEqual(screening['Lipinski_violations'],0)
        ionic=evaluator.calculate_properties('C[N+](C)(C)C.[Cl-]')
        warnings=evaluator.screening_properties(ionic,[])['ADMET_warnings']
        self.assertIn('Disconnected',warnings)

    def test_optional_model_enrichment_preserves_identity_and_cannot_overwrite_local_bbb(self):
        def predict(smiles):return {s:dict(ML_hERG=.8,ML_BBB_Martins=.3) for s in smiles},dict(provider='fixture',version='test')
        with patch('biomolexplorer.admet_predictors.predict_admet_ai',side_effect=predict):
            _,result=self.run_admet(admet_profile='custom',prediction_backend='admet_ai',admet_mode='filter',property_limits={'ML_hERG':[None,.2]})
        self.assertTrue(result.empty)
        evaluated=pd.read_csv(self.root/'output/compounds_alerts.csv')
        self.assertEqual(len(evaluated),3);self.assertEqual(evaluated['ML_hERG'].tolist(),[.8]*3)
        self.assertTrue(evaluated['BBB'].str.startswith('BBB').all())

    def test_operation_outputs_are_registered_and_primary_csv_respects_hypothesis(self):
        output=self.root/'worker'
        operation=execute_operation('admet',dict(base_input_path=str(self.source),input_file='compounds.csv',
            admet_profile='peripheral',admet_mode='filter'),output)
        names={Path(p).name for p in operation.artifacts}
        self.assertTrue({'compounds_BBB+.csv','compounds_BBB-.csv','compounds_HIA+.csv','compounds_HIA-.csv','compounds_alerts.csv','compounds_admet_report.json'}<=names)
        self.assertTrue((pd.read_csv(output/'compounds_BBB-.csv')['BBB']=='BBB-').all())

    def test_egg_resize_uses_full_viewport_and_preserves_hit_testing_and_outliers(self):
        from biomolexplorer.ui.results_viewer import ResultsViewer
        ui=SimpleNamespace(token='t',page=SimpleNamespace(update=lambda:None))
        model=egg_view(pd.DataFrame([dict(molecule_chembl_id='A',TPSA=20,WLOGP=1,BBB='BBB+'),
            dict(molecule_chembl_id='OUT',TPSA=420,WLOGP=-8,BBB='BBB-')]),'EGG')
        viewer=ResultsViewer(ui,'p',model,popup=True)
        viewer.resize_egg(SimpleNamespace(width=1200,height=560))
        self.assertEqual(viewer.WIDTH,1200);self.assertIsNone(viewer.scene.width);self.assertTrue(viewer.viewer.constrained)
        self.assertLess(viewer.project(200,0)[0],viewer.points['OUT'][0]);self.assertLess(viewer.points['OUT'][0],1160)
        self.assertGreater(viewer.points['OUT'][1],35);self.assertLess(viewer.points['OUT'][1],505)
        for code,point in viewer.points.items():self.assertEqual(viewer.index.nearest(*point),code)
        viewer.resize_egg(SimpleNamespace(width=400,height=560))
        self.assertTrue(all(65<x<360 for x,y in viewer.points.values()))
        self.assertTrue(egg_png(model).startswith(b'\x89PNG'))

    def test_guided_profiles_ranges_and_viewer_permissions(self):
        from biomolexplorer.ui.admet_settings import AdmetSettings
        form=SimpleNamespace(ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None)),writable=True,
            stage=new_stage('admet'),field_controls={},readers={})
        settings=AdmetSettings(form);settings.profile.value='oral';settings.apply_profile(None)
        settings.add_limit('QED',(.4,None));config=settings.read()
        self.assertEqual(config['rule_set'],'oral');self.assertEqual(config['admet_mode'],'filter')
        self.assertEqual(config['property_limits'],{'QED':[.4,None]})
        settings.add_limit('QED',(.6,None))
        with self.assertRaises(ValueError):settings.read()
        form.writable=False;readonly=AdmetSettings(form)
        self.assertTrue(readonly.profile.disabled);self.assertTrue(all(c.disabled for c in readonly.fields.values()))

    def test_all_supported_structural_catalogs_are_configurable_and_explained(self):
        import flet as ft
        from biomolexplorer.admet_config import ALERT_CATALOGS
        from biomolexplorer.ui.admet_settings import AdmetSettings
        form=SimpleNamespace(ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None)),writable=True,
            stage=new_stage('admet'),field_controls={},readers={})
        settings=AdmetSettings(form)
        self.assertEqual(set(settings.catalogs),{'PAINS','BRENK','NIH','ZINC','CHEMBL'})
        icons=[c for c in settings.control.content.controls[2].controls[0].controls
               for c in c.content.controls if isinstance(c,ft.IconButton)]
        self.assertEqual(len(icons),len(ALERT_CATALOGS));self.assertTrue(all(i.tooltip for i in icons))
        for control in settings.catalogs.values():control.value=True
        self.assertEqual(set(settings.read()['alert_catalogs']),set(ALERT_CATALOGS))
