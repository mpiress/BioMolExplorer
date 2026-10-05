"""Independent graph modes, compound provenance, MCC fragments and exploration."""
import asyncio
import base64
import csv
import io
import json
import shutil
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from uuid import uuid4

import pandas as pd
from rdkit import Chem

from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import compatible,stage_issues
from biomolexplorer.graph_inputs import GraphInputs
from biomolexplorer.operations import execute_operation,validate_operation
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.stage_cache import stage_key
from biomolexplorer.visualizations import SUFFIX,load_view
from biomolexplorer.workspace import AccessDenied,WorkspaceStore
from caad.graph_results import analyze_inputs,common_fragment,exact_edges

try:
    import flet as ft
    from biomolexplorer.ui.guided import GuidedForm
    from biomolexplorer.ui.stage_results import StageResults
    from biomolexplorer.ui.results_viewer import degree_color
except ImportError:
    GuidedForm=StageResults=None


class GraphExplorerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.compounds=self.root/'compounds.csv'
        self.compounds.write_text('molecule_chembl_id,canonical_smiles,source\nA,CCO,Own\nB,CCCO,Own\nC,c1ccccc1,Own\n')
        self.fp=self.root/'morgan_compounds.csv'
        self.fp.write_text('molecule_chembl_id,canonical_smiles,fingerprint,source\nA,CCO,"[1, 1, 0, 0]",Own\nB,CCCO,"[1, 1, 1, 0]",Own\nC,c1ccccc1,"[0, 0, 0, 1]",Own\n')
        self.ready=self.root/'arbitrary_filename.csv'
        self.ready.write_text('source,target,value\nA,B,0.6666666666666666\nB,A,0.6666666666666666\n')

    def tearDown(self):self.temp.cleanup()

    def entries(self):
        return [{'kind':'fingerprints','file':str(self.fp),'compound_files':[],'label':'Fingerprints propios'},
            {'kind':'similarity','file':str(self.ready),'compound_files':[str(self.compounds)],'label':'Similaridade pronta'}]

    def test_both_modes_generate_separate_graphs_and_preserve_isolates(self):
        output=self.root/'out';analyze_inputs(self.entries(),output,threshold=60)
        views=sorted(output.glob('plots/*'+SUFFIX));self.assertEqual(len(views),2)
        models=[load_view(p.read_bytes()) for p in views]
        for model in models:
            self.assertEqual({n['id'] for n in model['nodes']},{'A','B','C'})
            self.assertEqual(model['mcc'],['A','B']);self.assertEqual(len(model['edges']),1)
            self.assertEqual(model['fragment']['status'],'complete')
            query=Chem.MolFromSmarts(model['fragment']['smarts'])
            for smiles in ('CCO','CCCO'):self.assertTrue(Chem.MolFromSmiles(smiles).HasSubstructMatch(query))
            self.assertTrue(base64.b64decode(model['fragment']['image']).startswith(b'\x89PNG'))
            self.assertEqual(set(model['mcc_positions']),{'A','B'})
            image=output/'plots'/(model['analysis_id']+'.png');self.assertTrue(image.read_bytes().startswith(b'\x89PNG'))
            frame=pd.read_csv(output/'Molecules'/(model['analysis_id']+'.csv'))
            self.assertEqual(set(frame['molecule_chembl_id']),{'A','B'})
        self.assertEqual(models[0]['origin']['kind'],'fingerprints');self.assertEqual(models[1]['origin']['kind'],'similarity')
        self.assertAlmostEqual(models[0]['edges'][0]['value'],models[1]['edges'][0]['value'])

    def test_selected_metric_threshold_only_changes_fingerprint_experiment(self):
        tanimoto,_=exact_edges(self.fp,'Tanimoto',75);dice,_=exact_edges(self.fp,'Dice',75)
        self.assertEqual(len(tanimoto),0);self.assertEqual(len(dice),1)
        mcconnaughey,_=exact_edges(self.fp,'McConnaughey',50)
        self.assertEqual(len(mcconnaughey),1)
        self.assertGreaterEqual(mcconnaughey.iloc[0]['value'],.5)
        output=self.root/'ready';analyze_inputs(self.entries()[1:],output,threshold=100,metric='Dice')
        model=load_view(next(output.glob('plots/*'+SUFFIX)).read_bytes())
        self.assertEqual(len(model['edges']),1);self.assertIsNone(model['origin']['threshold'])

    def test_missing_smiles_unknown_codes_and_conflicting_fingerprints_are_excluded(self):
        self.fp.write_text('molecule_chembl_id,fingerprint\nA,"[1, 0]"\n')
        analyze_inputs(self.entries()[:1],self.root/'out')
        report=json.loads((self.root/'out/molecule_exclusions.json').read_text())
        self.assertIn('sem molécula',report['records'][0]['reason'])
        self.fp.write_text('name,fingerprint\nA,"[1, 0]"\nA,"[0, 1]"\n')
        edges,ids=exact_edges(self.fp,'Tanimoto',0)
        self.assertEqual(ids,[]);self.assertTrue(edges.empty)
        self.ready.write_text('source,target,value\nA,OUTSIDE,0.9\n')
        analyze_inputs(self.entries()[1:],self.root/'unknown')
        model=load_view(next((self.root/'unknown').glob('plots/*'+SUFFIX)).read_bytes())
        self.assertEqual(model['edges'],[])
        self.assertTrue((self.root/'unknown/molecule_exclusions.json').exists())

    def test_aliases_and_leading_zero_identifiers_are_preserved(self):
        self.fp.write_text('name,smiles,fingerprint\n001,CCO,"[1, 0]"\n002,CCCO,"[1, 1]"\n')
        analyze_inputs(self.entries()[:1],self.root/'alias',threshold=40)
        model=load_view(next((self.root/'alias').glob('plots/*'+SUFFIX)).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'001','002'})

    def test_common_fragment_single_empty_and_timeout_are_labeled_truthfully(self):
        self.assertEqual(common_fragment([])['status'],'empty')
        self.assertEqual(common_fragment(['CCO'])['atoms'],3)
        self.assertEqual(common_fragment([None])['status'],'missing_structures')
        result=SimpleNamespace(smartsString='[#6]-[#6]',canceled=True)
        with patch('caad.graph_results.rdFMCS.FindMCS',return_value=result) as find:
            fragment=common_fragment(['CCO','CCCO'],timeout=1)
        self.assertEqual(fragment['status'],'partial');self.assertEqual(find.call_args.kwargs['threshold'],1.0)
        self.assertEqual(find.call_args.kwargs['timeout'],1)
        with patch('caad.graph_results.rdFMCS.FindMCS',return_value=SimpleNamespace(smartsString='',canceled=False)):
            self.assertEqual(common_fragment(['C','O'])['status'],'no_common_fragment')

    def test_independent_inputs_may_reuse_codes_without_losing_structures_in_union(self):
        other=self.root/'other.csv';other.write_text('molecule_chembl_id,canonical_smiles,fingerprint\nA,CCC,"[1, 0]"\nB,CCCC,"[1, 0]"\n')
        entries=[self.entries()[0],{'kind':'fingerprints','file':str(other),'compound_files':[]}]
        analyze_inputs(entries,self.root/'conflicts',threshold=60)
        models=[load_view(p.read_bytes()) for p in (self.root/'conflicts').glob('plots/*'+SUFFIX)]
        self.assertTrue(all(set(m['mcc'])=={'A','B'} for m in models))
        union=pd.read_csv(self.root/'conflicts/Molecules/molecules.csv')
        self.assertEqual(len(union),4);self.assertEqual(union['molecule_chembl_id'].nunique(),4)
        self.assertEqual(set(union['original_code']),{'A','B'})

    def test_operation_rejects_direct_fingerprints_and_accepts_ready_similarity(self):
        source=self.root/'fp';source.mkdir();shutil.copy(self.fp,source/self.fp.name)
        with self.assertRaises(ValueError):
            execute_operation('graphs',{'fingerprints_path':str(source),'threshold':60},self.root/'operation')
        result=execute_operation('graphs',{'graph_inputs':self.entries()[1:]},self.root/'operation')
        self.assertTrue(any(p.endswith(SUFFIX) for p in result.artifacts))
        from wrappers.molecular_analyzer import analyze_graphs
        legacy=self.root/'datasets/input';legacy.mkdir(parents=True)
        (legacy/'Similarity').mkdir();shutil.copy(self.compounds,legacy/'compounds.csv')
        shutil.copy(self.ready,legacy/'Similarity'/self.ready.name)
        with patch.dict('os.environ',{'BIOMOL_WORKSPACE':str(self.root)}):
            analyze_graphs(base_input_path='/datasets/input',base_output_path='/datasets/output')
        self.assertTrue(list((self.root/'datasets/output/plots').glob('*'+SUFFIX)))
        for params in ({'mcs_timeout':0},{'mcs_ring_matches_ring_only':'True'},{'threshold':101},{'graph_inputs':[{}]}):
            with self.subTest(params=params),self.assertRaises(ValueError):validate_operation('graphs',params)


class GraphProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@graphs.org','owner-password')
        self.pid=self.store.create_project(self.token,'Graphs')['id'];self.root=self.store.project_dir(self.pid)
        self.retrieval=new_stage('retrieve_compounds');self.fp=new_stage('fingerprints');self.sim=new_stage('similarity');self.graph=new_stage('graphs')
        self.fp['bindings']['base_input_path']={'stage':self.retrieval['id'],'selector':'CHEMBL220_MOLS.csv'}
        self.sim['bindings']['base_input_path']={'stage':self.fp['id'],'selector':'morgan_CHEMBL220_MOLS.csv'}
        self.sim2=new_stage('similarity');self.sim2['name']='Similaridade alternativa'
        self.sim2['bindings']['base_input_path']={'stage':self.fp['id'],'selector':'auto'}
        self.graph['bindings']={'similarity_path':{'sources':[{'stage':self.sim['id'],'selector':'auto'},{'stage':self.sim2['id'],'selector':'auto'}]}}
        self.pipeline=[self.retrieval,self.fp,self.sim,self.sim2,self.graph]
        self.store.save_pipeline(self.token,self.pid,self.pipeline,0)
        self.c=self.root/'CHEMBL220_MOLS.csv';self.c.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\nB,CCCO\nC,c1ccccc1\n')
        self.f=self.root/'morgan_CHEMBL220_MOLS.csv';self.f.write_text('molecule_chembl_id,fingerprint\nA,"[1, 0]"\nB,"[1, 1]"\nC,"[0, 0]"\n')
        self.s=self.root/'Tanimoto_morgan_CHEMBL220_MOLS.csv';self.s.write_text('source,target,value\nA,B,0.5\n')
        self.s2=self.root/'Dice_morgan_CHEMBL220_MOLS.csv';self.s2.write_text('source,target,value\nA,B,0.8\n')
        self.results={self.retrieval['id']:[str(self.c)],self.fp['id']:[str(self.f)],self.sim['id']:[str(self.s)],self.sim2['id']:[str(self.s2)]}

    def tearDown(self):self.temp.cleanup()

    def resolve(self):return GraphInputs(self.store,self.pid,self.pipeline,self.results).resolve(self.graph)

    def test_both_generated_sources_infer_smiles_through_producer_chain(self):
        params=self.resolve();self.assertNotIn('base_input_path',params)
        self.assertEqual(len(params['graph_inputs']),2)
        self.assertTrue(all(e['compound_files']==[str(self.c)] for e in params['graph_inputs']))
        self.assertFalse(compatible(self.fp,self.graph,'fingerprints_path'));self.assertTrue(compatible(self.sim,self.graph,'similarity_path'))
        self.assertEqual(stage_issues(self.graph),[])
        service=PipelineService(self.store)
        try:
            resolved=service._resolve(self.pid,self.store.user(self.token)['id'],self.graph,self.results)
            self.assertEqual(resolved,params)
        finally:service.close()

    def test_multiple_similarity_files_remain_independent(self):
        second=self.root/'maccs_second.csv';second.write_text('molecule_chembl_id,canonical_smiles,fingerprint\nOWN,CCC,"[1, 1, 0]"\n')
        edges=self.root/'Dice_maccs_second.csv';edges.write_text('source,target,value\nOWN,OWN,1\n')
        self.results[self.fp['id']].append(str(second))
        self.results[self.sim2['id']].append(str(edges))
        params=self.resolve();self.assertEqual(len(params['graph_inputs']),3)
        separate=next(e for e in params['graph_inputs'] if e['file']==str(edges))
        self.assertEqual(separate['compound_files'],[str(second)])
        self.assertEqual(separate['fingerprint'],'maccs')
        original=next(e for e in params['graph_inputs'] if e['file']==str(self.s))
        self.assertEqual(original['compound_files'],[str(self.c)])

    def test_merged_similarity_uses_fingerprint_metadata_instead_of_incomplete_retrieval(self):
        merged=self.root/'morgan_selected_compounds.csv'
        merged.write_text('molecule_chembl_id,canonical_smiles,fingerprint\nA,CCO,"[1,0]"\nALIAS,CCCO,"[1,1]"\n')
        ready=self.root/'Tanimoto_selected_fingerprints.csv'
        ready.write_text('source,target,value\nA,ALIAS,0.8\n')
        self.sim['bindings']['base_input_path']={'stage':self.fp['id'],'selector':'auto'}
        self.results[self.fp['id']]=[str(merged)]
        self.results[self.sim['id']]=[str(ready)]
        self.graph['bindings']['similarity_path']={'stage':self.sim['id'],'selector':'auto'}
        params=self.resolve()
        self.assertEqual(params['graph_inputs'][0]['compound_files'],[str(merged)])
        result=execute_operation('graphs',params,self.root/'recovered')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'A','ALIAS'})
        self.assertEqual(result.details['excluded_records'],0)

    def test_cache_ignores_run_paths_but_tracks_inferred_compound_contents(self):
        params=self.resolve();original=stage_key(self.graph,params,{},'software')
        destination=self.root/'another_run';destination.mkdir()
        rewritten=json.loads(json.dumps(params))
        for entry in rewritten['graph_inputs']:
            for filename in [entry['file']]+entry['compound_files']:
                shutil.copy(filename,destination/Path(filename).name)
            entry['file']=str(destination/Path(entry['file']).name)
            entry['compound_files']=[str(destination/Path(f).name) for f in entry['compound_files']]
        self.assertEqual(original,stage_key(self.graph,rewritten,{},'software'))
        self.c.write_text(self.c.read_text().replace('CCO','CCN'))
        self.assertNotEqual(original,stage_key(self.graph,params,{},'software'))

    def test_raw_paths_and_assets_cannot_escape_project(self):
        other=self.store.create_project(self.token,'Other')['id'];file=self.store.project_dir(other)/'compound.csv';file.write_text('molecule_chembl_id,canonical_smiles\nX,CCO\n')
        self.graph['parameters']['graph_inputs']=[{'kind':'similarity','file':str(file),'compound_files':[]}];self.graph['bindings']={}
        with self.assertRaises(AccessDenied):self.resolve()
        self.graph['parameters']={};self.graph['bindings']={'similarity_path':{'asset':uuid4().hex}}
        with self.assertRaises(AccessDenied):self.resolve()

    def test_relative_paths_resolve_inside_project(self):
        folder=self.root/'similarity';folder.mkdir();shutil.copy(self.s,folder/self.s.name)
        self.graph['bindings']={};self.graph['parameters'].update(similarity_path='similarity',base_input_path='.')
        params=self.resolve();self.assertEqual(Path(params['graph_inputs'][0]['file']).parent,folder)
        similarity=self.root/'Similarity';similarity.mkdir();shutil.copy(self.s,similarity/self.s.name)
        self.graph['parameters']['similarity_path']='Similarity'
        legacy=self.resolve();self.assertEqual(legacy['graph_inputs'][0]['file'],str(similarity/self.s.name))

    def test_running_snapshot_keeps_compound_provenance_after_project_edits(self):
        snapshot=json.loads(json.dumps(self.pipeline))
        self.store.save_pipeline(self.token,self.pid,[],1)
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.pid,self.store.user(self.token)['id'],self.graph,self.results,pipeline=snapshot)
            self.assertTrue(all(e['compound_files']==[str(self.c)] for e in params['graph_inputs']))
        finally:service.close()

    def test_real_pipeline_with_uploads_both_modes_and_incremental_admet(self):
        def upload(path,kind):
            ticket=self.store.prepare_upload(self.token,self.pid,path.name,kind)
            (self.store.staging/ticket).write_bytes(path.read_bytes())
            return self.store.finish_upload(self.token,ticket)
        compounds=upload(self.c,'compounds');similarity=upload(self.s,'similarity');second=upload(self.s2,'similarity')
        stage=new_stage('graphs')
        stage['bindings']={'similarity_path':{'sources':[{'asset':similarity},{'asset':second}]},'base_input_path':{'asset':compounds}}
        self.store.save_pipeline(self.token,self.pid,[stage],1)
        def wait(run):
            deadline=time.monotonic()+35
            while time.monotonic()<deadline:
                run=self.store.get_run(self.token,run['id'])
                if run['status'] not in ('queued','running'):return run
                time.sleep(.03)
            self.fail('Scientific worker did not finish in 35 seconds')
        service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
        try:
            run=wait(service.submit(self.token,self.pid))
            self.assertEqual(run['status'],'succeeded',run['error'])
            choices=ResultFiles(self.store).graph_datasets(self.token,self.pid,run['id'],stage['id'])
            self.assertEqual(len(choices),2)
            artifacts=run['stages'][0]['artifacts']
            selected=next(p for p in artifacts if '/Molecules/001_' in p)
            next_stage=new_stage('admet');next_stage['bindings']['base_input_path']={'stage':stage['id'],'selector':'Molecules/'+Path(selected).name}
            self.store.save_pipeline(self.token,self.pid,[stage,next_stage],2)
            subsequent=wait(service.submit(self.token,self.pid))
            self.assertEqual(subsequent['status'],'awaiting_input')
            subsequent=wait(service.resume(self.token,subsequent['id'],next_stage))
            self.assertEqual(subsequent['status'],'succeeded',subsequent['error'])
            self.assertTrue(subsequent['stages'][0]['reused'])
            self.assertEqual(subsequent['stages'][0]['artifacts'],artifacts)
            self.assertTrue(any(p.endswith('_egg'+SUFFIX) for p in subsequent['stages'][1]['artifacts']))
        finally:service.close()

    def fixture_run(self):
        from caad.graph_results import analyze_inputs
        params=self.resolve();output=self.root/'runs/artifacts'
        analyze_inputs(params['graph_inputs'],output,threshold=40)
        artifacts=[str(p) for p in output.rglob('*') if p.is_file()]
        self.rid=uuid4().hex
        stage=dict(id=self.graph['id'],name=self.graph['name'],operation='graphs',status='succeeded',artifacts=artifacts)
        with self.store.connect() as db:db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
            (self.rid,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps([stage]),time.time(),time.time(),None))
        return stage

    def ui(self):
        dialogs=[];saved=[]
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):await action()
        async def save_file(**kwargs):saved.append(kwargs)
        return SimpleNamespace(token=self.token,current=self.store.project(self.token,self.pid),store=self.store,
            page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append,width=1440),call=call,guard=guard,
            picker=SimpleNamespace(save_file=save_file),dialogs=dialogs,saved=saved,notify=lambda s:None,
            service=SimpleNamespace(dock6_path=None))

    def test_authorized_graph_choices_exports_and_viewer_access(self):
        self.fixture_run();service=ResultFiles(self.store)
        choices=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id']);self.assertEqual(len(choices),2)
        guest=self.store.register('Viewer','reader@graphs.org','reader-password')
        with self.assertRaises(AccessDenied):service.graph_model(guest,self.pid,self.rid,self.graph['id'],choices[0]['path'])
        self.store.invite(self.token,self.pid,'reader@graphs.org','viewer');self.store.accept_invitation(guest,self.pid)
        args=(guest,self.pid,self.rid,self.graph['id'],choices[0]['path'])
        model=service.graph_model(*args);rows=list(csv.DictReader(io.StringIO(service.graph_mcc_csv(*args).decode())))
        self.assertEqual({r['molecule_chembl_id'] for r in rows},set(model['mcc']))
        self.assertTrue(service.graph_png(*args).startswith(b'\x89PNG'))
        with self.assertRaises(AccessDenied):service.graph_model(guest,self.pid,self.rid,self.graph['id'],str(self.c))
        second=Path(choices[1]['path']);model=json.loads(second.read_text());model['title']=choices[0]['name'];second.write_text(json.dumps(model))
        unique=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id'])
        self.assertEqual(len({c['name'] for c in unique}),2)

    @unittest.skipIf(GuidedForm is None,'Install ui extra')
    def test_guided_graph_mode_hides_irrelevant_calculation_parameters(self):
        ui=self.ui();form=GuidedForm(ui,self.graph,[],True);form.layout()
        self.assertNotIn('graph_inputs',form.readers)
        self.assertNotIn('metric',form.field_controls);self.assertNotIn('fingerprint',form.field_controls)
        self.assertNotIn('threshold',form.field_controls);self.assertNotIn('fingerprints_path',form.input_editors)
        self.assertEqual(len(form.input_editors['similarity_path'].rows),2)
        options=form.input_editors['similarity_path'].source_options()
        self.assertEqual({o.key for o in options if o.key.startswith('stage:')},{'stage:'+self.sim['id'],'stage:'+self.sim2['id']})
        self.assertFalse(any(o.key.startswith('stage:') for o in form.input_editors['base_input_path'].source_options()))
        self.assertTrue(form.field_controls['mcs_timeout'].visible)

    @unittest.skipIf(StageResults is None,'Install ui extra')
    def test_inline_graph_selector_mcc_fragment_color_controls_click_and_download(self):
        stage=self.fixture_run();ui=self.ui();results=StageResults(ui,self.pid,self.rid,stage,True)
        asyncio.run(results.load());self.assertTrue(results.loaded);self.assertEqual(len(results.dataset.options),2)
        self.assertTrue(results.graph_files.files)
        self.assertTrue(all(Path(f['name']).stem==Path(results.dataset.value).name.removesuffix(SUFFIX) for f in results.graph_files.files))
        self.assertTrue(all(f['description'] for f in results.graph_files.files))
        viewer=results.child;self.assertEqual(set(viewer.points),{'A','B','C'})
        self.assertFalse(viewer.fragment_panel.visible);self.assertEqual(len(viewer.color_choice.options),2)
        async def reset():return None
        viewer.viewer.reset=reset
        asyncio.run(viewer.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc'))))
        self.assertEqual(set(viewer.points),{'A','B'});self.assertTrue(viewer.fragment_panel.visible)
        self.assertTrue(any(isinstance(c,ft.Image) for c in viewer.fragment_panel.controls))
        x,y=viewer.points['A'];viewer.hover(SimpleNamespace(local_position=SimpleNamespace(x=x,y=y)))
        self.assertEqual(viewer.tip_text.value,'A');asyncio.run(viewer.select('A'))
        self.assertIn('molécula 2D',ui.dialogs[-1].title.value)
        asyncio.run(viewer.actions[1].on_click(None));self.assertTrue(ui.saved[-1]['file_name'].endswith('_mcc.csv'))
        results.dataset.value=results.dataset.options[1].key
        asyncio.run(results.change_graph(None));self.assertFalse(viewer.active)
        self.assertTrue(all(Path(f['name']).stem==Path(results.dataset.value).name.removesuffix(SUFFIX) for f in results.graph_files.files))
        self.assertNotEqual(degree_color(0,4),degree_color(4,4));self.assertEqual(degree_color(4,4),'#FDE725')
        results.close();self.assertFalse(results.child.active);self.assertFalse(results.graph_files.active)
