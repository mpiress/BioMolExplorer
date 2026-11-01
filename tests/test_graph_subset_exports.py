"""Graph subsets retain viewer numbering and can be selected as docking inputs."""
import copy
from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd

from biomolexplorer.artifact_choices import choices
from biomolexplorer.bindings import sources
from biomolexplorer.catalog import new_stage
from biomolexplorer.graph_exports import subset_label
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_origins import ResultOrigins
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.visualizations import SUFFIX,load_view
from caad.graph_results import analyze_inputs
import test_pipeline_execution as fixtures


class GraphSubsetExportTests(unittest.TestCase):
    setUp=fixtures.PipelineExecutionTests.setUp
    tearDown=fixtures.PipelineExecutionTests.tearDown

    def generate(self, second=False):
        root=self.store.project_dir(self.project_id)
        compounds=root/'input.csv'
        pd.DataFrame(dict(molecule_chembl_id=['001','B','C','D','E','F'],
            canonical_smiles=['CCO','CCCO','CCCCO','CCN','CCCN','c1ccccc1'],
            activity=[1,2,3,4,5,6])).to_csv(compounds,index=False)
        edges=root/'Tanimoto_morgan_target.csv'
        edges.write_text('source,target,value\n001,B,0.9\nB,C,0.9\nD,E,0.8\n')
        entries=[dict(kind='similarity',file=str(edges),compound_files=[str(compounds)],label='Target · Morgan')]
        if second:
            other=root/'Dice_maccs_other.csv';other.write_text('source,target,value\nD,E,0.9\nE,F,0.9\n')
            entries.append(dict(kind='similarity',file=str(other),compound_files=[str(compounds)],label='Other · MACCS'))
        output=root/'graphs'
        # Figures are unrelated to file membership; existing graph tests render them.
        with patch('caad.graph_results.report_png',return_value=b'figure'):
            analyze_inputs(entries,output,mcs_timeout=1)
        models=[load_view(path.read_bytes()) for path in sorted(output.glob('plots/*'+SUFFIX))]
        return output,models

    def test_full_graph_is_the_single_indexed_component_source(self):
        output,models=self.generate(second=True)
        for model in models:
            membership={n['id']:n['properties']['component'] for n in model['nodes']}
            exports=model['compound_exports']
            self.assertEqual(len(exports),1)
            for item in exports:
                frame=pd.read_csv(output/item['file'],dtype={'molecule_chembl_id':str})
                expected=set(membership)
                self.assertEqual(set(frame['molecule_chembl_id']),expected)
                self.assertEqual(set(frame['graph_scope']),{item['scope']})
                self.assertEqual(set(frame['graph_analysis_id']),{model['analysis_id']})
                self.assertEqual(frame['graph_component'].tolist(),[membership[n] for n in frame['molecule_chembl_id']])
                self.assertTrue(frame['canonical_smiles'].notna().all())
                self.assertIn('activity',frame)
                self.assertEqual(len(frame),item['nodes'])
            self.assertEqual(len(model['component_index']),model['statistics']['components'])
            self.assertEqual({r['component'] for r in model['component_index']},set(membership.values()))
            self.assertFalse(list((output/'Molecules').glob(model['analysis_id']+'__component_*.csv')))
            legacy=pd.read_csv(output/'Molecules'/(model['analysis_id']+'.csv'),dtype={'molecule_chembl_id':str})
            self.assertEqual(set(legacy['molecule_chembl_id']),set(model['mcc']))
        self.assertEqual(set(pd.read_csv(output/'Molecules/molecules.csv',dtype={'molecule_chembl_id':str})['molecule_chembl_id']),{'001','B','C','D','E','F'})

    def test_docking_picker_identifies_scopes_and_materializes_only_selected_subset(self):
        output,models=self.generate()
        model=models[0];producer=new_stage('graphs');vina=new_stage('docking_vina')
        producer['name']='Graph strategy · Target'
        artifacts=[str(p) for p in output.rglob('*') if p.is_file()]
        item=dict(producer,status='succeeded',configuration=copy.deepcopy(producer),artifacts=artifacts)
        available=choices(item)
        self.assertFalse(any('maxcomp' in name or 'data/similarity' in name for name in available))
        vina['bindings']['base_selected_mols']={'stage':producer['id'],'selector':'auto'}
        run=dict(id='run',stages=[item])
        service=PipelineService(self.store)
        self.addCleanup(service.close)
        for exported in model['compound_exports']:
            pending=dict(name=vina['name'],configuration=copy.deepcopy(vina))
            form=FileSelection(run,pending)
            picker=form.pickers['base_selected_mols']
            target=Path(exported['file']).name
            index=next(i for i,(_,ref) in enumerate(picker.rows) if Path(ref['selector']).name==target)
            label=picker.label(picker.rows[index][0])
            from biomolexplorer.artifact_naming import parse_filename
            self.assertEqual(parse_filename(label)['result_type'],'Graph-Full')
            self.assertIn('Graph strategy',picker.rows[index][0].origin_label)
            picker.set_selected([index],True)
            selected=form.read()
            revision=self.store.project(self.token,self.project_id)['revision']
            self.store.save_pipeline(self.token,self.project_id,[producer,selected],revision)
            selected=self.store.project(self.token,self.project_id)['pipeline'][1]
            refs=sources(selected['bindings']['base_selected_mols'])
            self.assertEqual(len(refs),1)
            self.assertIn('artifact_id',refs[0])
            folder,name,paths=service._materialize_inputs(self.project_id,selected,'base_selected_mols',refs,{producer['id']:artifacts})
            frame=pd.read_csv(folder/(name+'.csv'),dtype={'molecule_chembl_id':str})
            expected=pd.read_csv(output/exported['file'],dtype={'molecule_chembl_id':str})
            self.assertEqual(set(frame['molecule_chembl_id']),set(expected['molecule_chembl_id']))
            self.assertEqual(paths,[output/exported['file']])
        # A saved legacy MCC reference remains materializable, including its old path.
        legacy=output/'Molecules'/(model['analysis_id']+'.csv')
        folder,name,_=service._materialize_inputs(self.project_id,vina,'base_selected_mols',
            [{'stage':producer['id'],'selector':'Molecules/'+legacy.name}],{producer['id']:artifacts})
        self.assertEqual(set(pd.read_csv(folder/(name+'.csv'),dtype={'molecule_chembl_id':str})['molecule_chembl_id']),set(model['mcc']))

    def test_result_labels_preserve_origin_and_expose_component_number(self):
        producer=new_stage('graphs')
        path='/graphs/Molecules/001_similarity_morgan__component_0002.csv'
        producer.update(status='succeeded',artifacts=[path])
        details=ResultOrigins([producer],'r').describe(path,producer['id'],'r')
        from biomolexplorer.artifact_naming import parse_filename
        self.assertEqual(parse_filename(details['display_name'])['result_type'],'Molecules-C2')
        self.assertIn(producer['name'],details['origin'])
        self.assertIn('União dos MCCs',subset_label('/graphs/Molecules/molecules.csv'))
        from biomolexplorer.ui.localization import Translator
        self.assertIn('Component 2',subset_label(path,Translator('en')))
        self.assertIn('Full graph',subset_label('/x/001__full_graph.csv',Translator('en')))


if __name__=='__main__':unittest.main()
