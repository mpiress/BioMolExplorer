"""Provider-independent molecular investigation survives every graph paint mode."""
import asyncio
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

import networkx as nx
import pandas as pd
import flet as ft

from biomolexplorer.compound_properties import compound_properties
from biomolexplorer.ui.results_viewer import ResultsViewer
from biomolexplorer.visualizations import graph_view,load_view,SUFFIX,molecule_sdf
from biomolexplorer.operations import execute_operation
from crawlers.pubchem import PubChemSimilarMols
from wrappers.crawlers import _save_compound_dataset


class CompoundInspectionTests(unittest.TestCase):
    def test_chembl_handoff_has_only_disjoint_mol_sim_and_keeps_annotations(self):
        from kernel.utilities import fileHandling
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);mols=root/'mols';sims=root/'sims';mols.mkdir();sims.mkdir()
            def row(code,smiles):
                return dict(molecule_chembl_id=code,pref_name='Name '+code,
                    molecule_structures={'canonical_smiles':smiles,'standard_inchi_key':'key-'+code},
                    molecule_properties={'full_mwt':'46'})
            pd.DataFrame([row('M1','CCO')]).to_csv(mols/'provider.csv',index=False)
            pd.DataFrame([row('M1','CCO'),row('DUP','CCO'),row('S1','CCC')]).to_csv(sims/'provider.csv',index=False)
            output=root/'output'
            fileHandling(output_path=str(output)).prepare_datamols('target',str(mols),str(sims))
            self.assertEqual({p.name for p in output.glob('*.csv')},{'target_MOLS.csv','target_SIMS.csv'})
            first=pd.read_csv(output/'target_MOLS.csv');second=pd.read_csv(output/'target_SIMS.csv')
            self.assertEqual(first['molecule_chembl_id'].tolist(),['M1'])
            self.assertEqual(second['molecule_chembl_id'].tolist(),['S1'])
            self.assertEqual(first['pref_name'].tolist(),['Name M1'])
            self.assertEqual(first['standard_inchi_key'].tolist(),['key-M1'])

    def viewer(self,graph):
        graph=nx.relabel_nodes(graph,str)
        rows=pd.DataFrame([dict(molecule_chembl_id=code,canonical_smiles='CCO',
            source='ChEMBL' if i%2==0 else 'PubChem') for i,code in enumerate(graph)])
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):return await action()
        self.page=SimpleNamespace(update=Mock(),show_dialog=Mock(),width=1440,height=1000)
        self.page.show_dialog.side_effect=lambda dialog:setattr(dialog,'open',True)
        ui=SimpleNamespace(token='token',current={'id':'project'},page=self.page,
            store=SimpleNamespace(project=Mock()),call=call,guard=guard,notify=Mock())
        view=ResultsViewer(ui,'project',graph_view(graph,rows,'Molecular network'));view.build()
        async def noop(*args,**kwargs):pass
        view.viewer.reset=noop;view.viewer.pan=noop;view.viewer.zoom=noop
        return view

    def test_node_click_and_details_have_same_structures_and_properties_for_both_providers(self):
        view=self.viewer(nx.path_graph(2))
        for identifier,source in [('0','ChEMBL'),('1','PubChem')]:
            x,y=view.points[identifier]
            asyncio.run(view.click(SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))))
            inspector=view.explorer.inspector
            self.assertEqual(view.selected,identifier)
            self.assertTrue(inspector.structure['image'].startswith(b'\x89PNG'))
            self.assertEqual(inspector.properties['source'],source)
            for key in ('MW','TPSA','WLOGP','HBD','HBA','RotatableBonds','HeavyAtomCount','FormalCharge'):
                self.assertIsNotNone(inspector.properties.get(key),key)
                self.assertIn(key,inspector.properties['property_sources'])
            asyncio.run(inspector.details())
            self.assertTrue(inspector.dialog.open is not False)
            self.assertTrue(inspector.dialog.content)
            self.assertIn('V2000',molecule_sdf(inspector.properties['canonical_smiles']))
            inspector.close_details()
        self.assertEqual(view.explorer.molecular.neighbor.value,'0')
        asyncio.run(view.follow_relation('0'));self.assertEqual(view.selected,'0')
        asyncio.run(view.explorer.molecular.analyze())
        fields={option.key for option in view.explorer.molecular.property.options}
        self.assertTrue({'MW','TPSA','WLOGP','HBD','HBA','RotatableBonds'}<=fields)

    def test_dense_graph_paints_every_compound_and_retains_click_inspection_and_svg(self):
        view=self.viewer(nx.complete_graph(160))
        self.assertEqual(view.mode,'full');self.assertFalse(view.explorer.aggregate_mode)
        self.assertTrue(view.raster.visible);self.assertTrue(view.raster.src.startswith(b'\x89PNG'))
        self.assertEqual(len(view.points),160)
        self.assertEqual(len(view.displayed()),160)
        self.assertEqual(len(view.data['edges']),12720)
        x,y=view.points['0']
        asyncio.run(view.click(SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))))
        self.assertEqual(view.selected,'0')
        self.assertTrue(view.explorer.inspector.structure['image'])
        self.assertTrue(view.drawing.shapes) # Exact selection/neighbor overlay.
        self.assertIn(b'<circle',view.explorer.svg())
        asyncio.run(view.follow_relation('1'));self.assertEqual(view.selected,'1')

    def test_community_click_opens_its_molecules_using_the_same_inspector(self):
        view=self.viewer(nx.path_graph(8))
        view.explorer.aggregate_choice.value=True;asyncio.run(view.explorer.change())
        mapping=view.analysis.community_cache['full'][0]
        community=next(iter(mapping.values()))
        x,y=next((x,y) for cell in view.index.cells.values() for code,x,y in cell if code==community)
        asyncio.run(view.click(SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))))
        self.assertFalse(view.explorer.aggregate_mode)
        self.assertEqual(set(view.points),{code for code,group in mapping.items() if group==community})
        identifier=next(iter(view.points));x,y=view.points[identifier]
        asyncio.run(view.click(SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))))
        self.assertEqual(view.selected,identifier)
        self.assertTrue(view.explorer.inspector.structure['image'])

    def test_selection_in_expanded_view_opens_molecular_details(self):
        view=self.viewer(nx.path_graph(2));view.explorer.focus_mode=True
        asyncio.run(view.select('0'))
        self.assertIsNotNone(view.explorer.inspector.dialog)
        self.page.show_dialog.assert_called()


class CompoundExtractionTests(unittest.TestCase):
    def test_pubchem_fetches_complementary_properties_in_the_same_cached_batch(self):
        with tempfile.TemporaryDirectory() as temporary:
            crawler=PubChemSimilarMols(temporary)
            def request(endpoint,data,params=None):
                if endpoint.startswith('smiles/cids'):return {'IdentifierList':{'CID':[1]}}
                if endpoint.startswith('fastsimilarity'):return {'IdentifierList':{'CID':[3]}}
                self.assertEqual(data,{'cid':'3'})
                for field in ('InChI','Title','XLogP','TPSA','HBondDonorCount','HBondAcceptorCount','RotatableBondCount'):
                    self.assertIn(field,endpoint)
                return {'PropertyTable':{'Properties':[dict(CID=3,SMILES='CCC',MolecularWeight='44.1',
                    Title='Propane',IUPACName='propane',InChI='InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',
                    InChIKey='ATUOYWHBWRKTHZ-UHFFFAOYSA-N',XLogP=2.3,TPSA=0,
                    HBondDonorCount=0,HBondAcceptorCount=0,RotatableBondCount=0,HeavyAtomCount=3,Charge=0)]}}
            crawler._request=Mock(side_effect=request)
            data=crawler.search(pd.DataFrame([dict(molecule_chembl_id='CHEMBL1',canonical_smiles='CCO')]))
            crawler.session.close()
            row=data.iloc[0].to_dict()
            self.assertEqual(row['standard_inchi'],row['InChI'])
            self.assertEqual(row['pref_name'],'Propane');self.assertEqual(row['TPSA'],0)
            self.assertEqual(row['XLogP'],2.3)
            self.assertNotEqual(row['WLOGP'],row['XLogP'])
            self.assertIn('RDKit',row['property_sources']['WLOGP'])
            self.assertIn('PubChem',row['property_sources']['TPSA'])
            self.assertEqual(crawler._request.call_count,3)
            stored=pd.read_csv(Path(temporary)/'compounds.csv')
            self.assertIn('standard_inchi',stored.columns)

    def test_chembl_export_keeps_annotations_structures_and_identifier_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            original=pd.DataFrame([dict(molecule_chembl_id='CHEMBL1',canonical_smiles='CCO',
                pref_name='Ethanol',max_phase=4,molecule_properties={'full_mwt':46.07,'alogp':-0.3},
                molecule_structures={'canonical_smiles':'CCO','standard_inchi':'reference-inchi',
                                     'standard_inchi_key':'reference-key'})])
            result=_save_compound_dataset(original,None,Path(temporary),'source')
            row=result.iloc[0]
            self.assertEqual(row['pref_name'],'Ethanol');self.assertEqual(row['max_phase'],4)
            self.assertEqual(row['standard_inchi'],'reference-inchi')
            self.assertEqual(row['molecule_properties']['alogp'],-0.3)
            self.assertIn('RDKit',row['property_sources']['WLOGP'])

    def test_graph_exports_keep_common_and_original_properties_for_each_database(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);compounds=root/'compounds.csv';edges=root/'edges.csv'
            pd.DataFrame([dict(molecule_chembl_id='CHEMBL1',canonical_smiles='CCO',source='ChEMBL',pref_name='ethanol'),
                          dict(molecule_chembl_id='PUBCHEM3',canonical_smiles='CCCO',source='PubChem',InChIKey='pubchem-key')]).to_csv(compounds,index=False)
            edges.write_text('source,target,value\nCHEMBL1,PUBCHEM3,0.8\n')
            result=execute_operation('graphs',{'graph_inputs':[dict(kind='similarity',file=str(edges),compound_files=[str(compounds)])]},root/'out')
            model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
            for node in model['nodes']:
                for field in ('MW','TPSA','WLOGP','HBD','HBA','RotatableBonds'):
                    self.assertIsNotNone(node['properties'].get(field))
            pubchem=next(node['properties'] for node in model['nodes'] if node['id']=='PUBCHEM3')
            self.assertEqual(pubchem['standard_inchi_key'],'pubchem-key')

    def test_logp_methods_are_never_mislabeled_and_missing_data_is_not_invented(self):
        row=compound_properties(dict(molecule_chembl_id='PUBCHEM3',XLogP=4.5,canonical_smiles=None),True)
        self.assertEqual(row['XLogP'],4.5)
        self.assertNotIn('WLOGP',row)
        self.assertNotIn('TPSA',row)
        self.assertNotIn('standard_inchi',row)

    def test_identifiers_names_and_descriptors_survive_real_fingerprints_similarity_and_graphs(self):
        from biomolexplorer.compound_properties import compound_frame
        with tempfile.TemporaryDirectory() as temporary,patch.dict(os.environ,{'BIOMOL_CPU_WORKERS':'1'}):
            root=Path(temporary);inputs=root/'inputs';inputs.mkdir()
            rows=[dict(molecule_chembl_id='CHEMBL1',canonical_smiles='CCO',source='ChEMBL',
                       pref_name='Ethanol',standard_inchi_key='chembl-key'),
                  dict(molecule_chembl_id='PUBCHEM3',canonical_smiles='CCCO',source='PubChem',
                       Title='Propanol',InChIKey='pubchem-key',XLogP=0.8)]
            compound_frame(pd.DataFrame(rows)).to_csv(inputs/'compounds.csv',index=False)
            fp=execute_operation('fingerprints',{'base_input_path':str(inputs),'morgan':False,
                'maccs':True,'pharmacophore':False},root/'fp')
            fp_path=Path(next(p for p in fp.artifacts if Path(p).name=='maccs_compounds.csv'))
            sim=execute_operation('similarity',{'base_input_path':str(fp_path.parent),'filename':fp_path.name,
                'fingerprint':'maccs','threshold':1,'approximate':False},root/'sim')
            sim_path=next(p for p in sim.artifacts if p.endswith('.csv'))
            graph=execute_operation('graphs',{'graph_inputs':[dict(kind='similarity',file=sim_path,
                compound_files=[str(fp_path)])]},root/'graph')
            model=load_view(Path(next(p for p in graph.artifacts if p.endswith(SUFFIX))).read_bytes())
            nodes={node['id']:node['properties'] for node in model['nodes']}
            self.assertEqual(nodes['CHEMBL1']['pref_name'],'Ethanol')
            self.assertEqual(nodes['PUBCHEM3']['pref_name'],'Propanol')
            self.assertEqual(nodes['PUBCHEM3']['standard_inchi_key'],'pubchem-key')
            self.assertIn('RDKit',nodes['PUBCHEM3']['property_sources']['WLOGP'])
            self.assertEqual(nodes['PUBCHEM3']['XLogP'],0.8)
