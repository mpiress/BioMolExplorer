"""Readable compound cards, scoped MCS restoration and safe molecular details."""
import asyncio
import unittest

import flet as ft
import test_mcs_ranking_refinements as fixtures
from biomolexplorer.ui.molecular_cards import property_sections
from biomolexplorer.ui.localization import LocalizedPage


def walk(control):
    yield control
    for key in ('content','controls','title','actions'):
        child=getattr(control,key,None)
        for item in child if isinstance(child,list) else [child]:
            if item is not None and not isinstance(item,(str,int,float,bool)):
                yield from walk(item)


class MolecularInspectorPresentationTests(unittest.IsolatedAsyncioTestCase):
    def viewer(self):
        view=fixtures.MCSConfigurationTests().viewer()
        view.nodes['0']['properties'].update(pref_name='Example compound',MW=60.095,TPSA=20.23,BBB='BBB+',HIA='HIA+',
            standard_inchi='InChI=1S/C3H8O',custom_measurement='value from source')
        return view

    async def test_mcs_is_after_selected_compound_and_individual_scaffold(self):
        view=self.viewer();inspector=view.explorer.inspector
        content=view.explorer.inspector_card.content.controls[-1].controls
        self.assertLess(content.index(inspector.summary),content.index(inspector.mcs_tile))
        self.assertLess(content.index(inspector.scaffold),content.index(inspector.mcs_tile))
        await view.select('0',show_popup=False)
        text=[c.value for c in walk(inspector.summary) if isinstance(c,ft.Text)]
        for value in ('0','Example compound','BBB+','HIA+','Massa molecular','Conexões no universo'):
            self.assertIn(value,text)

    async def test_details_group_properties_and_retain_identifiers_values_and_provenance(self):
        view=self.viewer();await view.select('0',show_popup=False)
        await view.explorer.inspector.details()
        dialog=view.explorer.inspector.dialog
        labels=[c.title.value for c in walk(dialog.content) if isinstance(c,ft.ExpansionTile)]
        self.assertIn('Identificadores químicos',labels)
        self.assertIn('Propriedades físico-químicas',labels)
        self.assertIn('Predições ADMET',labels)
        self.assertIn('Outros atributos',labels)
        fields={c.label:c.value for c in walk(dialog.content) if isinstance(c,ft.TextField)}
        self.assertEqual(fields['SMILES'],'CCCO')
        self.assertEqual(fields['InChI'],'InChI=1S/C3H8O')
        texts=[c.value for c in walk(dialog.content) if isinstance(c,ft.Text)]
        self.assertIn('value from source',texts)
        self.assertIn('Origem dos dados',texts)
        self.assertIn('Posição na rede',texts)
        self.assertIn('Compostos relacionados',texts)
        self.assertIs(dialog.content.content.scroll,ft.ScrollMode.AUTO)

    async def test_detached_detail_request_cannot_mix_compounds_during_authorization(self):
        view=self.viewer();inspector=view.explorer.inspector;await view.select('0',show_popup=False)
        original=view.ui.call;entered,release=asyncio.Event(),asyncio.Event()
        async def call(function,*args,**kwargs):
            if function is view.ui.store.project:
                entered.set();await release.wait()
            return await original(function,*args,**kwargs)
        view.ui.call=call
        pending=asyncio.create_task(inspector.details());await entered.wait()
        view.selected='1';view.selection_version+=1
        release.set();await pending
        self.assertIsNone(inspector.dialog)
        view.ui.page.show_dialog.assert_not_called()

    async def test_narrow_details_are_bounded_and_sections_stack(self):
        view=self.viewer();view.ui.page.width=360;view.ui.page.height=640
        await view.select('0',show_popup=False);await view.explorer.inspector.details()
        dialog=view.explorer.inspector.dialog
        self.assertLessEqual(dialog.content.width,360)
        self.assertLess(dialog.content.height,640)
        columns=[c for c in walk(dialog.content) if isinstance(c,ft.Column) and isinstance(c.col,dict) and 'lg' in c.col]
        self.assertTrue(columns)
        self.assertTrue(all(c.col['xs']==12 for c in columns))

    async def test_missing_structure_retains_other_information(self):
        sections=property_sections({'canonical_smiles':None,'BBB':'BBB-','note':'preserved'}, {})
        texts=[c.value for s in sections for c in walk(s) if isinstance(c,ft.Text)]
        self.assertIn('BBB-',texts);self.assertIn('preserved',texts)

    async def test_open_details_resize_with_the_window(self):
        view=self.viewer();await view.select('0',show_popup=False)
        await view.explorer.inspector.details()
        view.ui.page.width=1440;view.ui.page.height=1100
        await view.explorer.resize_handler(None)
        content=view.explorer.inspector.dialog.content
        self.assertEqual((content.width,content.height),(1120,820))
        view.ui.page.width=360;view.ui.page.height=640
        await view.explorer.resize_handler(None)
        self.assertEqual((content.width,content.height),(296,480))

    async def test_english_details_translate_sections_without_translating_molecular_data(self):
        view=self.viewer();view.ui.page=LocalizedPage(view.ui.page,'en')
        await view.select('0',show_popup=False);await view.explorer.inspector.details()
        text=[c.value for c in walk(view.explorer.inspector.dialog) if isinstance(c,ft.Text)]
        self.assertIn('Chemical identifiers',text)
        self.assertIn('Physicochemical properties',text)
        self.assertIn('ADMET predictions',text)
        self.assertIn('Molecular details · 2D',text)
        self.assertIn('Example compound',text)

    async def test_scaffold_depiction_can_be_rendered_without_background_or_border(self):
        from biomolexplorer.ui.molecular_cards import image_panel
        image=image_panel(b'test image',170,framed=False)
        self.assertIsInstance(image,ft.Image)
        self.assertEqual(image.height,170)
