"""Every block inherits the settings layout without changing scientific requests."""
import copy
from types import SimpleNamespace
import unittest
import flet as ft
from biomolexplorer.catalog import TITLES,new_stage
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.settings_layout import dialog_size,fit_dialog
from test_frontend_navigation import walk


class SharedSettingsTests(unittest.TestCase):
    def test_all_forms_preserve_requests_after_layout_and_disable_read_only_controls(self):
        stages=[new_stage(op) for op in TITLES]
        ui=SimpleNamespace(current={'pipeline':stages},service=SimpleNamespace(dock6_path=None),
            page=SimpleNamespace(update=lambda:None))
        for stage in stages:
            if stage['operation']=='redocking':stage['parameters']['pdb_codes']=[['1LGU','AAH',163,'A']]
            for writable in (True,False):
                with self.subTest(operation=stage['operation'],writable=writable):
                    original=copy.deepcopy(stage)
                    form=GuidedForm(ui,stage,[],writable)
                    before=form.read()
                    layout=form.layout()
                    self.assertEqual(form.read(),before)
                    self.assertEqual(stage,original)
                    self.assertEqual(layout.scroll,ft.ScrollMode.AUTO)
                    if not writable:
                        for control in form.field_controls.values():
                            if isinstance(control,(ft.TextField,ft.Dropdown,ft.Switch,ft.Checkbox)):
                                self.assertTrue(control.disabled,control.label)

    def test_vina_ligand_preparation_uses_wrapped_responsive_cells_without_offsets(self):
        stage=new_stage('docking_vina')
        ui=SimpleNamespace(current={'pipeline':[stage]},service=SimpleNamespace(dock6_path=None),
            page=SimpleNamespace(update=lambda:None))
        form=GuidedForm(ui,stage,[],True);form.layout()
        group=form.preparation_settings.groups['ligand']
        controls=list(walk(group));rows=[c for c in controls if isinstance(c,ft.ResponsiveRow)]
        self.assertTrue(rows)
        for row in rows:
            for cell in row.controls:
                self.assertIsNone(cell.height);self.assertIsNone(cell.left);self.assertIsNone(cell.top)
                self.assertTrue(cell.col==12 or cell.col.get('xs')==12)
        for key,control in form.preparation_settings.settings['ligand'].items():
            self.assertIn(control,controls)
            self.assertIsNone(control.height)

    def test_dialogs_resize_in_both_directions_and_remain_scrollable(self):
        page=SimpleNamespace(width=1440,height=1000)
        body=ft.Container(width=960,height=660,content=ft.Column([ft.Text('Settings')]))
        dialog=ft.AlertDialog(content=body)
        fit_dialog(dialog,page)
        for width,height in ((1024,800),(560,760),(360,640),(1440,1000)):
            page.width,page.height=width,height;fit_dialog(dialog,page)
            self.assertLessEqual(body.width,width-64)
            self.assertLessEqual(body.height,height-240)
            self.assertEqual(body.content.scroll,ft.ScrollMode.AUTO)
        self.assertEqual(body.width,960)

    def test_short_confirmation_keeps_natural_height_with_bounded_dialog_scroll(self):
        page=SimpleNamespace(width=900,height=700)
        body=ft.Container(content=ft.Column([ft.Text('Short confirmation')],tight=True))
        dialog=ft.AlertDialog(content=body)
        fit_dialog(dialog,page)
        self.assertIsNone(body.height);self.assertIsNone(body.width)
        self.assertTrue(dialog.scrollable)
        self.assertIsNone(body.content.scroll)

    def test_settings_dialog_uses_responsive_width_but_natural_height(self):
        page=SimpleNamespace(width=1200,height=800)
        body=ft.Container(width=960,height=660,content=ft.Column([ft.Text('Settings')]))
        body._settings_auto_size=True
        fit_dialog(ft.AlertDialog(content=body),page)
        self.assertEqual(body.width,960);self.assertIsNone(body.height)


if __name__=='__main__':unittest.main()
