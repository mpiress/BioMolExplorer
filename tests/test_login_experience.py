"""Login presentation must retain authentication and finite motion lifecycles."""
import asyncio
from types import SimpleNamespace
import unittest
from unittest.mock import AsyncMock, patch

import flet as ft

from biomolexplorer.ui.app import WorkspaceUI
from test_frontend_navigation import Page, walk


class LoginExperienceTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.page = Page()
        self.ui = WorkspaceUI(self.page, SimpleNamespace(login=lambda *args: None), SimpleNamespace(), language='en')
        self.ui.show_login()

    async def test_pending_authentication_is_single_flight_and_restores_error_controls(self):
        entered = asyncio.Event()
        release = asyncio.Event()

        async def reject(function, *args):
            entered.set()
            await release.wait()
            raise ValueError('E-mail ou senha inválidos.')

        self.ui.call = AsyncMock(side_effect=reject)
        experience = self.ui._login_experience
        first = asyncio.create_task(experience.submit.on_click(None))
        await entered.wait()
        self.assertTrue(experience.submit.disabled)
        self.assertIsInstance(experience.submit.icon, ft.ProgressRing)
        await experience.submit.on_click(None)
        self.assertEqual(self.ui.call.await_count, 1)
        release.set()
        await first
        self.assertFalse(experience.submit.disabled)
        self.assertTrue(experience.error_box.visible)
        self.assertEqual(experience.error.value, 'Invalid email or password.')
        self.assertEqual(experience.submit.icon, ft.Icons.ARROW_FORWARD)
        self.assertIsNone(self.ui.token)

    async def test_registration_mismatch_keeps_validation_before_backend(self):
        self.ui.show_login(True)
        fields = {control.label: control for control in walk(self.page.controls) if isinstance(control, ft.TextField)}
        fields['Password'].value = 'original-password'
        fields['Confirm password'].value = 'different-password'
        self.ui.call = AsyncMock()
        experience = self.ui._login_experience
        await experience.submit.on_click(None)
        self.ui.call.assert_not_awaited()
        self.assertTrue(experience.error_box.visible)
        self.assertFalse(experience.submit.disabled)
        self.assertEqual(fields['Password'].value, 'original-password')
        self.assertIsNone(self.ui.token)

    async def test_replacing_view_cancels_reveal_and_stale_pointer_updates(self):
        old = self.ui._login_experience
        old.entry_task = asyncio.create_task(asyncio.sleep(30))
        task = old.entry_task
        self.ui.show_login(True)
        await asyncio.sleep(0)
        self.assertTrue(task.cancelled())
        self.assertFalse(old.active)
        self.assertFalse(old.valid())
        old.pointer(None)
        old.drift()
        old.leave(None)
        current = self.ui._login_experience
        current.toggle_motion(None)
        self.assertTrue(all(group.opacity == 1 for group in current.groups))
        self.assertTrue(all(group.animate_offset is None for group in current.groups))
        self.assertTrue(current.valid())

    async def test_official_identity_visual_and_site_link_are_real_flet_controls(self):
        experience = self.ui._login_experience
        texts=[control.value for control in walk(experience.control) if isinstance(control,ft.Text)]
        self.assertIn('Molecular Data Exploration for',texts)
        self.assertIn('Intelligent Drug Discovery',texts)
        self.assertIn('An integrated computational platform for molecular exploration and information retrieval, '
            'designed to support computer-aided drug discovery through molecular similarity analysis, virtual '
            'screening, redocking validation, and interactive 3D inspection of docking poses and receptor–ligand '
            'interactions.', texts)
        for capability in ('PDB', 'ChEMBL', 'PubChem', 'ZINC', 'Similarity', 'Docking', 'ADMET', 'Visualization'):
            self.assertIn(capability, texts)
        self.assertTrue(any(isinstance(control,ft.Semantics) and 'Carbon connected' in str(control.label)
            for control in walk(experience.control)))
        link=next(control for control in walk(experience.control) if isinstance(control,ft.TextButton)
            and control.content=='Visit the official BioMolExplorer website')
        with patch('flet.UrlLauncher.launch_url',new_callable=AsyncMock) as launch:
            await link.on_click(None)
            self.assertEqual(launch.await_args.args[0],'https://mpiress.github.io/BioMolExplorer/')
            self.assertEqual(launch.await_args.kwargs['web_only_window_name'],'_blank')
