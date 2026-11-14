"""The V5 presentation must retain authentication and responsive lifecycles."""
import asyncio
from types import SimpleNamespace
import unittest
from unittest.mock import AsyncMock

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
        self.assertEqual(experience.submit.icon, ft.Icons.LOGIN)
        self.assertIsNone(self.ui.token)

    async def test_registration_mismatch_keeps_validation_before_backend(self):
        self.ui.show_login(True)
        fields = {control.hint_text: control for control in walk(self.page.controls) if isinstance(control, ft.TextField)}
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

    async def test_v5_submit_uses_existing_authentication_and_opens_workspace(self):
        self.ui.call = AsyncMock(return_value='existing-session')
        self.ui.show_workspace = AsyncMock()
        experience = self.ui._login_experience
        experience.email.value = 'researcher@example.org'
        experience.password.value = 'existing-password'
        await experience.submit.on_click(None)
        self.ui.call.assert_awaited_once_with(
            self.ui.store.login, 'researcher@example.org', 'existing-password')
        self.ui.show_workspace.assert_awaited_once_with()
        self.assertEqual(self.ui.token, 'existing-session')

    async def test_password_recovery_reuses_card_for_otp_and_password_reset(self):
        experience = self.ui._login_experience
        experience.email.value='researcher@example.org'
        experience.forgot.on_click(None)
        self.assertEqual(experience.recovery_email_field.value,'researcher@example.org')
        self.ui.recovery_service.request=AsyncMock(return_value='12345678')
        self.ui.recovery_service.validate=AsyncMock(return_value='reset-grant')
        self.ui.recovery_service.reset=AsyncMock()
        async def direct(function,*args):return await function(*args)
        self.ui.call=direct
        await experience.send_recovery_code()
        self.assertEqual(experience.code_field.max_length,8)
        experience.code_field.value='12345678'
        await experience.validate_recovery_code()
        self.assertEqual(experience.reset_grant,'reset-grant')
        experience.new_password.value='new-password'
        experience.new_confirmation.value='different-password'
        await experience.save_new_password()
        self.ui.recovery_service.reset.assert_not_awaited()
        experience.new_confirmation.value='new-password'
        await experience.save_new_password()
        self.ui.recovery_service.reset.assert_awaited_once_with(
            'researcher@example.org','reset-grant','new-password')
        self.assertEqual(self.ui._login_experience.email.value,'researcher@example.org')

    async def test_failed_email_delivery_keeps_recovery_email_step_visible(self):
        experience=self.ui._login_experience
        experience.email.value='researcher@example.org';experience.begin_recovery()
        self.ui.recovery_service.request=AsyncMock(side_effect=RuntimeError(
            'Não foi possível enviar o código por e-mail. Verifique a configuração SMTP e tente novamente.'))
        async def direct(function,*args):return await function(*args)
        self.ui.call=direct
        await experience.send_recovery_code()
        self.assertTrue(hasattr(experience,'recovery_email_field'))
        self.assertFalse(hasattr(experience,'code_field'))
        self.assertIn('SMTP',experience.recovery_error.value)

    async def test_replacing_view_cancels_stale_tasks_and_mobile_hides_artwork(self):
        old = self.ui._login_experience
        old.entry_task = asyncio.create_task(asyncio.sleep(30))
        task = old.entry_task
        self.ui.show_login(True)
        await asyncio.sleep(0)
        self.assertTrue(task.cancelled())
        self.assertFalse(old.active)
        self.assertFalse(old.valid())
        current = self.ui._login_experience
        self.page.width=360;self.page.height=640
        current.resize()
        self.assertFalse(current.left.visible)
        self.assertEqual(current.card.width,332)
        self.assertLessEqual(current.card.height,564)
        self.assertTrue(current.valid())

    async def test_v5_uses_bundled_artwork_without_card_logo_and_real_controls(self):
        experience = self.ui._login_experience
        images=[control for control in walk(experience.control) if isinstance(control,ft.Image)]
        self.assertEqual(len(images),2)
        self.assertTrue(all(isinstance(image.src,bytes) for image in images))
        self.assertEqual(experience.title.value,'Sign In')
        self.assertEqual(experience.email.hint_text,'Username or email')
        self.assertEqual(experience.password.hint_text,'Password')
        self.assertEqual(experience.submit.content,'Sign in to workspace')
        self.assertFalse(experience.remember.disabled)
        self.assertTrue(experience.left.visible)
