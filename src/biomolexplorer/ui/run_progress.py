"""Accessible run feedback, independent of the pipeline editor's redraws."""
import time
import flet as ft
from .localization import verbatim, stage_control, Translator, default_stage_names
from .feedback import readable_error
from biomolexplorer.i18n import RuntimeMessages

ACTIVE = {'queued', 'running', 'awaiting_input'}


def elapsed(start, end):
    seconds = max(0, int(end - start))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours:02}:{minutes:02}:{seconds:02}'


class RunProgress:
    def __init__(self, ui, run):
        self.ui, self.run = ui, run
        self.language = getattr(getattr(ui.page, 'translator', None), 'language', getattr(ui, 'language', 'pt'))
        self.messages = RuntimeMessages(self.language)
        self.translate = Translator(self.language)
        self.token, self.project_id = ui.token, run['project_id']
        self.is_open = False
        self.closing = False
        self.reopen_requested = False
        self.cancel_pending = False
        self.heading = ft.Text(size=22, weight=ft.FontWeight.W_600)
        self.message = ft.Text(size=15, selectable=True)
        self.message_semantics = ft.Semantics(container=True, exclude_semantics=True, content=self.message)
        self.timing = ft.Text(size=13, color='#64748B')
        self.count = ft.Text(size=13, color='#64748B')
        self.bar = ft.ProgressBar(color='#2563EB')
        self.rows = ft.Column(spacing=12)
        self.error = ft.Text(color='#B91C1C', selectable=True)
        self.error_semantics = ft.Semantics(container=True, exclude_semantics=True, content=self.error)
        self.hint = ft.Text(size=13, color='#64748B')
        self.badge_text = ft.Text(size=14, weight=ft.FontWeight.W_600)
        self.badge = ft.Container(bgcolor='#EAF1FF', padding=14, border_radius=12,
            content=ft.Row([ft.Icon(ft.Icons.TIMELAPSE, color='#2563EB'),
                ft.Container(self.badge_text, expand=True),
                ft.TextButton(self.messages('run.action.follow'), on_click=lambda e: self.show())]))
        self.minimize = ft.TextButton(self.messages('run.action.minimize'), on_click=self.dismiss)
        self.cancel = ft.TextButton(self.messages('run.action.cancel'), on_click=self.request_cancel)
        self.log = ft.TextButton(self.messages('run.action.log'), on_click=self.open_log)
        self.results = ft.TextButton(self.messages('run.action.results'), on_click=self.open_results)
        self.configure = ft.TextButton(self.messages('run.action.configure'), on_click=self.configure_pending)
        self.dialog = ft.AlertDialog(modal=True, title=self.heading,
            content=ft.Container(width=600, height=max(220, min(420, 260 + len(run['stages']) * 34, (ui.page.height or 900) - 260)),
                content=ft.Column([self.message_semantics, self.bar, self.timing, self.count, self.rows,
                                   self.error_semantics, self.hint], spacing=18, scroll=ft.ScrollMode.AUTO)),
            actions=[self.configure, self.log, self.cancel, self.minimize, self.results],
            on_dismiss=self.dismissed)
        self.update(run)

    def valid(self):
        return (self.ui.token == self.token and self.ui.current is not None
                and self.ui.current['id'] == self.project_id
                and getattr(self.ui, 'run_progress', None) is self)

    def show(self):
        if not self.valid():
            return
        if self.closing:
            self.reopen_requested = True
        elif not self.is_open:
            self.ui.page.show_dialog(self.dialog)
            self.is_open = True

    def dismissed(self,e):
        self.is_open = False
        self.closing = False
        if self.reopen_requested:
            self.reopen_requested = False
            self.show()

    async def dismiss(self, e=None):
        if not self.valid():
            return
        self.dialog.open = False
        self.is_open = False
        self.closing = True
        if self.run['status'] not in ACTIVE:
            self.ui.run_progress = None
        await self.ui.guard(self.ui.draw_project)

    async def request_cancel(self, e):
        if not self.valid() or self.ui.current['role'] == 'viewer' or self.run['status'] not in ACTIVE or self.cancel_pending:
            return
        self.cancel_pending = True
        self.update(self.run)
        self.ui.page.update()
        async def action():
            try:
                await self.ui.cancel_run(self.run['id'])
            except Exception:
                self.cancel_pending = False
                self.update(self.run)
                raise
        await self.ui.guard(action)

    async def open_log(self, e):
        if not self.valid():
            return
        stage = next((s for s in self.run['stages'] if s['status'] in ('running', 'failed', 'cancelled') and s.get('log_path')), None)
        if stage is None:
            stage = next((s for s in reversed(self.run['stages']) if s.get('log_path')), None)
        if stage:
            await self.ui.guard(lambda: self.ui.open_stage_log(self.project_id, stage['log_path']))

    async def open_results(self, e):
        if self.valid():
            self.ui.tab = 'Execuções'
            await self.dismiss()

    async def configure_pending(self, e):
        if self.valid():await self.ui.guard(lambda:self.ui.configure_pending(self.run['id']))

    def update(self, run, now=None):
        self.run = run
        m = self.messages
        active = run['status'] in ACTIVE
        waiting = run['status'] == 'awaiting_input'
        current = next((s for s in run['stages'] if s['status'] == 'running'), None) if active else None
        heading_key = 'run.heading.' + run['status']
        self.heading.value = m(heading_key if heading_key in m.catalog else 'run.heading.default')
        phase = (current or {}).get('progress') or {}
        key = phase.get('message_key')
        if key in m.catalog:
            arguments = dict(phase.get('message_args') or {})
            if key == 'preparation.step' and 'step' in arguments:
                arguments['step'] = self.translate(arguments['step'])
            phase_message = m(key, **arguments)
        else:
            phase_message = self.translate(phase.get('message') or '')
        name = (current or {}).get('name', '')
        if name in default_stage_names():name = self.translate(name)
        self.message.value = (self.translate(run.get('error')) or m('run.message.configure') if waiting else
            m('run.message.current', stage=name, message=phase_message or m('run.message.processing')) if current else
            m('run.message.queued') if active else
            m('run.message.ready') if run['status'] == 'succeeded' else m('run.message.error'))
        if active and isinstance(phase.get('completed'), int) and isinstance(phase.get('total'), int):
            self.message.value += m('run.message.records', completed=phase['completed'], total=phase['total'])
        self.message_semantics.label = self.message.value
        end = (now if now is not None else time.time()) if active else run['updated']
        self.timing.value = m('run.time.elapsed', time=elapsed(run['created'], end))
        if current and current.get('started_at'):
            self.timing.value += m('run.time.stage', time=elapsed(current['started_at'], end))
        completed = sum(s['status'] == 'succeeded' for s in run['stages'])
        total = sum(s.get('configuration', {}).get('enabled', True) for s in run['stages'])
        self.count.value = m('run.count', completed=completed, total=total, run=run['id'][:8])
        self.bar.visible = active and not waiting
        self.error.value = '' if waiting else self.translate(readable_error(run.get('error'),self.language))
        self.error.visible = bool(self.error.value)
        self.error_semantics.label = self.error.value
        self.error_semantics.visible = self.error.visible
        hint = ('cancelling' if self.cancel_pending and active else 'waiting' if waiting else
            'external' if current and current['operation'] in ('retrieve_compounds', 'retrieve_pubchem',
                'expand_similar_compounds', 'retrieve_structures', 'retrieve_zinc') else
            'active' if active else 'done')
        self.hint.value = m('run.hint.' + hint)
        self.rows.controls = [ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE if s['status'] == 'succeeded' else
                ft.Icons.ERROR_OUTLINE if s['status'] == 'failed' else ft.Icons.TIMELAPSE,
                color='#15803D' if s['status'] == 'succeeded' else '#B91C1C' if s['status'] == 'failed' else '#64748B', size=19),
            stage_control(ft.Text(s['name'], expand=True), s),
            verbatim(ft.Text(m('run.status.' + ('reused' if s.get('reused') else 'provided' if s.get('provided') else s['status'])), size=12))], spacing=10) for s in run['stages']]
        self.log.visible = any(s.get('log_path') for s in run['stages'])
        self.cancel.visible = active and self.ui.current is not None and self.ui.current['role'] != 'viewer'
        self.cancel.disabled = self.cancel_pending
        self.configure.visible = waiting and self.ui.current is not None and self.ui.current['role'] != 'viewer'
        self.minimize.content = m('run.action.minimize' if active else 'run.action.close')
        self.results.visible = not active
        self.badge_text.value = f"{m('run.status.' + run['status'])} · {name if current else self.count.value} · {elapsed(run['created'], end)}"
        # Already rendered from this run's session catalog: never reverse-translate
        # the composite strings, their filenames or user-selected stage names.
        for control in (self.heading, self.message, self.message_semantics, self.timing, self.count,
                self.error, self.error_semantics, self.hint, self.badge_text, self.minimize,
                self.cancel, self.log, self.results, self.configure):
            verbatim(control)
