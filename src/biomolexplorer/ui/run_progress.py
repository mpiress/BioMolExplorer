"""Accessible run feedback, independent of the pipeline editor's redraws."""
import time
import flet as ft
from .localization import verbatim
from .feedback import readable_error

ACTIVE = {'queued', 'running', 'awaiting_input'}
LABELS = {'queued': 'Aguardando', 'running': 'Executando', 'succeeded': 'Concluído',
          'failed': 'Falhou', 'cancelled': 'Cancelado', 'interrupted': 'Interrompido', 'skipped': 'Não executado',
          'awaiting_input': 'Aguardando seleção de arquivos'}


def elapsed(start, end):
    seconds = max(0, int(end - start))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours:02}:{minutes:02}:{seconds:02}'


class RunProgress:
    def __init__(self, ui, run):
        self.ui, self.run = ui, run
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
                ft.TextButton('Acompanhar execução', on_click=lambda e: self.show())]))
        self.minimize = ft.TextButton('Minimizar', on_click=self.dismiss)
        self.cancel = ft.TextButton('Cancelar execução', on_click=self.request_cancel)
        self.log = ft.TextButton('Ver log da etapa', on_click=self.open_log)
        self.results = ft.TextButton('Ver resultados', on_click=self.open_results)
        self.configure = ft.TextButton('Configurar entradas e continuar', on_click=self.configure_pending)
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
        active = run['status'] in ACTIVE
        waiting = run['status']=='awaiting_input'
        current = next((s for s in run['stages'] if s['status'] == 'running'), None) if active else None
        self.heading.value = {'queued': 'Pipeline na fila', 'running': 'Seu pipeline está em execução',
            'awaiting_input': 'Escolha os arquivos para continuar',
            'succeeded': 'Pipeline concluído', 'failed': 'A execução encontrou um erro',
            'cancelled': 'Execução cancelada', 'interrupted': 'Execução interrompida'}.get(run['status'], 'Acompanhamento do pipeline')
        phase = (current or {}).get('progress') or {}
        self.message.value = (run.get('error') or 'Configure as entradas do próximo bloco.' if waiting else
            f"Etapa atual: {current['name']}\n{phase.get('message') or 'Processando os dados…'}" if current else
            'Aguardando um executor disponível…' if active else
            'Os arquivos gerados estão disponíveis na aba Execuções.' if run['status'] == 'succeeded' else
            'Consulte o log da etapa para entender o ocorrido antes de executar novamente.')
        if active and isinstance(phase.get('completed'), int) and isinstance(phase.get('total'), int):
            self.message.value += f"\nRegistros processados: {phase['completed']} de {phase['total']}"
        self.message_semantics.label = self.message.value
        end = (now if now is not None else time.time()) if active else run['updated']
        self.timing.value = 'Tempo decorrido: ' + elapsed(run['created'], end)
        if current and current.get('started_at'):
            self.timing.value += ' · Etapa: ' + elapsed(current['started_at'], end)
        completed = sum(s['status'] == 'succeeded' for s in run['stages'])
        total = sum(s.get('configuration', {}).get('enabled', True) for s in run['stages'])
        self.count.value = f'{completed} de {total} etapas concluídas · Execução {run["id"][:8]}'
        self.bar.visible = active and not waiting
        self.error.value = '' if waiting else readable_error(run.get('error'))
        self.error.visible = bool(self.error.value)
        self.error_semantics.label = self.error.value
        self.error_semantics.visible = self.error.visible
        self.hint.value = ('Cancelamento solicitado. Aguardando a interrupção do processo…' if self.cancel_pending and active else
            'O pipeline está pausado. Selecione os arquivos do bloco indicado e aplique a configuração para continuar.' if waiting else
            'Consultas externas podem demorar e repetir tentativas automaticamente. Você pode minimizar esta janela; a execução continua.' if current and current['operation'] in ('retrieve_compounds', 'expand_similar_compounds', 'retrieve_complex', 'retrieve_zinc') else
            'Você pode minimizar esta janela e continuar usando o projeto.' if active else
            'Os detalhes técnicos e os resultados ficam na aba Execuções.')
        self.rows.controls = [ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE if s['status'] == 'succeeded' else
                ft.Icons.ERROR_OUTLINE if s['status'] == 'failed' else ft.Icons.TIMELAPSE,
                color='#15803D' if s['status'] == 'succeeded' else '#B91C1C' if s['status'] == 'failed' else '#64748B', size=19),
            verbatim(ft.Text(s['name'], expand=True)), ft.Text('Reaproveitado' if s.get('reused') else 'Resultados fornecidos' if s.get('provided') else LABELS[s['status']], size=12)], spacing=10) for s in run['stages']]
        self.log.visible = any(s.get('log_path') for s in run['stages'])
        self.cancel.visible = active and self.ui.current is not None and self.ui.current['role'] != 'viewer'
        self.cancel.disabled = self.cancel_pending
        self.configure.visible = waiting and self.ui.current is not None and self.ui.current['role']!='viewer'
        self.minimize.content = 'Minimizar' if active else 'Fechar'
        self.results.visible = not active
        self.badge_text.value = f"{LABELS[run['status']]} · {current['name'] if current else self.count.value} · {elapsed(run['created'], end)}"
