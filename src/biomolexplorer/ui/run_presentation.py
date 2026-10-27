"""Presentation of run outcomes; result loading remains lazy and authorized."""
import flet as ft
from .run_progress import LABELS,elapsed,ACTIVE
from .localization import verbatim

STATES={
    'removed':('#64748B','#F1F5F9',ft.Icons.DELETE_OUTLINE),
    'succeeded':('#15803D','#F0FDF4',ft.Icons.CHECK_CIRCLE_OUTLINE),
    'failed':('#B91C1C','#FEF2F2',ft.Icons.ERROR_OUTLINE),
    'running':('#2563EB','#EFF6FF',ft.Icons.PLAY_CIRCLE_OUTLINE),
    'awaiting_input':('#B45309','#FFFBEB',ft.Icons.EDIT_OUTLINED),
    'queued':('#2563EB','#EFF6FF',ft.Icons.SCHEDULE),
    'skipped':('#B45309','#FFFBEB',ft.Icons.SKIP_NEXT),
    'cancelled':('#64748B','#F1F5F9',ft.Icons.STOP_CIRCLE_OUTLINED),
    'interrupted':('#64748B','#F1F5F9',ft.Icons.PAUSE_CIRCLE_OUTLINE),
}


def result_state(stage):
    """Execution history is distinct from currently available results."""
    if stage.get('outputs_removed') and not stage.get('artifacts'):
        return 'removed','Resultados excluídos · execute o bloco para gerar novamente'
    if stage.get('results_stale'):
        return 'removed','Resultados desatualizados · execute o bloco para atualizar'
    return stage['status'],('Reaproveitado' if stage.get('reused') else
        'Resultados fornecidos' if stage.get('provided') else LABELS[stage['status']])


def state_icon(status):
    color,background,icon=STATES.get(status,STATES['interrupted'])
    return ft.Container(ft.Icon(icon,size=22,color=color),width=38,height=38,
        alignment=ft.Alignment.CENTER,bgcolor=background,border_radius=12)


def outcome_overview(run,actions):
    color,background,_=STATES.get(run['status'],STATES['interrupted'])
    stages=run['stages']
    computed=[s for s in stages if s.get('operation')!='import_results' and not s.get('provided')]
    removed=bool(computed) and any(s.get('outputs_removed') for s in computed) and not any(s.get('artifacts') for s in computed)
    display_status='removed' if removed else run['status']
    if removed:color,background,_=STATES['removed']
    completed=sum(s['status']=='succeeded' for s in stages)
    total=sum(s.get('configuration',{}).get('enabled',True) for s in stages)
    files=len({str(p) for s in stages for p in s.get('artifacts',[])})
    reused=sum(s['status']=='succeeded' and bool(s.get('reused')) and bool(s.get('artifacts')) for s in stages)
    description=('Explore os resultados de cada etapa abaixo.' if run['status']=='succeeded' else
        'Acompanhe as etapas e consulte seus resultados disponíveis.' if run['status'] in ACTIVE else
        'Consulte os resultados disponíveis e o diagnóstico de cada etapa.')
    if removed:description='Os resultados calculados foram excluídos. Execute os blocos para gerar novamente. O histórico foi preservado.'
    metrics=[('Etapas concluídas (histórico)' if removed else 'Etapas concluídas',f'{completed} / {total}'),('Arquivos registrados',str(files)),
        ('Etapas reaproveitadas',str(reused))]
    if run['status'] not in ACTIVE:
        metrics.append(('Duração',elapsed(run['created'],run['updated'])))
    cards=[ft.Container(col={'xs':6,'lg':3},padding=12,bgcolor='#FFFFFF',border_radius=10,
        border=ft.Border.all(1,'#E2E8F0'),content=ft.Column([
            verbatim(ft.Text(value,size=21,weight=ft.FontWeight.W_600,color='#172B4D')),
            ft.Text(label,size=12,color='#64748B')],spacing=3)) for label,value in metrics]
    return ft.Container(padding=18,bgcolor=background,border_radius=14,content=ft.Column([
        ft.Row([state_icon(display_status),ft.Column([
            ft.Text('Resultados excluídos' if removed else LABELS[run['status']],size=22,weight=ft.FontWeight.W_600,color=color),
            ft.Text(description,size=13,color='#64748B')],expand=True,spacing=4)],spacing=12),
        ft.ResponsiveRow(cards,spacing=10,run_spacing=10),
        *([ft.Row(actions,wrap=True,spacing=10)] if actions else [])],spacing=16))
