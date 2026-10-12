"""Lightweight responsive charts; resizing never recalculates scientific metrics."""
import flet as ft
import flet.canvas as cv


def card(title, content, col=None):
    return ft.Container(col=col or {'xs':12,'md':6},padding=14,bgcolor='#FFFFFF',
        border=ft.Border.all(1,'#E2E8F0'),border_radius=12,
        content=ft.Column([ft.Text(title,size=14,weight=ft.FontWeight.W_600),content],spacing=12,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH))


def distribution(title, bins, color='#2563EB', x_label='Valor', y_label='Compostos'):
    """At most 25 bars, with actual axes, full values on hover and explicit truncation."""
    items=list(bins); shown=items[:25]
    chart=cv.Canvas(height=220,expand=True)
    tip=ft.Text(size=11,color='#64748B',value=f'{x_label} · {y_label}')
    regions=[]
    def draw(width):
        width=max(160,width);left=42;top=15;bottom=185;right=width-10
        peak=max((count for _,count in shown),default=1) or 1
        shapes=[];regions.clear()
        for fraction in (0,.5,1):
            y=bottom-(bottom-top)*fraction
            shapes.extend([cv.Line(left,y,right,y,paint=ft.Paint(color='#E2E8F0',stroke_width=1)),
                cv.Text(2,y-7,f'{peak*fraction:.3g}',style=ft.TextStyle(size=10,color='#64748B'))])
        stride=(right-left)/max(1,len(shown))
        for i,(label,count) in enumerate(shown):
            x=left+i*stride+2;w=max(1,stride-4);h=(bottom-top)*count/peak
            shapes.append(cv.Rect(x,bottom-h,w,h,paint=ft.Paint(color=color)))
            regions.append((x,x+w,str(label),count))
            if len(shown)<=8 or i%max(1,len(shown)//5)==0:
                shapes.append(cv.Text(x,bottom+6,str(label)[:12],style=ft.TextStyle(size=9,color='#64748B')))
        if not shown:shapes.append(cv.Text(left,80,'Sem dados disponíveis',style=ft.TextStyle(size=12)))
        chart.shapes=shapes
    def resize(e):draw(e.width);chart.update()
    def hover(e):
        if e.local_position is None:return
        tip.value=next((f'{x_label}: {label} · {y_label}: {count}' for a,b,label,count in regions if a<=e.local_position.x<=b),f'{x_label} · {y_label}')
        tip.update()
    chart.on_resize=resize;draw(360)
    result=card(title,ft.Column([ft.Container(ft.GestureDetector(chart,on_hover=hover,hover_interval=80),height=220),tip,
        *([ft.Text('Primeiras 25 classes; distribuição completa no relatório JSON.',size=11)] if len(items)>25 else [])]))
    return result
