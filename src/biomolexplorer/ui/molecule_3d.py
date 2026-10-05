"""Offline molecular ball-and-stick view with native drag rotation and zoom."""
import math
import flet as ft
import flet.canvas as canvas
from .zoom import MIN_SCALE, MAX_SCALE

COLORS={'C':'#475569','H':'#E2E8F0','N':'#2563EB','O':'#DC2626','S':'#EAB308','P':'#F97316',
        'F':'#16A34A','Cl':'#16A34A','Br':'#B45309','I':'#7C3AED'}


class Molecule3D:
    def __init__(self,page,model,width=620,height=360):
        self.page,self.model,self.width,self.height=page,model,width,height
        self.yaw,self.pitch,self.scale=.45,.25,1.0
        self.drawing=canvas.Canvas(width=width,height=height)
        # Logarithmic positions make both very small and very large scales
        # reachable without sacrificing control around the normal size.
        self.zoom_label=ft.Text('100%',width=75)
        self.zoom=ft.Slider(min=math.log10(MIN_SCALE),max=math.log10(MAX_SCALE),
            value=0,on_change=self.change_zoom,label='100%')
        self.redraw()

    def redraw(self):
        atoms=self.model['atoms']
        center=[sum(a[axis] for a in atoms)/len(atoms) for axis in ('x','y','z')]
        cy,sy,cp,sp=math.cos(self.yaw),math.sin(self.yaw),math.cos(self.pitch),math.sin(self.pitch)
        radius=max((math.sqrt(sum((a[k]-center[i])**2 for i,k in enumerate(('x','y','z')))) for a in atoms),default=1)
        factor=min(self.width,self.height)*.38/max(radius,1)*self.scale
        points=[]
        for atom in atoms:
            x,y,z=[atom[k]-center[i] for i,k in enumerate(('x','y','z'))]
            rx,rz=x*cy+z*sy,-x*sy+z*cy
            ry,rz=y*cp-rz*sp,y*sp+rz*cp
            points.append((self.width/2+rx*factor,self.height/2-ry*factor,rz))
        shapes=[]
        for bond in sorted(self.model['bonds'],key=lambda b:(points[b['a']][2]+points[b['b']][2])/2):
            a,b=points[bond['a']],points[bond['b']]
            shapes.append(canvas.Line(a[0],a[1],b[0],b[1],paint=ft.Paint(color='#94A3B8',stroke_width=(3 if bond['order']<2 else 5)*self.scale)))
        for index in sorted(range(len(atoms)),key=lambda i:points[i][2]):
            atom,point=atoms[index],points[index]
            radius=(5 if atom['element']=='H' else 11)*self.scale
            shapes.append(canvas.Circle(point[0],point[1],radius=radius,paint=ft.Paint(color=COLORS.get(atom['element'],'#DB2777'))))
            if atom['element']!='H':
                shapes.append(canvas.Text(point[0]-4*self.scale,point[1]-7*self.scale,atom['element'],style=ft.TextStyle(color='#FFFFFF',size=11*self.scale)))
        self.drawing.shapes=shapes

    def rotate(self,e):
        if e.local_delta is None:return
        self.yaw+=e.local_delta.x*.01
        self.pitch+=e.local_delta.y*.01
        self.redraw();self.page.update()

    def change_zoom(self,e):
        self.scale=max(MIN_SCALE,min(MAX_SCALE,10**float(e.control.value)))
        self.zoom_label.value=self.zoom.label=f'{self.scale*100:g}%'
        self.redraw();self.page.update()

    def reset(self,e):
        self.yaw,self.pitch,self.scale=.45,.25,1.0
        self.zoom.value=0
        self.zoom_label.value=self.zoom.label='100%'
        self.redraw();self.page.update()

    def build(self):
        return ft.Column([ft.Container(ft.GestureDetector(content=self.drawing,on_pan_update=self.rotate,drag_interval=25),
            bgcolor='#F8FAFC',border_radius=16,clip_behavior=ft.ClipBehavior.HARD_EDGE),
            ft.Row([ft.Text('Zoom'),ft.Container(self.zoom,expand=True),self.zoom_label,ft.TextButton('Recentrar 3D',on_click=self.reset)]),
            ft.Text('Arraste a molécula para girar. Conformero gerado localmente a partir do SMILES; não representa uma pose de docking.',size=12,color='#64748B')],spacing=12)
