"""Shared EGG geometry and publication/export presentation.

The regions reproduce BioMolExplorer's local legacy approximation, not the
published SwissADME classifier. Do not describe these as validated BOILED-Egg.
"""
import io
import math

REGIONS=((75,2,75,3,'#FFFFFF'),(42,2.3,47,2.2,'#FFE39A'))
COLORS={'BBB+':'#D55E00','BBB-':'#0072B2','unknown':'#64748B'}


def bounds(nodes):
    xs=[n['x'] for n in nodes];ys=[n['y'] for n in nodes]
    xmin,xmax=min([0]+xs),max([200]+xs);ymin,ymax=min([-2]+ys),max([7]+ys)
    return xmin-0.025*(xmax-xmin),xmax+0.025*(xmax-xmin),ymin-0.04*(ymax-ymin),ymax+0.04*(ymax-ymin)


def tick_values(low,high,count=5):
    raw=(high-low)/count;power=10**math.floor(math.log10(raw))
    step=next(v for v in (1,2,2.5,5,10) if v*power>=raw)*power
    start=math.ceil(low/step)*step
    return [start+i*step for i in range(int((high-start)/step)+1)]


def egg_png(model):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.patches import Ellipse,Patch
    from matplotlib.ticker import MaxNLocator
    fig=Figure(figsize=(12,7),layout='constrained');FigureCanvasAgg(fig)
    ax=fig.subplots();ax.set_facecolor('#EEF2F7');ax.set_axisbelow(True)
    ax.grid(color='#CBD5E1',linewidth=.6,alpha=.45)
    for cx,cy,rx,ry,color in REGIONS:
        ax.add_patch(Ellipse((cx,cy),2*rx,2*ry,facecolor=color,edgecolor='#94A3B8',linewidth=1,zorder=1))
    markers={'BBB+':'o','BBB-':'s','unknown':'^'}
    for label in ('BBB+','BBB-','unknown'):
        nodes=[n for n in model['nodes'] if (n['properties'].get('BBB') if n['properties'].get('BBB') in COLORS else 'unknown')==label]
        if nodes:ax.scatter([n['x'] for n in nodes],[n['y'] for n in nodes],c=COLORS[label],marker=markers[label],
            edgecolors='white',linewidths=.6,alpha=.8,s=45,label=label,zorder=3)
    xmin,xmax,ymin,ymax=bounds(model['nodes'])
    ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax),xlabel='TPSA (Å²)',ylabel='WLOGP (Wildman–Crippen)',title=model['title'])
    ax.xaxis.set_major_locator(MaxNLocator(6));ax.yaxis.set_major_locator(MaxNLocator(6))
    for spine in ('top','right'):ax.spines[spine].set_visible(False)
    handles,labels=ax.get_legend_handles_labels()
    handles+=[Patch(facecolor='white',edgecolor='#94A3B8'),Patch(facecolor='#FFE39A',edgecolor='#94A3B8')]
    labels+=['HIA · local region','BBB · local region']
    ax.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1),ncol=3,fontsize=9,framealpha=.9)
    fig.supxlabel('Local heuristic approximation of passive permeability · not an experimental measure or validated SwissADME output',fontsize=9,color='#64748B')
    stream=io.BytesIO();fig.savefig(stream,format='png',dpi=200);return stream.getvalue()
