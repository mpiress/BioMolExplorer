"""Bounded graph paint layer; hit testing stays on the exact molecular nodes."""
import io
from functools import lru_cache

import flet as ft
import flet.canvas as canvas


@lru_cache(maxsize=16)
def _font(size):
    from PIL import ImageFont
    try:return ImageFont.truetype('DejaVuSans.ttf',size)
    except OSError:return ImageFont.load_default()


def paint_graph(shapes,width,height):
    """Paint every primitive at the same scene coordinates as the point index."""
    from PIL import Image,ImageDraw
    image=Image.new('RGBA',(max(1,round(width)),max(1,round(height))),(248,250,253,255))
    draw=ImageDraw.Draw(image)
    for shape in shapes:
        paint=getattr(shape,'paint',None)
        color=getattr(paint,'color',None) or '#64748B'
        stroke=max(1,round(getattr(paint,'stroke_width',None) or 1))
        outline=getattr(paint,'style',None)==ft.PaintingStyle.STROKE
        if isinstance(shape,canvas.Line):
            draw.line((shape.x1,shape.y1,shape.x2,shape.y2),fill=color,width=stroke)
        elif isinstance(shape,canvas.Circle):
            bounds=(shape.x-shape.radius,shape.y-shape.radius,shape.x+shape.radius,shape.y+shape.radius)
            draw.ellipse(bounds,outline=color if outline else None,fill=None if outline else color,width=stroke)
        elif isinstance(shape,canvas.Rect):
            draw.rectangle((shape.x,shape.y,shape.x+shape.width,shape.y+shape.height),
                           outline=color if outline else None,fill=None if outline else color,width=stroke)
        elif isinstance(shape,canvas.Text):
            style=shape.style
            draw.text((shape.x,shape.y),shape.value or '',fill=getattr(style,'color',None) or '#172B4D',
                      font=_font(max(8,round(getattr(style,'size',None) or 10))))
    output=io.BytesIO();image.save(output,format='PNG')
    return output.getvalue()
