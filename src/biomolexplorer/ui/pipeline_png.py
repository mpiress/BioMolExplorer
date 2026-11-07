"""Render the complete canvas model to a bounded, high resolution PNG."""
import io
import math
from collections import Counter

from biomolexplorer import flow
from biomolexplorer.bindings import sources
from biomolexplorer.catalog import TITLES
from .graph_raster import _font
from .localization import default_stage_names

from .pipeline_palette import COLORS
MAX_PIXELS = 32_000_000


def render_pipeline(stages, translate=lambda value: value):
    """Keep canvas positions, ports and curves, independently of viewport zoom."""
    from PIL import Image, ImageDraw
    positions = {s['id']: flow.position(s, i) for i, s in enumerate(stages)}
    by_id = {s['id']: s for s in stages}
    block_counts=Counter(stage.get('block_id') for stage in stages if stage.get('block_id'))
    pad = 130  # Also covers the control points of backwards connections.
    left = min((p['x'] for p in positions.values()), default=0) - pad
    top = min((p['y'] for p in positions.values()), default=0) - pad
    right = max((p['x'] + flow.NODE_WIDTH for p in positions.values()), default=400) + pad
    bottom = max((positions[s['id']]['y'] + flow.node_height(s) for s in stages), default=200) + pad
    width, height = right - left, bottom - top
    scale = min(2.0, math.sqrt(MAX_PIXELS / (width * height)))
    image = Image.new('RGB', (max(1, int(width * scale)), max(1, int(height * scale))), '#F5F7FC')
    draw = ImageDraw.Draw(image)

    def point(x, y):
        return ((x - left) * scale, (y - top) * scale)

    def box(x, y, w, h):
        return (*point(x, y), *point(x + w, y + h))

    def text(x, y, label, size, color, max_width=None):
        font = _font(max(1, round(size * scale)))
        # Match the canvas's one-line labels without letting them escape a block.
        if max_width is not None:
            full = label
            while label and draw.textlength(label, font=font) > max_width * scale:
                label = label[:-1]
            if label != full:
                label = label[:-1] + '…'
        draw.text(point(x, y), label, font=font, fill=color)

    stroke = max(1, round(3 * scale));connection_labels=[];connection_index=0
    for index,stage in enumerate(stages):
        pos = positions[stage['id']]
        for j, port in enumerate(flow.input_ports(stage)):
            binding = stage.get('bindings', {}).get(port['field'])
            for reference in sources(binding) if binding else []:
                source = by_id.get(reference.get('stage'))
                if source is None:
                    continue
                connection_index+=1
                start = positions[source['id']]
                x, y = start['x'] + flow.NODE_WIDTH - 30, start['y'] + 83
                ex, ey = pos['x'] + 28, pos['y'] + 117 + j * 34
                points = []
                for step in range(101):
                    t = step / 100; u = 1 - t
                    points.append(point(u**3*x + 3*u*u*t*(x+100) + 3*u*t*t*(ex-100) + t**3*ex,
                                        u**3*y + 3*u*u*t*y + 3*u*t*t*ey + t**3*ey))
                color=COLORS[TITLES[source['operation']][2]]
                # A light halo keeps crossings distinguishable on dense
                # canvases; the arrowhead terminates at the exact input port.
                draw.line(points, fill='#F5F7FC', width=stroke+max(2,round(4*scale)))
                draw.line(points, fill=color, width=stroke)
                tip=points[-1];previous=points[-5]
                angle=math.atan2(tip[1]-previous[1],tip[0]-previous[0]);length=max(7,12*scale);spread=.55
                arrow=[tip,(tip[0]-length*math.cos(angle-spread),tip[1]-length*math.sin(angle-spread)),
                    (tip[0]-length*math.cos(angle+spread),tip[1]-length*math.sin(angle+spread))]
                draw.polygon(arrow,fill=color)
                connection_labels.append((f'C{connection_index:02d}',points[58],color))

    for label,(lx,ly),color in connection_labels:
        font=_font(max(8,round(10*scale)));bounds=draw.textbbox((0,0),label,font=font)
        tw,th=bounds[2]-bounds[0],bounds[3]-bounds[1];padding=max(2,round(3*scale))
        draw.rounded_rectangle((lx-tw/2-padding,ly-th/2-padding,lx+tw/2+padding,ly+th/2+padding),
            radius=max(2,round(3*scale)),fill='white',outline=color,width=max(1,round(scale)))
        draw.text((lx-tw/2,ly-th/2-bounds[1]),label,font=font,fill=color)

    for stage in stages:
        pos = positions[stage['id']]; x, y = pos['x'], pos['y']
        category = TITLES[stage['operation']][2]; color = COLORS[category]
        draw.rounded_rectangle(box(x, y, flow.NODE_WIDTH, flow.node_height(stage)),
                               radius=round(14*scale), fill='white', outline='#DDE5F0', width=max(1, round(2*scale)))
        draw.rounded_rectangle(box(x+2, y+2, flow.NODE_WIDTH-4, 64), radius=round(12*scale), fill=color)
        stored_block=stage.get('block_id')
        block_id=stored_block if stored_block and block_counts[stored_block]==1 else f'B{index+1:02d}'
        text(x+12, y+10, str(block_id), 10, 'white', 72)
        category_label=translate(category).upper()
        category_width=draw.textlength(category_label,font=_font(max(1,round(10*scale))))/scale
        text(x+flow.NODE_WIDTH-12-category_width, y+10, category_label, 10, 'white', 120)
        name = translate(stage['name']) if stage['name'] in default_stage_names() else stage['name']
        text(x+12, y+34, name, 13, 'white', flow.NODE_WIDTH-24)
        text(x+16, y+76, translate('Saída'), 12, color)
        draw.ellipse(box(x+flow.NODE_WIDTH-35, y+78, 10, 10), fill=color)
        for j, port in enumerate(flow.input_ports(stage)):
            py = y+117+j*34
            draw.ellipse(box(x+23, py-5, 10, 10), fill=color if port['field'] in stage.get('bindings', {}) else 'white',
                         outline=color, width=max(1, round(scale)))
            text(x+45, py-8, translate(port['label']), 12, '#334155', flow.NODE_WIDTH-58)
        footer = 'Resultados fornecidos · concluído' if stage.get('provided_results') else 'Pronto para executar'
        text(x+16, y+flow.node_height(stage)-30, translate(footer), 10, '#64748B', flow.NODE_WIDTH-32)
    output = io.BytesIO()
    image.save(output, format='PNG', dpi=(300, 300))
    return output.getvalue()
