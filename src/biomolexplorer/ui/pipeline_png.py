"""Publication-quality, off-screen renderers for the complete pipeline."""
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
EDGE_LABEL_SIZE = 14  # px at 96 DPI: 10.5 pt


def _connections(stages):
    by_id = {stage['id']: stage for stage in stages}
    result = []
    for target in stages:
        for port_index, port in enumerate(flow.input_ports(target)):
            binding = target.get('bindings', {}).get(port['field'])
            for reference in sources(binding) if binding else []:
                source = by_id.get(reference.get('stage'))
                if source is not None:
                    result.append((source, target, port_index))
    return result


def render_pipeline(stages, translate=lambda value: value, mode='structural'):
    """Render ``structural`` (canvas-faithful) or compact ``scientific`` PNG."""
    if mode not in ('structural', 'scientific'):
        raise ValueError('Pipeline image mode must be structural or scientific.')
    from PIL import Image, ImageDraw

    connections = _connections(stages)
    scientific = mode == 'scientific'
    node_width = 230 if scientific else flow.NODE_WIDTH
    node_height = 80 if scientific else None
    # Both exports honor the user's saved spatial organization.  The scientific
    # mode changes only node presentation and resolution, so a carefully
    # arranged canvas never gains new crossings or loses its visual grouping.
    positions = {stage['id']: flow.position(stage, index) for index, stage in enumerate(stages)}
    pad = 75 if scientific else 130
    left = min((point['x'] for point in positions.values()), default=0) - pad
    top = min((point['y'] for point in positions.values()), default=0) - pad
    right = max((point['x'] + node_width for point in positions.values()), default=400) + pad
    bottom = max((positions[stage['id']]['y'] + (node_height or flow.node_height(stage))
                  for stage in stages), default=200) + pad
    width, height = right - left, bottom - top
    preferred_scale = 3.0 if scientific else 2.0
    scale = min(preferred_scale, math.sqrt(MAX_PIXELS / max(1, width * height)))
    image = Image.new('RGBA', (max(1, round(width * scale)), max(1, round(height * scale))),
                      (255, 255, 255, 0))
    draw = ImageDraw.Draw(image)

    def point(x, y):
        return ((x - left) * scale, (y - top) * scale)

    def box(x, y, box_width, box_height):
        return (*point(x, y), *point(x + box_width, y + box_height))

    def text(x, y, label, size, color, max_width=None, weight=False):
        font = _font(max(1, round(size * scale)))
        label = str(label)
        if max_width is not None:
            full = label
            while label and draw.textlength(label, font=font) > max_width * scale:
                label = label[:-1]
            if label != full:
                label = label[:-1] + '…'
        draw.text(point(x, y), label, font=font, fill=color,
                  stroke_width=round(.25 * scale) if weight else 0)

    stroke = max(1, round((2.2 if scientific else 3) * scale))
    labels = []
    for connection_index, (source, target, port_index) in enumerate(connections, 1):
        start, end = positions[source['id']], positions[target['id']]
        if scientific:
            x, y = start['x'] + node_width, start['y'] + node_height / 2
            ex, ey = end['x'], end['y'] + node_height / 2
        else:
            x, y = start['x'] + flow.NODE_WIDTH - 30, start['y'] + 83
            ex, ey = end['x'] + 28, end['y'] + 117 + port_index * 34
        control = max(55, abs(ex - x) * .32)
        points = []
        for step in range(101):
            t = step / 100; u = 1 - t
            points.append(point(u**3*x + 3*u*u*t*(x+control) + 3*u*t*t*(ex-control) + t**3*ex,
                                u**3*y + 3*u*u*t*y + 3*u*t*t*ey + t**3*ey))
        color = COLORS[TITLES[source['operation']][2]]
        draw.line(points, fill=(255, 255, 255, 230), width=stroke + max(2, round(4 * scale)))
        draw.line(points, fill=color, width=stroke)
        tip, previous = points[-1], points[-5]
        angle = math.atan2(tip[1]-previous[1], tip[0]-previous[0])
        length, spread = max(7, 10*scale), .55
        arrow = [tip, (tip[0]-length*math.cos(angle-spread), tip[1]-length*math.sin(angle-spread)),
                 (tip[0]-length*math.cos(angle+spread), tip[1]-length*math.sin(angle+spread))]
        draw.polygon(arrow, fill=color)
        labels.append((f'C{connection_index:02d}', points[55], color))

    for label, (label_x, label_y), color in labels:
        font = _font(max(8, round(EDGE_LABEL_SIZE * scale)))
        bounds = draw.textbbox((0, 0), label, font=font, stroke_width=max(0, round(.2*scale)))
        text_width, text_height = bounds[2]-bounds[0], bounds[3]-bounds[1]
        padding = max(3, round(3 * scale))
        draw.rounded_rectangle((label_x-text_width/2-padding, label_y-text_height/2-padding,
                                label_x+text_width/2+padding, label_y+text_height/2+padding),
                               radius=max(2, round(4*scale)), fill='white', outline=color,
                               width=max(1, round(scale)))
        draw.text((label_x-text_width/2, label_y-text_height/2-bounds[1]), label, font=font,
                  fill=color, stroke_width=max(0, round(.2*scale)))

    block_counts = Counter(stage.get('block_id') for stage in stages if stage.get('block_id'))
    for index, stage in enumerate(stages):
        pos = positions[stage['id']]; x, y = pos['x'], pos['y']
        category = TITLES[stage['operation']][2]; color = COLORS[category]
        stored = stage.get('block_id')
        block_id = stored if stored and block_counts[stored] == 1 else f'B{index+1:02d}'
        name = translate(stage['name']) if stage['name'] in default_stage_names() else stage['name']
        if scientific:
            draw.rounded_rectangle(box(x, y, node_width, node_height), radius=round(7*scale),
                                   fill='white', outline=color, width=max(2, round(3*scale)))
            draw.rectangle(box(x, y, 9, node_height), fill=color)
            text(x+22, y+13, name, 15, '#172B4D', node_width-34, True)
            text(x+22, y+46, f'{block_id}  ·  {translate(category)}', 11, color, node_width-34, True)
            continue
        height = flow.node_height(stage)
        draw.rounded_rectangle(box(x, y, flow.NODE_WIDTH, height), radius=round(14*scale),
                               fill='white', outline='#DDE5F0', width=max(1, round(2*scale)))
        draw.rounded_rectangle(box(x+2, y+2, flow.NODE_WIDTH-4, 64), radius=round(12*scale), fill=color)
        text(x+12, y+10, block_id, 10, 'white', 72)
        category_label = translate(category).upper()
        category_width = draw.textlength(category_label, font=_font(max(1, round(10*scale)))) / scale
        text(x+flow.NODE_WIDTH-12-category_width, y+10, category_label, 10, 'white', 120)
        text(x+12, y+34, name, 13, 'white', flow.NODE_WIDTH-24)
        text(x+16, y+76, translate('Saída'), 12, color)
        draw.ellipse(box(x+flow.NODE_WIDTH-35, y+78, 10, 10), fill=color)
        for port_index, port in enumerate(flow.input_ports(stage)):
            port_y = y + 117 + port_index * 34
            draw.ellipse(box(x+23, port_y-5, 10, 10),
                         fill=color if port['field'] in stage.get('bindings', {}) else 'white',
                         outline=color, width=max(1, round(scale)))
            text(x+45, port_y-8, translate(port['label']), 12, '#334155', flow.NODE_WIDTH-58)
        footer = 'Resultados fornecidos · concluído' if stage.get('provided_results') else 'Pronto para executar'
        text(x+16, y+height-30, translate(footer), 10, '#64748B', flow.NODE_WIDTH-32)

    output = io.BytesIO()
    image.save(output, format='PNG', dpi=(300, 300))
    return output.getvalue()
