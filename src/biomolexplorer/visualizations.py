"""Versioned, offline visualization artifacts shared by workers and the UI."""
import json
import math
from pathlib import Path
from .storage import write_json

SUFFIX = '.biomol-view.json'
MAX_VIEW_BYTES = 32 * 1024 * 1024


def compound_records(frame):
    # Pandas converts missing and numpy values to JSON-compatible scalar values.
    return json.loads(frame.to_json(orient='records', date_format='iso'))


def write_view(path, data):
    path = Path(path)
    write_json(data, path)
    return path


def separated_graph_layout(graph, circular=False):
    """Lay out each component with unit node spacing, then pack disjoint boxes.

    Coordinates have a physical spacing instead of being squeezed into a fixed
    square. The viewer can enlarge its canvas without merging adjacent nodes.
    """
    import networkx as nx
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), min(map(str, c))))
    groups = []
    for component in components:
        subgraph = nx.Graph()
        subgraph.add_nodes_from(sorted(component, key=str))
        subgraph.add_edges_from(sorted(graph.subgraph(component).edges(), key=lambda pair: tuple(map(str, pair))))
        count = len(subgraph)
        if count == 1:
            positions = {next(iter(subgraph)): (0.0, 0.0)}
        elif circular or count > 1000:
            # Concentric rings retain linear cost and readable spacing for large
            # graphs, where a force calculation would require a dense matrix.
            positions = {}
            remaining = iter(sorted(subgraph, key=lambda node: (-graph.degree(node), str(node))))
            positions[next(remaining)] = (0.0, 0.0)
            radius = 1.0
            placed = 1
            while placed < count:
                slots = min(count-placed, max(6, int(2*math.pi*radius)))
                for index in range(slots):
                    positions[next(remaining)] = (radius*math.cos(2*math.pi*index/slots),
                                                  radius*math.sin(2*math.pi*index/slots))
                placed += slots
                radius += 1.25
        else:
            raw = nx.spring_layout(subgraph, seed=42, iterations=80)
            scale = max(1.0, math.sqrt(count))
            positions = _separate_positions({node: (float(point[0])*scale, float(point[1])*scale)
                                              for node, point in raw.items()}, subgraph)
        xmin = min(point[0] for point in positions.values())
        ymin = min(point[1] for point in positions.values())
        width = max(point[0] for point in positions.values())-xmin+3
        height = max(point[1] for point in positions.values())-ymin+3
        groups.append((positions, xmin, ymin, width, height))
    row_limit = max([width for _, _, _, width, _ in groups] +
                    [math.sqrt(sum(width*height for _, _, _, width, height in groups)*1.6)])
    layout = {}
    x = y = row_height = 0.0
    for positions, xmin, ymin, width, height in groups:
        if x and x+width > row_limit:
            x = 0.0
            y += row_height
            row_height = 0.0
        for node, point in positions.items():
            layout[node] = (point[0]-xmin+x+1.5, point[1]-ymin+y+1.5)
        x += width
        row_height = max(row_height, height)
    return layout


def _separate_positions(positions, graph):
    """Resolve close force-layout points using an inexpensive spatial grid."""
    cells = {}
    result = {}
    for node in sorted(positions, key=lambda value: (-graph.degree(value), str(value))):
        ox, oy = positions[node]
        attempt = 0
        while True:
            radius = .55*math.sqrt(attempt)
            angle = attempt*math.pi*(3-math.sqrt(5))
            x, y = ox+radius*math.cos(angle), oy+radius*math.sin(angle)
            cx, cy = math.floor(x), math.floor(y)
            occupied = (point for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                        for point in cells.get((cx+dx, cy+dy), []))
            if all((x-px)**2+(y-py)**2 >= 1 for px, py in occupied):
                result[node] = (x, y)
                cells.setdefault((cx, cy), []).append((x, y))
                break
            attempt += 1
    return result


def graph_view(graph, dataset, title):
    import networkx as nx
    records = {str(row['molecule_chembl_id']): row for row in compound_records(dataset)}
    layout = separated_graph_layout(graph)
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), min(map(str, c))))
    largest = set(components[0]) if components else set()
    mcc = graph.subgraph(largest)
    mcc_layout = {node: layout[node] for node in mcc}
    component_index = {node: index for index, group in enumerate(components) for node in group}
    nodes = []
    for node in graph:
        info = dict(records.get(str(node), {}))
        info.update(molecule_chembl_id=str(node), degree=graph.degree(node),
                    degree_density=graph.degree(node)/max(1,len(graph)-1), component=component_index[node] + 1)
        x, y = layout[node]
        nodes.append({'id': str(node), 'x': float(x), 'y': float(y), 'properties': info})
    return {'version': 1, 'kind': 'graph', 'title': title, 'nodes': nodes,
            'edges': [{'source': str(a), 'target': str(b), 'value': float(attrs.get('value', 0))}
                      for a, b, attrs in graph.edges(data=True)],
            'mcc': sorted(map(str, largest)),
            'mcc_positions': {str(n):[float(x),float(y)] for n,(x,y) in mcc_layout.items()},
            'layout': 'components'}


def egg_view(dataset, title):
    nodes = []
    seen = {}
    for row in compound_records(dataset):
        x, y = row.get('TPSA'), row.get('WLOGP')
        if x is None or y is None or not math.isfinite(x) or not math.isfinite(y):
            continue
        identifier = str(row['molecule_chembl_id'])
        if identifier in seen:
            if seen[identifier] != row:
                raise ValueError('O mesmo código identifica compostos ou propriedades diferentes: '+identifier)
            continue
        seen[identifier] = row
        nodes.append({'id': identifier, 'x': x, 'y': y, 'properties': row})
    return {'version': 1, 'kind': 'egg', 'title': title, 'nodes': nodes, 'edges': [], 'mcc': [],
            'description': 'Descritores e classificações heurísticas calculados pela plataforma.'}


def load_view(content):
    """Validate uploaded/local artifacts before constructing native controls."""
    if len(content) > MAX_VIEW_BYTES:
        raise ValueError('Visualização maior que o limite de 32 MB.')
    data = json.loads(content)
    if not isinstance(data, dict) or type(data.get('version')) is not int or data['version'] != 1 or data.get('kind') not in ('graph', 'egg'):
        raise ValueError('Formato de visualização não reconhecido.')
    if not isinstance(data.get('title'), str):
        raise ValueError('Título da visualização inválido.')
    nodes, edges, mcc = data.get('nodes'), data.get('edges'), data.get('mcc')
    if not all(isinstance(items, list) for items in (nodes, edges, mcc)):
        raise ValueError('Dados da visualização inválidos.')
    ids = set()
    for node in nodes:
        if not isinstance(node, dict) or not isinstance(node.get('id'), str) or not node['id'] or node['id'] in ids:
            raise ValueError('Identificador molecular vazio ou duplicado na visualização.')
        ids.add(node['id'])
        if any(type(node.get(axis)) not in (int, float) or not math.isfinite(node[axis]) for axis in ('x', 'y')):
            raise ValueError('Coordenadas inválidas na visualização.')
        if not isinstance(node.get('properties'), dict):
            raise ValueError('Informações do composto inválidas.')
        if data['kind'] == 'graph':
            props = node['properties']
            if 'component' in props and (type(props['component']) is not int or props['component'] < 1):
                raise ValueError('Componente conectado inválido.')
            if 'degree' in props and (type(props['degree']) is not int or props['degree'] < 0):
                raise ValueError('Grau do nó inválido.')
    for edge in edges:
        if not isinstance(edge, dict) or any(not isinstance(edge.get(key), str) or edge[key] not in ids for key in ('source','target')):
            raise ValueError('Aresta sem composto correspondente.')
        if 'value' in edge and (type(edge['value']) not in (int,float) or not math.isfinite(edge['value'])):
            raise ValueError('Valor de similaridade inválido.')
    if any(not isinstance(node, str) or node not in ids for node in mcc):
        raise ValueError('MCC contém compostos inexistentes.')
    if len(set(mcc)) != len(mcc):
        raise ValueError('MCC contém compostos duplicados.')
    if data['kind'] == 'egg' and (edges or mcc):
        raise ValueError('O EGG contém pontos, sem arestas ou MCC.')
    if 'mcc_positions' in data:
        positions=data['mcc_positions']
        if not isinstance(positions,dict) or set(positions)!=set(mcc):
            raise ValueError('Layout do MCC inválido.')
        for point in positions.values():
            if not isinstance(point,list) or len(point)!=2 or any(type(v) not in (int,float) or not math.isfinite(v) for v in point):
                raise ValueError('Coordenadas do MCC inválidas.')
    if 'fragment' in data:
        fragment=data['fragment']
        if not isinstance(fragment,dict) or fragment.get('status') not in ('empty','missing_structures','partial','no_common_fragment','complete'):
            raise ValueError('Fragmento molecular inválido.')
        if any(not isinstance(fragment.get(k),str) for k in ('smarts','image')) or any(type(fragment.get(k)) is not int or fragment[k]<0 for k in ('atoms','bonds','compounds','timeout')):
            raise ValueError('Dados do fragmento molecular inválidos.')
        if 'smiles' in fragment and not isinstance(fragment['smiles'], str):
            raise ValueError('SMILES do fragmento molecular inválido.')
        if fragment['image']:
            import base64
            import io
            from PIL import Image
            try:
                raw=base64.b64decode(fragment['image'],validate=True)
                if len(raw)>2*1024*1024:raise ValueError('Imagem do fragmento muito grande.')
                with Image.open(io.BytesIO(raw)) as picture:
                    if picture.format!='PNG' or picture.width>2048 or picture.height>2048:raise ValueError('Imagem do fragmento inválida.')
                    picture.verify()
            except (ValueError,OSError) as exc:
                raise ValueError('Imagem do fragmento molecular inválida.') from exc
    return data


class PointIndex:
    """Spatial hit testing avoids scanning every compound for every hover event."""
    def __init__(self, points, cell_size=24):
        self.cell_size = cell_size
        self.cells = {}
        for identifier, (x, y) in points.items():
            key = (math.floor(x / cell_size), math.floor(y / cell_size))
            self.cells.setdefault(key, []).append((identifier, x, y))

    def nearest(self, x, y, radius=12):
        best, distance = None, radius * radius
        for cx in range(math.floor((x-radius)/self.cell_size), math.floor((x+radius)/self.cell_size)+1):
            for cy in range(math.floor((y-radius)/self.cell_size), math.floor((y+radius)/self.cell_size)+1):
                for identifier, px, py in self.cells.get((cx, cy), []):
                    squared = (px-x)**2 + (py-y)**2
                    if squared <= distance:
                        best, distance = identifier, squared
        return best

    def nearby(self, x, y, radius=1):
        found = []
        for cx in range(math.floor((x-radius)/self.cell_size), math.floor((x+radius)/self.cell_size)+1):
            for cy in range(math.floor((y-radius)/self.cell_size), math.floor((y+radius)/self.cell_size)+1):
                for identifier, px, py in self.cells.get((cx, cy), []):
                    if (px-x)**2+(py-y)**2 <= radius*radius:
                        found.append(identifier)
        return sorted(found)


def molecule_image(smiles):
    """Render a structure locally; never execute HTML or fetch external assets."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    from io import BytesIO
    from rdkit import Chem
    from rdkit.Chem import Draw
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    output = BytesIO()
    Draw.MolToImage(molecule, size=(320, 230)).save(output, format='PNG')
    return output.getvalue()


def _conformer_molecule(smiles):
    """Generate a reproducible local 3D conformer; this is not a docking pose."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule=Chem.MolFromSmiles(smiles) if isinstance(smiles,str) else None
    if molecule is None or molecule.GetNumAtoms()==0:
        raise ValueError('Não foi possível interpretar o SMILES deste composto.')
    molecule=Chem.AddHs(molecule)
    if molecule.GetNumAtoms()>512:
        raise ValueError('A visualização 3D está limitada a moléculas com até 512 átomos, incluindo hidrogênios.')
    if AllChem.EmbedMolecule(molecule,randomSeed=42,maxAttempts=20)!=0:
        raise ValueError('Não foi possível gerar um conformero 3D para este SMILES. A estrutura 2D continua disponível.')
    if AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule,maxIters=200)
    return molecule


def molecule_sdf(smiles):
    """Preserve coordinates, elements and bond orders for the shared 3D viewer."""
    from rdkit import Chem
    return Chem.MolToMolBlock(_conformer_molecule(smiles))+'\n$$$$\n'


def molecule_conformer(smiles):
    """Return coordinates and bonds for local molecular analyses."""
    molecule=_conformer_molecule(smiles)
    conformer=molecule.GetConformer()
    return {'atoms':[{'element':atom.GetSymbol(),'x':float(conformer.GetAtomPosition(atom.GetIdx()).x),
        'y':float(conformer.GetAtomPosition(atom.GetIdx()).y),'z':float(conformer.GetAtomPosition(atom.GetIdx()).z)} for atom in molecule.GetAtoms()],
        'bonds':[{'a':bond.GetBeginAtomIdx(),'b':bond.GetEndAtomIdx(),'order':float(bond.GetBondTypeAsDouble())} for bond in molecule.GetBonds()]}
