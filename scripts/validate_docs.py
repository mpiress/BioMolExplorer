"""Check local documentation links, HTML IDs and target anchors."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import ast
import re
class Page(HTMLParser):
    def __init__(self,path):
        super().__init__(convert_charrefs=True);self.path=path;self.ids=set();self.links=[];self.operations=[]
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if 'data-operation' in attrs:self.operations.append(attrs['data-operation'])
        if 'id' in attrs:
            if attrs['id'] in self.ids:raise ValueError(f'Duplicate id {self.path}: {attrs["id"]}')
            self.ids.add(attrs['id'])
        for attr in ('href','src'):
            if attr in attrs:self.links.append(attrs[attr])
root=Path(__file__).resolve().parents[1];pages={}
for p in [root/'index.html',root/'index_doc.html',*(root/'docs').rglob('*.html')]:
    page=Page(p);page.feed(p.read_text());pages[p.resolve()]=page
errors=[]
for p,page in pages.items():
    for link in page.links:
        u=urlsplit(link)
        if u.scheme or u.netloc:continue
        target=(p.parent/unquote(u.path)).resolve() if u.path else p
        if not target.exists():errors.append(f'{p.relative_to(root)}: missing {link}')
        elif u.fragment and target in pages and unquote(u.fragment) not in pages[target].ids:
            errors.append(f'{p.relative_to(root)}: unknown anchor {link}')
# Check documentation against source without importing scientific dependencies.
registry=ast.parse((root/'src/biomolexplorer/operations.py').read_text())
operations=next(n.value for n in registry.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='OPERATIONS' for t in n.targets))
active={ast.literal_eval(key) for key in operations.keys}
for key,spec in zip(operations.keys,operations.values):
    module,function=(ast.literal_eval(v) for v in spec.args[:2])
    source=root/'src'/Path(*module.split('.')).with_suffix('.py')
    if not source.is_file():errors.append(f'{key.value}: missing implementation {module}')
    elif not any(isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==function for n in ast.parse(source.read_text()).body):
        errors.append(f'{key.value}: missing entry point {module}.{function}')
expected=active|{'import_results'}
features=pages[(root/'index.html').resolve()].operations
if set(features)!=expected or len(features)!=len(expected):
    errors.append('index.html: feature cards must match each active operation and the import block exactly once')
for path in (root/'docs/user_manual.md',root/'docs/en/user_manual.md'):
    examples=re.findall(r'^### [^\n]+\(`([a-z_0-9]+)`\)\n\n```text\n(.*?)\n```',path.read_text(),re.M|re.S)
    if {op for op,_ in examples}!=expected or len(examples)!=len(expected):
        errors.append(f'{path.relative_to(root)}: missing or duplicate operation examples')
    for op,example in examples:
        if not 3<=len(example.splitlines())<=5:errors.append(f'{path.relative_to(root)}: {op} example must have 3 to 5 lines')
for path in [root/'README.md',root/'docs/README.md',*(root/'docs').glob('*.md'),*(root/'docs/en').glob('*.md')]:
    for link in re.findall(r'\]\(([^\s)]+)\)',path.read_text()):
        u=urlsplit(link)
        if u.scheme or u.netloc or not u.path:continue
        if not (path.parent/unquote(u.path)).exists():errors.append(f'{path.relative_to(root)}: missing {link}')
print('\n'.join(errors) or f'{len(pages)} HTML pages: local files, anchors and IDs valid; {len(active)} scientific operations and import documented in both manuals and index')
raise SystemExit(bool(errors))
