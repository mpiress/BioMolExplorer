"""Check local documentation links, HTML IDs and target anchors."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
class Page(HTMLParser):
    def __init__(self,path):
        super().__init__(convert_charrefs=True);self.path=path;self.ids=set();self.links=[]
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
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
print('\n'.join(errors) or f'{len(pages)} HTML pages: local files, anchors and IDs valid')
raise SystemExit(bool(errors))
