"""Verify real WebGL rendering and mouse interaction in a local headless browser.

Run with PYTHONPATH=src and the UI environment. No browser download or external
service is needed; use --chrome to select an installed Chromium/Chrome binary.
"""
import argparse,json,time,tempfile,subprocess,base64,shutil
from pathlib import Path
import requests
from PIL import Image
from websockets.sync.client import connect
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.pdb_view import StructureViewers
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--chrome',default=shutil.which('google-chrome') or shutil.which('chromium'))
parser.add_argument('--pdb','--structure',dest='pdb',type=Path,default=Path(__file__).resolve().parents[1]/'tests/fixtures/dms/1CRN.pdb')
parser.add_argument('--redocking-root',type=Path,help='Temporary validation root with input/<target> and output/<target>.')
parser.add_argument('--docking-root',type=Path,help='Native Vina or DOCK6 output directory with docking_results.csv.')
parser.add_argument('--smiles',help='Validate a generated compound in the shared viewer instead of a structure file.')
parser.add_argument('--screenshot',type=Path,default=Path('/tmp/biomol-pdb-modern.png'))
args=parser.parse_args()
if not args.chrome:parser.error('Select an installed browser with --chrome.')
with tempfile.TemporaryDirectory(prefix='biomol-webgl-') as tmp:
 store=WorkspaceStore(Path(tmp)/'workspace');token=store.register('Test','test@example.org','viewer-test-password')
 pid=store.create_project(token,'Browser verification')['id'];views=StructureViewers(store)
 if args.redocking_root or args.docking_root:
  from biomolexplorer.catalog import new_stage
  from biomolexplorer.redocking_results import RedockingResults
  from biomolexplorer.docking_results import DockingResults
  if args.redocking_root:
   target=next(p for p in (args.redocking_root/'input').iterdir() if p.is_dir())
   root=store.project_dir(pid)/'artifacts';shutil.copytree(target,root/'structures'/target.name)
   shutil.copytree(args.redocking_root/'output'/target.name,root/target.name)
   stage=new_stage('redocking');kind='redocking'
  else:
   root=store.project_dir(pid)/'artifacts';shutil.copytree(args.docking_root,root)
   table=root/'docking_results.csv'
   import csv
   with table.open() as stream:engine=next(csv.DictReader(stream))['engine']
   stage=new_stage('docking_'+engine);kind='docking'
  sid=stage['id'];rid='browser-scene'
  result=dict(id=sid,operation=stage['operation'],status='succeeded',configuration=stage,artifacts=[str(p) for p in root.rglob('*') if p.is_file()])
  with store.connect() as db:db.execute('INSERT INTO runs(id,project_id,user_id,status,stages,created,updated) VALUES(?,?,?,?,?,?,?)',
   (rid,pid,store.user(token)['id'],'succeeded',json.dumps([result]),1,1))
  if kind=='redocking':selection=RedockingResults(store).simulations(token,pid,rid,sid)[0]['id']
  else:selection=dict(table=str(table),index=0,version=DockingResults(store).page(token,pid,rid,sid,str(table))['version'],column='conformer_file')
  key=views.issue_result(token,pid,rid,sid,kind,selection,'en');name='docking.pdb'
 elif args.smiles:
  name='compound.sdf';key=views.issue_compound(token,pid,args.smiles,'compound')
 else:
  path=store.project_dir(pid)/args.pdb.name;path.write_bytes(args.pdb.read_bytes())
  name=path.name;key=views.issue(token,pid,str(path))
 desktop_url=views.url(key)
 # Flet exposes a WebSocket page URL: exercise its conversion in a real browser.
 page_url=desktop_url.split('/molecular-viewer')[0].replace('http://','ws://')+'/ws'
 url=views.url(key,web=True,page_url=page_url)
 assert url==desktop_url,'Flet WebSocket URL was not converted to HTTP'
 profile=Path(tmp)/'chrome';profile.mkdir()
 error=(Path(tmp)/'chrome.log').open('w')
 proc=subprocess.Popen([args.chrome,'--headless=new','--no-sandbox','--disable-dev-shm-usage','--no-first-run','--disable-background-networking','--use-angle=swiftshader','--enable-unsafe-swiftshader','--remote-debugging-port=0','--user-data-dir='+str(profile),'about:blank'],stdout=error,stderr=error)
 try:
  deadline=time.monotonic()+15
  while not (profile/'DevToolsActivePort').exists():
   if time.monotonic()>deadline:raise RuntimeError('Chrome startup failed')
   time.sleep(.1)
  port=(profile/'DevToolsActivePort').read_text().splitlines()[0]
  targets=requests.get('http://127.0.0.1:'+port+'/json',timeout=5).json();target=next(x for x in targets if x['type']=='page')
  with connect(target['webSocketDebuggerUrl'],max_size=None,open_timeout=5) as ws:
   count=0;events=[]
   def command(method,params={}):
    global count
    count+=1;ws.send(json.dumps({'id':count,'method':method,'params':params}))
    while True:
     message=json.loads(ws.recv())
     if message.get('id')==count:
      if 'error' in message:raise RuntimeError(message['error'])
      return message.get('result',{})
     events.append(message)
   def evaluate(expression):
    data=command('Runtime.evaluate',{'expression':expression,'returnByValue':True,'awaitPromise':True})
    if 'exceptionDetails' in data:raise RuntimeError(data['exceptionDetails'])
    return data.get('result',{}).get('value')
   command('Page.enable');command('Runtime.enable');command('Network.enable')
   command('Emulation.setDeviceMetricsOverride',{'width':1280,'height':800,'deviceScaleFactor':1,'mobile':False})
   command('Page.navigate',{'url':url})
   deadline=time.monotonic()+20
   while evaluate('document.body && document.body.dataset.ready')!='true':
    if evaluate('document.body && document.body.dataset.ready')=='error' or time.monotonic()>deadline:
     raise RuntimeError({'page':evaluate('document.body.innerText'), 'errors':[e for e in events if e.get('method') in ('Runtime.exceptionThrown','Runtime.consoleAPICalled')]})
    time.sleep(.1)
   first=evaluate('biomolViewer.getView()')
   command('Input.dispatchMouseEvent',{'type':'mousePressed','x':800,'y':400,'button':'left','clickCount':1})
   command('Input.dispatchMouseEvent',{'type':'mouseMoved','x':945,'y':460,'button':'left','buttons':1})
   command('Input.dispatchMouseEvent',{'type':'mouseReleased','x':945,'y':460,'button':'left','clickCount':1})
   time.sleep(.2);rotated=evaluate('biomolViewer.getView()');assert first!=rotated,'Rotation failed'
   command('Input.dispatchMouseEvent',{'type':'mouseWheel','x':800,'y':400,'deltaX':0,'deltaY':-240})
   time.sleep(.2);zoomed=evaluate('biomolViewer.getView()');assert rotated!=zoomed,'Zoom failed'
   representations=('sticks','spheres','lines') if args.smiles else ('sticks','spheres','lines','cartoon')
   for value in representations:
    control_id='ligand-representation' if args.smiles else 'representation'
    evaluate(f'document.getElementById("{control_id}").value="{value}";document.getElementById("{control_id}").dispatchEvent(new Event("change"));true')
    if args.smiles:
     styles=evaluate('biomolViewer.getModel().selectedAtoms({}).filter(atom => !["H","D"].includes(atom.elem)).map(atom => Object.keys(atom.style))')
     expected={'sticks':'stick','spheres':'sphere','lines':'line'}[value]
     assert all(expected in style for style in styles),'Compound representation did not change'
   if args.smiles:
    evaluate('document.getElementById("ligand-representation").value="sticks";document.getElementById("ligand-representation").dispatchEvent(new Event("change"));true')
    assert evaluate('!document.getElementById("conformer").hidden'),'Conformer provenance missing'
   evaluate('document.getElementById("center").click();true')
   downloads=Path(tmp)/'downloads';downloads.mkdir()
   command('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(downloads)})
   evaluate('document.getElementById("snapshot").click();true')
   deadline=time.monotonic()+5
   snapshot=downloads/(Path(name).stem+'.png')
   while not snapshot.exists():
    if time.monotonic()>deadline:raise RuntimeError('PNG export failed')
    time.sleep(.1)
   assert snapshot.read_bytes().startswith(b'\x89PNG'),'Invalid PNG export'
   if not args.smiles and not args.redocking_root and not args.docking_root and args.pdb.suffix.lower()=='.pdbqt':
    first_model=[]
    for line in args.pdb.read_text().splitlines():
     if line.startswith('ENDMDL'):break
     if line.startswith(('ATOM  ','HETATM')):first_model.append(line)
    assert evaluate('biomolViewer.getModel().selectedAtoms({}).length')==len(first_model),'PDBQT atoms lost at torsion boundaries'
   # Save an actual rendered frame for visual inspection.
   time.sleep(.4)
   shot=command('Page.captureScreenshot',{'format':'png'})
   args.screenshot.write_bytes(base64.b64decode(shot['data']))
   # A ready WebGL canvas can still contain no visible geometry (for example,
   # when a ligand is incorrectly assigned the polymer cartoon style).
   image=Image.open(args.screenshot).convert('RGB')
   pixels=image.crop((420,120,1200,650)).tobytes()
   colored=sum(max(pixel)>45 for pixel in zip(pixels[0::3],pixels[1::3],pixels[2::3]))
   assert colored>100,'No molecular geometry visible in the rendered scene'
   requests_seen=[event['params']['request']['url'] for event in events if event.get('method')=='Network.requestWillBeSent']
   external=[request for request in requests_seen if request.startswith(('http://','https://')) and not request.startswith(url.split('/molecular-viewer')[0])]
   assert not external,external
   errors=[event['params'] for event in events if event.get('method')=='Runtime.exceptionThrown']
   assert not errors,errors
   if args.redocking_root or args.docking_root:
    if evaluate('biomolDockingScene.interactions.length'):
     assert evaluate('biomolDockingScene.layers.find(l=>l.role==="receptor").model.selectedAtoms({}).some(a=>a.style.stick)'), 'Interacting residues not visible in cartoon mode'
    assert evaluate('biomolDockingScene.layers.length')>=2,'Receptor/pose overlay missing'
    assert evaluate('biomolDockingScene.layers.every(layer => layer.model.selectedAtoms({}).length>0)'),'Empty layer'
    evaluate('document.getElementById("representation").value="lines";document.getElementById("representation").dispatchEvent(new Event("change"));document.getElementById("ligand-representation").value="spheres";document.getElementById("ligand-representation").dispatchEvent(new Event("change"));true')
    assert evaluate('biomolDockingScene.layers.find(l=>l.role==="pose").model.selectedAtoms({}).filter(a=>!["H","D"].includes(a.elem)).every(a=>a.style.sphere && !a.style.line)'), 'Ligand style is not independent'
    assert evaluate('biomolDockingScene.layers.find(l=>l.role==="receptor").model.selectedAtoms({}).filter(a=>a.resn==="PHE").every(a=>a.style.line && !a.style.sphere)'), 'Receptor style changed with ligand'
    assert evaluate('biomolDockingScene.layers.find(l=>l.role==="pose").model.selectedAtoms({elem:"H"}).every(a=>Object.keys(a.style).length===0)'), 'Ligand H not hidden'
    evaluate('document.getElementById("hydrogens").click();true')
    assert evaluate('biomolDockingScene.layers.find(l=>l.role==="pose").model.selectedAtoms({elem:"H"}).every(a=>a.style.sphere)'), 'Ligand H not shown'
    evaluate('document.getElementById("hydrogens").click();true')
    if evaluate('biomolDockingScene.interactions.length'):
     assert evaluate('biomolDockingScene.shapes.length===biomolDockingScene.interactions.length'), '3D interaction lines missing'
     assert evaluate('biomolDockingScene.shapes.every(s=>s.hoverable && s.intersectionShape.cylinder.length>1)'), 'Dashed interactive geometry missing'
     kinds=evaluate('[...new Set(biomolDockingScene.interactions.map(r=>r.kind))]')
     for kind in kinds:
      # Use real projected dash centres: callbacks alone cannot verify picking.
      candidates=evaluate(f'biomolDockingScene.shapes.filter(s=>s.interaction.kind==={json.dumps(kind)}).flatMap(s=>s.intersectionShape.cylinder.map(c=>biomolViewer.modelToScreen({{x:(c.c1.x+c.c2.x)/2,y:(c.c1.y+c.c2.y)/2,z:(c.c1.z+c.c2.z)/2}})))')
      found=False
      for point in candidates:
       command('Input.dispatchMouseEvent',{'type':'mouseMoved','x':point['x'],'y':point['y']})
       time.sleep(.2)
       if evaluate(f'biomolDockingScene.hoveredInteraction?.kind==={json.dumps(kind)} && !document.getElementById("interaction-tooltip").hidden'):
        found=True;break
      assert found, 'Mouse hover did not identify '+kind
      assert evaluate('document.getElementById("interaction-tooltip").textContent.includes("Å")'), 'Hover distance missing'
      command('Input.dispatchMouseEvent',{'type':'mouseMoved','x':1200,'y':100});time.sleep(.2)
      assert evaluate('document.getElementById("interaction-tooltip").hidden'), 'Tooltip did not clear on mouse exit'
     evaluate('document.getElementById("show-interactions").click();true')
     assert evaluate('document.getElementById("interaction-tooltip").hidden'), 'Tooltip survived hidden interactions'
     assert evaluate('biomolDockingScene.shapes.length===0'), 'Interaction visibility failed'
     evaluate('document.getElementById("show-interactions").click();true')
     evaluate('document.querySelector("#interaction-types input").click();true')
     assert evaluate('biomolDockingScene.shapes.length===biomolDockingScene.interactions.filter(r=>biomolDockingScene.enabledTypes.has(r.kind)).length'), 'Interaction type filter failed'
     evaluate('document.querySelector("#interaction-types input").click();true')
    before=evaluate('document.querySelectorAll(".contact-residue").length')
    evaluate('document.getElementById("cutoff").value="2";document.getElementById("cutoff").dispatchEvent(new Event("change"));true')
    assert evaluate('document.querySelectorAll(".contact-residue").length')<=before,'Contact cutoff did not filter residues'
    evaluate('document.getElementById("cutoff").value="4";document.getElementById("cutoff").dispatchEvent(new Event("change"));true')
    if before:
     evaluate('document.querySelector(".contact-residue").click();true')
     assert evaluate('biomolDockingScene.selected!==null'),'Residue highlighting missing'
    if args.redocking_root:
     assert evaluate('biomolDockingScene.layers.some(layer => layer.role==="reference")'),'Crystal reference missing'
    evaluate('document.querySelector(".scene-layer.pose input").click();true')
    assert not evaluate('biomolDockingScene.layers.find(layer => layer.role==="pose").visible'),'Pose visibility toggle failed'
    evaluate('document.querySelector(".scene-layer.pose input").click();true')
   errors=[event['params'] for event in events if event.get('method')=='Runtime.exceptionThrown' or (event.get('method')=='Runtime.consoleAPICalled' and event['params']['type']=='error')]
   assert not errors,errors
   print(json.dumps({'websocket_url_converted':True,'rendered':True,'atoms':evaluate('biomolViewer.getModel().selectedAtoms({}).length'),'mouse_rotation':first!=rotated,'wheel_zoom':rotated!=zoomed,'representations':len(representations),'png_export':True,'external_requests':len(external),'javascript_errors':len(errors),'screenshot':str(args.screenshot), 'interaction_hover':bool(args.docking_root or args.redocking_root) and bool(evaluate('window.biomolDockingScene?.interactions.length')), 'interactions_3d':evaluate('window.biomolDockingScene ? biomolDockingScene.interactions.length : 0'), 'layers':evaluate('window.biomolDockingScene ? biomolDockingScene.layers.map(l => ({role:l.role,atoms:l.model.selectedAtoms({}).length})) : []'), 'residue_contacts':evaluate('document.querySelectorAll(".contact-residue").length')}))
 finally:
  proc.terminate();proc.wait(timeout=15);error.close();views.close()
