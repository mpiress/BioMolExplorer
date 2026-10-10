/* Shared local docking scene; every layer keeps its original coordinates. */
window.BioMolDockingScene = async function (config) {
  'use strict';
  const el = id => document.getElementById(id), en = config.language === 'en';
  const t = (pt, english) => en ? english : pt;
  let viewer;
  const state = {layers: [], cutoff: 4, selected: null};
  const protein = new Set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP CYX SEC PYL'.split(' '));
  const water = atom => ['HOH', 'WAT', 'DOD'].includes(atom.resn);
  function pdbqt(data) {
    const elements = {A:'C', OA:'O', NA:'N', SA:'S', HD:'H', HS:'H'};
    return data.split(/\r?\n/).filter(line => /^(ATOM  |HETATM|CONECT|TER   )/.test(line)).map(line => {
      if (!/^(ATOM  |HETATM)/.test(line)) return line;
      const type = line.trim().split(/\s+/).pop();
      return line.slice(0,66).padEnd(76) + (elements[type] || type).padStart(2);
    }).join('\n');
  }
  function residueKey(row) {return `${row.chain}|${row.resi}|${row.icode || ''}|${row.resn}`;}
  function draw() {
    viewer.removeAllLabels();
    const receptor = state.layers.find(layer => layer.role === 'receptor');
    const chain = el('chain').value;
    const selected = chain ? {chain:chain === '__blank__' ? '' : chain} : {};
    state.layers.forEach(layer => {
      if (!layer.visible) {layer.model.hide(); return;}
      layer.model.show(); layer.model.setStyle({}, {});
      if (!['receptor','complex'].includes(layer.role)) {
        const color = layer.role === 'pose' ? '#22d3ee' : '#fbbf24';
        const rep=el('representation').value;
        const style=rep==='spheres' ? {sphere:{color,scale:.65}} : rep==='lines' ? {line:{color}} : {stick:{color,radius:.22},sphere:{color,scale:.2}};
        if (el('ligands').checked) layer.model.setStyle({},style);
      } else {
        const coloring = el('color').value === 'element' ? {colorscheme:'Jmol'} :
          el('color').value === 'spectrum' ? {color:'spectrum'} : {color:'#94a3b8'};
        const representation = el('representation').value;
        const style = representation === 'cartoon' && layer.format !== 'mol2' ? {cartoon:{...coloring,arrows:true}} :
          representation === 'spheres' ? {sphere:{...coloring,scale:.6}} :
          representation === 'lines' ? {line:{...coloring}} : {stick:{...coloring,radius:.14}};
        layer.model.setStyle({...selected,predicate:atom => protein.has(atom.resn)},style);
        if (el('color').value==='chain') {
          const palette=['#94a3b8','#a5b4fc','#60a5fa','#c084fc'];
          [...new Set(layer.model.selectedAtoms({}).map(a => a.chain || ''))].sort().forEach((value,index) => {
            if (chain && value !== (chain==='__blank__' ? '' : chain)) return;
            const tinted=Object.fromEntries(Object.entries(style).map(([key,val]) => [key,{...val,color:palette[index%palette.length]}]));
            layer.model.setStyle({chain:value,predicate:atom => protein.has(atom.resn)},tinted);
          });
        }
        if (el('ligands').checked) layer.model.setStyle({...selected,predicate:atom => !protein.has(atom.resn) && !water(atom)}, {stick:{colorscheme:'Jmol',radius:.15}});
        if (el('water').checked) layer.model.setStyle({...selected,predicate:water},{sphere:{color:'#60a5fa',scale:.18}});
        const keys = new Set((layer.role==='complex' ? state.referenceContacts : state.contacts).filter(row => row.distance <= state.cutoff).map(residueKey));
        layer.model.setStyle({...selected,predicate:atom => keys.has(residueKey(atom))},{stick:{color:'#86efac',radius:.15}});
      }
      layer.model.setClickable({},true,atom => {
        viewer.removeAllLabels();
        const label = `${atom.resn} ${atom.resi}${atom.icode || ''} · ${atom.chain || '—'} · ${atom.atom} (${atom.elem})`;
        el('atom').textContent = label;
        viewer.addLabel(label,{position:atom,backgroundColor:'#111b2b',fontColor:'#e5edf9',fontSize:12});
        viewer.render();
      });
    });
    if (state.selected && receptor.visible) {
      const row = state.selected, selection = {chain:row.chain,resi:row.resi,icode:row.icode || ''};
      receptor.model.addStyle(selection,{stick:{color:'#fb7185',radius:.24}});
      const atom = receptor.model.selectedAtoms(selection).find(a => a.atom === 'CA') || receptor.model.selectedAtoms(selection)[0];
      if (atom) viewer.addLabel(`${row.resn} ${row.resi}${row.icode || ''} · ${row.chain || '—'}`,{position:atom,backgroundColor:'#111b2b',fontColor:'#ffffff',fontSize:12});
    }
    viewer.render();
  }
  function contactTable() {
    const rows = new Map();
    state.contacts.filter(r => r.distance <= state.cutoff).forEach(r => rows.set(residueKey(r),{...r,pose:r.distance}));
    state.referenceContacts.filter(r => r.distance <= state.cutoff).forEach(r => {
      const key = residueKey(r); rows.set(key,{...(rows.get(key) || r),reference:r.distance});
    });
    const container = el('contacts');container.replaceChildren();
    const heading = document.createElement('strong');heading.textContent = t('Resíduos próximos','Nearby residues');container.append(heading);
    const header = document.createElement('p');header.textContent = state.referenceContacts.length ? t('Resíduo · Pose / Referência (Å)','Residue · Pose / Reference (Å)') : t('Resíduo · Distância mínima (Å)','Residue · Minimum distance (Å)');container.append(header);
    if (!rows.size) {const empty=document.createElement('p');empty.textContent=t('Nenhum contato nesta distância.','No contacts at this distance.');container.append(empty);}
    [...rows.values()].sort((a,b) => (a.pose ?? a.reference)-(b.pose ?? b.reference)).forEach(row => {
      const button=document.createElement('button');button.className='contact-residue';
      button.textContent=`${row.resn} ${row.resi}${row.icode || ''} · ${row.chain || '—'} · ${row.pose?.toFixed(3) || '—'}${state.referenceContacts.length ? ' / '+(row.reference?.toFixed(3) || '—') : ''}`;
      button.addEventListener('click',() => {
        state.selected=row;draw();
        const receptor=state.layers.find(layer => layer.role==='receptor');
        viewer.zoomTo({model:receptor.model.getID(),chain:row.chain,resi:row.resi,icode:row.icode || ''});viewer.render();
      });container.append(button);
    });
  }
  try {
    viewer=$3Dmol.createViewer(el('scene'),{backgroundColor:'#0b1220',antialias:true});viewer.setProjection('perspective');
    const response=await fetch(location.pathname.replace(/\/$/,'')+'/scene',{cache:'no-store',credentials:'same-origin'});
    if (!response.ok) throw new Error('Scene unavailable');
    const data=await response.json();state.contacts=data.contacts;state.referenceContacts=data.reference_contacts;
    data.layers.forEach(layer => {
      let format=layer.format,structure=layer.data;
      if (format==='pdbqt') {structure=pdbqt(structure);format='pdb';}
      const model=viewer.addModel(structure,format,{keepH:true});
      if (!model.selectedAtoms({}).length) throw new Error('Empty layer');
      state.layers.push({...layer,data:undefined,model,visible:layer.role!=='complex'});
      const label=document.createElement('label');label.className='scene-layer '+layer.role;
      const check=document.createElement('input');check.type='checkbox';check.checked=layer.role!=='complex';
      const name=layer.role==='receptor' ? t('Receptor utilizado','Docking receptor') : layer.role==='complex' ? t('Complexo cristalográfico','Crystallographic complex') : layer.role==='reference' ?
        data.reference_kind==='crystal' ? t('Ligante cristalográfico','Crystallographic ligand') : t('Referência preparada','Prepared reference') :
        layer.selection==='lowest_score' ? t('Melhor pose','Best pose') : t('Primeira pose sem score','First unscored pose');
      label.append(check,document.createTextNode(name+(layer.role==='pose' ? ` · ${layer.score ?? '—'} · ${t('modelo','model')} ${layer.model || 1}` : '')));
      check.addEventListener('change',() => {state.layers.find(l => l.role===layer.role).visible=check.checked;draw();});el('layers').append(label);
    });
    const receptor=state.layers.find(layer => layer.role==='receptor');
    [...new Set(receptor.model.selectedAtoms({}).map(a => a.chain || ''))].sort().forEach(chain => el('chain').add(new Option(chain || '—',chain || '__blank__')));
    el('model').parentElement.querySelector('label[for="model"]').remove();el('model').hidden=true;
    document.body.dataset.scene='true';el('docking-controls').hidden=false;el('cutoff-label').textContent=t('Distância de contato (Å)','Contact distance (Å)');
    el('contact-note').textContent=t('Proximidade entre átomos pesados; não identifica tipos de ligação. Coordenadas originais, sem realinhamento.','Heavy-atom proximity; does not identify bond types. Original coordinates, without realignment.');
    el('summary').textContent=data.name;
    el('representation').value=receptor.format==='mol2' ? 'sticks' : 'cartoon';
    contactTable();draw();viewer.zoomTo({model:state.layers.find(layer => layer.role==='pose').model.getID()});viewer.render();
    for (const id of ['representation','color','ligands','water','chain']) el(id).addEventListener('change',draw);
    el('cutoff').addEventListener('change',() => {const value=Number(el('cutoff').value);if (!Number.isFinite(value) || value<2 || value>8) {el('cutoff').value=state.cutoff;return;}state.cutoff=value;contactTable();draw();});
    el('center').addEventListener('click',() => {viewer.zoomTo();viewer.render();});
    el('snapshot').addEventListener('click',() => {const link=document.createElement('a');link.href=viewer.pngURI();link.download='docking.png';link.click();});
    el('fullscreen').addEventListener('click',async () => {if (document.fullscreenElement) await document.exitFullscreen();else await document.documentElement.requestFullscreen();});
    new ResizeObserver(() => {viewer.resize();viewer.render();}).observe(el('scene'));
    el('status').hidden=true;window.biomolViewer=viewer;window.biomolDockingScene=state;document.body.dataset.ready='true';
  } catch (error) {el('status').textContent=config.labels.error;document.body.dataset.ready='error';}
};
