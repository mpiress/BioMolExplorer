/* Shared local docking scene; every layer keeps its original coordinates. */
window.BioMolDockingScene = async function (config) {
  'use strict';
  const el = id => document.getElementById(id), en = config.language === 'en';
  const t = (pt, english) => en ? english : pt;
  let viewer;
  const state = {layers: [], cutoff: 4, selected: null, interactions: [], enabledTypes: new Set(), shapes: [], hoveredInteraction: null};
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
  function residueKey(row) {return `${(row.chain || '').trim()}|${row.resi}|${(row.icode || '').trim()}|${row.resn}`;}
  const residueSelection = row => ({chain:row.chain,resi:row.resi,predicate:atom => residueKey(atom)===residueKey(row)});
  const hydrogenVisible = atom => el('hydrogens').checked || !['H','D'].includes(atom.elem.toUpperCase());
  function ligandStyle(coloring) {
    const rep=el('ligand-representation').value;
    return rep==='spheres' ? {sphere:{...coloring,scale:.65}} : rep==='lines' ? {line:{...coloring}} : {stick:{...coloring,radius:.22},sphere:{...coloring,scale:.2}};
  }
  const interactionNames = {
    pi_parallel:t('π–π paralelo','Parallel π–π'), pi_t:t('π–π em T','T-shaped π–π'),
    hydrogen_bond:t('Ligação de hidrogênio','Hydrogen bond'), hydrophobic:t('Contato hidrofóbico','Hydrophobic contact'),
    van_der_waals:t('Contato van der Waals (geométrico)','Van der Waals contact (geometric)')
  };
  const position = xyz => ({x:xyz[0],y:xyz[1],z:xyz[2]});
  function interactionText(row) {
    return `${interactionNames[row.kind]} · ${row.resn} ${row.resi}${row.icode || ''} · ${row.chain || '—'} · ${row.distance.toFixed(2)} Å`;
  }
  function clearInteractionHover() {
    state.hoveredInteraction=null;el('interaction-tooltip').hidden=true;
  }
  function showInteractionHover(row,event) {
    state.hoveredInteraction=row;
    const tooltip=el('interaction-tooltip');tooltip.textContent=interactionText(row);
    tooltip.style.borderColor=row.color;tooltip.hidden=false;
    const box=tooltip.getBoundingClientRect();
    const x=event?.clientX ?? event?.pageX ?? window.innerWidth/2;
    const y=event?.clientY ?? event?.pageY ?? window.innerHeight/2;
    tooltip.style.left=`${Math.max(8,Math.min(x+14,window.innerWidth-box.width-8))}px`;
    tooltip.style.top=`${Math.max(8,Math.min(y+14,window.innerHeight-box.height-8))}px`;
  }
  function interactionPath(row,index,count) {
    const start=position(row.ligand_position),end=position(row.receptor_position);
    if (count===1) return [start,end];
    // Separate coincident relations while keeping both chemical endpoints fixed.
    const delta=[end.x-start.x,end.y-start.y,end.z-start.z];
    const axis=Math.abs(delta[0])<Math.abs(delta[1]) ? [1,0,0] : [0,1,0];
    const normal=[delta[1]*axis[2]-delta[2]*axis[1],delta[2]*axis[0]-delta[0]*axis[2],delta[0]*axis[1]-delta[1]*axis[0]];
    const length=Math.hypot(...normal) || 1,offset=(index-(count-1)/2)*.3;
    const middle={x:(start.x+end.x)/2+offset*normal[0]/length,y:(start.y+end.y)/2+offset*normal[1]/length,z:(start.z+end.z)/2+offset*normal[2]/length};
    return [start,middle,end];
  }
  function drawInteractions() {
    viewer.setHover(null);clearInteractionHover();
    state.shapes.forEach(shape => viewer.removeShape(shape));state.shapes=[];
    const receptor=state.layers.find(l => l.role==='receptor'),pose=state.layers.find(l => l.role==='pose');
    const container=el('interaction-list');container.replaceChildren();
    if (!el('show-interactions').checked || !el('ligands').checked || !receptor.visible || !pose.visible) return;
    const chain=el('chain').value;
    const labels=new Map();
    const visible=state.interactions.filter(row => state.enabledTypes.has(row.kind) && (!chain || row.chain===(chain==='__blank__' ? '' : chain)));
    const paths=new Map();
    const pathKey=row => [...row.ligand_position,...row.receptor_position].map(v => v.toFixed(3)).join('|');
    visible.forEach(row => {const key=pathKey(row);if(!paths.has(key))paths.set(key,[]);paths.get(key).push(row);});
    visible.forEach(row => {
      const group=paths.get(pathKey(row)),points=interactionPath(row,group.indexOf(row),group.length),end=points[points.length-1];
      const shape=viewer.addShape({color:row.color,hoverable:true,
        hover_callback:(_shape,_viewer,event) => showInteractionHover(row,event),
        unhover_callback:clearInteractionHover});
      for (let i=1;i<points.length;i++) shape.addDashedCylinder({start:points[i-1],end:points[i],radius:.055,color:row.color,dashLength:.22,gapLength:.12,fromCap:1,toCap:1});
      shape.finalize();shape.interaction=row;state.shapes.push(shape);
      const key=residueKey(row);
      if (!labels.has(key)) labels.set(key,{position:end,rows:[]});labels.get(key).rows.push(row);
      if (el('representation').value==='cartoon') receptor.model.setStyle(residueSelection(row),{stick:{colorscheme:'Jmol',radius:.14}},true);
      const button=document.createElement('button');button.className='interaction-row';button.textContent=interactionText(row);button.style.borderLeftColor=row.color;
      button.addEventListener('click',() => {state.selected=row;draw();viewer.zoomTo({model:receptor.model.getID(),...residueSelection(row)});viewer.render();});container.append(button);
    });
    const shortNames={pi_parallel:'π–π',pi_t:'π–π (T)',hydrogen_bond:t('lig. H','H-bond'),hydrophobic:t('hidrofóbico','hydrophobic'),van_der_waals:'vdW'};
    labels.forEach(({position,rows}) => viewer.addLabel(`${rows[0].resn} ${rows[0].resi}${rows[0].icode || ''} · ${rows[0].chain || '—'} · ${rows.map(row => shortNames[row.kind]).join(' / ')}`,
      {position,backgroundColor:'#111b2b',fontColor:'#e5edf9',fontSize:10,showBackground:true}));
  }
  function setupInteractions(capability) {
    el('interaction-controls').hidden=false;
    el('interaction-toggle-label').textContent=t('Mostrar interações em 3D','Show 3D interactions');
    el('show-interactions').disabled=!capability?.available;
    el('interaction-status').textContent=!capability?.available ?
      t('Classificação indisponível: faltam topologia química ou identificação dos resíduos. Os contatos por distância continuam disponíveis.','Classification unavailable: chemical topology or residue identities are missing. Distance contacts remain available.') :
      `${capability.interactions.length} ${t('interações geométricas. Ligações de hidrogênio exigem H explícitos; ausência não comprova inexistência. van der Waals indica proximidade pelos raios atômicos, sem cálculo de energia.','geometric interactions. Hydrogen bonds require explicit H; absence does not prove no bonds exist. van der Waals indicates proximity by atomic radii, without energy calculation.')}`;
    if (!capability?.available) return;
    state.interactions=capability.interactions;
    capability.types.forEach(kind => {
      state.enabledTypes.add(kind);
      const label=document.createElement('label'),check=document.createElement('input');check.type='checkbox';check.checked=true;check.dataset.kind=kind;
      const sample=state.interactions.find(row => row.kind===kind),color=capability.colors?.[kind] || sample?.color || '#94a3b8';
      const swatch=document.createElement('span');swatch.className='interaction-swatch';swatch.style.borderColor=color;swatch.setAttribute('aria-hidden','true');
      label.append(check,swatch,document.createTextNode(interactionNames[kind]));
      check.addEventListener('change',() => {if(check.checked)state.enabledTypes.add(kind);else state.enabledTypes.delete(kind);draw();});el('interaction-types').append(label);
    });
  }
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
        if (el('ligands').checked) layer.model.setStyle({predicate:hydrogenVisible},ligandStyle({color}));
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
        if (el('ligands').checked) layer.model.setStyle({...selected,predicate:atom => !protein.has(atom.resn) && !water(atom) && hydrogenVisible(atom)}, ligandStyle({colorscheme:'Jmol'}));
        if (el('water').checked) layer.model.setStyle({...selected,predicate:water},{sphere:{color:'#60a5fa',scale:.18}});
        const keys = new Set((layer.role==='complex' ? state.referenceContacts : state.contacts).filter(row => row.distance <= state.cutoff).map(residueKey));
        if (representation==='cartoon') layer.model.setStyle({...selected,predicate:atom => keys.has(residueKey(atom))},{stick:{color:'#86efac',radius:.15}},true);
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
      const row = state.selected, selection = residueSelection(row);
      const rep=el('representation').value;
      receptor.model.setStyle(selection,rep==='spheres' ? {sphere:{color:'#fb7185',scale:.6}} : rep==='lines' ? {line:{color:'#fb7185'}} : {stick:{color:'#fb7185',radius:.24}});
      const atom = receptor.model.selectedAtoms(selection).find(a => a.atom === 'CA') || receptor.model.selectedAtoms(selection)[0];
      if (atom) viewer.addLabel(`${row.resn} ${row.resi}${row.icode || ''} · ${row.chain || '—'}`,{position:atom,backgroundColor:'#111b2b',fontColor:'#ffffff',fontSize:12});
    }
    drawInteractions();
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
        viewer.zoomTo({model:receptor.model.getID(),...residueSelection(row)});viewer.render();
      });container.append(button);
    });
  }
  try {
    viewer=$3Dmol.createViewer(el('scene'),{backgroundColor:'#0b1220',antialias:true});viewer.setProjection('perspective');viewer.setHoverDuration(120);
    const response=await fetch(location.pathname.replace(/\/$/,'')+'/scene',{cache:'no-store',credentials:'same-origin'});
    if (!response.ok) throw new Error('Scene unavailable');
    const data=await response.json();state.contacts=data.contacts;state.referenceContacts=data.reference_contacts;
    data.layers.forEach(layer => {
      let format=layer.format,structure=layer.data;
      if (format==='pdbqt') {structure=pdbqt(structure);format='pdb';}
      const model=viewer.addModel(structure,format,{keepH:true});
      if (!model.selectedAtoms({}).length) throw new Error('Empty layer');
      state.layers.push({...layer,poseModel:layer.model,data:undefined,model,visible:layer.role!=='complex'});
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
    setupInteractions(data.interactions);contactTable();draw();viewer.zoomTo({model:state.layers.find(layer => layer.role==='pose').model.getID()});viewer.render();
    for (const id of ['representation','ligand-representation','hydrogens','show-interactions','color','ligands','water','chain']) el(id).addEventListener('change',draw);
    el('cutoff').addEventListener('change',() => {const value=Number(el('cutoff').value);if (!Number.isFinite(value) || value<2 || value>8) {el('cutoff').value=state.cutoff;return;}state.cutoff=value;contactTable();draw();});
    el('center').addEventListener('click',() => {viewer.zoomTo();viewer.render();});
    el('snapshot').addEventListener('click',() => {const link=document.createElement('a');link.href=viewer.pngURI();link.download='docking.png';link.click();});
    el('fullscreen').addEventListener('click',async () => {if (document.fullscreenElement) await document.exitFullscreen();else await document.documentElement.requestFullscreen();});
    el('scene').addEventListener('mouseleave',() => {viewer.setHover(null);clearInteractionHover();});
    new ResizeObserver(() => {viewer.resize();viewer.render();}).observe(el('scene'));
    el('status').hidden=true;window.biomolViewer=viewer;window.biomolDockingScene=state;document.body.dataset.ready='true';
  } catch (error) {console.error(error);el('status').textContent=config.labels.error;document.body.dataset.ready='error';}
};
