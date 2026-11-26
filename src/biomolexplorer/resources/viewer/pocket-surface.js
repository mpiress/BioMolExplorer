/* Binding-site surface shared by the structure and the docking viewers.

   Solvent-excluded surface colored by the Kyte & Doolittle (1982)
   hydrophobicity of each residue, the convention current structural papers
   use for ligand-pocket complementarity: dodger blue for hydrophilic, white
   for neutral and orange red for hydrophobic.  The residues lining the ligand
   stay solid while the rest of the receptor is ghosted, so the site reads as a
   cavity without hiding the fold around it.  No provider or CDN requests. */
window.BioMolPocket = (function () {
  'use strict';
  const SCALE = {ILE:4.5, VAL:4.2, LEU:3.8, PHE:2.8, CYS:2.5, MET:1.9, ALA:1.8, GLY:-0.4, THR:-0.7,
    SER:-0.8, TRP:-0.9, TYR:-1.3, PRO:-1.6, HIS:-3.2, GLU:-3.5, GLN:-3.5, ASP:-3.5, ASN:-3.5,
    LYS:-3.9, ARG:-4.5};
  // Protonation and force-field variants carry the hydrophobicity of the residue they name.
  const ALIAS = {HID:'HIS', HIE:'HIS', HIP:'HIS', HSD:'HIS', HSE:'HIS', HSP:'HIS', CYX:'CYS',
    CYM:'CYS', ASH:'ASP', GLH:'GLU', LYN:'LYS', MSE:'MET'};
  const HYDROPHILIC = [30,144,255], NEUTRAL = [255,255,255], HYDROPHOBIC = [255,69,0];
  const OTHER = '#94a3b8';              // cofactors, ions and anything off the scale
  const EXTREME = 4.5;                  // arginine and isoleucine bound the scale
  const CONTEXT = 0.55;                 // the ghost keeps the fold readable behind the site
  const LARGE = 12000;                  // heavy atoms above which the exact surface is too slow
  const GRADIENT = 'linear-gradient(90deg,rgb(30,144,255),rgb(255,255,255),rgb(255,69,0))';
  const cache = {};

  function hydrophobicity(resn) {
    const name = String(resn || '').trim().toUpperCase();
    const value = SCALE[name] ?? SCALE[ALIAS[name]];
    return value === undefined ? null : value;
  }
  function paint(value) {
    const ratio = Math.max(-1, Math.min(1, value / EXTREME));
    const [from, to] = ratio >= 0 ? [NEUTRAL, HYDROPHOBIC] : [NEUTRAL, HYDROPHILIC];
    const weight = Math.abs(ratio);
    return '#' + [0,1,2].map(index => Math.round(from[index] + (to[index] - from[index]) * weight)
      .toString(16).padStart(2, '0')).join('');
  }
  function color(atom) {
    const name = String(atom.resn || '').trim().toUpperCase();
    if (cache[name] === undefined) {
      const value = hydrophobicity(name);
      cache[name] = value === null ? OTHER : paint(value);
    }
    return cache[name];
  }

  function controller(viewer) {
    let signature = null, ids = [], weights = new Map(), token = 0, level = 0.8;
    const clamp = value => Math.max(0.05, Math.min(1, value));
    const style = id => ({colorfunc: color, opacity: clamp(weights.get(id) * level)});
    function repaint() {
      ids.forEach(id => {try {viewer.setSurfaceMaterialStyle(id, style(id));} catch (error) {}});
    }
    return {
      clear() {
        ids.forEach(id => {try {viewer.removeSurface(id);} catch (error) {}});
        ids = []; weights.clear(); signature = null;
      },
      get drawn() {return ids.length > 0;},
      // Opacity is a material change: it never recomputes the mesh.
      setLevel(value) {
        level = clamp(value);
        if (ids.length) {repaint(); viewer.render();}
      },
      /* Surfaces are the expensive part of the scene: rebuild one only when the
         atoms behind it change, and keep every identifier so a redraw started
         while a mesh is still being generated can still dispose of it. */
      build(key, specs, atoms) {
        if (key === signature) return null;
        this.clear(); signature = key;
        const mine = ++token;
        const kind = atoms > LARGE ? 'VDW' : 'SES';
        const tasks = specs.map(spec => {
          let task;
          try {
            task = viewer.addSurface(kind, {colorfunc: color, opacity: clamp(spec.weight * level)},
                                     spec.atoms, spec.context);
          } catch (error) {return Promise.resolve(null);}
          if (task && task.surfid !== undefined) {ids.push(task.surfid); weights.set(task.surfid, spec.weight);}
          return Promise.resolve(task).catch(() => null);
        });
        return Promise.all(tasks).then(() => {
          if (mine !== token) return false;
          viewer.render();
          return true;
        });
      }
    };
  }

  return {hydrophobicity, color, controller, CONTEXT, GRADIENT, OTHER};
})();
