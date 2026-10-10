/* BioMolExplorer local structure explorer. No provider or CDN requests. */
(async function () {
  "use strict";
  const config = JSON.parse(document.getElementById("configuration").textContent);
  if (config.scene) {await window.BioMolDockingScene(config); return;}
  const labels = config.labels;
  const control = id => document.getElementById(id);
  const status = control("status");
  let viewer, models, modelIndex = 0;
  const waterNames = ["HOH", "WAT", "DOD"];
  const proteinNames = new Set(["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "SEC", "PYL"]);
  const palette = ["#2dd4bf", "#60a5fa", "#c084fc", "#fbbf24", "#fb7185", "#34d399", "#f97316", "#a5b4fc"];
  const currentModel = () => models[modelIndex];
  const selection = () => control("chain").value ? {chain: control("chain").value === "__blank__" ? "" : control("chain").value} : {};
  const isWater = atom => waterNames.includes(atom.resn);
  // MOL2 lacks PDB's HETATM distinction. Keep small molecules visible in
  // the ribbons+ligands representation and honor the ligand visibility switch.
  const isLigand = atom => (config.ligand || atom.hetflag || (config.format === "mol2" && !proteinNames.has(atom.resn))) && !isWater(atom);
  const isPolymer = atom => !isLigand(atom) && !isWater(atom);
  const download = (url, name) => {
    const link = document.createElement("a"); link.href = url; link.download = name; link.click();
  };
  function chains() {
    const values = [...new Set(currentModel().selectedAtoms({}).map(atom => atom.chain || ""))].sort();
    control("chain").replaceChildren(new Option(labels.all, ""));
    values.forEach(chain => control("chain").add(new Option(chain || "—", chain || "__blank__")));
    const atoms = currentModel().selectedAtoms({}).length;
    control("summary").textContent = `${atoms.toLocaleString()} ${labels.atoms} · ${values.length} ${values.length === 1 ? labels.chain.toLocaleLowerCase() : labels.chains}`;
    return values;
  }
  function draw() {
    models.forEach((model, index) => index === modelIndex ? model.show() : model.hide());
    const model = currentModel(), chainValue = control("chain").value;
    const selected = chainValue ? {chain: chainValue === "__blank__" ? "" : chainValue} : {};
    model.setStyle({}, {}); viewer.removeAllLabels();
    const colorBy = control("color").value;
    const coloring = colorBy === "element" ? {colorscheme: "Jmol"} : colorBy === "spectrum" ? {color: "spectrum"} : {};
    const representation = control("representation").value;
    let style;
    if (representation === "sticks") style = {stick: {...coloring, radius: .16}};
    else if (representation === "spheres") style = {sphere: {...coloring, scale: .65}};
    else if (representation === "lines") style = {line: {...coloring, linewidth: 1.4}};
    else style = {cartoon: {...coloring, thickness: .35, arrows: true}};
    // Style the polymer separately, so ligand/water switches work in every representation.
    const polymer = {...selected, predicate: isPolymer};
    model.setStyle(polymer, style);
    if (colorBy === "chain") {
      const values = [...new Set(model.selectedAtoms({}).map(atom => atom.chain || ""))].sort();
      values.forEach((chain, index) => {
        if (chainValue && chain !== (chainValue === "__blank__" ? "" : chainValue)) return;
        const tinted = {};
        Object.entries(style).forEach(([key, value]) => tinted[key] = {...value, color: palette[index % palette.length]});
        model.setStyle({chain, predicate: isPolymer}, tinted);
      });
    }
    const ligandStyle = control("ligand-representation").value;
    const ligandVisible = atom => isLigand(atom) && (control("hydrogens").checked || !["H", "D"].includes(atom.elem.toUpperCase()));
    if (control("ligands").checked) model.setStyle({...selected, predicate: ligandVisible}, ligandStyle === "spheres" ? {sphere: {colorscheme: "Jmol", scale: .65}} :
      ligandStyle === "lines" ? {line: {colorscheme: "Jmol"}} : {stick: {colorscheme: "Jmol", radius: .2}, sphere: {colorscheme: "Jmol", scale: .22}});
    if (control("water").checked) model.setStyle({...selected, predicate: isWater}, {
      sphere: {color: "#60a5fa", scale: .18, opacity: .7}});
    model.setClickable(selected, true, atom => {
      viewer.removeAllLabels();
      const info = config.ligand ? `${atom.atom || atom.serial} (${atom.elem})` : `${atom.resn} ${atom.resi}${atom.icode || ""} · ${labels.chain} ${atom.chain || "—"} · ${atom.atom} (${atom.elem})`;
      control("atom").textContent = info;
      viewer.addLabel(info, {position: atom, backgroundColor: "#111b2b", fontColor: "#e5edf9", fontSize: 12, backgroundOpacity: .9});
      viewer.render();
    });
    viewer.render();
  }
  try {
    if (!window.$3Dmol) throw new Error("Viewer library unavailable");
    try {
      viewer = $3Dmol.createViewer(control("scene"), {backgroundColor: "#0b1220", antialias: true});
    } catch (error) {status.textContent = labels.webgl; return;}
    viewer.setProjection("perspective");
    const response = await fetch(location.pathname.replace(/\/$/, "") + "/structure", {cache: "no-store", credentials: "same-origin"});
    if (!response.ok) throw new Error("Structure unavailable");
    let structure = await response.text();
    let format = config.format || "pdb";
    if (format === "pdbqt") {
      // PDB parsers treat ENDROOT/ENDBRANCH as model boundaries. Strip torsion
      // records while retaining every atom and the actual Vina model boundaries.
      if (structure.includes("REMARK VINA RESULT:")) {
        config.ligand = true;
        config.representation = "sticks";
      }
      const elements = {A: "C", OA: "O", NA: "N", SA: "S", HD: "H", HS: "H"};
      structure = structure.split(/\r?\n/).filter(line => /^(ATOM  |HETATM|MODEL |ENDMDL|CONECT|TER   )/.test(line)).map(line => {
        if (!/^(ATOM  |HETATM)/.test(line)) return line;
        const type = line.trim().split(/\s+/).pop();
        const element = elements[type] || type;
        return line.slice(0, 66).padEnd(76) + element.padStart(2);
      }).join("\n");
      format = "pdb";
    }
    const loaded = viewer.addModels(structure, format, {keepH: true});
    models = loaded.filter(model => model.selectedAtoms({}).length);
    loaded.filter(model => !model.selectedAtoms({}).length).forEach(model => viewer.removeModel(model));
    if (!models.length || !models[0].selectedAtoms({}).length) throw new Error("No atoms in structure");
    models.forEach((model, index) => control("model").add(new Option(String(index + 1), String(index))));
    control("model").disabled = models.length === 1;
    control("representation").value = config.representation || "cartoon";
    if (config.ligand) {
      control("color").value = "element";
      control("color").querySelector('option[value="spectrum"]').disabled = true;
      control("representation").disabled = true;
      control("representation").querySelector('option[value="cartoon"]').disabled = true;
      control("chain").disabled = true;
      control("water").disabled = true;
    }
    if (config.generated) {
      control("conformer").textContent = labels.conformer;
      control("conformer").hidden = false;
    }
    chains(); draw(); viewer.zoomTo({model: modelIndex}); viewer.render();
    for (const id of ["representation", "ligand-representation", "hydrogens", "color", "ligands", "water"]) control(id).addEventListener("change", draw);
    control("chain").addEventListener("change", () => {draw(); viewer.zoomTo({model: modelIndex, ...selection()}); viewer.render();});
    control("model").addEventListener("change", () => {modelIndex = Number(control("model").value); chains(); draw(); viewer.zoomTo({model: modelIndex}); viewer.render();});
    control("center").addEventListener("click", () => {viewer.zoomTo({model: modelIndex, ...selection()}); viewer.render();});
    control("snapshot").addEventListener("click", () => download(viewer.pngURI(), config.name.replace(/\.(pdb|pdbqt|mol2|sdf)$/i, "") + ".png"));
    control("fullscreen").addEventListener("click", async () => {
      try {if (document.fullscreenElement) await document.exitFullscreen(); else await document.documentElement.requestFullscreen();} catch (_) {}
    });
    new ResizeObserver(() => {viewer.resize(); viewer.render();}).observe(control("scene"));
    status.hidden = true;
    // A public viewer handle supports inspection and future local integrations.
    window.biomolViewer = viewer;
    document.body.dataset.ready = "true";
  } catch (error) {status.textContent = labels.error; document.body.dataset.ready = "error";}
})();
