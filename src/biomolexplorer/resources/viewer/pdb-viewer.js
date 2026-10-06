/* BioMolExplorer local structure explorer. No provider or CDN requests. */
(async function () {
  "use strict";
  const config = JSON.parse(document.getElementById("configuration").textContent);
  const labels = config.labels;
  const control = id => document.getElementById(id);
  const status = control("status");
  let viewer, models, modelIndex = 0;
  const waterNames = ["HOH", "WAT", "DOD"];
  const palette = ["#2dd4bf", "#60a5fa", "#c084fc", "#fbbf24", "#fb7185", "#34d399", "#f97316", "#a5b4fc"];
  const currentModel = () => models[modelIndex];
  const selection = () => control("chain").value ? {chain: control("chain").value === "__blank__" ? "" : control("chain").value} : {};
  const isWater = atom => waterNames.includes(atom.resn);
  const isLigand = atom => atom.hetflag && !isWater(atom);
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
    const polymer = {...selected, hetflag: false};
    model.setStyle(polymer, style);
    if (colorBy === "chain") {
      const values = [...new Set(model.selectedAtoms({}).map(atom => atom.chain || ""))].sort();
      values.forEach((chain, index) => {
        if (chainValue && chain !== (chainValue === "__blank__" ? "" : chainValue)) return;
        const tinted = {};
        Object.entries(style).forEach(([key, value]) => tinted[key] = {...value, color: palette[index % palette.length]});
        model.setStyle({chain, hetflag: false}, tinted);
      });
    }
    if (control("ligands").checked) model.setStyle({...selected, predicate: isLigand}, {
      stick: {colorscheme: "Jmol", radius: .2}, sphere: {colorscheme: "Jmol", scale: .22}});
    if (control("water").checked) model.setStyle({...selected, predicate: isWater}, {
      sphere: {color: "#60a5fa", scale: .18, opacity: .7}});
    model.setClickable(selected, true, atom => {
      viewer.removeAllLabels();
      const info = `${atom.resn} ${atom.resi}${atom.icode || ""} · ${labels.chain} ${atom.chain || "—"} · ${atom.atom} (${atom.elem})`;
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
    models = viewer.addModels(await response.text(), "pdb");
    if (!models.length || !models[0].selectedAtoms({}).length) throw new Error("No atoms in structure");
    models.forEach((model, index) => control("model").add(new Option(String(index + 1), String(index))));
    control("model").disabled = models.length === 1;
    chains(); draw(); viewer.zoomTo({model: modelIndex}); viewer.render();
    for (const id of ["representation", "color", "ligands", "water"]) control(id).addEventListener("change", draw);
    control("chain").addEventListener("change", () => {draw(); viewer.zoomTo({model: modelIndex, ...selection()}); viewer.render();});
    control("model").addEventListener("change", () => {modelIndex = Number(control("model").value); chains(); draw(); viewer.zoomTo({model: modelIndex}); viewer.render();});
    control("center").addEventListener("click", () => {viewer.zoomTo({model: modelIndex, ...selection()}); viewer.render();});
    control("snapshot").addEventListener("click", () => download(viewer.pngURI(), config.name.replace(/\.pdb$/i, "") + ".png"));
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
