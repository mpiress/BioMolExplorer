# Molecular viewer: poses and 3D interactions

[Documentation](../README.md) · English · [Português](../molecular_viewer.md)

The local WebGL explorer inspects structures, compares calculated poses with
their receptors and identifies geometric interactions. Molecular files remain
within the application; viewing does not send structures to external services.

## Open a structure or result

| Source | Action and contents |
| --- | --- |
| Retrieved PDB or molecular file | **View 3D structure** opens the selected file |
| Compound table | **3D** generates a local SMILES conformer; it is not a docking pose |
| Completed redocking | **3D** overlays the best Vina pose, receptor and available reference |
| Vina or DOCK6 docking | **3D** overlays the best pose and its receptor; a prepared reference may accompany the result |
| Consensus | **3D Vina** and **3D DOCK6** open each engine's poses when their files and receptors are available |

The structure component accepts PDB, PDBQT, MOL2 and SDF. When opening a file
with multiple models, **Model** selects a model. Result scenes automatically
choose the lowest-scoring pose; without a score, they show the first pose,
explicitly labeled as such. These scenes do not offer **Model** selection. To
inspect other poses, open the corresponding file under **View simulation** or
download it.

Consensus can compute scores without structural attachments. Viewing a pose with
its receptor requires both to remain available and associated with the result.
Missing or ambiguous files are reported; another receptor is not substituted.

## Navigate and choose styles

Drag to rotate, use the mouse wheel to zoom and right-drag to pan. **Recenter**
restores the overall view, **Fullscreen** expands the scene and **Save image**
exports it as PNG.

1. Under **Receptor style**, choose ribbons, sticks, spheres or lines. Docking scenes start MOL2 receptors in sticks.
2. Under **Ligand style**, choose sticks, spheres or lines independently of the receptor. This also controls visible reference ligands.
3. Use **Color by** and **Chain** to explore the receptor. The pose appears in cyan and the reference in gold.
4. Enable **Ligand hydrogens** to show H atoms present in the file; disable it to hide them. This does not generate atoms or change preparation or scores.
5. Use **Ligands**, **Water** and the layer checkboxes to choose visible elements. Click an atom for its details.

Vina and DOCK6 may retain different hydrogen sets. H atoms absent from Vina output
cannot appear when the field is enabled. H visibility is a presentation choice:
classification continues to use atoms available in the original file.

## Compare pose, receptor and reference

Layers identify the **Docking receptor**, **Best pose** or **First unscored pose**,
and the available reference. In redocking, the crystallographic ligand appears
in gold; the rest of the crystallographic complex can be enabled as a comparison
layer. Without the original PDB, an available prepared reference is labeled as
such.

Original coordinates are preserved without automatic realignment. The
**Residues** table compares minimum heavy-atom distances for the pose and
reference and can be downloaded as CSV. In the viewer, **Nearby residues** lets
you select and highlight a residue. The initial cutoff is 4 Å, adjustable from
2 to 8 Å. This filters the proximity list, not the chemical criteria used to
classify interactions.

## Identify interactions

When chemical topology and residue identities allow classification, **Show 3D
interactions** provides filters by type. Each relation receives a dashed trace
between its participating atoms, or between ring centers for π–π. Traces sharing
endpoints receive a small visual offset so that each remains identifiable;
chemical endpoints and reported distances remain unchanged.

| Legend and trace color | Identified type |
| --- | --- |
| Green | Hydrogen bond |
| Pink | Parallel π–π |
| Coral | T-shaped π–π |
| Lilac | Hydrophobic contact |
| Lime green | Geometric van der Waals contact |

Hover over a dash: the tooltip shows **type, residue, chain and distance in Å**,
with a border matching the relation's color. The legend keeps these colors even
when a type is absent from the pose. The interaction list displays the same
data; click a row to highlight and zoom to the residue.

Filters hide individual types. Hiding the pose, receptor, ligands or all
interactions also removes their traces and tooltip. In ribbons mode, involved
residues additionally appear in sticks. Classified interactions describe the
selected pose, not the comparison reference.

## Interpret the results

Classification is geometric and limited to the five legend types. It uses
ligand chemistry and original positions; it does not calculate energies or
affinity or establish experimental binding. Hydrogen bonds require explicit H
and suitable orientation; aromatic rings require valid topology. PDBQT poses
rely on the associated SMILES to recover bond orders. Without this information,
the interface explains unavailability and retains distance contacts. A relation
missing from the scene does not prove its absence.

Hydrophobic, hydrogen-bond and π–π thresholds follow the
[PLIP parameters](https://github.com/pharmai/plip/blob/master/plip/basic/config.py),
but the tool does not run the complete PLIP classifier. Van der Waals contacts
use proximity relative to RDKit atomic radii. Hydrophobic and van der Waals
contacts retain the closest pair per residue.

**Footprint** is a separate analysis: DOCK6 results and consensus with an
associated PDF show per-residue van der Waals and electrostatic energies,
comparing the reference with the final pose. A van der Waals contact drawn in
3D is not this energy. See the [docking and consensus manual](user_manual.md#vina-dock6-and-consensus).

## Access, requirements and updates

Readers can inspect and download results. Authorization is checked when opening
structures and fetching the scene. The viewer requires a browser with WebGL and
graphics acceleration; the opening limit is 32 MB. Desktop mode opens a local
page; web mode uses the application's origin.

After updating the code, restart the application and open a new viewer window.
Existing windows retain previous scripts. Interface asset URLs are versioned by
content. If opening fails, return to the project and open the viewer again.

Geometry tests and actual Chrome mouse movements are described in the
[viewer validation report](../validation/viewer_controls_2026-10-10.md) (Portuguese).
