# Projects, files and versions

[Documentation](../README.md) · English · [Português](../projects.md)

## Choose the project folder

**New project** asks for a name, description, tags, color and **Project folder**.
Choose a new or empty folder using the button inside the field, which cannot be typed into. Desktop provides a system picker; web mode provides folder navigation and folder creation. In a browser, the
path belongs to the computer running the backend, rather than the browser's
computer. The application needs write access. Overlapping project folders and
paths containing the account database are rejected.

The folder contains `project.json` (portable configuration), `assets/` (original
inputs), `runs/` (experiments and results) and `.history/` (versions).
`project.json` follows saved edits; internal references use `@project` so that
import can move the folder. **Export** transfers the complete project; the JSON
alone does not contain the scientific data. Accounts, sessions, current permissions
and indexes stay in the workspace SQLite database. Logs remain centralized in the
application's configured log directory.

## Move to another computer

1. Finish or cancel any active execution.
2. Click **Export** on the project card and save the `.bme.zip` package.
3. Install the scientific environment on the destination computer and sign in.
4. Click **Import project**, select the package and choose a new or empty folder.
5. Open the imported project and run its pipeline normally.

The package includes configuration, uploaded files, results, run metadata and
history. Import verifies hashes and paths before creating the project, rebases
references and handles identifier collisions. Compatible completed stages are
reused. Missing artifacts or changed inputs and parameters require the affected
stages to run again. Researchers should assess differences in installed scientific
tools; choose **No, run again** in the reuse prompt when recalculation is needed.

Installation paths are machine settings. New DOCK6 calculations use the backend
installation; reusing a completed result does not require the previous installation
path.

The importing user becomes the owner. Accounts, passwords and sessions are never
exported. Existing collaborators on the destination workspace must accept a new
invitation; historical names and emails remain attributed to their original authors.
Restoring imported versions also preserves the acceptance requirement.

Browser uploads and downloads through the UI are limited to 200 MB. Larger exports
remain in `.exports/` for direct copying; use desktop to import them. Import accepts
up to 4 GB of uncompressed contents and 100,000 files. Readers may export projects;
import creates a private project owned by the authenticated user.

## Share and collaborate

The project owner's **Share** button invites registered accounts as **Editor** or
**Viewer**. Recipients accept from their workspace. Everyone accesses the same
persisted project, rather than a separate pipeline copy.

Applying block configuration saves immediately. Adding or removing blocks,
connections and block movements are saved after a short editing pause. **Save**
remains available. Open projects check collaborator updates every two seconds;
open configuration forms and local drafts are preserved.

Independent changes are merged when saved. Concurrent edits to the same field
produce a conflict message and preserve the local draft. Review the shared version
before reapplying that field. Backend authorization also checks revocations and
role changes; viewers cannot edit or execute projects.

## Review history and restore a version

**History** on the card or inside the project opens a newest-first table with the
date, time, author and description. Changes include configuration, metadata,
invitations, uploads, curation and completed executions. Times use the time zone
of the computer running the Python interface.

Owners may choose **Restore before**. Confirmation restores the state preceding
the selected change, including configuration, files, results and permissions.
The internal revision increases to invalidate stale drafts; ownership and the
audit history remain intact. Restoration itself is recorded. Later versions remain
preserved and can be recovered through the state preceding a subsequent restoration.
Active executions prevent rollback.

Versions use content-addressed copies in `.history/blobs/`: identical files occupy
one stored copy and modified files create new contents. Preserve this folder in
backups. This release does not automatically expire or purge versions. Project
deletion remains logical; deleted projects are outside normal access and cannot
be restored through the history button.

## Select multiple block inputs

Open block configuration and choose a source under **Inputs**. For completed stages,
**Result used for this input** lists files from the latest successful execution,
using relative paths to distinguish identical filenames. For CHEMBL220 retrieval,
choose `CHEMBL220_FULL.csv`, `CHEMBL220_MOLS.csv`, `CHEMBL220_SIMS.csv` or the
integrated `compounds.csv`.

**Add input** includes files from the same stage or other sources. You can mix stage outputs and uploaded files; a new connection adds a source. During execution, the popup offers **Process individually** or **Merge files**. Individual retains per-file results; merge creates a normalized CSV union. Duplicate codes for the same molecule are deduplicated; conflicting codes are reported and excluded by molecular cleaning. Other columns are preserved, and the first selected row wins for equivalent duplicates. Originals remain unchanged. Structural files with identical names and different contents must be renamed before combining.

Blocks show inputs and runs record resolved files. Individual forms combinations between distinct ports; DOCK6 keeps only corresponding combinations and consensus pairs poses/scores by identifier. The popup preserves choices when opening configuration, and confirmed selections are saved in the block. Projects remain limited to 100 blocks. When intact earlier results exist, execution asks whether to reuse them; compatible stages do not repeat selection. See the [manual](user_manual.md).

## Own inputs and completed results

**Upload my files** supplies inputs for a block's calculation. Upload validates
columns, SMILES, fingerprints, scores or coordinates according to the selected type.
An error explains the expected format and prevents the file from entering the
library. Valid originals are preserved; execution normalizes a separate copy.

For calculations already completed elsewhere, choose **How should this block be
used? → Use my completed results**, select or upload results and apply configuration.
The block is marked **Provided results**. Execution only organizes those files
in its result directory; the block's scientific calculation is bypassed.
Downstream stages consume them, and previous inputs and dependencies are inactive
in this mode.

| Result | Expected format |
| --- | --- |
| Compounds | UTF-8 CSV `molecule_chembl_id,canonical_smiles`; aliases `name,smiles` accepted; missing identifiers generated |
| Completed ADMET | CSV with identifier, SMILES and finite `TPSA,WLOGP`; other properties preserved |
| Fingerprints | CSV identifier and `fingerprint`, a consistent-length list of 0/1 bits |
| Similarity | CSV `source,target,value`, with finite score between 0 and 1 |
| PDB | ATOM/HETATM records with valid coordinates; optional `pdb_codes.csv` metadata |
| Prepared receptors | PDBQT, `pdb_codes.csv` and `centers.csv` with three coordinates per complex |
| Vina / DOCK6 | PDBQT containing `REMARK VINA RESULT` or `*_scored.mol2` with molecule, atom and Grid_Score records |

Completed ADMET creates an interactive EGG using supplied descriptors. Expand the stage in **Runs** and choose the result dataset; hover to identify compounds and click to view their 2D structures
and properties. The viewer does not recalculate descriptors.

The graph block receives ready relationships from similarity blocks or external files through the form; fingerprints are not a direct input. SMILES come from lineage or corresponding external tables. Individual retains independent analyses; merge combines selected relationships into one analysis. Select a specific MCC to feed the next stage. Interactive models, fragment images and degree presentations accompany project export. The fragment is displayed as SMILES extracted from the reference.
