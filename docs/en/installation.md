# Installing and configuring BioMolExplorer

> Source synchronization: 2026-11-13. See [current source reference](current_reference.md) for exact registered parameters, defaults, compatibility and unsupported behaviors. This describes the checkout, not an installation or global release acceptance.


[Documentation](../README.md) · English · [Português](../installation.md)

Follow this order: download the project → install the external tools → prepare the Conda environment → configure and check the executables → start BioMolExplorer. Terminal examples use Linux and Bash, the current backend execution environment. Install Git and Anaconda or Miniconda beforehand.

**Before starting the application, install and configure UCSF Chimera 1.17 and DOCK6 6.11 on the computer that will run the calculations.** Installing the Python package or interface does not install these two tools. The interface may open without them, but preparation, redocking and docking stages that depend on them will fail during execution.

## 1. Download the project from GitHub

Open the [official BioMolExplorer repository](https://github.com/mpiress/BioMolExplorer). With Git installed, choose where to keep the source code and run:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

To download without Git, click **Code → Download ZIP** in the repository, extract the archive and open a terminal inside the extracted folder, usually `BioMolExplorer-master`. Run subsequent commands there, where `requirements.yml` and `pyproject.toml` should be present; do not use the `src` folder or run commands inside the ZIP.

The source code folder is separate from each research project's folder. You will select where inputs and results are stored when creating a project in the interface.

## 2. Install the required external tools

| Tool | Version required by this guide | Use in BioMolExplorer | Download and installation |
| --- | --- | --- | --- |
| UCSF Chimera | 1.17 | Active receptor and ligand preparation | [Previous releases from UCSF](https://www.cgl.ucsf.edu/chimera/olddownload.html) and [installation instructions](https://www.cgl.ucsf.edu/chimera/download.html) |
| DOCK6 | 6.11 | Docking, refinement and scoring, with its accessory programs | [Official DOCK 6 page](https://dock.compbio.ucsf.edu/DOCK_6/index.htm), [6.11 release notes](https://dock.compbio.ucsf.edu/DOCK_6/new_in_6.11.txt) and the manual included in the distribution |

Install UCSF Chimera 1.17, Open Babel and Vina for preparation/redocking; install DOCK6 6.11 for its docking stages. No standalone AmberTools environment is required. ChimeraX does not interpret Chimera scripts.

Install Chimera using the installer for your system. For DOCK6, follow the 6.11 distribution instructions, build the accessory programs as well, and retain the complete installation, including `bin/` and `parameters/`. Run the tests supplied by DOCK6 as described in its manual.

The molecular surface is computed natively in Python using NumPy and SciPy. No DMS installation is needed; see the [port and validation report](dms_migration.md).

Install the tools required by the selected protocol on the worker computer, including in web deployments. Chimera is required for preparation; other analysis branches do not require it.

## 3. Create the scientific environment and install the interface

From the downloaded source code root:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
python -m pip check
command -v biomolexplorer-ui
biomolexplorer-ui --help
```

Run one line at a time and wait for it to finish without errors. The `pip install` command is required: installing libraries through Conda does not register BioMolExplorer commands. The `-e` option links the installation to the source folder; do not move or delete it afterwards. `.[ui]` includes the web and desktop interface. The `--help` command should display options without opening the application.

If the environment already exists, update it with `conda env update -f requirements.yml --prune`, then activate it. Repeat `python -m pip install -e '.[ui]'` and the checks above. The environment includes Open Babel, Vina and the UI. It does not replace Chimera and DOCK6, which remain external installations. Compatible ranges — notably NumPy 1.x and pandas 2.x — prevent machine-learning extras from silently replacing the scientific worker stack. The package requires Python 3.12.

## 4. Configure paths and check the installation

With the Conda environment active, add the executable directories to `PATH`. Replace every example path with your actual installation path:

```bash
export BIOMOL_DOCK6_ROOT="/path/to/dock6-6.11"
export PATH="/path/to/chimera-1.17/bin:$BIOMOL_DOCK6_ROOT/bin:$PATH"
```

`BIOMOL_DOCK6_ROOT` is a convenience variable used in this guide. The application receives that path through `--dock6-path`, shown below. Supply the DOCK6 root containing `bin/` and `parameters/`, rather than only the executable or `bin/` directory. Start the application from the same terminal so workers inherit `PATH`. To retain these settings between sessions, add the exports with actual paths to your shell startup file and open a new terminal.

Check executables before starting:

```bash
for tool in chimera dock6 sphgen showbox grid obabel vina; do
    if command -v "$tool"; then
        printf 'OK: %s available\n' "$tool"
    else
        printf 'MISSING: %s — install it or fix PATH before starting\n' "$tool"
    fi
done
if test -d "$BIOMOL_DOCK6_ROOT/parameters"; then
    echo 'OK: DOCK6 parameters directory found'
else
    echo 'MISSING: parameters directory — check the DOCK6 installation root'
fi
```

**Do not execute a protocol while a required tool or directory is missing.** Confirm the paths of the tools actually used; Chimera 1.17 is required for preparation and DOCK6 6.11 only for its protocols, especially if multiple versions are installed. Locating an executable confirms its presence on `PATH`; it does not validate its version or the protocol's behavior. Check versions in the tools/distributions and run a reference case before starting your study.

## 5. Start the configured application

In the same terminal, with the environment active and checks completed:

```bash
biomolexplorer-ui --web --language en --dock6-path "$BIOMOL_DOCK6_ROOT"
```

Open `http://127.0.0.1:8550`. For desktop mode, omit `--web`. To serve without automatically opening a browser, use `--web --no-browser`. If you select another interpreter through `--worker-python`, prepare that scientific environment too and ensure workers can access the external tools.

Continue with the [user manual](user_manual.md) to create an account, select a project folder and configure the pipeline. For CLI execution and technical parameters, see [backend usage](backend_usage.md).

## Checks before execution

Redocking with preparation checks Chimera, Open Babel and Vina. Prepared inputs require only Vina. Historical automatic/Python settings now use Chimera.

**Prepare for docking** uses Chimera for receptors, ligands, `gas`/`am1` charges and minimization. Open Babel converts Vina/DOCK6 formats. These blocks do not require an AmberTools environment. See [molecular preparation](molecular_preparation.md).

**Amber Score · DOCK6** requires a complete DOCK6 **6.11** installation, including its bundled AMBER helpers and parameter library. No standalone AmberTools environment, sander or cpptraj is used. Follow the native Amber Score section in the [manual](user_manual.md).

See [Redocking configuration](redocking_configuration.md) for pair selection, cofactors, prepared inputs, diagnosis and results.

## Optional ADMET-AI predictor

The **Local RDKit** backend uses the main scientific environment. Selecting **ADMET-AI** also requires its predictor on the worker computer. It intentionally stays in an isolated environment so it cannot replace pandas, RDKit, NumPy or docking-worker dependencies.

From the repository root, with compatible Python 3.12:

```bash
python scripts/setup_admet_ai.py
```

The script creates `.admet-env`, installs CPU PyTorch and `admet-ai==2.0.1`, runs an inference and writes `.ready` only after success. Installation needs access to package distributions and models. Workers discover the repository-local environment automatically. To use another directory:

```bash
python scripts/setup_admet_ai.py --directory /path/to/admet-env
export BIOMOL_ADMET_PYTHON="/path/to/admet-env/bin/python"
```

Start the application from that terminal so workers inherit the setting. This variable selects the predictor Python independently of `--worker-python`. When moving computers, recreate the environment and update paths; exporting a project does not install scientific libraries.

For **ADMET-AI is not installed in the worker environment**, check the availability shown in ADMET settings, the interpreter path, installation and worker access to models. Restart after configuring the environment. Alternatively select **Local RDKit** explicitly: it provides local descriptors, estimates and rules without `ML_` endpoints. There is no silent model substitution. Use the isolated environment described above.

See [hypothesis-driven ADMET](admet_research.md) for units, profiles, alerts, outputs and model limits.


## Reopen the application and troubleshoot

In a new terminal, run `conda activate BioMolExplorer`, repeat the paths from step 4 if you did not save them in your shell, then run the command from step 5. Reinstallation is not needed each time.

| Message or situation | Solution |
| --- | --- |
| `conda: command not found` | Install Miniconda/Anaconda and open a new Bash terminal. Initialize Bash following the installer instructions if needed. |
| `CondaError: Run conda init before conda activate` | Run `conda init bash`, close the terminal and open another one. |
| `biomolexplorer-ui: command not found` | Activate `BioMolExplorer`, return to the folder containing `pyproject.toml` and repeat `python -m pip install -e '.[ui]'`. Check `python -m pip --version` and `command -v biomolexplorer-ui`: both should point into the active environment. |
| Missing `pyproject.toml` | Enter the extracted/cloned source folder, rather than `src`. |
| Installation or `pip check` failure | Stop and record the full message. Do not use `sudo pip` or the system Python. For old environments, prefer a fresh one using `conda env create -n BioMolExplorer-clean -f requirements.yml`, then activate it before editable installation. |
| Desktop interface does not open | Try `--web --no-browser` and open the address shown in the terminal. |

## Scope and scientific environment

The recommended route uses `requirements.yml` for scientific dependencies and `python -m pip install -e '.[ui]'` to register the application. Chimera, DOCK6 and Conda are external tools, not pip dependencies. The `ui` extra alone does not install the full scientific environment. ADMET-AI remains optional and isolated; do not install `.[admet-ml]` in the docking environment.

This guide uses Linux/Bash and Python 3.12. No Linux lock is distributed in this checkout; do not use a nonexistent lock file. A clean installation requires network access and availability of the declared packages. `pip check` checks Python metadata, not scientific protocols. Validate external tools and a reference case before using results in research.

## Native AMBER structural support

Native Amber Score includes target-independent sanitation, 32 ff03/PTM models, CCD-backed isolated-organic-component parameters and execution-local DOCK6 torsion compatibility when required. It preserves the original DOCK6 installation. The default `residue_parameter_policy` is `auto_curated`; `native_ff94` and `ff03_ffptm` are alternatives. Compiler/make/flex/zlib requirements and the noncommercial **CC BY-NC 4.0** license of bundled PTM data are separate from the MIT code license. Unknown covalent PTMs and unsupported metals require validated parameters. See [structural support](amber_score_structural_support.md) for supported states, caches, limits and reproduction.
