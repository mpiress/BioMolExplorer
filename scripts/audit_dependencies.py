"""Audit project Python imports and environment declarations using the stdlib.

Run with Python 3.12: python scripts/audit_dependencies.py --output audit.json
Optional ADMET-AI imports are checked against their isolated environment. Native
executables and build/runtime infrastructure are retained explicitly.
"""
import argparse
import ast
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
DISTRIBUTIONS = {
    'Bio':'biopython', 'PIL':'pillow', 'admet_ai':'admet-ai', 'datasketch':'datasketch',
    'fastapi':'fastapi', 'flet':'flet', 'flet_web':'flet-web', 'matplotlib':'matplotlib',
    'networkx':'networkx', 'numpy':'numpy', 'pandas':'pandas', 'pymol':'pymol-open-source',
    'rcsbsearchapi':'rcsbsearchapi', 'rdkit':'rdkit', 'requests':'requests', 'scipy':'scipy',
    'seaborn':'seaborn', 'setuptools':'setuptools', 'urllib3':'urllib3', 'uvicorn':'uvicorn',
    'gemmi':'gemmi','meeko':'meeko','openbabel':'openbabel',
    'vina':'vina', 'websockets':'websockets',
}
INFRASTRUCTURE = {'python':'interpreter', 'pip':'environment installer',
    'openbabel':'obabel CLI: docking conversion/preparation',
    'torch':'CPU tensor runtime used by ADMET-AI (pinned to avoid CUDA installs)',
    'poppler':'pdftoppm CLI: PDF preview', 'setuptools':'setup.py/build backend'}
SKIP = {'.git','.conda','.admet-env','.venv','venv','__pycache__','.pytest_cache','build','dist','node_modules'}


def python_files(root):
    for folder, dirs, files in os.walk(root):
        dirs[:] = sorted(d for d in dirs if d not in SKIP and not d.endswith('.egg-info')
            and not (Path(folder)/d/'pyvenv.cfg').is_file() and not (Path(folder)/d/'conda-meta').is_dir())
        for name in sorted(files):
            if name.endswith('.py'):yield Path(folder)/name


def environment(path):
    """Read direct declarations in this project's simple Conda YAML manifests."""
    result = {}
    for line in path.read_text().splitlines():
        match = re.match(r'^\s+-\s+([a-zA-Z0-9_-]+)(?:\[[^\]]+\])?([=<>!].*)?\s*$',line)
        if match and match[1] not in ('conda-forge','defaults','nodefaults','pip:'):
            result[match[1].lower().replace('_','-')] = (match[2] or '').strip()
    return result


def audit(root=ROOT):
    stdlib = sys.stdlib_module_names
    local = {p.name for p in (root/'src').iterdir() if p.is_dir()}
    local |= {p.stem for p in (root/'src').glob('*.py')} | {p.stem for p in (root/'tests').glob('*.py')}
    found = defaultdict(set);dynamic=[];files=list(python_files(root))
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'),filename=str(path))):
            modules=[]
            if isinstance(node,ast.Import):modules=[a.name.split('.')[0] for a in node.names]
            elif isinstance(node,ast.ImportFrom) and node.module and not node.level:modules=[node.module.split('.')[0]]
            elif isinstance(node,ast.Call) and isinstance(node.func,(ast.Name,ast.Attribute)):
                name=getattr(node.func,'id','') or getattr(node.func,'attr','')
                if name in ('__import__','import_module') and node.args and isinstance(node.args[0],ast.Constant) and isinstance(node.args[0].value,str):
                    modules=[node.args[0].value.split('.')[0]]
                    dynamic.append(f'{path.relative_to(root)}:{node.lineno}: {node.args[0].value}')
            for module in modules:
                if module not in stdlib and module not in local:
                    found[module].add(f'{path.relative_to(root)}:{node.lineno}')
    main=environment(root/'requirements.yml');optional=environment(root/'requirements-admet.yml') if (root/'requirements-admet.yml').is_file() else main
    missing=[]
    for module in found:
        package=DISTRIBUTIONS.get(module)
        declared=optional if module=='admet_ai' else main
        if package not in declared and module!='playwright':missing.append(module)
    used={DISTRIBUTIONS.get(module) for module in found}
    orphans=sorted(set(main)-used-set(INFRASTRUCTURE))
    return dict(optional_development={'playwright':'Browser audit script only; not installed by this change'} if 'playwright' in found else {},python_files=len(files),imports={k:dict(package=DISTRIBUTIONS.get(k),locations=sorted(v),
        environment='optional ADMET' if k=='admet_ai' else 'main') for k,v in sorted(found.items())},
        dynamic_imports=dynamic,infrastructure=INFRASTRUCTURE,main=main,optional_admet=optional,
        missing=sorted(missing),orphans=orphans)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=audit()
    if args.output:args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('python_files','missing','orphans')},ensure_ascii=False))
    return bool(result['missing'] or result['orphans'])


if __name__=='__main__':raise SystemExit(main())
