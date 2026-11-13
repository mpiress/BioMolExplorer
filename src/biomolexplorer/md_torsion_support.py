"""Execution-local DOCK6 6.11 support for general AMBER periodic torsions."""
import math
import os
from pathlib import Path
import shlex
import shutil
import tempfile
import sys

from .calculation_cache import CalculationCache,identity,copy_artifact
from .processes import run_command,ScientificToolError
from .stage_cache import file_digest

HEADER=Path(__file__).with_name('resources')/'amber_support/periodic_torsion.h'
PATCH_VERSION='periodic-torsion-forward-ad-v2'
HOOK='''
        /* BioMolExplorer full periodic torsion compatibility, local build. */
        if (bme_extended_torsion(Pn[atyp],Phase[atyp])) {
          double positions[12], gradient[12];
          int atom_offsets[4], bi, bj;
          if (dim != 3) { fprintf(stderr,"Extended torsions require 3D coordinates\\n"); exit(1); }
          atom_offsets[0]=at1; atom_offsets[1]=at2; atom_offsets[2]=at3; atom_offsets[3]=at4;
          for(bi=0;bi<4;bi++) for(bj=0;bj<3;bj++) positions[3*bi+bj]=x[atom_offsets[bi]+bj];
          e=bme_periodic_torsion(positions,Pn[atyp],Pk[atyp],Phase[atyp],gradient);
          e_tors+=e;
          for(bi=0;bi<4;bi++) for(bj=0;bj<3;bj++) f[foff+atom_offsets[bi]+bj]+=gradient[3*bi+bj];
          if(Pn[atyp]<0.) { atyp++; goto multi_term; }
          continue;
        }
'''


def read_torsion_parameters(path):
    values={};active=None;width=16
    for line in Path(path).read_text().splitlines():
        if line.startswith('%FLAG '):
            name=line.split()[1]
            active=name if name in ('DIHEDRAL_PERIODICITY','DIHEDRAL_PHASE') else None
            if active:values[active]=[]
        elif line.startswith('%FORMAT'):
            import re
            match=re.search(r'[eEfF](\d+)',line)
            if active and not match:raise ValueError('Unsupported AMBER torsion format')
            if match:width=int(match[1])
        elif active and not line.startswith('%'):
            values[active].extend(float(line[i:i+width].replace('D','E')) for i in range(0,len(line),width) if line[i:i+width].strip())
    pn=values.get('DIHEDRAL_PERIODICITY',[]);phase=values.get('DIHEDRAL_PHASE',[])
    if len(pn)!=len(phase):raise ValueError('Incomplete AMBER torsion parameters')
    unsupported=[]
    for index,(number,angle) in enumerate(zip(pn,phase)):
        if not math.isfinite(number) or not math.isfinite(angle) or abs(number)<1 or abs(abs(number)-round(abs(number)))>1e-8:
            raise ValueError('Invalid AMBER periodic torsion parameter')
        reduced=abs(math.remainder(angle,2*math.pi))
        if round(abs(number)) not in (1,2,3,4,6) or min(reduced,abs(reduced-math.pi))>1e-5:
            unsupported.append(dict(index=index,periodicity=number,phase_radians=angle))
    return unsupported


def compiler_command(name,fallback):
    command=shlex.split(os.environ.get(name,'') or fallback)
    executable=shutil.which(command[0])
    if not executable and not os.environ.get(name):
        candidates=sorted((Path(sys.prefix)/'bin').glob('*-'+fallback))
        executable=next((str(path) for path in candidates if path.is_file() and os.access(path,os.X_OK)),None)
    if not executable:raise ScientificToolError('Compiler unavailable for DOCK6 torsion support: '+command[0],
        'TOOL_NOT_FOUND','Install the native compiler dependencies declared in requirements.yml.')
    return [executable,*command[1:]]


def compatible_binary(root,directory):
    """Compile only in a local copy; immutable installation files stay intact."""
    root=Path(root);directory=Path(directory)
    source=root/'src/dock';eff=source/'nab/eff.c'
    required=[eff,source/'nab/sff.c',root/'install/config.h',root/'install/rules.h']
    if not all(p.is_file() for p in required):
        raise ScientificToolError('DOCK6 6.11 sources are required for full periodic torsion support.',
            'DOCK6_SOURCE_REQUIRED','Keep src/dock and install/config.h/rules.h from the licensed DOCK6 6.11 distribution.')
    cc=compiler_command('CC','gcc');cxx=compiler_command('CXX','g++')
    make=shutil.which('make')
    if not make:raise ScientificToolError('make is required for DOCK6 torsion support.','TOOL_NOT_FOUND','Install make from requirements.yml.')
    text=eff.read_text();anchor='      multi_term:\n'
    if text.count(anchor)!=1 or 'bad value for Pn:' not in text:
        raise ValueError('Unsupported DOCK6 torsion source; refusing an unverified patch.')
    source_files=[p for p in source.rglob('*') if p.is_file() and p.suffix in ('.c','.cpp','.h','.o')]
    recipe=dict(version=PATCH_VERSION,header=file_digest(HEADER),hook=identity(HOOK),
        native=file_digest(root/'bin/dock6'),sources={str(p.relative_to(root)):file_digest(p) for p in [*required,*source_files]},
        cc=cc,cxx=cxx,compiler_files={command[0]:file_digest(command[0]) for command in (cc,cxx)})
    cache=CalculationCache('dock6-periodic-kernel',directory);key=identity(recipe)
    target=directory/'.amber-support/dock6';target.parent.mkdir(parents=True,exist_ok=True)
    with cache.locked(key):
        if not cache.restore(key,{'dock6':target}):
            with tempfile.TemporaryDirectory(prefix='dock6-build-',dir=cache.root) as folder:
                build=Path(folder);(build/'install').mkdir()
                for name in ('config.h','rules.h'):shutil.copy2(root/'install'/name,build/'install'/name)
                shutil.copytree(source,build/'src/dock')
                local=build/'src/dock/nab'
                (local/'periodic_torsion.h').write_bytes(HEADER.read_bytes())
                (local/'eff.c').write_text('#include "periodic_torsion.h"\n'+text.replace(anchor,anchor+HOOK))
                (local/'sff.o').unlink(missing_ok=True)
                run_command([make,'-j1','dock','CC='+shlex.join(cc),'CXX='+shlex.join(cxx),'LOAD='+shlex.join(cxx)],
                    cwd=build/'src/dock',timeout=600,tool_name='DOCK6 periodic torsion compatibility build')
                binary=build/'src/dock/dock6'
                if not binary.is_file():raise ValueError('DOCK6 compatibility build produced no executable')
                copy_artifact(binary,target);target.chmod(0o700)
                run_command([str(target),'-h'],timeout=10)
                cache.remember(key,{'dock6':target},dict(recipe=recipe))
    target.chmod(0o700)
    return target,dict(version=PATCH_VERSION,recipe_sha256=key,native_sha256=recipe['native'],binary_sha256=file_digest(target))
