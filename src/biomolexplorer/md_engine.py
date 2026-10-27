"""DOCK6 AMBER-score adapter with isolated preparation and optional AmberTools MD."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil

from .md_config import DEFAULTS,validate_settings
from .md_inputs import validate_structure
from .processes import run_command,ScientificToolError
from .stage_cache import file_digest
from .storage import write_json,write_text


def tools(dock6_app_path,ambertools_path,settings,ambertools_env=''):
    from .docking_tools import dock6_root
    root=dock6_root(dock6_app_path)
    if not root:raise ValueError('Configure a instalação DOCK6 antes da dinâmica.')
    needed=['dock6','prepare_amber.pl','amberize_receptor','amberize_ligand','amberize_complex','tleap','teLeap','antechamber','parmchk']
    found={name:root/'bin'/name for name in needed}
    if settings['amber_score_type']=='dock6_sander':
        from .preparation_tools import amber_tool
        for name in ['sander']+(['cpptraj'] if settings['analyze_trajectory'] else []):
            found[name]=amber_tool(name,ambertools_path,environment=ambertools_env)
    for name,path in found.items():
        if not path.is_file() or name!='prepare_amber.pl' and not os.access(path,os.X_OK):
            raise ScientificToolError('Ferramenta DOCK6 indisponível: '+name,'TOOL_NOT_FOUND','Verifique a instalação DOCK6 e suas ferramentas AMBER.')
    return root,found


def tool_identity(root,found):
    library=root/'parameters/leap'
    return dict(tools={name:file_digest(path) for name,path in found.items()},
        parameters={p.relative_to(library).as_posix():file_digest(p) for p in sorted(library.rglob('*')) if p.is_file()},
        adapter={name:file_digest(Path(__file__).with_name(name)) for name in ('md_engine.py','md_analysis.py','md_config.py','md_cache.py','md_inputs.py','processes.py')})


def prepare(directory,root,found,existing_charges):
    """Relocate copies of installation scripts; never change installed tools.

    Some DOCK distributions embed their original build directory. Helpers are
    invoked by relative names, avoiding Perl interpolation of workspace paths.
    """
    directory=Path(directory)
    for name in ('prepare_amber.pl','amberize_receptor','amberize_ligand','amberize_complex'):
        text=found[name].read_text()
        if name=='prepare_amber.pl':
            text=re.sub(r"'[^']*'/bin/(amberize_(?:receptor|ligand|complex))",r'./\1',text)
            if existing_charges:
                # Old bundled antechamber overflows its charge-type label
                # when no -c is supplied. Read the exact input charges instead.
                text=text.replace('$charge_method = "-s 2"','$charge_method = "-s 2 -c rc -cf existing.charges"')
        else:
            text=re.sub(r"export DOCK_HOME='[^']*'",lambda _: 'export DOCK_HOME='+shlex.quote(str(root)),text)
            text=re.sub(r"'[^']*'/bin/([A-Za-z0-9_]+)",lambda m:shlex.quote(str(root/'bin'/m[1])),text)
            text=text.replace(shlex.quote(str(root/'bin/tleap')),'./tleap')
        target=directory/name;write_text(text,target);target.chmod(0o700)
    # tleap is another relocated shell launcher in several DOCK distributions.
    arguments=[str(found['teLeap'])]+['-I'+str(root/'parameters/leap'/name) for name in ('prep','lib','parm','cmd')]
    write_text('#!/bin/sh\nexec '+shlex.join(arguments)+' "$@"\n',directory/'tleap')
    (directory/'tleap').chmod(0o700)
    if existing_charges:
        atom_block=(directory/'lig.mol2').read_text().split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
        charges=[line.split()[8] for line in atom_block.splitlines() if line.strip()]
        total=sum(map(float,charges))
        if not math.isfinite(total) or abs(total-round(total))>.1:
            raise ValueError('As cargas parciais do ligante não somam uma carga molecular inteira; revise a protonação e as cargas.')
        write_text('\n'.join(' '.join(charges[i:i+8]) for i in range(0,len(charges),8))+'\n',directory/'existing.charges')
    env=dict(os.environ);env.pop('AMBERHOME',None);env['DOCK_HOME']=str(root)
    env['PATH']=str(root/'bin')+os.pathsep+env.get('PATH','')
    from .calculation_cache import CalculationCache,identity
    from .md_analysis import topology
    cache=CalculationCache('amber-preparation',directory)
    software=tool_identity(root,found)
    components=[]
    for prefix,helper,sources,extra in [('rec','amberize_receptor',['rec.pdb'],[]),
            ('lig.1','amberize_ligand',['lig.mol2'],['.gaff.mol2','.frcmod']),
            ('rec.lig.1','amberize_complex',['rec.pdb','lig.mol2'],[])]:
        key=identity(dict(component=prefix,inputs=[file_digest(directory/source) for source in sources],software=software,
                          existing_charges=existing_charges if prefix!='rec' else None))
        files={prefix+ext:directory/(prefix+ext) for ext in ['.prmtop','.inpcrd','.amber.pdb']+extra}
        with cache.locked(key):
            if cache.restore(key,files):
                write_text('#!/bin/sh\nexit 0\n',directory/helper)
                (directory/helper).chmod(0o700)
        components.append((prefix,key,files))
    try:
        run_command(['perl','prepare_amber.pl']+(['-c'] if existing_charges else [])+['lig.mol2','rec.pdb'],cwd=directory,env=env)
    except ScientificToolError as error:
        if error.error_code in ('TOOL_NOT_FOUND','TOOL_START_FAILED'):raise
        logs=[]
        for name in ('amberize_receptor.out','amberize_ligand.lig.1.out','amberize_complex.lig.1.out'):
            path=directory/name
            if path.is_file():logs.append(name+':\n'+path.read_text(errors='replace')[-3000:])
        raise ScientificToolError('Preparação AMBER falhou. Revise protonação, resíduos não padrão, cofatores, metais e parâmetros ausentes.\n'+'\n'.join(logs),
            'AMBER_PARAMETERIZATION_FAILED','Inspect the preserved Amber helper logs; provide compatible residue/cofactor/metal parameters and retry.') from error
    finally:
        # A successful receptor remains reusable even if ligand preparation fails.
        for prefix,key,files in components:
            if all(p.is_file() and p.stat().st_size for p in files.values()):
                try:topology(directory/(prefix+'.prmtop'))
                except ValueError:continue
                with cache.locked(key):cache.remember(key,files,dict(component=prefix))
    for prefix in ('rec','lig.1','rec.lig.1'):
        for extension in ('.prmtop','.inpcrd','.amber.pdb'):
            if not (directory/(prefix+extension)).is_file():raise ValueError('Preparação AMBER incompleta: '+prefix+extension)
    if not (directory/'lig.amber_score.mol2').is_file():raise ValueError('MOL2 sem identificador AMBER_SCORE_ID.')
    return env


def dock_input(settings):
    # Deliberately no orient_ligand: AMBER score refines a supplied pose.
    lines=dict(conformer_search_type='rigid',use_internal_energy='no',ligand_atom_file='lig.amber_score.mol2',
        limit_max_ligands='no',skip_molecule='no',read_mol_solvation='no',calculate_rmsd='no',
        use_database_filter='no',orient_ligand='no',bump_filter='no',score_molecules='yes')
    for name in ('contact','grid','multigrid','dock3.5','continuous','footprint_similarity','pharmacophore','hbond','internal_energy','descriptor','gbsa_zou','gbsa_hawkins','SASA'):
        lines[name+'_score_primary']='no'
    lines.update(amber_score_primary='yes',amber_score_receptor_file_prefix='rec')
    for name,value in settings.items():
        if name.startswith('amber_score_') and name!='amber_score_type':lines[name]=('yes' if value else 'no') if type(value) is bool else str(value)
    if settings['amber_score_movable_region']=='distance':lines['receptor_site_file']='site.sph'
    lines.update(ligand_outfile_prefix='amber',write_orientations='no',num_scored_conformers=1,rank_ligands='no')
    return '\n'.join(f'{name:<56} {value}' for name,value in lines.items())+'\n'


def energies(path):
    text=Path(path).read_text();result={}
    for source,target in [('Amber_Score','amber_score'),('Amber_complex_energy','complex_energy'),('Amber_receptor_energy','receptor_contribution'),('Amber_ligand_energy','ligand_contribution')]:
        match=re.search(r'\b'+source+r':\s*([-+0-9.eE]+)',text)
        if match is None or not math.isfinite(float(match[1])):raise ValueError('Energia AMBER ausente ou não finita: '+source)
        result[target]=float(match[1])
    if not math.isclose(result['amber_score'],sum(result[k] for k in ('complex_energy','receptor_contribution','ligand_contribution')),abs_tol=5e-5):raise ValueError('Termos AMBER incompatíveis com o score total.')
    return result


def plot(path,x,y,title,xlabel,ylabel):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    figure=Figure(figsize=(8,4.5));FigureCanvasAgg(figure);axis=figure.subplots()
    axis.plot(x,y);axis.set(title=title,xlabel=xlabel,ylabel=ylabel);axis.grid(alpha=.2)
    figure.tight_layout();figure.savefig(path,dpi=140)


def trajectory(directory,found,settings,env,ambertools_path='',ambertools_env=''):
    from .preparation_tools import run_amber_tool
    steps=settings['amber_score_md_steps'];stride=settings['trajectory_stride']
    mask=settings['md_movable_mask']
    movable=", ibelly=1, bellymask='"+mask+"'" if mask else ''
    cycles=settings['amber_score_before_md_minimization_cycles']
    source='rec.lig.1.final_pose.amber.restart'
    if not (directory/source).exists():source='rec.lig.1.inpcrd'
    if cycles:
        write_text(f"Complex minimization\n&cntrl imin=1, maxcyc={cycles}, ncyc={min(cycles,10)}, ntb=0, igb={settings['amber_score_gb_model']}, gbsa=1, cut={settings['amber_score_nonbonded_cutoff']}, ntpr=1{movable}, /\n",directory/'min.in')
        run_amber_tool('sander',['-O','-i','min.in','-p','rec.lig.1.prmtop','-c',source,'-r','min.rst7','-o','min.out'],cwd=directory,env=env,prefix=ambertools_path,environment=ambertools_env)
        source='min.rst7'
    common=(f"ntb=0, igb={settings['amber_score_gb_model']}, gbsa=1, "
            f"cut={settings['amber_score_nonbonded_cutoff']}, ntt=3, gamma_ln=2., "
            f"dt={settings['md_time_step']}, ntc=2, ntf=2{movable}")
    temperature=settings['amber_score_temperature']
    write_text(f"Implicit-solvent heating\n&cntrl imin=0, irest=0, ntx=1, {common}, "
               f"tempi=0., temp0={temperature}, nstlim={settings['md_heating_steps']}, "
               f"ig={settings['random_seed']}, ntpr={max(1,min(stride,settings['md_heating_steps']))}, ntwx=0, /\n",
               directory/'heat.in')
    run_amber_tool('sander',['-O','-i','heat.in','-p','rec.lig.1.prmtop','-c',source,
        '-r','heat.rst7','-o','heat.out'],cwd=directory,env=env,prefix=ambertools_path,environment=ambertools_env)
    write_text(f"Implicit-solvent equilibration\n&cntrl imin=0, irest=1, ntx=5, {common}, "
               f"tempi={temperature}, temp0={temperature}, nstlim={settings['md_equilibration_steps']}, "
               f"ig=-1, ntpr={max(1,min(stride,settings['md_equilibration_steps']))}, ntwx=0, /\n",
               directory/'equil.in')
    run_amber_tool('sander',['-O','-i','equil.in','-p','rec.lig.1.prmtop','-c','heat.rst7',
        '-r','equil.rst7','-o','equil.out'],cwd=directory,env=env,prefix=ambertools_path,environment=ambertools_env)
    write_text(f"Implicit-solvent production trajectory\n&cntrl imin=0, irest=1, ntx=5, {common}, "
               f"tempi={temperature}, temp0={temperature}, nstlim={steps}, ig=-1, "
               f"ntpr={stride}, ntwx={stride}, ioutfm=1, /\n",directory/'md.in')
    run_amber_tool('sander',['-O','-i','md.in','-p','rec.lig.1.prmtop','-c','equil.rst7',
        '-r','trajectory.rst7','-x','trajectory.nc','-o','md.out'],cwd=directory,env=env,prefix=ambertools_path,environment=ambertools_env)
    if not (directory/'trajectory.nc').is_file() or not (directory/'trajectory.nc').stat().st_size:raise ValueError('sander não produziu uma trajetória.')
    values=[]
    # Amber also prints averaged/RMS sections with NSTEP headers. They are
    # statistics, not additional frames or the final instantaneous energy.
    output=re.split(r'A\s+V\s+E\s+R\s+A\s+G\s+E\s+S',(directory/'md.out').read_text(),maxsplit=1)[0]
    for match in re.finditer(r'\bNSTEP\s*=\s*(\d+).*?TIME\(PS\)\s*=\s*(\S+).*?\bEPtot\s*=\s*(\S+)',output,re.S):
        step=int(match[1])
        time_value,energy=map(float,match.groups()[1:])
        if not all(math.isfinite(v) for v in (time_value,energy)):raise ValueError('Energias sander não finitas.')
        if step>0 and step%stride==0:values.append((time_value,energy))
    if not values:raise ValueError('sander não produziu energias potenciais válidas.')
    with (directory/'trajectory_energy.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['time_ps','potential_energy_kcal_mol']);writer.writerows(values)
    plot(directory/'trajectory_energy.png',[r[0] for r in values],[r[1] for r in values],'Energia potencial · trajetória sander','Tempo (ps)','kcal/mol')
    if settings['analyze_trajectory']:
        from .md_analysis import run
        run(directory,settings,found,env,ambertools_path,ambertools_env)
    return dict(sander_final_potential_energy=values[-1][1],sander_mean_potential_energy=sum(v[1] for v in values)/len(values))


def run_pair(pair,directory,root,found,settings,ambertools_path='',ambertools_env=''):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    def phase(name):
        write_json(dict(stage=name),directory/'progress.json')
        from .diagnostics import event,get_logger
        event(get_logger('scientific_tools'),'dynamics.phase','Etapa da dinâmica: '+name,
            compound_identifier=pair['ligand'].get('name'),configuration=pair.get('id'))
    phase('input_validation')
    receptor=validate_structure(pair['receptor']['path'],'receptor');shutil.copy2(receptor,directory/'rec.pdb')
    ligand=pair['ligand'];source=Path(ligand['path'])
    if source.suffix.lower()=='.csv':
        row=ligand.get('compound',{})
        pose=next((row.get(k) for k in ('pose_file','conformer_file','prepared_mol2','prepared_pdbqt') if row.get(k)),None)
        if not pose:raise ValueError('O composto não possui pose 3D vinculada ao receptor; execute o docking antes da dinâmica.')
        source=Path(pose)
    index=int(ligand.get('compound',{}).get('conformer_index') or 1)
    if index<1:raise ValueError('conformer_index deve ser um inteiro positivo.')
    if source.suffix.lower()!='.mol2' or index!=1 or source.read_text().count('@<TRIPOS>MOLECULE')>1:
        from .docking_data import convert_structure
        convert_structure(source,directory/'lig.mol2',index=index)
    else:shutil.copy2(validate_structure(source,'ligand'),directory/'lig.mol2')
    validate_structure(directory/'lig.mol2','ligand')
    phase('preparation_parameterization_complex')
    existing_charges=settings['use_existing_ligand_charges']
    charge_origin='input charges' if existing_charges else 'AM1-BCC'
    metadata=ligand.get('compound',{}).get('preparation_metadata')
    if metadata:
        try:metadata=json.loads(metadata)
        except (ValueError,TypeError):metadata={}
        if not isinstance(metadata,dict):metadata={}
        if metadata.get('charge_method')=='am1bcc' and isinstance(metadata.get('exports'),dict) and metadata['exports'].get('mol2')==file_digest(directory/'lig.mol2'):
            from .calculation_cache import CalculationCache
            checkpoint=CalculationCache('molecular-charges',directory)
            cached=checkpoint.restore(metadata.get('charge_cache_key',''),{'charges':directory/'prepared_charges.json'})
            if cached and cached.get('recipe',{}).get('method')=='am1bcc':
                import numpy as np
                block=(directory/'lig.mol2').read_text().split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
                actual=[float(line.split()[8]) for line in block.splitlines() if line.strip()]
                stored=json.loads((directory/'prepared_charges.json').read_text()).get('charges',[])
                if len(actual)==len(stored) and np.allclose(actual,stored,atol=1e-7,rtol=0):
                    existing_charges=True;charge_origin='verified prepared AM1-BCC'
    env=prepare(directory,root,found,existing_charges)
    from .md_analysis import selection
    topology,masks=selection(directory)
    write_json(dict(version=1,settings=settings,software=tool_identity(root,found),
        inputs={kind:dict(name=pair[kind].get('name',Path(pair[kind]['path']).name),identity=pair[kind].get('identity'),
            source=pair[kind]['path'],sha256=file_digest(directory/filename)) for kind,filename in [('receptor','rec.pdb'),('ligand','lig.mol2')]},
        preparation=dict(protein_force_field='ff94',ligand_force_field='GAFF',
            charge_method=charge_origin,
            solvent='implicit GB',water_model=None,pressure_control=None,complex_atom_order='ligand first, receptor second',
            atoms=topology['natom'],masks={k:v for k,v in masks.items() if isinstance(v,str)}),
        trajectory=dict(engine='sander' if settings['amber_score_type']=='dock6_sander' else None,
            equilibrium_claim=False,production_claim=False,thermostat='Langevin',gamma_ln=2.,
            constraints='SHAKE bonds involving hydrogen',duration_ps=settings['amber_score_md_steps']*settings['md_time_step'],
            phases=dict(minimization_cycles=settings['amber_score_before_md_minimization_cycles'],
                heating_steps=settings['md_heating_steps'],equilibration_steps=settings['md_equilibration_steps'],
                production_steps=settings['amber_score_md_steps']),
            warning='Short implicit-solvent refinement, not an equilibrated production simulation or binding free energy.')),
        directory/'protocol.json')
    if settings['amber_score_movable_region']=='distance':
        rows=[];active=False
        for line in (directory/'lig.mol2').read_text().splitlines():
            if line=='@<TRIPOS>ATOM':active=True;continue
            if line.startswith('@<TRIPOS>') and active:break
            if active and line.strip():
                row=line.split();rows.append([float(v) for v in row[2:5]])
        spheres='DOCK site from the selected ligand\ncluster     1   number of spheres in cluster '+str(len(rows))+'\n'
        spheres+=''.join(f'{i:5d}{x:10.5f}{y:10.5f}{z:10.5f}{1.:8.3f}{i:5d}{0:2d}{0:3d}\n' for i,(x,y,z) in enumerate(rows,1))
        write_text(spheres,directory/'site.sph')
    write_text(dock_input(settings),directory/'amber.in')
    phase('dock6_min_md_min')
    run_command([str(found['dock6']),'-i','amber.in','-o','amber.out','-v'],cwd=directory,env=env)
    result=energies(directory/'amber_scored.mol2')
    plot(directory/'energies.png',list(range(4)),list(result.values()),'AMBER Score · contribuições finais (DOCK6)','0: score · 1: complexo · 2: receptor · 3: ligante','kcal/mol')
    if settings['amber_score_type']=='dock6_sander':
        phase('sander_trajectory_analysis');result.update(trajectory(directory,found,settings,env,ambertools_path,ambertools_env))
    write_json(dict(units='kcal/mol',receptor_and_ligand_terms='Signed contributions reported by DOCK6',**result),directory/'energies.json')
    phase('completed')
    if os.environ.get('BIOMOL_MD_DEBUG','0')!='1':
        for name in ('prepare_amber.pl','amberize_receptor','amberize_ligand','amberize_complex','tleap','existing.charges','prepared_charges.json'):
            (directory/name).unlink(missing_ok=True)
    return result
