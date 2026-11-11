"""Explicit DOCK6 sampling policies; scoring never silently becomes a search."""
from pathlib import Path

MODES=('auto','guided','denovo','score_only','restrained','redock')


def mode_for_row(requested,row):
    coordinates=bool(row.get('conformer_file') or row.get('pose_file') or row.get('prepared_mol2'))
    bound=coordinates and (row.get('receptor_id') or row.get('engine') in ('vina','dock6') or
        row.get('prepared_origin')=='pose' or row.get('conformer_origin') in ('pose','crystal','control','docking'))
    return resolve_mode(requested,bool(bound) if requested=='auto' else coordinates)


def resolve_mode(requested, has_pose):
    if requested not in MODES:raise ValueError('Modo de operação DOCK6 inválido: '+str(requested))
    mode=('guided' if has_pose else 'denovo') if requested=='auto' else requested
    if mode in ('guided','score_only','restrained') and not has_pose:
        raise ValueError('O modo '+mode+' requer uma pose estrutural de entrada; SMILES não define coordenadas de referência.')
    return mode


def policy(mode):
    if mode not in MODES or mode=='auto':raise ValueError('Resolva o modo DOCK6 antes de gerar o arquivo de execução.')
    if mode in ('denovo','redock'):return {}
    return dict(conformer_search_type='rigid',orient_ligand='no',
        minimize_ligand='no' if mode=='score_only' else 'yes',
        use_advanced_simplex_parameters='no',num_scored_conformers='1',
        write_orientations='no',rank_ligands='no',
        calculate_rmsd='yes',use_rmsd_reference_mol='yes',
        use_internal_energy='no' if mode=='score_only' else 'yes',
        simplex_restraint_min='yes' if mode=='restrained' else 'no',
        simplex_coefficient_restraint='10.0')


def apply_policy(path, mode, *, seed=None):
    """Protect mode-defining switches, including customized legacy templates."""
    values=policy(mode)
    if not values:return
    if seed is not None:values['simplex_random_seed']=str(seed)
    path=Path(path);lines=[];seen=set()
    for line in path.read_text().splitlines():
        fields=line.split()
        if fields and fields[0] in values:
            name=fields[0]
            if name in seen:continue
            lines.append(f'{name:<58} {values[name]}');seen.add(name)
        else:lines.append(line)
    lines.extend(f'{name:<58} {value}' for name,value in values.items() if name not in seen)
    path.write_text('\n'.join(lines)+'\n')
