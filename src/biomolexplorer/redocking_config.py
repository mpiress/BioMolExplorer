"""Validated pair preparation and shared Chimera settings."""
import re
from pathlib import Path

DEFAULTS = {'remove_solvent': True, 'remove_hydrogens': True,
            'add_hydrogens': True, 'minimize': True, 'charge_type': 'gas'}


def pair_key(record):
    return '|'.join(str(record[i]) for i in range(4))


def validate_pairs(records, settings):
    if not isinstance(records, (list, tuple)) or not records:
        raise ValueError('Selecione pelo menos um par receptor / ligante e sua cadeia na aba Input Data.')
    if not isinstance(settings, dict):
        raise ValueError('Configurações de preparação inválidas.')
    seen = set(); receptors = {}
    for record in records:
        if not isinstance(record, (list, tuple)) or len(record) < 4 or any(
                not re.fullmatch(r'[A-Za-z0-9_+-]{1,40}', str(v)) for v in record[:4]):
            raise ValueError('Informe PDB, ligante, resíduo e cadeia válidos.')
        try: int(record[2])
        except (ValueError, TypeError): raise ValueError('Número do resíduo inválido.')
        key = pair_key(record)
        if key in seen: raise ValueError('O par receptor / ligante já está selecionado.')
        seen.add(key)
        config = settings.get(key, {})
        if not isinstance(config, dict): raise ValueError('Configuração do par inválida.')
        if config.keys() - {'cofactors', 'ligand_chain', 'receptor', 'ligand'}: raise ValueError('Configuração do par desconhecida.')
        cofactors = config.get('cofactors', [])
        if not isinstance(cofactors, list) or any(not isinstance(c, str) or not re.fullmatch(r'[A-Z0-9]{1,5}', c) for c in cofactors):
            raise ValueError('Informe códigos de cofatores válidos, como FAD, separados por vírgula.')
        if record[1] in cofactors: raise ValueError('O ligante de redocking não pode ser também um cofator.')
        # Older projects may contain two chain fields; CHAIN is now authoritative.
        for role in ('receptor', 'ligand'):
            options = config.get(role, {})
            if not isinstance(options, dict) or options.keys() - DEFAULTS.keys(): raise ValueError('Opções de preparação inválidas.')
            for option, value in options.items():
                if option == 'charge_type':
                    if value not in ('gas', 'am1'): raise ValueError('Método de cargas inválido.')
                elif type(value) is not bool: raise ValueError('Opção de preparação deve ser booleana.')
        # Legacy docking identifies receptors by PDB and chain; prevent silent overwrites.
        receptor_key = (record[0], record[3])
        receptor_config = (sorted(set(cofactors)), dict(DEFAULTS, **config.get('receptor', {})))
        if receptor_key in receptors and receptors[receptor_key] != receptor_config:
            raise ValueError('Pares do mesmo receptor e cadeia devem usar os mesmos cofatores e preparação do receptor.')
        receptors[receptor_key] = receptor_config
    if settings.keys() - seen: raise ValueError('Há configurações de pares que não foram selecionados.')


def configure_template(source, name, config, record=None):
    """Apply one ligand configuration to preparation and conformation alike."""
    role = 'ligand' if name in ('prepare_ligand.template', 'prepare_better_conform.template') else 'receptor'
    options = dict(DEFAULTS, **config.get(role, {}))
    cofactors = config.get('cofactors', [])
    suffix = ''.join(' | :' + c for c in cofactors)
    lines = []
    for line in source.splitlines():
        stripped = line.strip()
        if name == 'prepare_complex.template' and stripped.startswith('select #0:'):
            line = 'select #0:.{chain}' + suffix
        if name == 'prepare_receptor.template':
            if stripped == 'delete ligand': continue
            if stripped == 'select protein': line = 'select protein' + suffix + (' | solvent' if not options['remove_solvent'] else '')
        if name == 'prepare_ligand.template' and stripped.startswith('select :{resnum}') and not options['remove_solvent']:
            line += ' | solvent'
        if name == 'prepare_complex.template' and stripped in ('delete solvent', 'delete element.H'): continue
        if stripped == 'delete solvent' and not options['remove_solvent']: continue
        if stripped == 'delete element.H' and not options['remove_hydrogens'] and name != 'prepare_receptor.template': continue
        if stripped == 'addh' and not options['add_hydrogens']: continue
        if stripped.startswith('minimize') and not options['minimize']: continue
        if stripped.startswith('addcharge'):
            line = re.sub(r'method (gas|am1|\{charge_type\})', 'method ' + options['charge_type'], line)
        lines.append(line)
    if name != 'prepare_complex.template':
        index = next((i for i, line in enumerate(lines) if line.startswith('addcharge')), 1)
        if options['add_hydrogens'] and 'addh' not in lines: lines.insert(index, 'addh')
        if options['minimize'] and not any(line.startswith('minimize') for line in lines):
            index = next((i for i, line in enumerate(lines) if line.startswith('write')), len(lines))
            lines.insert(index, 'minimize spec #0')
    if name in ('prepare_ligand.template', 'prepare_receptor.template') and options['remove_hydrogens'] and (name == 'prepare_receptor.template' or 'delete element.H' not in lines):
        index = next((i for i, line in enumerate(lines) if line == 'addh' or line.startswith('addcharge')), 1)
        lines.insert(index, 'delete element.H')
    # Keep solvent removal local to the selected ligand/receptor.
    if name != 'prepare_complex.template' and options['remove_solvent'] and 'delete solvent' not in lines:
        index = next((i for i, line in enumerate(lines) if line == 'addh' or line.startswith('addcharge')), 1)
        lines.insert(index, 'delete solvent')
    return '\n'.join(lines) + '\n'


def metadata_records(paths):
    """Read available ligand tuples and retain resolution without a UI field."""
    import csv
    records = []
    for path in paths:
        path = Path(path)
        if path.name != 'pdb_codes.csv' or not path.is_file(): continue
        with path.open(encoding='utf-8-sig', newline='') as stream:
            for row in csv.DictReader(stream):
                try:
                    record = [row['PDB_CODE'], row['LIGAND'], int(row['RESNUM']), row['CHAIN']]
                    resolution = row.get('RESOLUTION', '')
                    if resolution and resolution.lower() != 'nan': record.append(float(resolution))
                except (KeyError, ValueError): continue
                if record not in records: records.append(record)
    return records


def validate_structure_pairs(folder, records, settings, prepared=False):
    """Validate actual selected residues before running any external engine."""
    import csv
    import math
    folder = Path(folder)
    validate_pairs(records, settings)
    if prepared:
        centers_file = folder / 'Prepared' / 'centers.csv'
        if not centers_file.is_file():
            raise ValueError('Os complexos preparados precisam do arquivo centers.csv.')
        with centers_file.open(encoding='utf-8-sig', newline='') as stream:
            centers = list(csv.DictReader(stream))
        for record in records:
            code, ligand, number, chain = record[:4]
            identity = f'{code}_{ligand}_{number}{chain}'
            for filename in (f'{code}_{chain}.dockprep.pdbqt', f'{identity}.lig.pdbqt'):
                path = folder / 'Prepared' / filename
                if not path.is_file() or not any(line.startswith(('ATOM  ', 'HETATM')) for line in path.read_text().splitlines()):
                    raise ValueError(f'Arquivo preparado ausente ou sem átomos: {filename}')
            column = identity if centers and identity in centers[0] else f'{code}_{ligand}_{number}_{chain}'
            try:
                coordinates = [float(row[column]) for row in centers]
                if len(coordinates) != 3 or not all(math.isfinite(v) for v in coordinates): raise ValueError()
            except (KeyError, ValueError, TypeError):
                raise ValueError(f'Centro do ligante ausente ou inválido: {identity}') from None
        return
    structures = {}
    for record in records:
        code, ligand, number, chain = record[:4]
        if code not in structures:
            path = folder / f'{code}.pdb'
            if not path.is_file(): raise ValueError(f'Estrutura PDB ausente: {code}')
            residues = set(); protein_chains = set()
            with path.open() as stream:
                for line in stream:
                    if line.startswith('ENDMDL'): break
                    if not line.startswith(('ATOM  ', 'HETATM')): continue
                    try:
                        residue = (line[17:20].strip(), int(line[22:26]), line[21:22].strip())
                        coordinates = [float(line[start:start+8]) for start in (30,38,46)]
                    except ValueError: raise ValueError(f'Estrutura PDB contém átomos inválidos: {code}') from None
                    if not all(math.isfinite(v) for v in coordinates): raise ValueError(f'Estrutura PDB contém átomos inválidos: {code}')
                    residues.add(residue)
                    if line.startswith('ATOM  ') and residue[0] != ligand: protein_chains.add(residue[2])
            structures[code] = residues, protein_chains
        residues, protein_chains = structures[code]
        if chain not in protein_chains:
            raise ValueError(f'Cadeia do receptor não encontrada no PDB {code}: {chain}')
        if (ligand, int(number), chain) not in residues:
            raise ValueError(f'Ligante ou resíduo não encontrado na cadeia selecionada: {code} / {ligand} / {number} / {chain}')
        for cofactor in settings.get(pair_key(record), {}).get('cofactors', []):
            if not any(name == cofactor for name, _, _ in residues):
                raise ValueError(f'Cofator não encontrado no PDB {code}: {cofactor}')


def validate_redocking_tools(prepare=True, search_path=None):
    """Use the same executable search path as the scientific worker."""
    import os
    import shutil
    import sys
    search_path = search_path or str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')
    for executable in (('chimera', 'obabel', 'vina') if prepare else ('vina',)):
        if not shutil.which(executable, path=search_path):
            raise ValueError(f'Ferramenta necessária ao redocking não encontrada: {executable}. Configure a instalação antes de executar.')
