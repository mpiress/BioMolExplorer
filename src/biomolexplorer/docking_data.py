"""Compound identities, pose-preserving preparation and interoperable scores."""
import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

RESULT_FILE = 'docking_results.csv'
NO_INTERSECTION = 'O consenso não foi executado: Vina e DOCK6 não avaliaram compostos em comum para o mesmo receptor. Selecione resultados com códigos de molécula e receptores correspondentes.'


def result_tables(paths):
    """Use each batch summary, excluding its redundant per-receptor tables.

    An empty curated summary remains authoritative: deleted rows must never
    reappear from a nested table or from the original pose files.
    """
    tables=list(dict.fromkeys(Path(p) for p in paths))
    return [p for p in tables if not any(other != p and other.parent in p.parent.parents for other in tables)]


def read_compounds(path):
    from .input_validation import ALIASES
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return [{ALIASES.get(k,k):v for k,v in row.items()} for row in csv.DictReader(stream)]


def convert_structure(source, destination, *, ph=None, charges=None, index=1):
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    command = ['obabel', '-i'+source.suffix[1:], str(source), '-o'+destination.suffix[1:],
               '-O', str(destination), '-f', str(index), '-l', str(index)]
    if ph is not None: command.extend(['-p', str(ph)])
    if charges: command.extend(['--partialcharge', charges])
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    if result.returncode or not destination.is_file() or not destination.stat().st_size:
        raise ValueError('Não foi possível converter a conformação molecular: '+result.stderr[-1000:])
    return destination


def structure_smiles(path):
    result = subprocess.run(['obabel', '-i'+Path(path).suffix[1:], str(path), '-ocan', '-f', '1', '-l', '1'],
                            capture_output=True, text=True, timeout=120)
    if result.returncode or not result.stdout.strip():
        raise ValueError('Não foi possível obter o SMILES da conformação molecular.')
    return result.stdout.split()[0]


def binding_center(prepared, record):
    import pandas as pd
    centers=pd.read_csv(Path(prepared)/'centers.csv')
    native=f'{record[0]}_{record[1]}_{record[2]}{record[3]}'
    legacy=f'{record[0]}_{record[1]}_{record[2]}_{record[3]}'
    values=centers[native if native in centers else legacy].tolist()
    if len(values)!=3 or not all(math.isfinite(float(v)) for v in values):
        raise ValueError('Centro do sítio de docking ausente ou inválido.')
    return list(map(float,values))


def center_mol2(path, center):
    """Place a newly generated ligand in the selected site, retaining charges/types."""
    lines=Path(path).read_text().splitlines(keepends=True);atoms=[];active=False
    for i,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):active=line.strip()=='@<TRIPOS>ATOM';continue
        if active and line.strip():atoms.append((i,line.split()))
    if not atoms:raise ValueError('O ligante MOL2 não contém átomos.')
    centroid=[sum(float(a[1][k]) for a in atoms)/len(atoms) for k in (2,3,4)]
    for i,fields in atoms:
        fields[2:5]=[f'{float(fields[k+2])+center[k]-centroid[k]:.4f}' for k in range(3)]
        lines[i]=' '.join(fields)+'\n'
    Path(path).write_text(''.join(lines))


def prepare_ligands(dataset, destination, format, ph=7.4, charge_type='gas', center=None, preparation_options=None):
    """Convert an existing pose, or generate 3D only when the input is SMILES."""
    from .visualizations import molecule_sdf
    dataset, destination = Path(dataset), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if preparation_options is not None:
        from .docking_preparation import prepare_candidates
        dataset=prepare_candidates(dataset,destination/'Prepared',preparation_options,ph,
                                   'vina' if format=='pdbqt' else 'dock6')
    outputs=[]
    for row in read_compounds(dataset):
        code=row['molecule_chembl_id']
        if not re.fullmatch(r'[A-Za-z0-9_.+-]{1,100}', code) or '..' in code:
            raise ValueError('Código molecular inválido para o docking.')
        prepared=row.get('prepared_'+format)
        if row.get('docking_engines') and not prepared:
            raise ValueError('A saída de preparação não inclui o formato solicitado pelo docking.')
        if prepared:
            source=Path(prepared)
            if not source.is_absolute():source=dataset.parent/source
            output=destination/(code+'.lig.'+format)
            shutil.copy2(source,output)
            if format=='mol2' and center is not None and row.get('prepared_origin') in ('smiles','library'):center_mol2(output,center)
            outputs.append(output)
            continue
        pose=row.get('conformer_file') or row.get('pose_file')
        if pose:
            source=Path(pose)
            if not source.is_absolute():source=dataset.parent/source
            if not source.is_file():raise ValueError('A conformação selecionada não está disponível: '+code)
        else:
            source=destination/(code+'.sdf')
            source.write_text(molecule_sdf(row['canonical_smiles']),encoding='utf-8')
        output=destination/(code+'.lig.'+format)
        convert_structure(source,output,ph=ph,charges='gasteiger' if format=='mol2' and charge_type=='gas' else None,
                          index=int(row.get('conformer_index') or 1))
        if (not pose or row.get('conformer_origin')=='library') and format=='mol2' and center is not None:center_mol2(output,center)
        if format=='mol2' and charge_type=='am1':
            script=destination/(code+'.charges.com')
            script.write_text(f'open {output.name}\naddcharge all spec #0 chargeModel 14sb method am1\nwrite format mol2 #0 {output.name}\nclose all\n',encoding='utf-8')
            result=subprocess.run(['chimera','--nogui','--silent',script.name],cwd=destination,capture_output=True,text=True,timeout=300)
            if result.returncode:raise ValueError('Não foi possível atribuir cargas AM1: '+result.stderr[-1000:])
        outputs.append(output)
    return outputs


def input_records(selected, available):
    """Keep identities from result tables when individual poses are selected."""
    allowed={Path(p).resolve() for p in available}
    metadata={}
    def table_rows(path):
        rows=read_compounds(path)
        for row in rows:
            for key in ('conformer_file','pose_file','prepared_pdbqt','prepared_mol2'):
                pose=row.get(key)
                if not pose:continue
                resolved=(Path(pose) if Path(pose).is_absolute() else Path(path).parent/pose).resolve()
                if resolved not in allowed:raise ValueError('A conformação não pertence aos arquivos autorizados da entrada.')
                row[key]=str(resolved)
        return rows
    from .input_validation import columns
    for path in available:
        if Path(path).suffix.lower()=='.csv' and {'molecule_chembl_id','canonical_smiles'}<=columns(path):
            rows=read_compounds(path)
            for row in rows:
                for key in ('conformer_file','pose_file','vina_pose','dock6_pose','prepared_pdbqt','prepared_mol2'):
                    pose=row.get(key)
                    if not pose:continue
                    resolved=(Path(pose) if Path(pose).is_absolute() else Path(path).parent/pose).resolve()
                    if resolved not in allowed:raise ValueError('A conformação não pertence aos arquivos autorizados da entrada.')
                    metadata[resolved]=dict(row,conformer_file=str(resolved))
    result=[]
    for path in map(Path,selected):
        if path.suffix=='.csv':result.extend(table_rows(path));continue
        if path.resolve() in metadata:
            result.append(dict(metadata[path.resolve()]));continue
        output=subprocess.run(['obabel','-i'+path.suffix[1:],str(path),'-ocan'],capture_output=True,text=True,timeout=120)
        if output.returncode or not output.stdout.strip():raise ValueError('Estrutura molecular inválida: '+path.name)
        lines=output.stdout.strip().splitlines()
        for index,line in enumerate(lines,1):
            smiles,*name=line.split()
            code=path.name.removesuffix('_scored.mol2').removesuffix(path.suffix).removesuffix('.lig')
            if len(lines)>1:code=(name[0] if name and re.fullmatch(r'[A-Za-z0-9_+-]{1,100}',name[0]) else code+'_'+str(index))
            result.append(dict(molecule_chembl_id=code,canonical_smiles=smiles,conformer_file=str(path.resolve()),conformer_index=index))
    return result


def score_input_records(selected, available, engine):
    """Scores define consensus inputs; visualization files are optional metadata."""
    allowed={Path(p).resolve() for p in available}
    rows=[]
    for path in map(Path,selected):
        if path.suffix.lower()=='.csv':
            records=read_compounds(path)
        else:
            code=path.name.removesuffix('_scored.mol2').removesuffix('.pdbqt').removesuffix('.lig')
            records=[dict(molecule_chembl_id=code,canonical_smiles=structure_smiles(path),
                          score=score(path,engine),conformer_file=str(path.resolve()))]
        for record in records:
            row=dict(record)
            if row.get('engine') and row['engine']!=engine:
                raise ValueError('Resultado de docking incompatível com o motor selecionado.')
            if not row.get('molecule_chembl_id'):raise ValueError('Código do composto ausente no resultado.')
            try:value=float(row['score'])
            except (KeyError,ValueError,TypeError):raise ValueError('Score ausente ou inválido no resultado: '+engine) from None
            if not math.isfinite(value):raise ValueError('Score não finito no resultado: '+engine)
            row.update(engine=engine,score=value)
            row.setdefault('canonical_smiles','');row.setdefault('receptor_id','')
            for field in ('conformer_file','pose_file','prepared_pdbqt','prepared_mol2','receptor_file','reference_file','footprint_file'):
                value=row.get(field)
                if not value:continue
                file=(Path(value) if Path(value).is_absolute() else path.parent/value).resolve()
                if file in allowed and file.is_file():row[field]=str(file)
                else:row.pop(field,None)
            rows.append(row)
    return rows


def copy_input_records(rows, destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    result=[]
    for index,row in enumerate(rows):
        row=dict(row)
        for key in ('conformer_file','pose_file','prepared_pdbqt','prepared_mol2'):
            if not row.get(key):continue
            source=Path(row[key]);pose=destination/'poses'/f'{index}_{key}{source.suffix}'
            pose.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,pose)
            # Absolute project-scoped paths survive quality-input snapshots.
            row[key]=str(pose.resolve())
        result.append(row)
    fields=list(dict.fromkeys(k for row in result for k in row)) or [
        'molecule_chembl_id','canonical_smiles','receptor_id','engine','score','conformer_file']
    write_csv(destination/'compounds.csv',result,fields)
    return destination/'compounds.csv'


def score(path, engine):
    label='REMARK VINA RESULT' if engine=='vina' else 'Grid_Score'
    match=re.search(re.escape(label)+r':\s*(\S+)',Path(path).read_text())
    if not match:raise ValueError('Score ausente no resultado: '+label)
    value=float(match.group(1))
    if not math.isfinite(value):raise ValueError('Score não finito no resultado: '+label)
    return value


def write_results(engine, poses, dataset, records, output, prepared=None):
    """Retain input SMILES/IDs and link to the actual docked pose, not a new conformer."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    compounds=read_compounds(dataset);rows=[]
    for pose in sorted(map(Path,poses)):
        name=pose.name.removesuffix('_scored.mol2').removesuffix('.pdbqt').removesuffix('.lig')
        candidates=[r for r in compounds if name==r['molecule_chembl_id'] or name.endswith('_'+r['molecule_chembl_id'])]
        if not candidates:continue
        row=max(candidates,key=lambda r:len(r['molecule_chembl_id']))
        matches=[r for r in records if name.startswith(f'{r[0]}_{r[1]}_{r[2]}{r[3]}_') or name.startswith(f'{r[0]}_{r[3]}_')]
        record=matches[0] if matches else records[0] if len(records)==1 else None
        if record is None:raise ValueError('A pose não identifica um receptor único: '+pose.name)
        rows.append(dict(molecule_chembl_id=row['molecule_chembl_id'],canonical_smiles=row['canonical_smiles'],
                         receptor_id=f'{record[0]}_{record[3]}',engine=engine,score=score(pose,engine),
                         conformer_file=pose.resolve().relative_to(output.resolve()).as_posix()))
    if prepared is not None:
        import shutil
        prepared=Path(prepared)
        context=output/'receptors';context.mkdir(exist_ok=True)
        for row in rows:
            rid=row['receptor_id']
            source=next((prepared/(rid+suffix) for suffix in ('.noH.pdb','.dockprep.pdbqt','.dockprep.mol2') if (prepared/(rid+suffix)).is_file()),None)
            if source:
                target=context/source.name;shutil.copy2(source,target)
                row['receptor_file']=target.relative_to(output).as_posix()
            record=next(r for r in records if f'{r[0]}_{r[3]}'==rid)
            reference=prepared/f'{record[0]}_{record[1]}_{record[2]}{record[3]}.lig.pdb'
            if reference.is_file():
                target=context/reference.name;shutil.copy2(reference,target)
                row['reference_file']=target.relative_to(output).as_posix()
            if engine=='dock6':
                pdf=output/'Dock6/footprint/plots'/(row['molecule_chembl_id']+'.pdf')
                if pdf.is_file():
                    row['footprint_file']=pdf.relative_to(output).as_posix()
                    row['footprint_origin']='docked_pose'
    if not rows:raise ValueError('O docking não produziu conformações com scores válidos.')
    write_csv(output/RESULT_FILE,rows,list(dict.fromkeys(k for row in rows for k in row)))
    return rows


def write_csv(path, rows, fields=None):
    fields=fields or list(rows[0])
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def read_results(directory, engine):
    directory=Path(directory);rows=[]
    tables=result_tables(directory.rglob(RESULT_FILE))
    for path in tables:
        for row in read_compounds(path):
            if row.get('engine') and row['engine']!=engine:continue
            row['engine']=engine
            if row.get('conformer_file'):
                pose=Path(row['conformer_file'])
                row['conformer_file']=str(pose if pose.is_absolute() else path.parent/pose)
            row.setdefault('canonical_smiles','');row.setdefault('receptor_id','')
            for field in ('receptor_file','reference_file','footprint_file'):
                if row.get(field):
                    target=Path(row[field]);row[field]=str(target if target.is_absolute() else path.parent/target)
            value=float(row['score'])
            if not math.isfinite(value):raise ValueError('Score não finito no resultado: '+engine)
            row['score']=value;rows.append(row)
    if tables:return rows
    # Compatibility with completed results supplied as individual pose files.
    pattern='*.pdbqt' if engine=='vina' else '*_scored.mol2'
    for path in sorted(directory.rglob(pattern)):
        code=path.name.removesuffix('_scored.mol2').removesuffix('.pdbqt').removesuffix('.lig')
        rows.append(dict(molecule_chembl_id=code,canonical_smiles=structure_smiles(path),receptor_id='',engine=engine,
                         score=score(path,engine),conformer_file=str(path)))
    return rows


def consensus_rows(vina, dock6, weight=1.):
    def best(rows):
        result={}
        for row in rows:
            key=(row['receptor_id'],row['molecule_chembl_id'])
            if key not in result or row['score']<result[key]['score']:result[key]=row
        return result
    first,second=best(vina),best(dock6);rows=[]
    for key in sorted(first.keys() & second.keys()):
        a,b=first[key],second[key]
        if a['canonical_smiles'] and b['canonical_smiles']:
            from rdkit import Chem
            if Chem.MolToSmiles(Chem.MolFromSmiles(a['canonical_smiles']))!=Chem.MolToSmiles(Chem.MolFromSmiles(b['canonical_smiles'])):
                raise ValueError('O mesmo código molecular identifica estruturas diferentes: '+key[1])
        rows.append(dict(molecule_chembl_id=key[1],molecule=key[1],canonical_smiles=a['canonical_smiles'] or b['canonical_smiles'],
                         receptor_id=key[0],vina=a['score'],dock6=b['score']))
        for field,source in (('vina_pose',a),('dock6_pose',b)):
            if source.get('conformer_file'):rows[-1][field]=source['conformer_file']
        for engine,source in (('vina',a),('dock6',b)):
            for field in ('receptor_file','reference_file'):
                if source.get(field):rows[-1][engine+'_'+field]=source[field]
        for field in ('footprint_file','footprint_origin'):
            if b.get(field):rows[-1][field]=b[field]
    return rows
