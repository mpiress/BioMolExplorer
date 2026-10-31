"""Exact, project-scoped receptor/ligand choices, with stable pair identities."""
import hashlib
import json
import math
from pathlib import Path

from .artifacts import records,validate_reference,locator
from .bindings import sources


def molecular_identity(source):
    """Hash selected molecular content, never unrelated rows or snapshot paths."""
    from .stage_cache import file_digest
    row=source.get('compound')
    if row is None:return file_digest(source['path'])
    pose=next((row.get(k) for k in ('pose_file','conformer_file','prepared_mol2','prepared_pdbqt') if row.get(k)),None)
    try:preparation=json.loads(row.get('preparation_metadata') or '{}')
    except (TypeError,ValueError):preparation={}
    if not isinstance(preparation,dict):preparation={}
    return dict(compound=row.get('molecule_chembl_id'),pose=file_digest(pose) if pose else None,index=row.get('conformer_index',1),
                preparation={k:preparation.get(k) for k in ('charge_method','charge_cache_key','exports')},
                smiles=row.get('canonical_smiles') if not pose else None)


def pair_id(receptor,ligand):
    return hashlib.sha256(json.dumps([receptor['identity'],ligand['identity']],separators=(',',':')).encode()).hexdigest()


def validate_structure(path,kind):
    path=Path(path)
    if not path.is_file():raise ValueError('Arquivo molecular não encontrado: '+path.name)
    if kind=='receptor':
        from .input_validation import validate_file
        if path.suffix.lower()!='.pdb':raise ValueError('O receptor da dinâmica deve ser um arquivo PDB, não PDBQT.')
        validate_file(path,'structures')
        from .molecular_roles import macromolecule_pdb
        if not macromolecule_pdb(path):raise ValueError('A entrada de receptores requer uma estrutura de macromolécula; este PDB não contém resíduos de proteína ou ácido nucleico.')
    else:
        if path.suffix.lower()!='.mol2':raise ValueError('O ligante da dinâmica deve ser uma pose MOL2 com coordenadas.')
        from .molecular_roles import small_molecule_file
        if not small_molecule_file(path):raise ValueError('A entrada de ligantes requer uma molécula pequena, não uma estrutura de macromolécula.')
        text=path.read_text();lines=text.splitlines()
        if text.count('@<TRIPOS>MOLECULE')!=1:raise ValueError('Selecione um MOL2 contendo exatamente um ligante por arquivo.')
        try:
            start=lines.index('@<TRIPOS>ATOM')+1;atoms=[]
            for line in lines[start:]:
                if line.startswith('@<TRIPOS>'):break
                if not line.strip():continue
                row=line.split()
                if len(row)<9 or any(not math.isfinite(float(row[i])) for i in (2,3,4,8)):raise ValueError()
                atoms.append(row)
            if not atoms or '@<TRIPOS>BOND' not in lines:raise ValueError()
            counts=lines[lines.index('@<TRIPOS>MOLECULE')+2].split()
            if int(counts[0])!=len(atoms):raise ValueError()
        except (ValueError,IndexError):raise ValueError('MOL2 incompleto: confira átomos, coordenadas, tipos e cargas finitas.') from None
    return path


class MDInputs:
    def __init__(self,store,pid,results):self.store,self.pid,self.results=store,pid,results

    def expand(self,candidates,kind,suffixes,compound_id=None):
        from .docking_data import read_compounds
        result=[]
        for candidate in candidates:
            path=self.store.scoped_path(self.pid,candidate['path'])
            if path.suffix.lower() not in suffixes:continue
            if kind!='ligand' or path.suffix.lower()!='.csv':
                result.append(candidate);continue
            for row in read_compounds(path):
                code=row.get('molecule_chembl_id')
                if not code or compound_id and code!=compound_id:continue
                for name in ('conformer_file','pose_file','prepared_mol2','prepared_pdbqt','dock6_pose','vina_pose'):
                    if row.get(name):
                        pose=Path(row[name]);row[name]=str(self.store.scoped_path(self.pid,pose if pose.is_absolute() else path.parent/pose))
                if not any(row.get(k) for k in ('conformer_file','pose_file','prepared_mol2','prepared_pdbqt')) and any(row.get(k) for k in ('dock6_pose','vina_pose')):
                    for engine in ('dock6','vina'):
                        pose=row.get(engine+'_pose')
                        if not pose:continue
                        origin=str(row.get('receptor_id',''));index=str(row.get('conformer_index',1))
                        result.append(dict(candidate,name=code+' · '+engine+(' · '+origin if origin else ''),
                            identity=candidate['identity']+':compound:'+code+':'+origin+':'+index+':'+engine,
                            compound=dict(row,conformer_file=pose,engine=engine)))
                    continue
                origin=str(row.get('receptor_id',''));index=str(row.get('conformer_index',1))
                result.append(dict(candidate,name=code+(' · '+origin if origin else ''),
                    identity=candidate['identity']+':compound:'+code+':'+origin+':'+index,compound=row))
        return result

    def files(self,stage,field):
        kind='receptor' if field=='base_input_path' else 'ligand'
        suffixes={'.pdb'} if kind=='receptor' else {'.mol2','.pdbqt','.sdf','.csv'}
        result=[]
        for ref in sources(stage.get('bindings',{}).get(field,{})):
            if 'asset' in ref:
                with self.store.connect() as db:
                    asset=db.execute('SELECT path,name FROM assets WHERE project_id=? AND id=?',(self.pid,ref['asset'])).fetchone()
                if asset is None:raise ValueError('Arquivo molecular não autorizado.')
                candidates=[dict(path=asset['path'],name=asset['name'],identity='asset:'+ref['asset'])]
            elif 'stage' in ref:
                paths=self.results.get(ref['stage'],[])
                if ref.get('selector','auto')!='auto':
                    selected=validate_reference(ref,paths)
                    paths=[selected]
                elif kind=='ligand':
                    # Automatic docking outputs use final poses, never receptor
                    # MOL2s or intermediate copies of an already scored ligand.
                    from .docking_inputs import ligand_input_file
                    paths=[p for p in paths if ligand_input_file(p)]
                    from .artifact_choices import consolidated_docking_tables
                    tables=consolidated_docking_tables(paths) or [p for p in paths if Path(p).name=='docking_results.csv']
                    scored=[p for p in paths if Path(p).name.endswith('_scored.mol2')]
                    if tables:paths=tables
                    elif scored:paths=scored
                else:
                    paths=receptor_input_files(paths)
                candidates=[dict(path=str(p),name=Path(p).name,
                    identity='stage:'+ref['stage']+':'+locator(ref['stage'],p)) for p in paths]
            else:
                # The editor keeps an empty reference while a source is being
                # selected.  Treat that transient placeholder as no input.
                continue
            result.extend(self.expand(candidates,kind,suffixes,ref.get('compound_id')))
        direct=stage['parameters'].get(field)
        if not stage.get('bindings',{}).get(field) and direct:
            folder=self.store.scoped_path(self.pid,direct)
            paths=[folder] if folder.is_file() else sorted(p for p in folder.rglob('*') if p.suffix.lower() in suffixes)
            from .artifact_choices import consolidated_docking_tables
            if kind=='receptor':paths=receptor_input_files(paths)
            tables=consolidated_docking_tables(paths) or [p for p in paths if Path(p).name=='docking_results.csv']
            if kind=='ligand' and tables:paths=tables
            candidates=[dict(path=str(p),name=p.name,identity='path:'+str(p.relative_to(self.store.project_dir(self.pid)))) for p in paths]
            result.extend(self.expand(candidates,kind,suffixes))
        unique={c['identity']:c for c in result}
        return [dict(c,path=str(self.store.scoped_path(self.pid,c['path']))) for c in unique.values()]

    def resolve(self,stage,item=None,check_files=True):
        from .md_config import validate_settings
        params=dict(stage['parameters']);validate_settings(params)
        receptors=self.files(stage,'base_input_path');ligands=self.files(stage,'base_selected_mols')
        if not receptors or not ligands:raise ValueError('Selecione ao menos um receptor PDB e um ligante MOL2 para a dinâmica.')
        if len(receptors)*len(ligands)>100000:raise ValueError('Selecione menos de 100000 combinações para a dinâmica.')
        selected=set(params.get('selected_pairs') or [])
        if not selected:raise ValueError('Selecione explicitamente os pares receptor–ligante antes de executar a dinâmica.')
        excluded=set(params.get('excluded_pairs') or [])
        pairs=[]
        for receptor in receptors:
            for ligand in ligands:
                key=pair_id(receptor,ligand)
                if key in excluded or key not in selected:continue
                # Scientific validation runs per pair in the worker so that a
                # malformed molecule cannot prevent unrelated pairs completing.
                pairs.append(dict(id=key,receptor=receptor,ligand=ligand))
        if not pairs:raise ValueError('Nenhum par ativo de receptor e ligante foi selecionado.')
        if item is not None:item['input_files']={field:[c['path'] for c in values] for field,values in (('base_input_path',receptors),('base_selected_mols',ligands))}
        params.update(md_inputs=pairs,base_input_path=str(Path(receptors[0]['path']).parent),base_selected_mols=str(Path(ligands[0]['path']).parent))
        return params


def normalize_receptor_variants(source,destination):
    """Disambiguate HIS only from explicit hydrogen names; preserve coordinates."""
    lines=Path(source).read_text().splitlines(keepends=True);atoms={}
    for line in lines:
        if line.startswith(('ATOM  ','HETATM')) and line[17:20]=='HIS':atoms.setdefault(line[21:27],set()).add(line[12:16].strip())
    variants={key:('HIP' if {'HD1','HE2'}<=names else 'HID' if 'HD1' in names else 'HIE') for key,names in atoms.items() if {'HD1','HE2'}&names}
    Path(destination).write_text(''.join(line[:17]+variants[line[21:27]]+line[20:] if line.startswith(('ATOM  ','HETATM')) and line[17:20]=='HIS' and line[21:27] in variants else line for line in lines))
    return [dict(residue=key.strip(),input='HIS',output=value,evidence='Explicit HD1/HE2 hydrogen names; coordinates preserved') for key,value in variants.items()]


def receptor_input_files(paths):
    """Prefer a prepared receptor per source and chain without hiding other targets."""
    from .molecular_roles import macromolecule_pdb
    paths=list(dict.fromkeys(Path(p) for p in paths));selected={};raw=[]
    endings=('.final_pose.amber.pdb','.prepared.pdb','.noh.pdb','.amber.pdb')
    for path in paths:
        name=path.name.casefold()
        if path.suffix.lower()!='.pdb' or '.lig.' in name or '.complex.' in name or name.startswith('lig.'):continue
        if not macromolecule_pdb(path):continue
        parent=path.parent.parent if path.parent.name=='Prepared' else path.parent
        match=next(((i,suffix) for i,suffix in enumerate(endings) if name.endswith(suffix)),None)
        if match:
            priority,suffix=match;prefix=name[:-len(suffix)];key=(parent,prefix)
            if key not in selected or priority<selected[key][0]:selected[key]=(priority,path)
        elif path.name.count('.')==1:raw.append(path)
    covered={(parent,'_'.join(prefix.split('_')[:i])) for parent,prefix in selected for i in range(1,len(prefix.split('_'))+1)}
    return sorted([value[1] for value in selected.values()]+[p for p in raw if (p.parent,p.stem.casefold()) not in covered])
