"""Verified, immutable per-pair reuse using the existing artifact hash contracts."""
import json
from pathlib import Path
import shutil

from .stage_cache import artifact_manifest,file_digest
from .storage import write_json
from .calculation_cache import CalculationCache

NAME='md_cache.json'
SUFFIXES={'.csv','.png','.json','.pdb','.pdbqt','.mol2','.sdf','.pdf','.prmtop','.inpcrd','.nc','.rst7','.restart','.out','.dat','.in','.log','.frcmod'}


def seed_runs(store,pid,sid,destination):
    """Seed only intact successful pairs, including runs with unrelated failures."""
    with store.connect() as db:
        rows=db.execute('SELECT stages FROM runs WHERE project_id=? ORDER BY created DESC',(pid,)).fetchall()
    count=0;seen=set();destination=Path(destination)
    for row in rows:
        for item in json.loads(row[0]):
            if item['id']!=sid or item['status']!='succeeded' or item.get('outputs_removed'):continue
            for name in item.get('artifacts',[]):
                if Path(name).name!=NAME:continue
                try:
                    index=store.scoped_path(pid,name)
                    manifest=item.get('artifact_manifest',{})
                    if file_digest(index)!=manifest.get(str(index.resolve())):continue
                    entries=[];files={}
                    for entry in json.loads(index.read_text()).get('entries',[]):
                        try:
                            if entry['key'] in seen:continue
                            valid={}
                            for relative,digest in entry['files'].items():
                                path=store.scoped_path(pid,index.parent/relative)
                                path.relative_to(index.parent.resolve())
                                if file_digest(path)!=digest or manifest.get(str(path.resolve()))!=digest:raise ValueError('Cache de dinâmica alterado.')
                                valid[relative]=path
                            if valid:entries.append(entry);files.update(valid)
                        except (OSError,ValueError,KeyError,TypeError):continue
                    if not entries:continue
                    target=destination/str(count)
                    for relative,path in files.items():
                        out=target/relative;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,out)
                    write_json(dict(version=1,entries=entries),target/NAME);seen.update(e['key'] for e in entries);count+=1
                except (OSError,ValueError,KeyError,TypeError):continue
    return count


class PairCache:
    def __init__(self,output,seeds=None):
        self.output=Path(output);self.entries=[];self.previous={}
        self.index=CalculationCache('molecular-dynamics',output)
        if seeds:
            for index in sorted(Path(seeds).glob('*/'+NAME)):
                try:
                    for entry in json.loads(index.read_text()).get('entries',[]):self.previous.setdefault(entry['key'],(index.parent,entry))
                except (OSError,ValueError,KeyError,TypeError):continue

    def restore(self,key,directory):
        if not self.index.reuse:return None
        metadata=self.index.metadata(key)
        if metadata:
            restored=self.index.restore(key,{name:Path(directory)/name for name in metadata['files']})
            if restored:return dict(restored['result'],reused=True)
        candidate=self.previous.get(key)
        if not candidate:return None
        root,entry=candidate;directory=Path(directory)
        try:
            for relative,digest in entry['files'].items():
                path=(root/relative).resolve();path.relative_to(root.resolve())
                if file_digest(path)!=digest:return None
            for relative in entry['files']:
                path=root/relative
                # Cache file names are rooted in the source pair directory.
                suffix=Path(relative).relative_to(entry['directory'])
                out=directory/suffix;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,out)
            return dict(entry['record'],reused=True)
        except (OSError,ValueError,KeyError,TypeError):return None

    def remember(self,key,directory,record):
        directory=Path(directory)
        manifest=artifact_manifest([p for p in directory.rglob('*') if p.is_file() and p.suffix in SUFFIXES])
        self.entries.append(dict(key=key,directory=directory.relative_to(self.output).as_posix(),record=record,
            files={Path(path).relative_to(self.output).as_posix():value for path,value in manifest.items()}))
        write_json(dict(version=1,entries=self.entries),self.output/NAME)
        files={p.relative_to(directory).as_posix():p for p in directory.rglob('*') if p.is_file() and p.suffix in SUFFIXES}
        self.index.remember(key,files,dict(files=list(files),result=record,
            receptor_identifier=record['receptor'],compound_identifier=record['ligand']))
