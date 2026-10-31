"""Authorized RMSD summaries and simulation-specific redocking downloads."""
import csv
import hashlib
import io
import math
import re
import zipfile
from pathlib import Path

from .result_files import ResultFiles
from .workspace import AccessDenied


def best_rmsd_by_target(rows):
    best={}
    for row in rows:
        value=row.get('rmsd')
        if row.get('status','succeeded')!='succeeded' or not isinstance(value,(int,float)) or not math.isfinite(value):continue
        for target in row.get('targets',[]):best[target]=min(best.get(target,math.inf),value)
    return best


class RedockingResults:
    def __init__(self, store):
        self.store = store
        self.results = ResultFiles(store)

    @staticmethod
    def root(metadata):
        # Native operations place structures/<target> and Vina <target> side by side.
        # Imported datasets keep their own metadata folder as the simulation boundary.
        return metadata.parent.parent.parent if metadata.parent.parent.name == 'structures' else metadata.parent

    def simulations(self, token, pid, rid, sid):
        stage = self.results.stage(token, pid, rid, sid)
        if stage['operation'] != 'redocking':
            raise AccessDenied('Esta etapa não contém resultados de redocking.')
        if stage['status'] != 'succeeded':
            return []
        simulations = []
        origins = self.results.origins(token,pid,rid)
        metadata_seen = set()
        for item in self.results.files(token, pid, rid, sid):
            metadata = Path(item['path'])
            if metadata.name != 'pdb_codes.csv' or metadata in metadata_seen:
                continue
            metadata_seen.add(metadata)
            content = self.store.read_file(token, pid, str(metadata)).decode('utf-8-sig')
            routes = origins.describe(str(metadata),sid,rid)['source_routes']
            seen = set()
            for record in csv.DictReader(io.StringIO(content)):
                try:
                    code, ligand, chain = (record[k].strip() for k in ('PDB_CODE', 'LIGAND', 'CHAIN'))
                    number = record['RESNUM'].strip()
                    int(number)
                    status=record.get('STATUS') or 'succeeded';failed=status!='succeeded'
                    rmsd = None if failed else float(record['RMSD'])
                    if not all(re.fullmatch(r'[A-Za-z0-9_+-]{1,40}', value) for value in (code, ligand, chain)) or (not failed and (not math.isfinite(rmsd) or rmsd < 0)):
                        continue
                except (KeyError, ValueError, TypeError):
                    continue
                identity = (code, ligand, number, chain)
                if identity in seen:
                    continue
                seen.add(identity)
                key = hashlib.sha256((str(metadata)+'\0'+'|'.join(identity)).encode()).hexdigest()
                from .dataset_identity import redocking_targets
                targets=redocking_targets(routes,code)
                if not targets:
                    explicit=(record.get('Target') or record.get('TARGET') or '').strip()
                    targets=[explicit] if explicit else []
                simulations.append(dict(id=key, pdb=code, ligand=ligand, residue=number,
                    target=' + '.join(targets),targets=targets,engine=(record.get('ENGINE') or stage.get('configuration',{}).get('parameters',{}).get('redocking_engine','vina')).lower(),
                    rmsd_method=record.get('RMSD_METHOD') or 'legacy_prepared',rmsd_reference=record.get('RMSD_REFERENCE') or 'prepared',
                    chain=chain, rmsd=rmsd, metadata=str(metadata),
                    status=status if failed else 'succeeded',error=record.get('ERROR',''),
                    failed_step=record.get('FAILED_STEP','')))
        best=best_rmsd_by_target(simulations)
        for row in simulations:
            row['best_rmsd_for']=[t for t in row['targets'] if row['rmsd'] is not None and
                math.isclose(row['rmsd'],best.get(t,math.inf),rel_tol=0,abs_tol=1e-9)]
        return simulations

    def simulation(self, token, pid, rid, sid, key):
        row = next((r for r in self.simulations(token, pid, rid, sid) if r['id'] == key), None)
        if row is None:
            raise AccessDenied('Simulação de redocking indisponível nesta etapa.')
        metadata = Path(row['metadata'])
        root = self.root(metadata)
        reference = f"{row['pdb']}_{row['ligand']}_{row['residue']}{row['chain']}"
        receptor = f"{row['pdb']}_{row['chain']}"
        files = []
        seen = set()
        for item in self.results.files(token, pid, rid, sid):
            path = Path(item['path'])
            if not path.is_relative_to(root) or path in seen:
                continue
            if path == metadata:
                description = 'Metadados e RMSD'
            elif path == metadata.parent/'Prepared'/'centers.csv':
                description = 'Centros dos ligantes'
            elif path.name == row['pdb']+'.pdb':
                description = 'Estrutura de entrada'
            elif path.name.startswith(receptor+'.'):
                description = 'Complexo preparado' if '.complex.' in path.name else 'Receptor preparado'
            elif path.name.startswith(reference+'.'):
                description = 'Ligante de referência' if path.parent == metadata.parent/'Prepared' else 'Poses do redocking'
            else:
                continue
            seen.add(path)
            files.append(dict(item, description=description, archive_name=path.relative_to(root).as_posix(),
                molecular=path.suffix.lower() in ('.pdb', '.pdbqt', '.mol2')))
        return dict(row, files=files)

    def rmsd_csv(self,token,pid,rid,sid,keys=None):
        rows=self.simulations(token,pid,rid,sid)
        if keys is not None:
            selected=set(keys)
            if selected-{row['id'] for row in rows}:raise AccessDenied('Simulação de redocking indisponível nesta etapa.')
            rows=[r for r in rows if r['id'] in selected]
        output=io.StringIO(newline='')
        writer=csv.DictWriter(output,fieldnames=['Target','PDB','Ligand','Residue','Chain','Engine','RMSD','RMSD_method','RMSD_reference','Status','Best_for_target'])
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(Target=row['target'],PDB=row['pdb'],Ligand=row['ligand'],Residue=row['residue'],
                Chain=row['chain'],Engine=row['engine'],RMSD=row['rmsd'],RMSD_method=row['rmsd_method'],RMSD_reference=row['rmsd_reference'],Status=row['status'],Best_for_target=' + '.join(row['best_rmsd_for'])))
        return output.getvalue().encode('utf-8')

    def archive(self, token, pid, rid, sid, key):
        simulation = self.simulation(token, pid, rid, sid, key)
        output = io.BytesIO()
        with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as archive:
            for item in simulation['files']:
                archive.writestr(item['archive_name'], self.store.read_file(token, pid, item['path']))
        return output.getvalue()

    def remove_failure(self,token,pid,rid,sid,key):
        from .redocking_curation import remove_failure
        return remove_failure(self,token,pid,rid,sid,key)

    def scene(self,token,pid,rid,sid,key):
        simulation=self.simulation(token,pid,rid,sid,key)
        files=[Path(item['path']) for item in simulation['files'] if item['molecular']]
        metadata=Path(simulation['metadata']);prepared=metadata.parent/'Prepared'
        identity=f"{simulation['pdb']}_{simulation['ligand']}_{simulation['residue']}{simulation['chain']}"
        receptor_id=f"{simulation['pdb']}_{simulation['chain']}"
        extension='mol2' if simulation['engine']=='dock6' else 'pdbqt'
        poses=[p for p in files if p.name==identity+'.lig.'+extension and p.parent!=prepared]
        if len(poses)!=1:raise ValueError('A pose do redocking está ausente ou é ambígua.')
        complex_file=next((p for p in files if p==metadata.parent/(simulation['pdb']+'.pdb')),None)
        if complex_file is None:complex_file=next((p for p in files if p.name==receptor_id+'.complex.pdb'),None)
        prepared_receptor=next((p for name in (receptor_id+'.noH.pdb',receptor_id+'.dockprep.pdbqt') for p in files if p.name==name),None)
        receptor=prepared_receptor or complex_file
        if receptor is None:raise ValueError('O receptor associado a esta pose não está disponível.')
        residue=[simulation['ligand'],simulation['residue'],simulation['chain']]
        layers=[dict(role='receptor',path=str(receptor),format=receptor.suffix[1:]),
                dict(role='pose',path=str(poses[0]),format=extension)]
        if complex_file:
            if receptor==complex_file:layers[0]['exclude_residue']=residue
            else:layers.append(dict(role='complex',path=str(complex_file),format='pdb',exclude_residue=residue))
            layers.append(dict(role='reference',path=str(complex_file),format='pdb',residue=residue))
        else:
            reference=next((p for p in files if p.name==identity+'.reference.mol2'),None)
            if reference is None:reference=next((p for p in files if p.parent==prepared and p.name==identity+'.lig.pdb'),None)
            if reference:layers.append(dict(role='reference',path=str(reference),format=reference.suffix[1:]))
        return dict(name=identity+' · Redocking · '+simulation['engine'].upper(),layers=layers,
            reference_kind='crystal' if complex_file else 'prepared' if len(layers)>2 else None)
