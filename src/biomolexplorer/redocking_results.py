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
        metadata_seen = set()
        for item in self.results.files(token, pid, rid, sid):
            metadata = Path(item['path'])
            if metadata.name != 'pdb_codes.csv' or metadata in metadata_seen:
                continue
            metadata_seen.add(metadata)
            content = self.store.read_file(token, pid, str(metadata)).decode('utf-8-sig')
            seen = set()
            for record in csv.DictReader(io.StringIO(content)):
                try:
                    code, ligand, chain = (record[k].strip() for k in ('PDB_CODE', 'LIGAND', 'CHAIN'))
                    number = record['RESNUM'].strip()
                    int(number)
                    rmsd = float(record['RMSD'])
                    if not all(re.fullmatch(r'[A-Za-z0-9_+-]{1,40}', value) for value in (code, ligand, chain)) or not math.isfinite(rmsd) or rmsd < 0:
                        continue
                except (KeyError, ValueError, TypeError):
                    continue
                identity = (code, ligand, number, chain)
                if identity in seen:
                    continue
                seen.add(identity)
                key = hashlib.sha256((str(metadata)+'\0'+'|'.join(identity)).encode()).hexdigest()
                simulations.append(dict(id=key, pdb=code, ligand=ligand, residue=number,
                    chain=chain, rmsd=rmsd, metadata=str(metadata)))
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

    def archive(self, token, pid, rid, sid, key):
        simulation = self.simulation(token, pid, rid, sid, key)
        output = io.BytesIO()
        with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as archive:
            for item in simulation['files']:
                archive.writestr(item['archive_name'], self.store.read_file(token, pid, item['path']))
        return output.getvalue()

    def scene(self,token,pid,rid,sid,key):
        simulation=self.simulation(token,pid,rid,sid,key)
        files=[Path(item['path']) for item in simulation['files'] if item['molecular']]
        metadata=Path(simulation['metadata']);prepared=metadata.parent/'Prepared'
        identity=f"{simulation['pdb']}_{simulation['ligand']}_{simulation['residue']}{simulation['chain']}"
        receptor_id=f"{simulation['pdb']}_{simulation['chain']}"
        poses=[p for p in files if p.name==identity+'.lig.pdbqt' and p.parent!=prepared]
        if len(poses)!=1:raise ValueError('A pose do redocking está ausente ou é ambígua.')
        complex_file=next((p for p in files if p==metadata.parent/(simulation['pdb']+'.pdb')),None)
        if complex_file is None:complex_file=next((p for p in files if p.name==receptor_id+'.complex.pdb'),None)
        prepared_receptor=next((p for name in (receptor_id+'.noH.pdb',receptor_id+'.dockprep.pdbqt') for p in files if p.name==name),None)
        receptor=prepared_receptor or complex_file
        if receptor is None:raise ValueError('O receptor associado a esta pose não está disponível.')
        residue=[simulation['ligand'],simulation['residue'],simulation['chain']]
        layers=[dict(role='receptor',path=str(receptor),format=receptor.suffix[1:]),
                dict(role='pose',path=str(poses[0]),format='pdbqt')]
        if complex_file:
            if receptor==complex_file:layers[0]['exclude_residue']=residue
            else:layers.append(dict(role='complex',path=str(complex_file),format='pdb',exclude_residue=residue))
            layers.append(dict(role='reference',path=str(complex_file),format='pdb',residue=residue))
        else:
            reference=next((p for p in files if p.parent==prepared and p.name==identity+'.lig.pdb'),None)
            if reference:layers.append(dict(role='reference',path=str(reference),format='pdb'))
        return dict(name=identity+' · Redocking',layers=layers,
            reference_kind='crystal' if complex_file else 'prepared' if len(layers)>2 else None)
