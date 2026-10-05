"""PubChem FastSimilarity 2D expansion of downloaded ChEMBL molecules."""
import ast
import hashlib
import json
import logging
import os
import time
from biomolexplorer.storage import write_dataframe, write_json
from biomolexplorer.rate_limit import wait_for_slot
from biomolexplorer.progress import report_progress
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class PubChemSimilarMols:
    BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    COLUMNS = ['molecule_chembl_id', 'canonical_smiles', 'molecule_properties',
               'source', 'PubChem_CID', 'InChIKey', 'Molecular_Formula',
               'Molecular_Weight', 'PubChem_Threshold_Percent']

    def __init__(self, output_path, threshold=75, max_records=1000):
        if type(threshold) is not int or not 1 <= threshold <= 100:
            raise ValueError('PubChem threshold must be an integer between 1 and 100')
        if type(max_records) is not int or max_records <= 0:
            raise ValueError('PubChem max_records must be a positive integer')
        self.output_path = Path(output_path)
        shared_cache = os.environ.get('BIOMOL_CACHE_DIR')
        self.cache_path = Path(shared_cache) / 'pubchem' if shared_cache else self.output_path / 'cache'
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.cache_path.mkdir(parents=True, exist_ok=True)
        self.threshold, self.max_records = threshold, max_records
        self.logger = logging.getLogger(__name__)
        self.session = requests.Session()
        retry = Retry(total=2, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
                      allowed_methods=['GET', 'POST'], raise_on_status=False)
        self.session.mount('https://', HTTPAdapter(max_retries=retry))
        self.session.headers['User-Agent'] = 'BioMolExplorer PubChem retrieval'
        self.last_request = 0

    def _request(self, endpoint, data, params=None):
        key = hashlib.sha256(json.dumps([endpoint, data, params], sort_keys=True).encode()).hexdigest()
        path = self.cache_path / (key + '.json')
        if path.exists():
            return json.loads(path.read_text())
        wait_for_slot(self.cache_path / '.request-rate.lock')
        try:
            response = self.session.post(f'{self.BASE_URL}/compound/{endpoint}',
                                         data=data, params=params, timeout=(5, 20))
            if response.status_code == 404:
                result = {}
            else:
                response.raise_for_status()
                result = response.json()
            if endpoint.startswith('cid/property/'):
                expected = set(map(int, data['cid'].split(',')))
                returned = {p['CID'] for p in result.get('PropertyTable', {}).get('Properties', [])}
                if expected != returned:
                    raise RuntimeError('Incomplete PubChem property response; response was not cached')
            write_json(result, path)
            return result
        except requests.RequestException as exc:
            response = getattr(exc, 'response', None)
            status = response.status_code if response is not None else None
            message = f'HTTP {status}' if status is not None else 'falha de conexão ou tempo limite'
            raise RuntimeError(f'Não foi possível concluir a consulta PubChem ({message}), após as tentativas automáticas. '
                               'Tente novamente em alguns minutos; os detalhes estão no log da etapa.') from exc
        finally:
            self.last_request = time.monotonic()

    @staticmethod
    def _identity(smiles):
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
        if mol is None:
            return None, None
        return Chem.MolToSmiles(mol, isomericSmiles=True), Chem.MolToInchiKey(mol)

    @staticmethod
    def read_chembl(paths):
        rows = []
        for directory in paths:
            for path in sorted(Path(directory).glob('*.csv')):
                for row in pd.read_csv(path).to_dict('records'):
                    structures = row.get('molecule_structures', {})
                    if isinstance(structures, str):
                        try:
                            structures = ast.literal_eval(structures)
                        except (ValueError, SyntaxError):
                            structures = {}
                    structures = structures if isinstance(structures, dict) else {}
                    row['canonical_smiles'] = structures.get('canonical_smiles') or row.get('canonical_smiles')
                    row['InChIKey'] = structures.get('standard_inchi_key') or row.get('InChIKey')
                    rows.append(row)
        return pd.DataFrame(rows)

    def search(self, chembl):
        """Exclude ChEMBL CIDs before property retrieval; preserve query provenance.

        Threshold is a search cutoff, not an individual similarity score. Errors
        propagate so an incomplete retrieval is not silently reported as complete.
        Successful API responses are cached for resumable runs.
        """
        if chembl.empty:
            raise ValueError('No downloaded ChEMBL molecules available for PubChem search')
        references, known_smiles, known_keys = [], set(), set()
        for row in chembl.to_dict('records'):
            smiles, key = self._identity(row.get('canonical_smiles'))
            if not smiles:
                self.logger.warning('Skipping invalid ChEMBL SMILES: %s', row.get('molecule_chembl_id'))
                continue
            known_keys.add(key)
            supplied_key = row.get('InChIKey')
            if isinstance(supplied_key, str) and supplied_key:
                known_keys.add(supplied_key)
            if smiles not in known_smiles:
                references.append((row['molecule_chembl_id'], smiles))
            known_smiles.add(smiles)
        if not references:
            raise ValueError('No valid ChEMBL structures available for PubChem search')

        # Resolve every ChEMBL structure first, including structures with no hits
        # in their own similarity search, so exclusions apply across all queries.
        excluded = set()
        report_progress('Identificando no PubChem os compostos já obtidos no ChEMBL…', 0, len(references))
        for index, (_, smiles) in enumerate(references, 1):
            result = self._request('smiles/cids/JSON', {'smiles': smiles})
            excluded.update(result.get('IdentifierList', {}).get('CID', []))
            report_progress('Identificando no PubChem os compostos já obtidos no ChEMBL…', index, len(references))
        matches = []
        report_progress('Consultando compostos similares no PubChem…', 0, len(references))
        for index, (name, smiles) in enumerate(references, 1):
            result = self._request('fastsimilarity_2d/smiles/cids/JSON', {'smiles': smiles},
                                   {'Threshold': self.threshold, 'MaxRecords': self.max_records})
            cids = result.get('IdentifierList', {}).get('CID', [])
            if len(cids) >= self.max_records:
                self.logger.warning('PubChem result limit reached for %s; results may be truncated', name)
            matches.extend((name, smiles, cid) for cid in dict.fromkeys(cids) if cid not in excluded)
            report_progress('Consultando compostos similares no PubChem…', index, len(references))
        candidates = sorted({cid for _, _, cid in matches})
        records, accepted, new_structures = [], {}, {}
        properties = 'SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight'
        report_progress('Baixando as propriedades dos novos compostos PubChem…', 0, len(candidates))
        for offset in range(0, len(candidates), 100):
            batch = candidates[offset:offset + 100]
            result = self._request(f'cid/property/{properties}/JSON', {'cid': ','.join(map(str, batch))})
            props = result.get('PropertyTable', {}).get('Properties', [])
            if {p['CID'] for p in props} != set(batch):
                raise RuntimeError('Incomplete PubChem property response; rerun to resume from cache')
            for prop in props:
                smiles, key = self._identity(prop.get('SMILES') or prop.get('IsomericSMILES')
                                             or prop.get('ConnectivitySMILES') or prop.get('CanonicalSMILES'))
                cid = prop['CID']
                existing = new_structures.get(smiles) or new_structures.get(key) or new_structures.get(prop.get('InChIKey'))
                if existing:
                    accepted[cid] = existing
                    continue
                if not smiles or smiles in known_smiles or key in known_keys or prop.get('InChIKey') in known_keys:
                    continue
                known_smiles.add(smiles)
                known_keys.update([key, prop.get('InChIKey')])
                accepted[cid] = cid
                for identity in (smiles, key, prop.get('InChIKey')):
                    if identity:
                        new_structures[identity] = cid
                records.append({'molecule_chembl_id': f'PUBCHEM{cid}', 'canonical_smiles': smiles,
                                'molecule_properties': {'full_mwt': float(prop['MolecularWeight'])},
                                'source': 'PubChem', 'PubChem_CID': cid, 'InChIKey': prop.get('InChIKey'),
                                'Molecular_Formula': prop.get('MolecularFormula'),
                                'Molecular_Weight': prop.get('MolecularWeight'),
                                'PubChem_Threshold_Percent': self.threshold})
            report_progress('Baixando as propriedades dos novos compostos PubChem…',
                            min(offset + 100, len(candidates)), len(candidates))
        compounds = pd.DataFrame(records, columns=self.COLUMNS)
        provenance = pd.DataFrame([(name, smiles, accepted[cid], self.threshold) for name, smiles, cid in matches
                                   if cid in accepted], columns=['Query', 'Query_SMILES', 'PubChem_CID',
                                                               'PubChem_Threshold_Percent']).drop_duplicates()
        write_dataframe(compounds, self.output_path / 'compounds.csv')
        write_dataframe(provenance, self.output_path / 'matches.csv')
        return compounds
