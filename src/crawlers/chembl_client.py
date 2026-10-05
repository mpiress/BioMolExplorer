"""Paginated ChEMBL client with conservative queries and local field selection.

An exact target ID identifies one record. Its organism/type can therefore be
checked locally without asking the API to perform extra joins. Field selection
is also local: a server-side ``only`` projection can fail on nested fields even
when the corresponding full record remains available.
"""
from urllib.parse import urlsplit
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL = 'https://www.ebi.ac.uk/chembl/api/data'
KEYS = {'target':'targets', 'activity':'activities', 'molecule':'molecules', 'similarity':'molecules'}


def _target_matches(row, filters):
    for key, expected in filters.items():
        field = key.removesuffix('__in')
        value = row.get(field)
        choices = expected if key.endswith('__in') else [expected]
        if isinstance(choices, str):
            choices = choices.split(',') if key.endswith('__in') else [choices]
        if value not in choices:
            return False
    return True

class ChEMBLQuery:
    def __init__(self,resource,parameters,columns=None,maximum=None):
        self.resource=resource
        self.parameters=dict(parameters)
        self.columns=columns
        self.maximum=maximum
        self._records=None
    def only(self,columns):
        return ChEMBLQuery(self.resource,self.parameters,columns,self.maximum)
    def take(self,maximum):
        if type(maximum) is not int or maximum<1: raise ValueError('O limite deve ser um inteiro positivo.')
        return ChEMBLQuery(self.resource,self.parameters,self.columns,maximum)
    def _load(self):
        if self._records is not None:
            return self._records
        query_parameters = dict(self.parameters)
        for key, value in tuple(query_parameters.items()):
            if key.endswith('__in') and isinstance(value, (list, tuple)):
                if not value:
                    self._records = []
                    return self._records
                if len(value) == 1 and key.removesuffix('__in') not in query_parameters:
                    query_parameters[key.removesuffix('__in')] = query_parameters.pop(key)[0]
        local_filters = {}
        if self.resource == 'target' and query_parameters.get('target_chembl_id'):
            for key in ('organism', 'organism__in', 'target_type', 'target_type__in'):
                if key in query_parameters:
                    local_filters[key] = query_parameters.pop(key)
        parameters={k:','.join(map(str,v)) if isinstance(v,(list,tuple)) else str(v).lower() if isinstance(v,bool) else v for k,v in query_parameters.items()}
        path=self.resource
        if path=='similarity':
            identifier=parameters.pop('chembl_id')
            threshold=parameters.pop('similarity',70)
            path=f'similarity/{identifier}/{threshold}'
        exact_identifier = (self.resource in ('target', 'molecule')
                            and set(parameters) == {f'{self.resource}_chembl_id'})
        if not exact_identifier:
            parameters.setdefault('limit',min(100,self.maximum) if self.maximum else 100)
            parameters.setdefault('offset',0)
        session=requests.Session()
        retry=Retry(total=2,backoff_factor=.5,status_forcelist=(429,500,502,503,504),allowed_methods=('GET',),respect_retry_after_header=True,raise_on_status=False)
        session.mount('https://',HTTPAdapter(max_retries=retry))
        records=[]
        seen=set()
        url=f'{BASE_URL}/{path}.json'
        try:
            while url:
                if url in seen:
                    raise RuntimeError('A paginação da ChEMBL repetiu uma página.')
                seen.add(url)
                response=session.get(url,params=parameters,timeout=(5,20))
                if not response.ok:
                    raise RuntimeError(f'ChEMBL indisponível (HTTP {response.status_code}) na consulta {self.resource}. Tente novamente em alguns minutos.')
                data=response.json()
                if not isinstance(data, dict):
                    raise RuntimeError('Resposta inválida da API ChEMBL: objeto de dados ausente.')
                batch=data.get(KEYS[self.resource])
                if not isinstance(batch,list):
                    raise RuntimeError('Resposta inválida da API ChEMBL: lista de registros ausente.')
                if not all(isinstance(row, dict) for row in batch):
                    raise RuntimeError('Resposta inválida da API ChEMBL: registro molecular ausente.')
                if local_filters:
                    batch = [row for row in batch if _target_matches(row, local_filters)]
                records.extend(batch)
                metadata = data.get('page_meta') or {}
                if not isinstance(metadata, dict):
                    raise RuntimeError('Resposta inválida da API ChEMBL: metadados de paginação ausentes.')
                next_url=metadata.get('next')
                if self.maximum is not None and len(records)>=self.maximum:
                    records=records[:self.maximum]; next_url=None
                if next_url:
                    if not isinstance(next_url, str):
                        raise RuntimeError('Endereço de paginação ChEMBL inválido.')
                    from urllib.parse import urljoin
                    url=urljoin(BASE_URL+'/',next_url)
                    parsed=urlsplit(url)
                    if parsed.scheme!='https' or parsed.netloc!='www.ebi.ac.uk' or not parsed.path.startswith('/chembl/api/data/'):
                        raise RuntimeError('Endereço de paginação ChEMBL inválido.')
                else: url=None
                parameters=None
        except requests.RequestException as exc:
            identifier = self.parameters.get('molecule_chembl_id') or self.parameters.get('target_chembl_id') or self.parameters.get('chembl_id')
            context = f' ({identifier})' if identifier else ''
            raise RuntimeError(f'Falha de conexão ou tempo limite na consulta ChEMBL {self.resource}{context}, após as tentativas automáticas. Tente novamente em alguns minutos; os detalhes estão no log da etapa.') from exc
        finally:
            session.close()
        if self.columns:
            columns = list(dict.fromkeys(column.split('__', 1)[0] for column in self.columns))
            records = [{column: row.get(column) for column in columns} for row in records]
        self._records=records
        return records
    def __getitem__(self,index): return self._load()[index]
    def __iter__(self): return iter(self._load())
    def __len__(self): return len(self._load())

class ChEMBLResource:
    def __init__(self,name): self.name=name
    def filter(self,**parameters): return ChEMBLQuery(self.name,parameters)

class ChEMBLClient:
    def __init__(self):
        for name in KEYS: setattr(self,name,ChEMBLResource(name))
