"""Independent PubChem retrieval with the application's compound CSV contract."""
from pathlib import Path


def reference_choices(paths):
    """Read normalized reference IDs from already authorized compound files."""
    import pandas as pd
    from .input_validation import columns
    from .molecule_quality import clean_dataframe
    result={}
    for path in paths:
        if Path(path).suffix.lower()!='.csv' or 'canonical_smiles' not in columns(path):continue
        frame=clean_dataframe(pd.read_csv(path,dtype=str).fillna(''),source=str(path))
        for row in frame.to_dict('records'):
            result.setdefault(row['molecule_chembl_id'],row['canonical_smiles'])
    return result


def retrieve_pubchem(base_output_path:str, target:str='PubChem', reference_source:str='compounds',
                     reference_type:str='smiles', reference:str='', selection_mode:str='all',
                     compound_id:str='', threshold:int=75, max_records:int=1000,
                     base_input_path:str=None, input_file:str=None):
    """Retrieve only new PubChem similars; never mix input compounds into output."""
    import pandas as pd
    from crawlers.pubchem import PubChemSimilarMols
    from .paths import resolve_path
    from .retrieval import collection_name
    from .storage import write_dataframe,write_json
    from .molecule_quality import clean_dataframe
    from .input_validation import columns,validate_file
    root=resolve_path(base_output_path);label=collection_name(target)
    crawler=PubChemSimilarMols(root/'PubChem/similars'/label,threshold,max_records)
    try:
        if reference_source=='manual':
            references=crawler.reference(reference_type,reference)
        elif reference_source=='compounds':
            if not base_input_path:raise ValueError('Conecte compostos ou envie um arquivo CSV para a busca PubChem.')
            folder=resolve_path(base_input_path)
            if input_file and Path(input_file).name!=input_file:raise ValueError('Selecione um arquivo da pasta de entrada.')
            paths=[folder/input_file] if input_file else [p for p in sorted(folder.glob('*.csv')) if 'canonical_smiles' in columns(p)]
            if not paths:raise ValueError('Conecte compostos ou envie um arquivo CSV para a busca PubChem.')
            frames=[]
            for path in paths:
                validate_file(path,'compounds',validate_rows=False)
                frames.append(pd.read_csv(path,dtype=str).fillna(''))
            references=clean_dataframe(pd.concat(frames,ignore_index=True),source='referências PubChem')
            if selection_mode=='single':
                if not compound_id:raise ValueError('Selecione um composto de referência para a busca PubChem.')
                references=references[references['molecule_chembl_id']==compound_id].copy()
                if references.empty:raise ValueError('O composto selecionado não está nos arquivos de entrada PubChem.')
            elif selection_mode!='all':raise ValueError('Selecione todas as referências ou um composto específico.')
        else:raise ValueError('Escolha uma origem válida para as referências PubChem.')
        compounds=crawler.search(references)
    finally:
        crawler.session.close()
    compounds=clean_dataframe(compounds,source='PubChem')
    output=root/'compounds'/label
    write_dataframe(compounds,output/'compounds.csv')
    write_json({'provider':'PubChem','reference_source':reference_source,'reference_type':reference_type if reference_source=='manual' else None,
                'selection_mode':selection_mode if reference_source=='compounds' else 'single',
                'references':references[['molecule_chembl_id','canonical_smiles']].to_dict('records'),
                'threshold':threshold,'max_records_per_reference':max_records,'selected':len(compounds)},root/'retrieval_report.json')
    return compounds
