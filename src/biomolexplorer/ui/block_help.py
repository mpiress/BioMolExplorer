"""Brief, typed block guidance shared by the library and canvas."""
from biomolexplorer.catalog import TITLES,field_label
from biomolexplorer.flow import INPUTS,EXTERNAL_INPUTS,input_types,output_types

TYPE_DESCRIPTIONS={
    'compounds':'Compostos: CSV com identificador e SMILES; estruturas moleculares quando disponíveis.',
    'chembl':'Dados ChEMBL: compostos e bioatividades recuperados.',
    'structures':'Estruturas PDB brutas do receptor/complexo.',
    'prepared_structures':'Receptores preparados e metadados do sítio de docking.',
    'fingerprints':'Fingerprints: CSV com identificador e descritores moleculares.',
    'similarity':'Similaridades: CSV com pares de compostos e valores de similaridade.',
    'vina':'Resultados Vina: scores e poses PDBQT por composto e receptor.',
    'dock6':'Resultados DOCK6: scores e poses MOL2 por composto e receptor.',
    'scores':'Consenso: tabela de scores combinados por composto e receptor.',
    'redocking':'Validação do redocking: poses, scores e RMSD dos ligantes de referência.',
    'zinc_urls':'Lista de links ZINC: TXT/URI/URLS ou script exportado das tranches.',
    'other':'Arquivos locais compatíveis com o padrão indicado na configuração.',
    'visualization':'Visualizações: gráficos e relatórios HTML.',
}
NOTES={
    'similarity':'Fluxo recomendado: ChEMBL, PubChem ou ZINC → Gerar fingerprints → Calcular similaridade.',
    'prepare_structures':'Selecione Vina, DOCK6 ou ambos. Receptores já preparados no redocking são reutilizados; os compostos externos são preparados.',
    'docking_vina':'Use receptor e compostos preparados para Vina. Na entrada de compostos, selecione um composto específico ou mantenha todos.',
    'docking_dock6':'Selecione o receptor e, em Compostos selecionados, compostos ou resultados Vina. Receptores preparados são reutilizados. Escolha um composto específico ou mantenha todos.',
    'retrieve_pubchem':'Também aceita uma referência manual: SMILES, CID ou nome.',
    'graphs':'Conecte Calcular similaridade; compostos auxiliares podem ser enviados na configuração.',
}


def block_information(stage):
    operation=stage['operation']
    inputs=[]
    for field in INPUTS.get(operation,{}):
        if operation=='retrieve_pubchem' and stage['parameters'].get('reference_source')=='manual':continue
        inputs.append(field_label(operation,field))
        inputs.extend(TYPE_DESCRIPTIONS[k] for k in sorted(input_types(stage,field)) if k!='other')
    for field,kinds in EXTERNAL_INPUTS.get(operation,{}).items():
        inputs.append('Arquivos auxiliares enviados na configuração:')
        inputs.extend(TYPE_DESCRIPTIONS[k] for k in sorted(kinds))
    if not inputs:
        inputs=['Arquivos do disco selecionados e validados por tipo.' if operation=='import_results'
                else 'Consulta e filtros definidos na configuração do bloco.']
    outputs=[TYPE_DESCRIPTIONS[k] for k in sorted(output_types(stage))]
    if operation=='import_results':outputs=['Arquivos selecionados, publicados com os tipos informados na importação.']
    return TITLES[operation][1],inputs,outputs,NOTES.get(operation,'')
