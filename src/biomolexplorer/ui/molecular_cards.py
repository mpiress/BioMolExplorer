"""Responsive, read-only presentation of molecular identity and properties."""
import json
import flet as ft
from .localization import verbatim

INK='#0F172A'
MUTED='#64748B'
BLUE='#2552E8'

IDENTITY={'molecule_chembl_id':'ChEMBL ID','pref_name':'Nome do composto',
    'canonical_smiles':'SMILES','standard_inchi':'InChI','standard_inchi_key':'InChIKey',
    'inchi':'InChI','inchikey':'InChIKey','molecule_properties.full_molformula':'Fórmula molecular'}
DESCRIPTORS={
    'MW':('Massa molecular','g/mol'),'molecule_properties.full_mwt':('Massa molecular','g/mol'),
    'molecule_properties.mw_freebase':('Massa da base livre','g/mol'),
    'molecule_properties.alogp':('ALogP',''),'WLOGP':('WLOGP',''),'XLogP':('XLogP',''),
    'TPSA':('TPSA','Å²'),'molecule_properties.psa':('PSA','Å²'),
    'HBD':('Doadores de hidrogênio',''),'molecule_properties.hbd':('Doadores de hidrogênio',''),
    'HBA':('Aceptores de hidrogênio',''),'molecule_properties.hba':('Aceptores de hidrogênio',''),
    'RotatableBonds':('Ligações rotacionáveis',''),'molecule_properties.rtb':('Ligações rotacionáveis',''),
    'HeavyAtomCount':('Átomos pesados',''),'FormalCharge':('Carga formal',''),
}


def rendered(value):
    return json.dumps(value,ensure_ascii=False) if isinstance(value,(dict,list)) else str(value)


def badge(label,value,color=BLUE,system=False):
    return ft.Container(padding=ft.Padding.symmetric(horizontal=10,vertical=6),border_radius=20,bgcolor='#EFF6FF',
        content=ft.Row([ft.Text(label,size=11,color=MUTED),ft.Text(str(value),size=11,color=color,weight=ft.FontWeight.W_600) if system else verbatim(ft.Text(str(value),size=11,color=color,weight=ft.FontWeight.W_600))],spacing=5,tight=True))


def metric(label,value,unit='',col=6):
    return ft.Container(col=col,padding=12,bgcolor='#F8FAFC',border_radius=10,border=ft.Border.all(1,'#E2E8F0'),
        content=ft.Column([ft.Text(label,size=11,color=MUTED),
            verbatim(ft.Text(str(value),size=18,weight=ft.FontWeight.W_600,color=INK,selectable=True)),
            *([ft.Text(unit,size=10,color=MUTED)] if unit else [])],spacing=3))


def quick_metrics(properties,col=6):
    controls=[]
    for label,keys,unit in [('Massa molecular',('MW','molecule_properties.full_mwt'),'g/mol'),
            ('LogP',('WLOGP','molecule_properties.alogp','XLogP'),''),('TPSA',('TPSA','molecule_properties.psa'),'Å²')]:
        value=next((properties[k] for k in keys if properties.get(k) is not None and properties[k]!=''),None)
        if value is not None:
            try:value=f'{float(value):.6g}'
            except (TypeError,ValueError):pass
            controls.append(metric(label,value,unit,col=col))
    return ft.ResponsiveRow(controls,spacing=8,run_spacing=8,visible=bool(controls))


def image_panel(image,height=190,notice='Estrutura indisponível',framed=True):
    if not framed:
        return ft.Image(src=image,height=height,fit=ft.BoxFit.CONTAIN) if image else ft.Text(notice,size=12,color=MUTED)
    return ft.Container(padding=12,bgcolor='#F8FAFC',border=ft.Border.all(1,'#E2E8F0'),border_radius=12,
        alignment=ft.Alignment.CENTER,content=ft.Image(src=image,height=height,fit=ft.BoxFit.CONTAIN) if image else
            ft.Column([ft.Icon(ft.Icons.SCIENCE_OUTLINED,size=32,color=MUTED),ft.Text(notice,size=12,color=MUTED)],
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,spacing=10))


def property_sections(properties,original):
    """Keep every value and its origin, grouping long identifiers separately."""
    groups={name:[] for name in ('Identificadores químicos','Propriedades físico-químicas','Predições ADMET','Outros atributos')}
    origins=properties.get('property_sources',{})
    if not isinstance(origins,dict):origins={}
    for key,value in properties.items():
        if value is None or value=='':continue
        if key in ('degree','degree_density','component','property_sources'):continue
        if key in ('molecule_properties','properties') and any(k.startswith(key+'.') for k in properties):continue
        label,unit=DESCRIPTORS.get(key,(IDENTITY.get(key,key),''))
        origin=origins.get(key) or ('grafo' if original.get(key.split('.')[0]) is not None else 'fonte selecionada')
        if key in IDENTITY:
            field=ft.TextField(label=label,value=rendered(value),read_only=True,multiline=True,min_lines=1,max_lines=3,text_size=12)
            group='Identificadores químicos';col=12
        else:
            field=ft.Column([ft.Text(label,size=11,color=MUTED) if key in DESCRIPTORS else verbatim(ft.Text(label,size=11,color=MUTED)),
                verbatim(ft.Text(rendered(value)+(' '+unit if unit else ''),size=13,color=INK,selectable=True))],spacing=5)
            group=('Propriedades físico-químicas' if key in DESCRIPTORS or key.startswith('molecule_properties.') else
                'Predições ADMET' if key in ('BBB','HIA') or key.startswith(('ML_','ADMET','admet','BBB','HIA')) else 'Outros atributos')
            col={'xs':12,'sm':6}
        groups[group].append(ft.Container(col=col,padding=12,border_radius=10,bgcolor='#F8FAFC',
            content=ft.Column([field,ft.Text('Origem: '+str(origin),size=10,color=MUTED)],spacing=6,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)))
    return [ft.ExpansionTile(title=ft.Text(name,size=14,weight=ft.FontWeight.W_600),expanded=name!='Outros atributos',
        controls=[ft.ResponsiveRow(rows,spacing=8,run_spacing=8)]) for name,rows in groups.items() if rows]
