"""Shared molecular property names with explicit provider/method provenance."""
import ast
import json
import math

COMMON_COLUMNS=['MW','TPSA','WLOGP','HBD','HBA','RotatableBonds','HeavyAtomCount','Molecular_Formula',
                'FormalCharge','standard_inchi','standard_inchi_key','pref_name','property_sources']


def mapping(value):
    if isinstance(value,dict):return dict(value)
    if isinstance(value,str):
        try:value=json.loads(value)
        except (ValueError,TypeError):
            try:value=ast.literal_eval(value)
            except (ValueError,SyntaxError,TypeError):return {}
        if isinstance(value,dict):return value
    return {}


def present(value):
    if value is None or value=='':return False
    if isinstance(value,float) and not math.isfinite(value):return False
    return True


def compound_properties(record,calculate=False,source=None):
    """Fill missing common fields; retain provider values and distinguish logP methods."""
    row=dict(record)
    structures=mapping(row.get('molecule_structures'))
    nested=mapping(row.get('molecule_properties'))
    origins=mapping(row.get('property_sources'))
    provider=source or row.get('source') or ('ChEMBL' if str(row.get('molecule_chembl_id','')).startswith('CHEMBL') else 'arquivo')
    if not present(row.get('source')):row['source']=provider
    aliases={
        'canonical_smiles':[(row,'SMILES'),(structures,'canonical_smiles')],
        'MW':[(row,'Molecular_Weight'),(row,'MolecularWeight'),(nested,'full_mwt')],
        'Molecular_Formula':[(row,'MolecularFormula'),(nested,'full_molformula')],
        'TPSA':[(row,'PubChem_TPSA')],
        'HBD':[(row,'HBondDonorCount'),(nested,'hbd')],
        'HBA':[(row,'HBondAcceptorCount'),(nested,'hba')],
        'RotatableBonds':[(row,'RotatableBondCount'),(nested,'rtb')],
        'HeavyAtomCount':[(nested,'heavy_atoms')],
        'FormalCharge':[(row,'Charge')],
        'standard_inchi':[(row,'InChI'),(row,'inchi'),(structures,'standard_inchi')],
        'standard_inchi_key':[(row,'InChIKey'),(row,'inchikey'),(structures,'standard_inchi_key')],
        'pref_name':[(row,'Title'),(row,'IUPACName')],
    }
    for field,candidates in aliases.items():
        if present(row.get(field)):
            origins.setdefault(field,str(provider)+': '+field)
            continue
        for values,key in candidates:
            if present(values.get(key)):
                row[field]=values[key]
                origins[field]=str(provider)+': '+('molecule_properties.' if values is nested else 'molecule_structures.' if values is structures else '')+key
                break
    # PSA and ALogP from ChEMBL and XLogP from PubChem remain available under
    # their actual method names. Never relabel ALogP or XLogP as RDKit WLOGP.
    if calculate and any(not present(row.get(field)) for field in
                         ('MW','TPSA','WLOGP','HBD','HBA','RotatableBonds','HeavyAtomCount','FormalCharge','Molecular_Formula')):
        from .molecular_parsing import mol_from_smiles
        from rdkit import Chem
        from rdkit.Chem import Descriptors,rdMolDescriptors,Crippen,Lipinski
        smiles=row.get('canonical_smiles')
        molecule=mol_from_smiles(smiles) if isinstance(smiles,str) and smiles else None
        if molecule is not None:
            functions={'MW':Descriptors.MolWt,'TPSA':rdMolDescriptors.CalcTPSA,'WLOGP':Crippen.MolLogP,
                       'HBD':Lipinski.NumHDonors,'HBA':Lipinski.NumHAcceptors,
                       'RotatableBonds':Lipinski.NumRotatableBonds,'HeavyAtomCount':Lipinski.HeavyAtomCount,
                       'FormalCharge':Chem.GetFormalCharge,'Molecular_Formula':rdMolDescriptors.CalcMolFormula}
            methods={'MW':'Descriptors.MolWt','TPSA':'rdMolDescriptors.CalcTPSA','WLOGP':'Crippen.MolLogP',
                     'HBD':'Lipinski.NumHDonors','HBA':'Lipinski.NumHAcceptors',
                     'RotatableBonds':'Lipinski.NumRotatableBonds','HeavyAtomCount':'Lipinski.HeavyAtomCount',
                     'FormalCharge':'Chem.GetFormalCharge','Molecular_Formula':'rdMolDescriptors.CalcMolFormula'}
            for field,function in functions.items():
                if present(row.get(field)):continue
                try:row[field]=function(molecule)
                except (ValueError,RuntimeError):continue
                origins[field]='RDKit: '+methods[field]
    if nested:row['molecule_properties']=nested
    row['property_sources']=origins
    return row


def compound_frame(frame,calculate=True,source=None):
    import pandas as pd
    if frame.empty:return frame
    records=frame.astype(object).where(pd.notna(frame),None).to_dict('records')
    return pd.DataFrame([compound_properties(row,calculate,source) for row in records])
