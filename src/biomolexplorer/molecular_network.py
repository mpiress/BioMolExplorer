"""Explicit, source-scoped molecular joins and exploratory numerical summaries."""
from biomolexplorer.molecular_parsing import mol_from_smiles
import csv
import ast
import io
import json
import math


def enrich(model,content,source):
    from .input_validation import ALIASES
    from rdkit import Chem
    nodes={n['id']:n for n in model['nodes']};groups={};ambiguous=set()
    reader=csv.DictReader(io.StringIO(content) if isinstance(content,str) else content)
    if not reader.fieldnames or not any(ALIASES.get(f,f)=='molecule_chembl_id' for f in reader.fieldnames):
        raise ValueError('A fonte precisa de identificadores molecule_chembl_id.')
    for row in reader:
        row={ALIASES.get(k,k):v for k,v in row.items() if k is not None and k!='fingerprint'}
        identifier=row.get('molecule_chembl_id')
        if identifier not in nodes or identifier in ambiguous:continue
        if identifier in groups and row!=groups[identifier]:
            ambiguous.add(identifier);groups.pop(identifier)
        else:groups[identifier]=row
    canonical={}
    def signature(value):
        if value not in canonical:
            molecule=mol_from_smiles(value) if value else None
            canonical[value]=Chem.MolToSmiles(molecule) if molecule else None
        return canonical[value]
    matches={};conflicts=sorted(ambiguous)
    for identifier,row in groups.items():
        original=nodes[identifier]['properties']
        left,right=original.get('canonical_smiles'),row.get('canonical_smiles')
        if left and right and (signature(left) is None or signature(left)!=signature(right)):
            conflicts.append(identifier);continue
        # Stored graph properties take precedence; conflicting values are reported.
        properties={k:v for k,v in row.items() if v not in (None,'')}
        collisions=[k for k,v in properties.items() if k in original and original[k] is not None and str(original[k])!=str(v)]
        for key in collisions:properties.pop(key)
        matches[identifier]={'properties':properties,'conflicting_fields':collisions}
    return dict(source=source,matches=matches,conflicts=conflicts,matched=len(matches),
        unmatched=len(nodes)-len(matches)-len(conflicts),rule='Exact compound ID; canonical SMILES must agree when available. Graph values take precedence. Ambiguous duplicate rows are excluded.')


def numeric_properties(nodes,enrichment=None,calculate=False):
    properties={n['id']:dict(n['properties']) for n in nodes}
    matches=(enrichment or {}).get('matches',{})
    for identifier in properties:
        row=matches.get(identifier,{'properties':{}})
        for key,value in row['properties'].items():
            if properties[identifier].get(key) is None:properties[identifier][key]=value
    for row in properties.values():
        for key in ('molecule_properties','properties'):
            nested=row.get(key)
            if isinstance(nested,str) and len(nested)<=32000:
                try:nested=json.loads(nested)
                except (ValueError,TypeError):
                    try:nested=ast.literal_eval(nested)
                    except (ValueError,SyntaxError,TypeError):nested=None
            if isinstance(nested,dict):
                for field,value in nested.items():
                    if isinstance(value,(str,int,float,bool)) or value is None:row.setdefault(key+'.'+str(field),value)
        from .compound_properties import compound_properties
        row.update(compound_properties(row,calculate=calculate))
    numeric={}
    for identifier,row in properties.items():
        for key,value in row.items():
            if key in ('molecule_chembl_id','canonical_smiles','degree','component','degree_density') or key.lower().endswith('_id') or key.lower() in ('pubchem_cid','cid') or isinstance(value,bool):continue
            try:number=float(value)
            except (TypeError,ValueError):continue
            if math.isfinite(number):numeric.setdefault(key,{})[identifier]=number
    return properties,numeric


def histogram(values,bins=12):
    values=list(values)
    if not values:return []
    low,high=min(values),max(values)
    if low==high:return [(low,high,len(values))]
    width=(high-low)/bins;counts=[0]*bins
    for value in values:counts[min(bins-1,int((value-low)/width))]+=1
    return [(low+i*width,low+(i+1)*width,count) for i,count in enumerate(counts)]


def association(x,y):
    common=sorted(set(x)&set(y));pairs=[(n,x[n],y[n]) for n in common]
    if len(pairs)<3:return pairs,None
    mx=sum(x[n] for n in common)/len(common);my=sum(y[n] for n in common)/len(common)
    covariance=sum((x[n]-mx)*(y[n]-my) for n in common)
    variance=math.sqrt(sum((x[n]-mx)**2 for n in common)*sum((y[n]-my)**2 for n in common))
    return pairs,covariance/variance if variance else None
