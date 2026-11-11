"""Scale-free ECR relevance with conservative therapeutic-target aggregation.

Raw energies are never compared across engines or receptors. The multi-target
geometric aggregation is a breadth policy, not a calibrated binding probability.
"""
import json
import math
from collections import defaultdict

MODEL = 'target-balanced-ecr-v1'
TAU = 0.1


def consensus_rows(vina, dock6, weight=1.):
    """One row per shared compound; weight remains a legacy no-op argument.

    Require at least one paired receptor observation, retain partial coverage
    explicitly, and penalize absent observations rather than imputing energies.
    """
    from .molecular_parsing import mol_from_smiles
    from rdkit import Chem

    best = {'vina': {}, 'dock6': {}}
    structures = {}
    canonical = {}
    parsed = {}
    targets = set()
    for engine, records in (('vina', vina), ('dock6', dock6)):
        for original in records:
            row = dict(original)
            code = str(row.get('molecule_chembl_id') or '').strip()
            target = str(row.get('receptor_id') or '').strip()
            if not code:
                raise ValueError('Consenso: código molecular vazio.')
            value = float(row['score'])
            if not math.isfinite(value):
                raise ValueError('Consenso: score não finito para '+code)
            smiles = row.get('canonical_smiles') or ''
            if smiles:
                if smiles not in parsed:
                    parsed[smiles] = Chem.MolToSmiles(mol_from_smiles(smiles), isomericSmiles=True)
                identity = parsed[smiles]
                if code in structures and structures[code] != identity:
                    raise ValueError('O mesmo código molecular identifica estruturas diferentes: '+code)
                structures[code] = identity
                canonical[code] = identity
            targets.add(target)
            row.update(score=value, receptor_id=target, molecule_chembl_id=code)
            key = (target, code)
            current = best[engine].get(key)
            # Ties select the same structural reference independently of input order.
            if current is None or (value, str(row.get('conformer_file', ''))) < (current['score'], str(current.get('conformer_file', ''))):
                best[engine][key] = row

    paired = best['vina'].keys() & best['dock6'].keys()
    codes = sorted({code for target, code in paired})
    if not codes:
        return []
    targets = sorted(targets)
    ranks = {'vina': {}, 'dock6': {}}
    for engine in best:
        groups = defaultdict(list)
        for (target, code), row in best[engine].items():
            groups[target].append((row['score'], code))
        for target, values in groups.items():
            values.sort()
            size = len(values)
            start = 0
            while start < size:
                stop = start+1
                while stop < size and values[stop][0] == values[start][0]:
                    stop += 1
                # Zero-based mean tie rank divided by N-1: best=0, worst=1.
                quantile = ((start+stop-1)/2)/max(size-1, 1)
                for _, code in values[start:stop]:
                    ranks[engine][target, code] = quantile
                start = stop

    rows = []
    for code in codes:
        details = []
        factors = {}
        for engine in best:
            observed = [ranks[engine].get((target, code)) for target in targets]
            coverage = sum(q is not None for q in observed)/len(targets)
            # Geometric mean of exp(-q/tau), in log space; absent ranks=worst.
            factors[engine] = coverage*math.exp(-math.fsum(q if q is not None else 1. for q in observed)/(len(targets)*TAU))
        for target in targets:
            detail = {'receptor_id': target}
            for engine in best:
                source = best[engine].get((target, code))
                q = ranks[engine].get((target, code))
                detail[engine+'_raw'] = source['score'] if source else None
                detail[engine+'_rank_fraction'] = q
                detail[engine+'_relevance'] = math.exp(-q/TAU) if q is not None else None
            details.append(detail)
        candidates = [target for target in targets if (target, code) in paired]
        representative = min(candidates, key=lambda t: (
            -sum(math.exp(-ranks[e][t, code]/TAU) for e in best), t))
        a, b = (best[e][representative, code] for e in ('vina', 'dock6'))
        consensus = (factors['vina']+factors['dock6'])/2
        row = dict(molecule_chembl_id=code, molecule=code, canonical_smiles=canonical.get(code, ''),
                   receptor_id=representative, vina=factors['vina'], dock6=factors['dock6'],
                   normalized_score=consensus, consensus=consensus, consensus_model=MODEL,
                   target_count=len(targets), paired_target_count=len(candidates),
                   vina_target_coverage=sum((t, code) in best['vina'] for t in targets)/len(targets),
                   dock6_target_coverage=sum((t, code) in best['dock6'] for t in targets)/len(targets),
                   vina_raw=a['score'], dock6_raw=b['score'],
                   target_scores=json.dumps(details, ensure_ascii=False, separators=(',', ':'), allow_nan=False))
        for engine, source in (('vina', a), ('dock6', b)):
            if source.get('conformer_file'):
                row[engine+'_pose'] = source['conformer_file']
            for field in ('receptor_file', 'reference_file'):
                if source.get(field):
                    row[engine+'_'+field] = source[field]
        for field in ('footprint_file', 'footprint_origin'):
            if b.get(field):
                row[field] = b[field]
        rows.append(row)
    return sorted(rows, key=lambda row: (-row['consensus'], row['molecule_chembl_id']))
