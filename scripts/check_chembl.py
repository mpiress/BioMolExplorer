"""Bounded live IC50 retrieval probe, with dataset export and central diagnostics.

Run with the installed package or PYTHONPATH=src. This does not query PubChem or
structural similarities: it verifies actual measured IC50 reference compounds.
"""
import argparse
import json
from datetime import datetime,timezone
from pathlib import Path
from biomolexplorer.diagnostics import configure_logging, log_directory
from biomolexplorer.storage import write_dataframe
from crawlers.chembl_client import ChEMBLClient


def main():
    parser=argparse.ArgumentParser(description='Testar recuperação ChEMBL por IC50')
    parser.add_argument('--target',default='CHEMBL220')
    parser.add_argument('--limit',type=int,default=5,help='Máximo de atividades e compostos na amostra')
    parser.add_argument('--output',type=Path,default=Path('datasets/chembl-check'))
    args=parser.parse_args()
    if args.limit<1: parser.error('--limit deve ser positivo')
    logger=configure_logging('chembl-check')
    report={'target':args.target.upper(),'standard_type':'IC50','requested_limit':args.limit,'timestamp':datetime.now(timezone.utc).isoformat()}
    try:
        import pandas as pd
        client=ChEMBLClient()
        targets=list(client.target.filter(target_chembl_id=report['target']).take(1))
        if not targets: raise ValueError('Alvo não encontrado.')
        activities=list(client.activity.filter(target_chembl_id=report['target'],standard_type='IC50',standard_units='nM',pchembl_value__isnull=False,standard_value__lte=5000).take(args.limit))
        if not activities: raise ValueError('Nenhuma atividade IC50 corresponde aos filtros (nM, até 5000, com pChEMBL).')
        ids=list(dict.fromkeys(a['molecule_chembl_id'] for a in activities))
        molecules=list(client.molecule.filter(molecule_chembl_id__in=ids))
        rows=[{'molecule_chembl_id':m['molecule_chembl_id'],'canonical_smiles':(m.get('molecule_structures') or {}).get('canonical_smiles'),'source':'ChEMBL'} for m in molecules]
        dataset=pd.DataFrame(rows).dropna(subset=['canonical_smiles'])
        destination=args.output/report['target'];destination.mkdir(parents=True,exist_ok=True)
        write_dataframe(dataset,destination/'compounds.csv')
        write_dataframe(pd.DataFrame(activities),destination/'activities_IC50.csv')
        report.update(status='succeeded',activities=len(activities),compounds=len(dataset),output=str(destination.resolve()))
        logger.info('Consulta concluída: %s',report)
        result=0
    except Exception as exc:
        logger.exception('Falha na consulta ao vivo: target=%s standard_type=IC50',report['target'])
        report.update(status='failed',error=f'{type(exc).__name__}: {exc}')
        result=1
    log_directory().mkdir(parents=True,exist_ok=True)
    (log_directory()/'chembl-check-report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    print(json.dumps(report,indent=2,ensure_ascii=False))
    return result

if __name__=='__main__':raise SystemExit(main())
