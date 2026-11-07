"""Audit presentation literals, dynamic templates and dictionary placeholders.

Run with PYTHONPATH=src python scripts/audit_localization.py [--output report.json].
Native tool output, editable values and user/scientific data are not translated.
"""
import argparse
import ast
import json
import re
from pathlib import Path
from string import Formatter

from biomolexplorer.i18n import ROOT,dictionary,locales
from biomolexplorer.ui.localization import Translator

PORTUGUESE=re.compile(r'[ãõçáéíóúêô]|\b(?:Analisando|Organizando|Criando|Restaurando|Excluindo|Entrando|Saindo|Importando|Exportando|Removendo|Concluindo|Identificando|Selecione|Aguarde|Carregando|Preparando|Salvando|Escolha|Nenhum|Nenhuma|Informe|Excluir|Gerando|Verificando|Processando|Calculando|Lendo|Não|Continuar|Fechar|Baixar|Cancelar)\b',re.I)


def fields(value):
    return {field for _,field,_,_ in Formatter().parse(value) if field is not None}


def template(node,slots):
    if isinstance(node,ast.Constant) and isinstance(node.value,str):return node.value
    if isinstance(node,ast.JoinedStr):return ''.join(template(child.value if isinstance(child,ast.FormattedValue) else child,slots) for child in node.values)
    if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Add):return template(node.left,slots)+template(node.right,slots)
    index=len(slots);slots.append(index);return '{'+str(index)+'}'


def source_messages(root):
    for path in sorted((root/'src/biomolexplorer').rglob('*.py')):
        tree=ast.parse(path.read_text(encoding='utf-8'))
        parents={id(child):node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}
        for node in ast.walk(tree):
            parent=parents.get(id(node))
            if isinstance(parent,(ast.JoinedStr,ast.BinOp,ast.FormattedValue)):continue
            if isinstance(node,ast.Constant) and isinstance(node.value,str):value=node.value
            elif isinstance(node,ast.JoinedStr) or isinstance(node,ast.BinOp) and isinstance(node.op,ast.Add):value=template(node,[])
            else:continue
            if not PORTUGUESE.search(value) or value.startswith(('SELECT ','UPDATE ','INSERT ','DELETE ')):continue
            if isinstance(parent,ast.Expr):continue  # Docstrings.
            owner=parent
            while owner is not None and not isinstance(owner,(ast.Call,ast.stmt)):owner=parents.get(id(owner))
            if isinstance(owner,ast.Call) and isinstance(owner.func,ast.Attribute) and owner.func.attr in ('exception','debug','warning','info','error'):
                continue  # Raw diagnostics remain unchanged; stage_log presents them.
            yield str(path.relative_to(root)),node.lineno,value


def audit(root):
    errors=[];counts={};translate=Translator('en');messages=list(source_messages(root))
    for path,line,source in messages:
        probe=re.sub(r'\{(\d+)\}',lambda m:'SLOT'+m[1]+'VALUE',source)
        rendered=translate(probe)
        if PORTUGUESE.search(rendered):errors.append({'path':path,'line':line,'source':source,'rendered':rendered})
    for code in locales():
        counts[code]={}
        for kind in ('phrases','patterns','runtime','viewer','scene','methods'):
            reference=dictionary('en',kind);target=dictionary(code,kind);counts[code][kind]=len(target)
            for key,value in reference.items():
                # Portuguese phrase/pattern dictionaries support legacy English
                # aliases; new dictionaries use the canonical source keys.
                translated=target.get(key,target.get(value) if kind in ('phrases','patterns') else None)
                if translated is None:
                    errors.append({'language':code,'kind':kind,'missing':key});continue
                try:
                    if fields(value)!=fields(translated):errors.append({'language':code,'kind':kind,'placeholders':key})
                except ValueError as exc:errors.append({'language':code,'kind':kind,'invalid_format':key,'error':str(exc)})
    return {'presentation_expressions':len(messages),'catalogs':counts,'errors':errors}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=audit(Path(__file__).resolve().parents[1]);data=json.dumps(result,ensure_ascii=False,indent=2)
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(data+'\n',encoding='utf-8')
    print(data);raise SystemExit(bool(result['errors']))


if __name__=='__main__':main()
