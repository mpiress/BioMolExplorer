# Auditoria das dependências remanescentes do Chimera — 23/10/2026

## Resultado

A chamada ao Chimera que apenas atribuía AM1-BCC na conversão legada de ligantes foi substituída pelo Antechamber direto. Nenhum protocolo de protonação, minimização ou carga de proteínas foi trocado. Não é seguro retirar o Chimera da instalação atual.

## Inventário de execução

| Caminho | Trabalho efetivo | Decisão |
| --- | --- | --- |
| `docking_data.prepare_ligands`, sem opções nativas, MOL2 + `am1` | Converter pose com Open Babel e atribuir AM1-BCC | **Substituído nesta revisão**: `charge_existing_mol2`, Antechamber direto, sem nova protonação/minimização |
| `docking_preparation.prepare_candidates` | Legado explícito; automático quando há minimização ou falta Meeko/Gemmi | Preservado. Casos compatíveis já seguem RDKit/Meeko ou Antechamber; não é uma chamada incondicional |
| `caad.docking.prepare_for_docking` / `prepare_complex.template` | Seleção de cadeia/cofatores, inversão/remoção, exportação PDB | Candidato à extração nativa, ainda não ativado: seletores personalizados, modelos, altlocs, conectividade inferida e conservação de resíduos precisam de validação própria |
| Mesmo método / `prepare_receptor.template` | Selecionar proteína/cofatores/águas, AddH, ff14SB nos resíduos padrão, método selecionado nos não padrão, minimização opcional | Preservado: não equivale a carregar todas as cargas com Gasteiger ou trocar ff14SB por ff94 do bundle DOCK6 |
| Mesmo método / `prepare_ligand.template` | Extrair resíduo cristalográfico, AddH, cargas, minimização opcional, PDB/MOL2 | Preservado: PDB isolado não garante ordens de ligação/cargas formais completas; H e geometria após minimização afetam centros e scores |
| `caad.docking.recover_better_conforms_of_vina` | Recuperação histórica de primeira pose, reatribuição de cargas/minimização | Sem chamador nos wrappers ativos pesquisados. Mantido para compatibilidade, não apresentado como caminho atual obrigatório Vina → DOCK6 |
| `caad.redocking.Docking` | Implementação histórica duplicada dos templates de preparação | Preservada; wrappers atuais importam `caad.docking`. Remover a API seria mudança de compatibilidade, não substituição científica |
| `prepare_md.template` | Recurso histórico editável com preparação PDB | Sem chamador científico identificado; **não** executa o bloco atual de dinâmica |
| `wrappers.docking.perform_consensus` | Mensagem para refinamento manual de loops quando falta PDB | Orientação, não subprocesso Chimera |
| `preparation_tools.amber_tool` | Localiza AmberTools no PATH, prefixo configurado ou bundle do Chimera | Não executa Chimera. Permite usar a mesma versão do Antechamber da referência |
| `ligand_cache`, `stage_cache`, validação e interface | Fingerprints de ferramentas, diagnóstico, seleção de templates e backend | Referências ao nome não significam execução. Pré-validação continua exigindo Chimera quando prepara receptores |

O bloco de dinâmica molecular usa `prepare_amber.pl`, Antechamber/LEaP, DOCK6 e, opcionalmente, sander/cpptraj. Não foi encontrada chamada ao Chimera em `md_engine.py` ou `wrappers/dynamics.py`. Superfícies usam a implementação SES nativa, não Chimera/DMS externo.

## Alteração implementada e integridade

`molecular_adapters.charge_existing_mol2` lê a química existente pelo núcleo compartilhado, exige uma molécula e usa os motores Antechamber existentes. Nomes duplicados importados recebem identificação inequívoca apenas no cálculo temporário, sem renomear os registros de saída. Os argumentos científicos continuam `-j 5 -dr n -at gaff`, com `qm_theory='AM1'` para AM1-BCC, como no AddCharge instalado. Não utiliza `Core.prepare`: isso impediria uma reprotonação acidental.

O adaptador valida quantidade de átomos, nomes, conectividade, coordenadas e soma/valores finitos das cargas antes de publicar. Substitui somente a coluna de cargas e o rótulo `USER_CHARGES`, preservando os outros registros MOL2. A publicação usa escrita atômica. Metais/química não suportada não recebem um modelo alternativo silencioso; falhas permanecem explícitas.

Reuso: SQLite/objetos SHA-256 já existentes, chave por conteúdo de entrada sem a carga antiga (saída a substituir), método, identidade dos executáveis/dados e implementação. Repetir a operação sobre o próprio resultado é idempotente. Mudar geometria ou método invalida somente o cálculo afetado. Falhas não publicam checkpoint; a retomada calcula o item pendente. Cópias com o mesmo conteúdo compartilham a preparação; não são criados scripts `.charges.com` nem temporários permanentes. O cache de conversão e o fingerprint do bloco também incluem o Antechamber efetivamente utilizado.

## Verificação científica real

Comando reproduzível, no ambiente científico do projeto:

```bash
PYTHONPATH=src python scripts/compare_antechamber_migration.py \
  --dock6 /caminho/para/dock6 --charge-only \
  --output /tmp/chimera-charge-only.json
```

Foram usados Chimera 1.17.3, seu Antechamber Amber18 e DOCK6 instalado. Oito moléculas cobrem neutro, cátion, ânion, heterociclo, aromático, halogênio, múltiplos heteroátomos e centro estereogênico; dois métodos, 16 comparações. Mesmos H explícitos, nomes, geometria, receptor e score contínuo de pose fixa. **16/16 passaram** nas tolerâncias predefinidas.

Máximos: diferença de carga atômica 0,0005; diferença de soma 0,001; RMSD de coordenadas 0; diferença de score total 0,000194 kcal/mol; eletrostática 0,000195; van der Waals 0. Tipos atômicos iguais nos 16 casos. Diferenças pequenas existem: não se alega igualdade numérica absoluta ou equivalência para todo espaço químico. Não é benchmark de qualidade de docking nem validação de minimização/receptores.

Dados completos: [JSON pareado](chimera_charge_only_equivalence_2026-10-23.json). Os testes automatizados cobrem preservação dos registros, cache, reuso desabilitado, mudança de geometria/método, rejeição de NaN, retomada após interrupção e integração da conversão legada. Nenhum ambiente principal foi modificado ou pacote novo instalado nesta revisão; usam-se as dependências AmberTools já declaradas.

Resultado de software: regressão de 74 testes aprovada (`test_chimera_charge_replacement`, `test_molecular_core`, `test_docking_interoperability`, `test_docking_configuration`, `test_docking_preparation`, `test_preparation_failures`, `test_stage_cache`, `test_redocking_incremental`), seguida de sete testes direcionados aprovados após o ajuste de nomes duplicados, com sobreposição entre execuções. Documentação regenerada; 34 páginas HTML com arquivos/âncoras/IDs válidos; compilação dos módulos alterados e `git diff --check` aprovados. Não é alegação de execução da suíte integral do repositório.

## Próximas substituições seguras em potencial

1. Extração/exportação de complexos: estabelecer fixtures pareadas para todos os seletores e casos estruturais, mantendo fallback para templates personalizados.
2. Recuperação de poses: transportar o grafo químico/mapeamento da molécula preparada até a pose; não reconstruir identidade apenas de PDBQT truncado.
3. Receptores: comparar estados de histidinas/terminais/pontes dissulfeto/cofatores, carga por átomo e parâmetros ff14SB em uma rota AmberTools/OpenMM explícita.
4. Minimização: definir protocolo e tolerâncias de energia/geometria/RMSD/score antes de habilitar outro minimizador. Desligá-la ou substituí-la por UFF/MMFF para eliminar dependência não seria migração sem perda.
