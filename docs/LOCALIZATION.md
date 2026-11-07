# Dicionários e idiomas

A seleção pertence à sessão. Alterar o idioma não modifica parâmetros, nomes
escritos pelo usuário, identificadores, arquivos, SMILES ou resultados científicos.
Os logs brutos e as saídas de ferramentas externas permanecem disponíveis sem
reescrita; sua apresentação traduz as mensagens reconhecidas do BioMolExplorer.

Os arquivos ficam em `src/biomolexplorer/resources/i18n/`:

| Arquivo | Conteúdo |
| --- | --- |
| `locales.json` | Idiomas disponíveis, nome nativo, locale e bandeira opcional |
| `<idioma>.json` | Frases, rótulos, ajuda, acessibilidade e mensagens da interface |
| `patterns_<idioma>.json` | Mensagens com parâmetros numéricos `{0}`, `{1}`, etc. |
| `runtime_<idioma>.json` | Mensagens com identificadores estáveis usadas por workers e telas de execução |
| `viewer_<idioma>.json` | Visualizador molecular 3D |
| `scene_<idioma>.json` | Contatos e interações das cenas de docking |
| `methods_<idioma>.json` | Métodos descritos nos relatórios científicos |

Por compatibilidade, os padrões ingleses continuam em `patterns.json`, e os
arquivos portugueses também reconhecem aliases ingleses. Novos dicionários devem
usar as mesmas chaves dos arquivos ingleses. As chaves de frases e padrões podem
ser frases portuguesas legadas; o valor é o texto traduzido. Não renomeie as chaves
para traduzir: elas identificam mensagens usadas pelo código.

## Acrescentar um idioma

1. Acrescente sua entrada em `locales.json`, por exemplo:

   ```json
   "es": {"name": "Español", "language_code": "es", "country_code": "ES"}
   ```

2. Copie os seis dicionários ingleses para os nomes correspondentes ao novo
   idioma e traduza seus valores. Use `patterns.json` como origem para
   `patterns_es.json`. Preserve as chaves e os parâmetros. A ordem dos parâmetros
   pode mudar: `{1}` pode aparecer antes de `{0}`. Os nomes de parâmetros de
   mensagens runtime, como `{stage}`, também devem ser preservados.
3. Opcionalmente, adicione uma bandeira em `resources/branding/` e informe seu
   nome no campo `flag`. Sem bandeira, o seletor mostra o código do idioma.
4. Execute a auditoria abaixo e os testes de localização. Reinicie a aplicação.
   O seletor e `--language` descobrem a entrada pelo manifesto, sem editar Python.

Idiomas incompletos usam inglês para chaves conhecidas ainda não traduzidas. A
auditoria exige cobertura completa antes da distribuição. Mensagens desconhecidas
não são traduzidas por heurísticas que poderiam corromper nomes ou dados; devem
ser acrescentadas ao dicionário.

## Novos textos e verificação

Frases estáticas nos controles são traduzidas por `LocalizedPage`. Textos gerados
fora da árvore Flet, como documentos e imagens, usam `Translator(language)` ou
`localized_dictionary(language, namespace)`. Para mensagens runtime, use
`RuntimeMessages(language)(key, **arguments)`.

Prefira traduzir a frase antes de interpolar valores do usuário. Para compatibilidade
com mensagens já construídas, os padrões reconhecem o texto do sistema e mantêm
os argumentos intactos. Apenas os parâmetros declarados como mensagens internas
em `SYSTEM_ARGUMENTS` são traduzidos recursivamente. Marque conteúdo do usuário
com `verbatim`; não traduza valores editáveis nem chaves de opções de seleção.

```bash
PYTHONPATH=src python scripts/audit_localization.py
PYTHONPATH=src:tests python -m unittest test_localization test_localization_catalogs test_runtime_messages
```

A auditoria percorre expressões Python do BioMolExplorer e verifica frases
portuguesas reconhecíveis, templates e parâmetros dos seis grupos de dicionários.
Os testes também percorrem todos os formulários e mensagens dos splashes, verificam
sessões independentes, troca de idioma, overlays, tooltips, conteúdo protegido e
carregamento de um idioma de teste. Essa análise estática não substitui inspeção
visual, nem reconhece automaticamente toda frase possível ou texto de ferramentas
externas. Textos com composição nova devem receber um caso de teste quando sua
tradução depender de argumentos internos.
