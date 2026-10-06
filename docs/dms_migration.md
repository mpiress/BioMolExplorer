# Porte nativo do algoritmo DMS

[Documentação](README.md) · Português · [English](en/dms_migration.md)

O BioMolExplorer gera agora a superfície molecular em Python. `Dock6.prepare_surface()` chama `biomolexplorer.molecular_surface.generate_surface()` diretamente, sem executar ou instalar `dms`/`dmsd`. O arquivo `{receptor}.dms`, os parâmetros de densidade e sonda e a sequência `sphgen → sphere_selector → showbox → grid → DOCK6` permanecem no fluxo.

## Código original estudado e licença

Foi baixada e compilada a [distribuição oficial UCSF](https://www.cgl.ucsf.edu/Overview/ftp/dms.zip), publicada na [página de software UCSF](https://www.cgl.ucsf.edu/Overview/software.html). SHA-256 do ZIP: `0699283fc4902d4073fba7ccb393e830c9d8b9e1a8b2c4a883b9f0dd102dc011`.

| Fonte original | Responsabilidade portada |
| --- | --- |
| `ms.c` | Parâmetros e intervalos aceitos no caminho utilizado |
| `input.c` e `radii.proto` | Seleção PDB, identificadores, raios por prefixo do nome do átomo |
| `compute.c` | Ordenação, sequência dos cálculos, vizinhos removidos e colapso de sondas coincidentes |
| `dmsd/server.c` | Vizinhos, centros de sonda, oclusão, amostragem esférica, superfícies de contato, toroidais e côncavas |
| `dmsd/viewat.c` e `lookat.c` | Transformações entre os referenciais geométricos |
| `output.c` | Formatação dos átomos/pontos, tipos de superfície, área e normais |

O original permite redistribuição e modificação com atribuição. A licença integral está em `src/biomolexplorer/resources/dms/LICENSE` e acompanha o pacote Python. O software foi originalmente desenvolvido pelo UCSF Computer Graphics Laboratory, com apoio do NIH National Center for Research Resources, grant P41-RR01081. O porte derivado conserva essa licença; a licença geral MIT do projeto não substitui os termos desse componente.

## Funcionalidades implementadas

O módulo calcula a superfície excluída ao solvente por rolamento de uma sonda, sem substituir a geometria por uma estimativa SASA. Inclui:

- Regiões convexas de contato com cada átomo (`SC0`).
- Regiões toroidais entre pares, inclusive tratamento de toros spindle (`SS0`).
- Regiões côncavas esféricas entre trios (`SR0`), com oclusão e colapso de sondas quase coincidentes.
- Associação dos pontos aos átomos segundo as regras do original, áreas por ponto e normais orientadas.
- Raios UCSF empacotados, regras por prefixo do nome atômico e sobreposição por arquivo `radii` no diretório corrente ou caminho explícito.
- Leitura dos registros PDB elegíveis, cadeias/códigos de inserção e inclusão opcional de HETATM. Como `dms -a`, a inclusão opcional também aceita resíduos não padronizados.
- Escrita com o mesmo formato consumido por `sphgen`, substituição atômica do arquivo após sucesso e resumo de pontos/área para o log.

Foram mantidas a constante histórica `PI=3.141592`, a amostragem por camadas, a densidade interna `2.75 × densidade`, as regras de arredondamento e a precisão de seis casas do protocolo C. Essas escolhas evitam alterar a superfície final. Busca espacial `scipy.spatial.cKDTree` e operações NumPy aceleram vizinhança e filtragem.

O porte cobre integralmente o caminho que a aplicação utilizava: `dms receptor.noH.pdb -d <densidade> -n -w <raio> -v -o receptor.dms`. A administração dos servidores DMS distribuídos, sua CLI completa e seletores de resíduos não usados pela aplicação não são reimplementados. A execução é local, dentro do worker Python.

## Uso e dependências

NumPy e SciPy já pertencem ao ambiente científico, portanto nenhuma biblioteca adicional é necessária. `requirements.yml` é o manifesto Conda disponível no checkout; o arquivo `environment.yml` já estava removido antes deste porte e não foi recriado. O instalador não baixa, compila nem instala mais DMS.

```python
from biomolexplorer.molecular_surface import generate_surface

summary = generate_surface(
    "receptor.noH.pdb", "receptor.dms",
    density=0.5, probe_radius=1.4,
)
print(summary)
```

Parâmetros: `density` de 0,1 a 10; `probe_radius` de 1 a 201 Å, como no original. Os padrões permanecem 0,5 e 1,4 Å. `normals=False` omite normais; `include_hetero=True` inclui todos os resíduos e HETATM; `radii_path` define uma tabela alternativa. A integração DOCK6 mantém normais ativadas e não inclui HETATM, assim como a chamada anterior.

## Evidência de equivalência

O C oficial foi compilado em diretório temporário. A única adaptação do fonte foi remover uma declaração antiga e não usada de `sbrk`, conflitante com os cabeçalhos atuais. Nenhuma equação ou regra geométrica foi modificada.

Os 31 pares PDB/superfície de referência estão em `tests/fixtures/dms`, com parâmetros e checksums em `manifest.json`. Eles cobrem esfera isolada, par, triângulo, tetraedro, toro spindle, átomos ocultos, 20 conjuntos de coordenadas com raios diferentes, quatro geometrias simétricas (quadrado, pentágono, cubo e octaedro) e a proteína real 1CRN (331 átomos). Os testes comparam a multiplicidade de cada registro e todos os seus campos. Apenas a ordem dos registros e o sinal de zero arredondado são normalizados: a ordem de conclusão dos servidores C não é uma garantia do algoritmo.

Para 1CRN nos padrões, ambos produziram 2.479 pontos de contato, 2.384 de sela e 1.390 côncavos (6.253 no total), com todos os campos equivalentes após a formatação. Comparações adicionais da proteína com `(density=0.1, probe_radius=1.0)` e `(1.0, 2.0)` também não apresentaram diferenças. No ambiente da revisão, o porte levou aproximadamente 8,6 segundos no caso padrão; isso é uma medição pontual, não uma garantia de desempenho para receptores grandes.

```bash
PYTHONPATH=src:tests python -m unittest test_molecular_surface test_docking_handoffs
PYTHONPATH=src python scripts/validate_native_surface.py --dms-executable /caminho/para/dms
```

A primeira verificação funciona sem DMS. A segunda serve somente para auditoria contra uma compilação independente. As instruções e a origem dos fixtures estão em `tests/fixtures/dms/README.md`. O teste de integração executa o gerador real e verifica que as chamadas externas restantes são `sphgen` e `sphere_selector`, com seus arquivos esperados.

## Limites da validação

Átomos colineares ou coincidentes recebiam tratamento indefinido no C, que pode escrever `NaN`. O porte trata explicitamente o bloqueio de círculos de sonda e não escreve coordenadas não finitas; testes separados cobrem esses casos. A equivalência com registros inválidos `NaN` do original não é uma meta.

Não foi executado `sphgen`/DOCK6 real nesta revisão porque seus executáveis não estão disponíveis no ambiente. A evidência demonstra equivalência dos dados da superfície nos casos testados e preservação do contrato de integração; não comprova identidade de scores, poses ou desempenho para qualquer receptor. A geração `.dms` foi substituída, mas Chimera, DOCK6 e seus auxiliares continuam sendo dependências dos demais passos existentes.
