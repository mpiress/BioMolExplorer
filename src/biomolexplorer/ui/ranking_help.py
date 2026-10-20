"""Explanations of the existing ranking metrics; no changes to their computation."""
from functools import lru_cache

HELP = {
    'core': ('É o maior k para o qual o composto pertence a um subgrafo com grau mínimo k.',
        'Identifica núcleos densamente conectados, mesmo quando o grau de um nó isoladamente é alto.',
        'Em A — B — C todos pertencem ao núcleo 1.',
        'Calculado sem laços próprios e sem pesos. Um núcleo não demonstra atividade biológica.'),
    'clustering': ('Fração dos pares de vizinhos que estão conectados entre si.',
        'Ajuda a distinguir vizinhanças coesas de conexões entre grupos separados.',
        'Em A — B — C não há triângulos: os três coeficientes são zero.',
        'Sem pesos. Interseções orientadas por grau, com orçamento de 20 milhões de comparações; reduza o escopo se excedido.'),
    'triangles': ('Conta os triângulos que passam por cada composto.',
        'Quantifica relações fechadas locais entre três compostos.',
        'Em A — B — C nenhum composto participa de um triângulo.',
        'Contagem exata sem pesos. Usa interseções orientadas por grau, até 20 milhões de comparações estimadas.'),
    'degree': (
        'Conta quantos vizinhos diretos um composto possui no universo analisado.',
        'Ajuda a encontrar compostos com muitas relações diretas e a distinguir nós pouco conectados.',
        'B tem duas conexões; A e C têm uma. B fica à frente no ranking decrescente.',
        'Um grau alto descreve conectividade; não confirma diversidade estrutural nem atividade biológica.'),
    'density': (
        'Divide o grau por n − 1, onde n é o número de compostos do universo. É centralidade de grau, não a densidade global do grafo.',
        'Mostra a fração dos outros compostos aos quais um nó está ligado. Considere o tamanho e o universo ao comparar análises.',
        'Há três compostos: B se conecta aos outros dois e recebe 1; A e C se conectam a um e recebem 0,5.',
        'O valor muda quando o escopo muda: Full, MCC ou grafo exibido com filtros. O denominador usa o número de nós desse escopo.'),
    'betweenness': (
        'Mede a participação de um nó nos menores caminhos entre outros pares de compostos.',
        'Ajuda a localizar pontes topológicas entre regiões ou famílias da rede para explorar transições entre grupos.',
        'O único caminho de A até C passa por B. B recebe 1; A e C recebem 0.',
        'Sem pesos. Árvores, cliques e componentes pequenos usam cálculo exato; os demais usam até 64 pivôs por componente, semente 42. O painel informa a aproximação. O orçamento é de 20 milhões de visitas estimadas; a normalização usa o escopo inteiro.'),
    'closeness': (
        'Mede a proximidade aos demais nós por caminhos topológicos. Caminhos mais curtos aumentam o valor; a correção Wasserman–Faust considera componentes desconectados.',
        'Ajuda a localizar compostos que alcançam os demais com poucas etapas na rede, considerando o universo escolhido.',
        'B alcança A e C em uma etapa: valor 1. A alcança B em uma e C em duas: valor 2/3.',
        'Cálculo exato por componente, com correção Wasserman–Faust pelo escopo inteiro. Árvores e cliques usam fórmulas diretas; outros componentes usam caminhos em lotes, até 200 milhões de visitas estimadas. Etapas topológicas não são distância química.'),
    'eigenvector': (
        'Dá mais importância às conexões com vizinhos que também são centrais. A plataforma normaliza o vetor para comprimento 1.',
        'Ajuda a explorar compostos integrados ao núcleo topológico de uma família, além de simplesmente contar vizinhos.',
        'B está ligado a A e C: seu valor é aproximadamente 0,707. A e C têm aproximadamente 0,500 cada.',
        'Exige um escopo conectado: escolha MCC ou filtre um componente e use Grafo exibido. O solver esparso tem tolerância 1e-8 e até 2.000 iterações. A falha de convergência é informada; não há bloqueio fixo em 50.000 arestas.'),
    'pagerank': (
        'Distribui importância pela rede como uma caminhada aleatória: em 85% das etapas segue uma conexão e nas demais reinicia uniformemente.',
        'Ajuda a priorizar compostos que recebem importância de seus vizinhos. A contribuição de cada vizinho é dividida entre suas conexões.',
        'A e C enviam sua contribuição a B; B a divide entre A e C. B recebe aproximadamente 0,486 e cada ponta 0,257.',
        'Usa matrizes esparsas, amortecimento 0,85, tolerância 1e-8 e até 1.000 iterações. Ignora pesos de escore; a falha de convergência é informada, sem bloqueio fixo por contagem de nós ou arestas.'),
}


@lru_cache(maxsize=6)
def toy_values(metric):
    """Use the same scientific implementation as the actual ranking."""
    from biomolexplorer.network_analysis import NetworkAnalysis
    data={'nodes':[{'id':n} for n in ('A','B','C')],
          'edges':[{'source':'A','target':'B'},{'source':'B','target':'C'}]}
    return NetworkAnalysis(data).metrics('full',metric)[0]
