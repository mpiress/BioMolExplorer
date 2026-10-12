"""Explanations of the existing ranking metrics; no changes to their computation."""
from functools import lru_cache

HELP = {
    'degree': (
        'Conta quantos vizinhos diretos um composto possui no universo analisado.',
        'Ajuda a encontrar compostos com muitas relações diretas e a distinguir nós pouco conectados.',
        'B tem duas conexões; A e C têm uma. B fica à frente no ranking decrescente.',
        'Um grau alto descreve conectividade; não confirma diversidade estrutural nem atividade biológica.'),
    'density': (
        'Divide o grau por n − 1, onde n é o número de compostos do universo. É centralidade de grau, não a densidade global do grafo.',
        'Mostra a fração dos outros compostos aos quais um nó está ligado. Considere o tamanho e o universo ao comparar análises.',
        'Há três compostos: B se conecta aos outros dois e recebe 1; A e C se conectam a um e recebem 0,5.',
        'O valor muda entre Full Graph e LCC porque o denominador depende do universo, mesmo quando o grau permanece igual.'),
    'betweenness': (
        'Mede a participação de um nó nos menores caminhos entre outros pares de compostos.',
        'Ajuda a localizar pontes topológicas entre regiões ou famílias da rede para explorar transições entre grupos.',
        'O único caminho de A até C passa por B. B recebe 1; A e C recebem 0.',
        'A plataforma ignora pesos de escore. Acima de 400 nós usa até 64 pivôs, com semente 42; a aproximação aparece no painel. Acima de 5.000 nós ou 50.000 arestas fica indisponível.'),
    'closeness': (
        'Mede a proximidade aos demais nós por caminhos topológicos. Caminhos mais curtos aumentam o valor; a correção Wasserman–Faust considera componentes desconectados.',
        'Ajuda a localizar compostos que alcançam os demais com poucas etapas na rede, considerando o universo escolhido.',
        'B alcança A e C em uma etapa: valor 1. A alcança B em uma e C em duas: valor 2/3.',
        'Etapas na rede não são distância química nem tempo de reação. Acima de 1.500 nós ou 20.000 arestas a métrica fica indisponível.'),
    'eigenvector': (
        'Dá mais importância às conexões com vizinhos que também são centrais. A plataforma normaliza o vetor para comprimento 1.',
        'Ajuda a explorar compostos integrados ao núcleo topológico de uma família, além de simplesmente contar vizinhos.',
        'B está ligado a A e C: seu valor é aproximadamente 0,707. A e C têm aproximadamente 0,500 cada.',
        'Exige um universo conectado: escolha LCC quando Full Graph tiver vários componentes. Fica indisponível acima de 50.000 arestas ou quando não converge em 300 iterações.'),
    'pagerank': (
        'Distribui importância pela rede como uma caminhada aleatória: em 85% das etapas segue uma conexão e nas demais reinicia uniformemente.',
        'Ajuda a priorizar compostos que recebem importância de seus vizinhos. A contribuição de cada vizinho é dividida entre suas conexões.',
        'A e C enviam sua contribuição a B; B a divide entre A e C. B recebe aproximadamente 0,486 e cada ponta 0,257.',
        'Usa amortecimento 0,85 e ignora pesos de escore. Fica indisponível acima de 20.000 nós ou 50.000 arestas, ou sem convergência em 200 iterações.'),
}


@lru_cache(maxsize=6)
def toy_values(metric):
    """Use the same scientific implementation as the actual ranking."""
    from biomolexplorer.network_analysis import NetworkAnalysis
    data={'nodes':[{'id':n} for n in ('A','B','C')],
          'edges':[{'source':'A','target':'B'},{'source':'B','target':'C'}]}
    return NetworkAnalysis(data).metrics('full',metric)[0]
