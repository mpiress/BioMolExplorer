from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
        "Module used to generate graph similarity relationship "
        "between molecules and similar extracted from ChEMBL."
    ),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
import warnings

import matplotlib
matplotlib.use('agg')
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
import networkx as nx
import re
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import seaborn as  sns
import random
import os

from typing import Optional
from pathlib import Path
from pandas import DataFrame, concat, read_csv
from statistics import mean
from collections import Counter
from networkx import Graph
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.utilities import fileHandling
from kernel.descriptors import Descriptors, similarityFunctions, fingerprints
from kernel.loggers import LoggerManager
from biomolexplorer.visualizations import graph_view, write_view, SUFFIX
#----------------------------------------------------------------------------------------------


class GraphAnalysis():

    def __init__(self, input_path:Optional[str]=None, output_path:Optional[str]=None) -> None:
        self.__path = str(Path.cwd())
        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/graph_analizer.log')
        self.set_inputpath(input_path)
        self.set_outputpath(output_path)


    def set_inputpath(self, inputpath:str) -> None:
        self.__inputpath = inputpath



    def set_outputpath(self, outputpath:str) -> None:
        self.__outputpath = outputpath  if outputpath != None else None
        if self.__outputpath != None and not os.path.exists(directory(self.__outputpath)):
           os.makedirs(directory(self.__outputpath), exist_ok=True)



    def count_nodes_by_degree(self, graph:Graph) -> dict:
        degrees = dict(graph.degree())
        degree_counts = Counter(degrees.values())
        return degree_counts



    def get_statisticals(self, G) -> dict:

        try:
            resp = {'vertices':0, 'arestas':0, 'Grau Medio Geral':0, 'Grau Medio Vizinhos':0,
                    'Densidade':0, 'Qtd Compoentes Conectados':0, 'Coeficiente de Cluster Medio':0, 'vertices isolados':0}

            resp['vertices'] = nx.number_of_nodes(G)
            resp['arestas'] = nx.number_of_edges(G)
            resp['vertices isolados'] = nx.number_of_isolates(G)
            resp['Grau Medio Vizinhos'] = nx.average_neighbor_degree(G)
            resp['Grau Medio Geral'] = mean(dict(G.degree()).values()) if G else 0
            resp['Densidade']  = nx.density(G)
            resp['Qtd Compoentes Conectados'] = nx.number_connected_components(G)
            resp['Coeficiente de Cluster Medio'] = nx.average_clustering(G) if G else 0

            return resp

        except Exception as e:
            self.logger.error(f'Error during to perform the get_statisticals function', exc_info=True)
            raise



    def plot_fruchterman_reingold(self, G, file_name:str, output_path:str):

        random.seed(42)

        try:
            componentes = list(nx.connected_components(G))
            cores = [plt.cm.viridis(random.random()) for _ in range(len(componentes))]

            mapeamento_cores = {}
            for i, componente in enumerate(componentes):
                for vertice in componente:
                    mapeamento_cores[vertice] = cores[i]

            posicoes = nx.spring_layout(G)

            for vertice, posicao in posicoes.items():
                cor = mapeamento_cores[vertice]
                nx.draw_networkx_nodes(G, posicoes, nodelist=[vertice], node_color=cor, node_size=10)

            nx.draw_networkx_edges(G, posicoes, alpha=0.5)

            plt.axis('off')
            os.makedirs(output_path) if not os.path.exists(output_path) else None

            plt.savefig(output_path+file_name)
            plt.clf()

        except Exception as e:
            self.logger.error(f'Error during to perform the plot_fruchterman_reingold function', exc_info=True)
            raise




    def plot_statistical_degree_analysis(self, G:Graph, MaxG:Graph, file_name:str, output_path:str):

        try:
            degree_sequence = sorted((d for _, d in MaxG.degree()), reverse=True)
            num_vertices = MaxG.number_of_nodes()
            num_arestas = MaxG.number_of_edges()

            fig = plt.figure("Degree Analysis", figsize=(10, 12))

            # Create a gridspec for adding subplots of different sizes
            sns.set_theme(style="whitegrid")
            axgrid = fig.add_gridspec(7, 4)

            ax0 = fig.add_subplot(axgrid[0:3, :])
            pos = nx.kamada_kawai_layout(MaxG)

            # Calcular os graus dos nós
            degrees = dict(MaxG.degree())
            node_colors = [degrees[node] for node in MaxG.nodes()]
            # Criar um objeto ScalarMappable para mapear os valores dos nós para as cores
            sm = ScalarMappable(cmap=plt.cm.viridis)
            sm.set_array(node_colors)

            nx.draw_networkx_nodes(MaxG, pos, node_color=node_colors, cmap=plt.cm.viridis, ax=ax0, node_size=20)
            nx.draw_networkx_edges(MaxG, pos, ax=ax0, alpha=0.1, edge_color='gray')

            # Adicione texto para representar o número de vértices e arestas
            plt.text(0.5, -0.1, f'Amount of Compounds: {num_vertices}\nNumber of Relationships: {num_arestas}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            ax0.set_title("Analysis and Overview by Maximum Connected Component (MaxComp Degree)")
            ax0.set_axis_off()

            # Adicionar uma barra de cores
            cbar = plt.colorbar(sm, ax=ax0)
            cbar.set_label('Node Degree')

            ax1 = fig.add_subplot(axgrid[3:5, :2])
            ax1.locator_params(axis='both', integer=True)
            sns.lineplot(data=degree_sequence, ax=ax1)
            ax1.set_title("Degree Rank Plot")
            ax1.set_ylabel("Degree")
            ax1.set_xlabel("Rank")
            ax1.xaxis.grid(False)
            ax1.yaxis.grid(False)


            ax2 = fig.add_subplot(axgrid[3:5, 2:])
            ax2.locator_params(axis='both', integer=True)
            sns.histplot(degree_sequence, ax=ax2, kde=False)
            ax2.set_title("Degree histogram")
            ax2.set_xlabel("Degree")
            ax2.set_ylabel("Number of Nodes")
            ax2.xaxis.grid(False)
            ax2.yaxis.grid(False)

            # Adicionando o gráfico desejado como ax3
            ax3 = fig.add_subplot(axgrid[5:, :])
            ax3.locator_params(axis='both', integer=True)

            degree_counts = self.count_nodes_by_degree(G)
            degrees, counts = zip(*degree_counts.items()) if degree_counts else ([], [])
            degrees = list(degrees)
            counts  = list(counts)

            sns.scatterplot(x=degrees, y=counts, color='gray', s=50, ax=ax3, label='Graph Degree')

            degree_counts = self.count_nodes_by_degree(MaxG)
            degrees, counts = zip(*degree_counts.items()) if degree_counts else ([], [])
            degrees = list(degrees)
            counts  = list(counts)

            sns.scatterplot(x=degrees, y=counts, color='#7995c4', s=50, ax=ax3, label='MaxComp Degree')

            plt.xlabel('Degree')
            plt.ylabel('Number of Nodes')
            num_vertices = G.number_of_nodes()
            num_arestas = G.number_of_edges()
            plt.title(f'Completely Graph Relationships Degree Distribution\nAmount of Compounds: {num_vertices}   Number of Relationships: {num_arestas}')
            plt.grid(True, linestyle='--', linewidth=0.5, color='lightgrey', which='both', axis='both')
            plt.locator_params(axis='both', integer=True)
            #plt.legend(loc='upper center', ncol=2)

            fig.tight_layout()

            os.makedirs(output_path) if not os.path.exists(output_path) else None

            plt.savefig(output_path+file_name)
            plt.close(fig)

        except Exception as e:
            self.logger.error(f'Error during to perform the plot_statistical_degree_analysis function', exc_info=True)
            raise



    def max_conected_component(self, df:DataFrame, nodes=None) -> tuple:

        try:

            G = nx.from_pandas_edgelist(df, source='source', target='target', edge_attr='value')
            if nodes is not None:
                G.add_nodes_from(nodes)

            component     = min(nx.connected_components(G), key=lambda group: (-len(group), min(map(str, group))), default=set())
            max_component = G.subgraph(component)

            edges     = DataFrame(max_component.edges(), columns=['source', 'target'])
            component = DataFrame(G.degree(component), columns=['molecule_chembl_id', 'degree'])

            return (max_component, component, edges, G)

        except Exception as e:
            self.logger.error(f'Error during to perform the max_conected_component function', exc_info=True)
            raise

        return tuple()




    def prepare_graph_analysis(self, metric:Optional[similarityFunctions]=similarityFunctions.TanimotoSimilarity,
                               fp:Optional[fingerprints]=fingerprints.Morgan, similarity_path:Optional[str]=None):
        """Keep the legacy entry point while using the shared graph exporter."""
        from caad.graph_results import export_analysis, read_compounds
        try:
            prefix = metric.value + '_' + fp.value + '_'
            folder = Path(similarity_path or self.__inputpath + 'Similarity/')
            files = sorted(folder.glob(prefix + '*.csv'))
            if not files:
                raise ValueError('No similarity CSV files found for the selected metric and fingerprint')
            selected = []
            fragments = []
            for path in files:
                dataset = read_compounds([Path(self.__inputpath)/(path.stem.removeprefix(prefix)+'.csv')])
                edges = read_csv(path, dtype={'source':str,'target':str})
                model, frame = export_analysis(dataset,edges,self.__outputpath,path.stem,path.stem)
                selected.append(frame)
                fragments.append({'id':path.stem,'smarts':model['fragment']['smarts'],
                                  'status':model['fragment']['status']})
            destination=Path(self.__outputpath)/'Molecules'/'molecules.csv'
            concat(selected,ignore_index=True).drop_duplicates('molecule_chembl_id').to_csv(destination,index=False)
            DataFrame(fragments).to_csv(Path(self.__outputpath)/'centroids'/'fragments.csv',index=False)
        except Exception:
            self.logger.error('during to perform the prepare_graph_analysis function',exc_info=True)
            raise
