import pandas as pd
import networkx as nx
from networkx.algorithms import community
import statistics
import math

import os


class module:
    def __init__(self, args) -> None:
        """
            this class is for computing modularity features of mutation modules

            features included:
                Var: Id of Variants
                CS: Conservation Score of Res (from original graph)
                CSS: Average css score of the subgraph
                Kw: Average Weighted Node Degree
                P: PageRank Centrality
                D: Degree Centrality
                B: Betweeness Centrality (from original graph, node_scale)
                C: Closeness Centrality
                Q: Modularity (subgraph-scale)
                QC: QC score defined in SAveRUNNER (subgraph-scale)
                P_Score: P Score (node-scale)
                Sub_Cen: Subgraph Centrality (node-scale)
                Num_cc: number of cc (subgraph-scale)
                Max_cc_pro: ratio of graph nodes in maximum cc
        """
        self.input_dir = r'{0}/procon'.format(args.pccn_dir)    # input_dir: the root_dir of input data
        self.out_dir = args.out_dir
        self.feat_names = [
            'Var','CS', 'CSS', 'Kw', 'P', 'D', 'B', 'C', 'Q', 'QC', 
            'P_Score', 'Sub_Cen', 'Num_cc', 'Max_cc_ratio'
        ]
        self.feats = pd.DataFrame(columns=self.feat_names)

    def aggregate(self):
        input_folder = sorted(os.listdir(self.input_dir))
        for sub_dir in input_folder:
            if sub_dir.endswith("_threshold_100"):
                var_name = list(sub_dir.split('_'))[0]
                nodes_file = r'{0}/{1}/{2}_node_info.csv'.format(self.input_dir, sub_dir, var_name)
                nodes_info_file = r'{0}/{1}/{2}_network_node_info.csv'.format(self.input_dir, sub_dir, var_name)
                edges_info_file = r'{0}/{1}/{2}_network_edge_info.csv'.format(self.input_dir, sub_dir, var_name)
                self.compute_module_feature(nodes_file, nodes_info_file, edges_info_file, var_name)

        out_file = r'{0}/module_feature.csv'.format(self.out_dir)
        self.feats.to_csv(out_file, index=False)

    def compute_module_feature(self, nodes_file, nodes_info_file, edges_info_file, var_name):
        feat_dict = {'Var': var_name.split('.fasta')[0]}
        nodes = pd.read_csv(nodes_file)
        network_nodes = pd.read_csv(nodes_info_file)
        network_edges = pd.read_csv(edges_info_file)

        all_nodes = set(network_nodes.iloc[:, 0])
        mut_nodes = set(network_nodes.iloc[list(network_nodes['is_mutation']), 0])
        non_mut_nodes = all_nodes - mut_nodes

        # average score for the mutation block
        # CS: conservation score; 
        mut_nodes_info = nodes.iloc[list(network_nodes['is_mutation']), :]
        mut_nodes_score = mut_nodes_info.mean(axis=0).to_dict()
        feat_dict.update(mut_nodes_score)

        # extract mutation edges
        is_source_mut = [(True if source in mut_nodes else False) for source in list(network_edges['source'])]
        is_target_mut = [(True if target in mut_nodes else False) for target in list(network_edges['target'])]
        mut_edges_bool = [is_source_mut[i] and is_target_mut[i] for i in range(len(is_source_mut))]
        mutation_edges = network_edges.iloc[mut_edges_bool, :]

        # create nx graph of ori & mut graph
        ori_graph = self.build_graph(network_edges)
        mut_graph = self.build_graph(mutation_edges)

        mut_max_cc, mut_cc_num, mut_max_cc_pro = self.extract_max_components(mut_graph)
        feat_dict.update({'Num_cc': mut_cc_num, 'Max_cc_ratio': mut_max_cc_pro})
        mut_max_subgraph = nx.subgraph(mut_graph, mut_max_cc)

        # compute modularity(Q) of max cc
        Q_value = nx.community.modularity(ori_graph, [mut_nodes, non_mut_nodes])
        feat_dict['Q'] = Q_value
        # non_mut_nodes = all_nodes - mut_max_cc
        # Q_value = nx.community.modularity(ori_graph, [mut_max_cc, non_mut_nodes])

        # compute QC
        QC_value = self.compute_qc(M=mut_max_subgraph, G=ori_graph)
        feat_dict['QC'] = QC_value

        # compute sub-graph centrality
        sub_centrality = self.average_subgraph_centrality(mut_max_subgraph)
        feat_dict['Sub_Cen'] = sub_centrality

        # compute P-score
        p_score, _ = self.average_p_score(mut_max_subgraph, ori_graph)
        feat_dict['P_Score'] = p_score

        # compute CSS
        css = self.average_CSS(mut_max_subgraph)
        feat_dict['CSS'] = css

        pre_feats = pd.DataFrame(feat_dict, index=[0])
        self.feats = pd.concat([self.feats, pre_feats], axis=0)


    def extract_max_components(self, G:nx.Graph):
        """
            calculate number of cc and extract the largest one
            input: G, a nx graph
            return:
                max_component: a set of nodes belonging to the maximum cc
                num_components: number of cc
                max_component_pro: ratio of graph nodes in the maximum cc
        """
        # Calculate the connected components
        connected_components = list(nx.connected_components(G))
        num_components = len(connected_components)
        num_total_nodes = G.number_of_nodes()

        # Print the number of nodes in each component
        component_id = 0
        component_max_node = 0
        for i, component in enumerate(connected_components):
            component_size = len(component)
            if component_size > component_max_node:
                component_id = i
                component_max_node = component_size

        max_component = connected_components[component_id]
        max_component_pro = component_max_node / num_total_nodes
        return max_component, num_components, max_component_pro


    def build_graph(self, edges:pd.DataFrame):
        """
            build nx graph from pandas data.frame
        """
        source = list(edges['source'])
        target = list(edges['target'])
        weight = list(edges['weight'])

        # construct mutation graph using networkX
        graph = nx.DiGraph()
        for source, target, weight in zip(source, target, weight):
            graph.add_edge(source, target, weight=weight)

        # convert to undirected graph
        graph = graph.to_undirected()
        return graph
    

    def average_subgraph_centrality(self, sub_graph:nx.Graph):
        """
            compute average subgraph centrality
            input: sub_graph, the target subgraph
            return: average subgraph centrality
        """
        sub_centrality = nx.centrality.subgraph_centrality(sub_graph)
        values = list(sub_centrality.values())
        average = math.log10(sum(values) / len(values))
        return average


    def average_z_score(self, sub_graph:nx.Graph):
        """
            compute average Z_score of the subgraph(based on degree)
            input:
                sub_graph: the target subgraph
            return:
                average_z_score: the average Z_score of the subgraph
                z_score_dict: a dict storing Z_score of every nodes in the subgraph
        """
        degrees = list(dict(sub_graph.degree).values())
        mean = statistics.mean(degrees)
        std = statistics.stdev(degrees)
        z_score = [(degree - mean) / std for degree in degrees]
        ave_z_score = statistics.mean(z_score)
        # also record zi
        keys = list(sub_graph.nodes())
        z_score_dict = dict(zip(keys, z_score))
        return ave_z_score, z_score_dict
    

    def average_p_score(self, sub_graph:nx.Graph, graph:nx.Graph):
        '''
            compute average P_score of the subgraph
            input:
                sub_graph: the target subgraph
                graph: the original graph
            return:
                average_p_score: the average P_score of the subgraph
                p_score_dict: a dict storing P_score of every nodes in the subgraph
        '''
        sub_degrees = sub_graph.degree()
        ori_degrees = graph.degree()

        p_score_dict = {}
        for node in sub_graph.nodes():
            in_degree = sub_degrees[node]
            out_degree = ori_degrees[node]
            p_score = 1 - (in_degree / out_degree)**2       # p-score = 1 - (Kin / K)^2
            p_score_dict[node] = p_score
        ave_p_score = statistics.mean(list(p_score_dict.values()))

        return ave_p_score, p_score_dict
    

    def compute_qc(self, M:nx.Graph, G:nx.Graph, lamb=1.0):
        '''
            compute QC score defined in SAveRUNNER
            input:
                M: the target subgraph
                G: original graph
                lamb: Hyper-parameters adjusting weight of Temp
            intermediate paramters:
                Win: total weights of edges within the subgraph
                Wout: total weights of edges connecting the subgraph to rest part of the graph
                Temp: Temparature factor, equals to the ratio of network nodes within the subgraph
            return:
                QC score, a float in (0,1)
        '''
        Win = sum(G.get_edge_data(u, v)['weight'] for u, v in M.edges())
    
        Wout = 0
        for node in M.nodes():
            neighbors = set(G.neighbors(node)) - set(M.nodes())
            for neighbor in neighbors:
                Wout += G.get_edge_data(node, neighbor)['weight']
        
        Temp = lamb * M.number_of_nodes() / G.number_of_nodes()
        ratio = Win / (Win + Wout + Temp)
        return ratio

    def average_CSS(self, graph:nx.Graph):
        total_weight = sum([graph.edges[edge]['weight'] for edge in graph.edges])
        num_edges = graph.number_of_edges()
        average_weight = total_weight / num_edges
        return average_weight



    

