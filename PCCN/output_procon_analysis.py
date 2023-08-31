import json
# 日志
import logging
import math
import os
import pickle
import re
import time
from collections import defaultdict
from functools import reduce
from itertools import combinations, permutations
from typing import List

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from networkx.algorithms.centrality import degree_centrality, betweenness_centrality, closeness_centrality, \
    edge_betweenness_centrality
from networkx.algorithms.shortest_paths import shortest_path_length
from pyecharts import options as opts
from pyecharts.charts import Graph
# 组合数
from scipy.special import comb
from scipy.stats import mannwhitneyu, ttest_1samp
from sklearn import preprocessing
from statannotations.Annotator import Annotator

import logconfig


logconfig.setup_logging()
log = logging.getLogger("cov2")

AA = ['-', 'C', 'F', 'I', 'L', 'M', 'V', 'W', 'Y', 'A', 'T', 'D', 'E', 'G', 'P', 'N', 'Q', 'S', 'H', 'R', 'K']

# sns.set()
sns.set_style("ticks")
plt.rcParams.update({'font.size': 16})
plt.rcParams["axes.titlesize"] = "medium"

# scaler
def scaler(d: dict):
    min_max_scaler = preprocessing.MinMaxScaler()
    keys = d.keys()
    # values = np.array(d.values(), dtype=float)
    values = np.array(list(d.values())).reshape(-1, 1)
    values = min_max_scaler.fit_transform(values)
    values = values.reshape(-1).tolist()
    res = dict(zip(keys, values))
    return res


class AnalysisMutationGroup:
    def __init__(self,floder_path ,analysis_res_path , file_name, seed=0, ):
        analysis = r"{0}/{1}_analysis.json".format(analysis_res_path,file_name)
        with open(analysis) as f:
            self.analysis: dict = json.load(f)
        self.file_name=file_name
        self.seed = seed
        self.fasta = next(SeqIO.parse('{0}/{1}'.format(floder_path,file_name), "fasta")).seq
        self.positions = [f"{i + 1}{aa.upper()}" for i, aa in enumerate(self.fasta)]

        self.non_duplicated_aas = self.get_non_duplicated_aas()
        self.non_duplicated_aas_positions = self.get_non_duplicated_aas_position()
        self.non_duplicated_aas_sample = self.resample_aas(self.non_duplicated_aas_positions)

        self.aa_groups, self.aa_groups_info = self.get_aa_groups()
        self.aa_groups_position = self.get_aa_groups_position()
        self.group_count_sample = self.resample_groups()

    @staticmethod
    def aa2position(aa: str):
        if aa[0] in AA and aa[-1] in aa and aa[1:-1].isdigit():
            return aa[1:-1] + aa[0]
        elif aa[:-1].isdigit() and aa[-1] in aa:
            return aa
        else:
            raise RuntimeError(f"无法转化aa={aa}")

    def get_aa_groups(self):
        aas = []
        names = []
        categroies = []
        for i, row in self.analysis.items():
            aas.append(row["aas"])
            categroies.append(row["WHO label"])
            # categroies.append(row["Category"])
            if type(row["WHO label"]) == str:
                # name = "{} ({})".format(row["Lineage + additional mutations"], row["WHO label"])
                name = "{}({})".format(row["Lineage + additional mutations"], row["WHO label"])
            else:
                name = "{}".format(row["Lineage + additional mutations"])
            names.append(name)

        info = pd.DataFrame({"name": names, "category": categroies})
        # color
        info["color"] = sns.hls_palette(n_colors=len(info), )
        return aas, info

    def get_aa_groups_position(self):
        groups, _ = self.get_aa_groups()
        groups = [[self.aa2position(aa) for aa in group] for group in groups]
        return groups

    def get_non_duplicated_aas(self):
        aas = []
        for i, row in self.analysis.items():
            aas += row["aas"]
        aas = list(set(aas))
        return aas

    def get_non_duplicated_aas_position(self):
        aas = self.get_non_duplicated_aas()
        positions = [self.aa2position(aa) for aa in aas]
        positions = list(set(positions))
        return positions

    def resample_aas(self, aas, N=1000):
        # resample
        positions = list(set(self.positions) - set(aas))
        positions = sorted(positions)
        positions = pd.Series(positions)
        log.debug("count %s", N)
        # sample_aas = [positions.sample(n=len(aas), random_state=self.seed + i).tolist() for i in range(N)]
        sample_aas = [positions.sample(n=len(aas), random_state=self.seed + i).tolist() for i in range(N)]
        return sample_aas

    def resample_groups(self, N=1000):
        """
        resample
        """
        groups = self.aa_groups
        positions = list(set(self.positions) - set(self.non_duplicated_aas_positions))
        positions = sorted(positions)
        positions = pd.Series(positions)

        sample_groups = {}
        group_counts = np.unique([len(group) for group in groups])
        log.debug("count %s", N)
        for count in group_counts:
            one_sample = [positions.sample(n=count, random_state=self.seed + j + count * 1000).tolist() for j in
                          range(N)]
            sample_groups[count] = one_sample
            log.debug(f"group len =  {count}")
        return sample_groups

    def display_seq_and_aa(self):
        aas = self.non_duplicated_aas
        aas = sorted(aas, key=lambda x: int(x[1:-1]))
        log.info(f"fasta({len(self.fasta)}): {self.fasta}")
        count = [len(group) for group in self.aa_groups_position]
        log.debug("pd.value_counts(count).sort_index() = %s", pd.value_counts(count).sort_index())

        all_aa = reduce(lambda x, y: x + y, self.aa_groups_position)
        log.debug("all aa count = %s", pd.value_counts(all_aa))

        log.info("number of variation: %s", len(self.aa_groups))
        log.info("number of aas: %s", len(self.non_duplicated_aas))
        log.info("number of site: %s", len(self.non_duplicated_aas_positions))

    def count_aa(self):
        group = self.aa_groups
        group_info = self.aa_groups_info
        names = group_info["name"].to_list()

        result = {}
        for i in range(len(group)):
            name = names[i]
            N = len(group[i])
            if N not in result:
                result[N] = {"Count": 0, "Variant": []}
            result[N]["Count"] += 1
            result[N]["Variant"].append(name)

        for k in result.keys():
            result[k]["Variant"] = ", ".join(result[k]["Variant"])

        result = pd.DataFrame(result).T
        result = result.sort_index()
        result.index.name = "Number of substitutions"
        result = result.rename(columns={"Count": "Number of variants"})
        result.to_excel(os.path.join("/PCCN/data/procon", "aa count.xlsx"))


class ProConNetwork:
    def __init__(self,file_name,type_file,data_dir,
                 analysis_mutation_groups: AnalysisMutationGroup,
                 threshold=100,
                 ):
        self.file_name=file_name
        parse1 = r"{2}/{0}_type/{1}_type1_parse.csv".format(file_name,file_name,type_file)
        parse2 = r"{2}/{0}_type/{1}_type2_parse.csv".format(file_name,file_name,type_file)
        self.analysis_mutation_group = analysis_mutation_groups

        log.info("构造图...")
        # read data
        self.type1 = pd.read_csv(parse1)
        self.type2 = pd.read_csv(parse2)
        self.type1_mms, self.type1["info_norm"] = self._normalize_info(self.type1["information"].values)
        self.type2_mms, self.type2["info_norm"] = self._normalize_info(self.type2["info"].values)

        if threshold > 1:
            self.threshold = threshold
        else:
            count_type2 = int(threshold * len(self.type2))
            threshold_score = self.type2["info"].sort_values(ascending=False)[count_type2]
            self.threshold = threshold_score
        log.debug("self.threshold = %s", self.threshold)

        # ouput dir
        self.data_dir = os.path.join(data_dir, "{0}_threshold_{1}".format(file_name,threshold))
        if not os.path.exists(self.data_dir):
            os.mkdir(self.data_dir)
        log.debug("self.data_dir = %s", self.data_dir)

        # fasta position
        self.fasta = self.analysis_mutation_group.fasta
        self.positions = self.analysis_mutation_group.positions

        # construct network
        nodes = self._get_nodes(self.type1)
        links = self._get_links(self.type2[self.type2["info"] >= self.threshold])
        self.G = self._get_G(links, nodes)
        self.G = self._add_neighbour_links(self.G)

        log.debug(f"\n\ttype1 {len(self.type1)}\n\tnode {self.G.number_of_nodes()}\n\t"
                  f"type2 {len(self.type2)}\n\tedge {self.G.number_of_edges()}\n\t"
                  f"neighbour edge {self.G.number_of_edges() - len(links)}")

        # centrality
        self.degree_c, self.betweenness_c, self.closeness_c, self.edge_betweenness_c = self._get_centralities()
        self.degree = self.calculate_degree()

        # page rank
        log.info("page rank ...")
        self.page_rank = nx.pagerank(self.G)
        log.info("shortest path length ...")
        self.shortest_path_length = self.get_shortest_path_length()

    @staticmethod
    def _normalize_info(info: np.ndarray):
        """
        normalize
        """
        info = info.reshape(-1, 1)
        mms = preprocessing.MinMaxScaler().fit(info)
        return mms, mms.transform(info)

    @staticmethod
    def _aa2position(aa: str):
        if aa[0] in AA and aa[-1] in aa and aa[1:-1].isdigit():
            return aa[1:-1] + aa[0]
        elif aa[:-1].isdigit() and aa[-1] in aa:
            return aa
        else:
            raise RuntimeError(f"error aa={aa}")

    def _get_nodes(self, type1):
        nodes = self.positions
        nodes = pd.DataFrame({"position": nodes})
        nodes = pd.merge(nodes, type1, how="left", left_on="position", right_on="position")
        nodes = nodes.loc[:, ["position", "info_norm"]]
        nodes.columns = ["name", "size"]
        nodes = nodes.fillna(0)
        return nodes

    def _get_links(self, type2, ):
        links = type2.loc[:, ["site1", "site2", "info_norm"]]
        links.columns = ["source", "target", "weight"]
        return links

    def _add_neighbour_links(self, G):
        neighbour_links = []
        for i in range(len(self.positions) - 1):
            node1 = self.positions[i]
            node2 = self.positions[i + 1]
            if not G.has_edge(node1, node2):
                neighbour_links.append([node1, node2, self.threshold])

        neighbour_links = pd.DataFrame(neighbour_links, columns=["source", "target", "weight"])
        neighbour_links["weight"] = self.type2_mms.transform(neighbour_links["weight"].values.reshape(-1, 1), )
        log.debug("neighbour_links = %s", neighbour_links)

        # add node
        for i, row in neighbour_links.iterrows():
            G.add_edge(row["source"], row["target"], weight=row["weight"])
        return G

    def _get_G(self, links, nodes):
        # draw relationship
        G = nx.Graph()
        for i, row in nodes.iterrows():
            if row["name"] in self.analysis_mutation_group.non_duplicated_aas_positions:
                G.add_node(row["name"], size=row["size"], is_mutation=True, pst=int(row["name"][:-1]))
            else:
                G.add_node(row["name"], size=row["size"], is_mutation=False, pst=int(row["name"][:-1]))

        for i, row in links.iterrows():
            G.add_edge(row["source"], row["target"], weight=row["weight"])
        return G

    def display(self):
        log.debug("self.type1 = %s", self.type1)
        log.debug("self.type2 = %s", self.type2)
        log.info("nodes = %s", self.nodes)
        log.info("links = %s", self.links)
        log.debug("dc = %s", self.degree_c)
        log.debug("cc = %s", self.closeness_c)
        log.debug("bc = %s", self.betweenness_c)

    @staticmethod
    def count_distribution(d: dict):
        values = list(d.values())
        counts = pd.value_counts(values)
        counts = counts.sort_index()
        counts = counts.to_dict()
        return counts.keys(), counts.values()

    def _get_centralities(self, is_weight=True):
        file_names = ["degree", "betweenness", "closeness", "e_betweenness"]
        if is_weight:
            weight = "weight"
            outpath = [os.path.join(self.data_dir, "cache", "{0}_{1}_centrality.json".format(self.file_name,i)) for i in
                       file_names]
        else:
            weight = None
            outpath = [os.path.join(self.data_dir, "cache", "{0}_{1}_centrality_no_weight.json".format(self.file_name,i)) for i in
                       file_names]
        for path in outpath:
            if not os.path.exists(os.path.dirname(path)):
                os.mkdir(os.path.dirname(path))
        self.centrality_cache_dir = outpath

        if not os.path.exists(outpath[0]):
            log.debug("no cache")
            log.info("degree")
            dc = degree_centrality(self.G)
            with open(outpath[0], "w") as f:
                f.write(json.dumps(dc))
        else:
            log.debug("with cache")
            log.info("degree")
            with open(outpath[0], ) as f:
                dc = json.load(f)

        if not os.path.exists(outpath[1]):
            log.debug("no cache")
            log.info("betweenness")
            bc = betweenness_centrality(self.G, weight=weight)
            with open(outpath[1], "w") as f:
                f.write(json.dumps(bc))
        else:
            log.debug("with cache")
            # betweenness
            log.info("betweenness")
            with open(outpath[1], ) as f:
                bc = json.load(f)

        if not os.path.exists(outpath[2]):
            log.debug("no cache")
            # closeness
            log.info("closeness")
            cc = closeness_centrality(self.G, distance=weight)
            with open(outpath[2], "w") as f:
                f.write(json.dumps(cc))
        else:
            log.debug("with cache")
            # closeness
            log.info("closeness")
            with open(outpath[2], ) as f:
                cc = json.load(f)

        if not os.path.exists(outpath[3]):
            log.debug("no cache")
            #  betweenness
            log.info("betweenness")
            e_bc = edge_betweenness_centrality(self.G, weight=weight, normalized=False)
            with open(outpath[3], "wb") as f:
                pickle.dump(e_bc, f)
        else:
            log.debug("with cache")
            # betweenness
            log.info("betweenness")
            with open(outpath[3], "rb") as f:
                e_bc = pickle.load(f)

        result = [dc, bc, cc, e_bc]
        # result = [scaler(i) for i in result]
        result = tuple(result)
        return result

    def analysis_aa(self, aa: str):
        aa = self._aa2position(aa)
        res = {
            "degree_centrality": self.degree_c[aa],
            "betweenness_centrality": self.betweenness_c[aa],
            "closeness_centrality": self.betweenness_c[aa]
        }
        return res

    def analysis_group(self, group: list):
        aas = [self._aa2position(i) for i in group]

    def get_degree(self):
        degrees = self.G
        return degrees

    def get_weighted_degree(self):
        degrees = {}
        for n, nbrs in self.G.adj.items():
            wts = []
            for nbr, eattr in nbrs.items():
                wt = eattr["weight"]
                wts.append(wt)
            if wts:
                avg_wt = np.sum(wts)
            else:
                avg_wt = 0.0
            # degrees.append(avg_wt)
            degrees[n] = avg_wt
        return degrees

    def get_avg_weighted_degree(self):
        degrees = {}
        for n, nbrs in self.G.adj.items():
            wts = []
            for nbr, eattr in nbrs.items():
                wt = eattr["weight"]
                wts.append(wt)
            if wts:
                avg_wt = np.mean(wts)
            else:
                avg_wt = 0.0
            # degrees.append(avg_wt)
            degrees[n] = avg_wt
        return degrees

    def apply_to_aas(self, aas, aas_sample, func):
        aas_score = func(aas)
        aas_sample_score = [func(a_sample) for a_sample in aas_sample]
        return aas_score, aas_sample_score

    def analysisG(self, ):
        # self._plot_origin_distribution()  # procon distribution
        # self._plot_mutations_relationship()  # mutation relationship: node-mutation site, size-occurrence count, edge-conservation
        self._collect_mutation_info()  # collection mutation info and create table
        # self._plot_2D()  # 2D figure

        # substitution
        # self._boxplot_for_all_kinds()
        # self._boxplot_for_all_kinds("BA.4(Omicron)")
        # self._boxplot_for_all_kinds("B.1.617.2(Delta)")

        # variant
        # self._group_plot_with_node()
        # self._group_plot_with_node_for_variant("BA.4(Omicron)")
        # self._group_plot_with_node_for_variant("B.1.617.2(Delta)")

    def output_for_gephi(self):
        # edge
        rows = []
        for source, target in self.G.edges:
            is_neighbour = abs(int(target[:-1]) - int(source[:-1])) == 1
            is_mutation = (source in self.analysis_mutation_group.non_duplicated_aas_positions) or (
                    target in self.analysis_mutation_group.non_duplicated_aas_positions)
            if is_mutation:
                tag = "mutation"
            else:
                tag = "normal"

            w = self.G.edges[source, target]['weight']
            if is_neighbour and w < 0.9:
                w = 0.9

            rows.append([source, target, is_neighbour, is_mutation, tag, w])
        rows = pd.DataFrame(rows, columns=["source", "target", "is_neighbour", "is_mutation", "tag", "weight"], )
        rows.to_csv(os.path.join(self.data_dir, "{}_network_edge_info.csv".format(self.file_name)), index=None)

        # nodes
        nodes = []
        for n in self.G.nodes:
            if n in self.analysis_mutation_group.non_duplicated_aas_positions:
                if self.G.nodes[n]!={}:
                    size=self.G.nodes[n]['size']
                    size = 0.001 if size == 0 else size
                    node = [n, n, True, 5 + 1 / size]
            else:
                node = [n, "", False, 1]
            nodes.append(node)
        nodes = pd.DataFrame(nodes, columns=["Id", "Label", "is_mutation", "new_size"])
        nodes.to_csv(os.path.join(self.data_dir, "{}_network_node_info.csv".format(self.file_name)), index=None)

    def get_functions(self):
        funcs = {
            "CS": self.calculate_conservation,
            "CCS": self.calculate_co_conservation,
            "Kw": self.calculate_avg_weighted_degree,
            "P": self.calculate_page_rank,
            "D": self.calculate_degree_centrality,
            "B": self.calculate_betweenness_centrality,
            "C": self.calculate_closeness_centrality,
            "L": self.calculate_weighted_shortest_path,
        }
        return funcs

    def calculate_degree_centrality(self, grp):
        dc = scaler(self.degree_c)
        return [dc[aa] for aa in grp]

    def calculate_betweenness_centrality(self, grp):
        bc = scaler(self.betweenness_c)
        return [bc[aa] for aa in grp]

    def calculate_closeness_centrality(self, grp):
        cc = scaler(self.closeness_c)
        return [cc[aa] for aa in grp]

    def calculate_degree(self):
        result = {}
        for node in self.G.nodes:
            weighted_degrees = []
            for nbr, datadict in self.G.adj[node].items():
                weighted_degrees.append(datadict.get("weight", 0))

            if weighted_degrees:
                avg_weighted_degree = np.mean(weighted_degrees)
            else:
                avg_weighted_degree = 0
            result[node] = avg_weighted_degree
        return result

    def calculate_avg_weighted_degree(self, grp):
        degree = scaler(self.degree)
        return [degree[aa] for aa in grp]

    def calculate_page_rank(self, grp):
        pr = scaler(self.page_rank)
        return [pr[aa] for aa in grp]

    def calculate_conservation(self, grp):
        return [self.G.nodes[aa]["size"] for aa in grp]

    def get_shortest_path_length(self):
        if not os.path.exists(os.path.join(self.data_dir, "{}_weighted shortest path length.json".format(self.file_name))):
            weighted_shortest_path_length = dict(nx.shortest_path_length(self.G, weight="weight"))
            with open(os.path.join(self.data_dir, "{}_weighted shortest path length.json".format(self.file_name)), "w") as f:
                log.info("save")
                f.write(json.dumps(weighted_shortest_path_length))
        else:
            with open(os.path.join(self.data_dir, "{}_weighted shortest path length.json".format(self.file_name)), ) as f:
                log.info("load")
                weighted_shortest_path_length = json.loads(f.read())
        return weighted_shortest_path_length

    def calculate_weighted_shortest_path(self, grp, is_dict=False):
        res = []
        names = []
        for n1, n2 in combinations(grp, 2):
            res.append(self.shortest_path_length[n1][n2])
            if is_dict:
                names.append(f"{n1}-{n2}")

        if is_dict:
            return dict(zip(names, res))
        else:
            return res

    def calculate_co_conservation(self, grp, is_dict=False):
        res = []
        names = []
        for n1, n2 in combinations(grp, 2):
            if self.G.has_edge(n1, n2):
                res.append(self.G.edges[n1, n2]["weight"])
                if is_dict:
                    names.append(f"{n1}-{n2}")
        if is_dict:
            return dict(zip(names, res))
        else:
            return res

    def calculate_edge_betweenness_centrality(self, grp):
        # if not hasattr(self, "edge_betweenness_centrality"):
        #     try:
        #         # 保存至 json 文件
        #         if not os.path.exists(os.path.join(self.data_dir, "edge betweenness centrality.json")):
        #             log.info("没有 self.edge_betweenness_centrality, 初始化")
        #             self.edge_betweenness_centrality = dict(nx.edge_betweenness_centrality(self.G, weight="weight"))
        #             # 重新建立索引
        #             info_map = defaultdict(dict)
        #             for (n1, n2), value in self.edge_betweenness_centrality.items():
        #                 info_map[n1][n2] = value
        #                 info_map[n2][n1] = value
        #
        #             with open(os.path.join(self.data_dir, "edge betweenness centrality.json"), "w") as f:
        #                 log.info("保存至文件")
        #                 f.write(json.dumps(info_map))
        #
        #         with open(os.path.join(self.data_dir, "edge betweenness centrality.json"), ) as f:
        #             log.info("加载文件")
        #             self.edge_betweenness_centrality = json.loads(f.read())
        #     except:
        #         log.error("保存或加载文件失败")
        #     log.debug("初始化完成")
        #
        # ebc = self.edge_betweenness_centrality
        # # ebc_values = []
        # # for k1, v1 in ebc.items():
        # #     for k2, v2 in v1.items():
        # #         ebc_values.append(v2)
        # # scaler = preprocessing.MinMaxScaler().fit(np.array(ebc_values))
        #
        # for n1, n2 in combinations(grp, 2):
        #     if n1 in ebc and n2 in ebc[n1]:
        #         res.append(ebc[n1][n2])
        #     elif n2 in ebc and n1 in ebc[n2]:
        #         res.append(ebc[n2][n1])
        #     else:
        #         res.append(0)
        if not hasattr(self, "edge_betweenness_centrality_scaler"):
            self.edge_betweenness_centrality_scaler = scaler(self.edge_betweenness_c)

        ebc = self.edge_betweenness_centrality_scaler
        res = []
        for n1, n2 in combinations(grp, 2):
            if (n1, n2) in self.edge_betweenness_c:
                res.append(ebc[(n1, n2)])
            elif (n2, n1) in self.edge_betweenness_c:
                res.append(ebc[(n2, n1)])
            # else:
            #     res.append(0)
        return res

    def calculate_co_conservation_rate(self):
        # find pairwise, create index
        # type2 = self.type2[self.type2["info"] >= self.threshold]
        type2 = self.type2
        idx = defaultdict(list)
        for i, row in type2.iterrows():
            idx[row.site1].append(row.site2)
            idx[row.site2].append(row.site1)

        def _cal(grp):
            rns = len(grp)
            rnp = comb(rns, 2)
            rpc = 0
            for p1, p2 in combinations(grp, 2):
                if p2 in idx[p1]:
                    rpc += 1
            res = [rpc / rnp, ] * len(grp)
            return res

        return _cal

    def generate_ebc(self):
        if not hasattr(self, "edge_betweenness_centrality"):
            try:
                # save
                if not os.path.exists(os.path.join(self.data_dir, "edge betweenness centrality.json")):
                    log.info("no self.edge_betweenness_centrality, init")
                    self.edge_betweenness_centrality = dict(nx.edge_betweenness_centrality(self.G, weight="weight"))
                    info_map = defaultdict(dict)
                    for (n1, n2), value in self.edge_betweenness_centrality.items():
                        info_map[n1][n2] = value
                        info_map[n2][n1] = value

                    with open(os.path.join(self.data_dir, "edge betweenness centrality.json"), "w") as f:
                        log.info("save")
                        f.write(json.dumps(info_map))

                with open(os.path.join(self.data_dir, "edge betweenness centrality.json"), ) as f:
                    log.info("load")
                    self.edge_betweenness_centrality = json.loads(f.read())
            except:
                log.error("error")
            log.debug("success")

    def generate_all_node_top_info(self, top=10):
        seq = self.fasta
        aas = []
        for i, aa in enumerate(seq):
            i = i + 1
            aa = f"{i}{aa}"
            aas.append(aa)
        print(len(aas))
        # cal
        funcs = self.get_functions()
        node_info = {"aas": aas}
        pair_info = {}
        for fname, func in funcs.items():
            if fname in ["CCS", "L"]:
                pair_info[fname] = func(aas, is_dict=True)
            else:
                node_info[fname] = func(aas)

        # node info
        node_info = pd.DataFrame(node_info)
        node_info = node_info.set_index("aas")
        node_info.to_csv(os.path.join(self.data_dir, "{}_node_info.csv".format(self.file_name)))
        # find top
        node_info_top = {}
        for label, content in node_info.items():
            content_sorted = content.sort_values(ascending=False)
            node_info_top[label] = content_sorted.index.to_list()[:top]
        node_info_top = pd.DataFrame(node_info_top)
        node_info_top.to_csv(os.path.join(self.data_dir, "{}_node_info_top.csv".format(self.file_name)))

        # pair info
        pair_info_top = {}
        for lable, content in pair_info.items():
            content = pd.Series(content)
            content_sorted = content.sort_values(ascending=False)
            pair_info_top[lable] = content_sorted.index.to_list()[:top]
        pair_info_top = pd.DataFrame(pair_info_top)
        pair_info_top.to_csv(os.path.join(self.data_dir, "{}_pair_info_top.csv".format(self.file_name)))

        '''
        # variation node
        mutations = self.analysis_mutation_group.non_duplicated_aas_positions
        mutations = pd.Series(mutations)
        mutations.to_csv(os.path.join(self.data_dir, "mutation positions.csv"), header=None)
        '''

if __name__ == '__main__':
    start_time = time.time()
    # conservation network
    # important variant
    mutation_groups = AnalysisMutationGroup('AY.33.txt')
    mutation_groups.display_seq_and_aa()
    mutation_groups.count_aa()
    pcn = ProConNetwork('AY.33.txt',mutation_groups,threshold=100)
    #pcn.analysisG()  # draw

    #pcn.generate_all_node_top_info()

    #pcn._collect_mutation_info()  # save info

    # generate data for other tools
    #pcn.generate_ebc()
    #pcn.output_for_gephi()
    #pcn.output_for_DynaMut2()
    #pcn.output_for_topnettree()
    end_time = time.time()
    log.info(f"run time: {end_time - start_time}")
    # print(mutation_groups.non_duplicated_aas)

    # thresholds = [50, 100, 150, 200, 250, 300]
    # for t in thresholds:
    #     if t != 100:
    #         continue
    #     pcn = ProConNetwork(mutation_groups, threshold=t)
    #
    #     pcn.analysisG()
    #     # pcn.random_sample_analysis(aas, groups.get_aa_groups())
    #
    #     end_time = time.time()
    #     log.info(f"run time: {end_time - start_time}")
