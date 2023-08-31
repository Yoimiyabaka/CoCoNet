#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@ author: Ziyun Zhou
@ Email: zhouziyun1900@hotmail.com
@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.
"""

from enm.Enm import * 
from HitCommute import *
import os
import argparse

# Get all network files in procon folder
parser = argparse.ArgumentParser()
parser.add_argument('--dir', help='pccn working directory, default:./data/procon')
args = parser.parse_args()

dirpath = args.dir
network_file_list = list()

for root,dirs,files in os.walk(dirpath):
    for file in files:
        network_file_list.append(os.path.join(root,file)) 
network_file_path = [s for s in network_file_list if '.fasta_network_edge_info.csv' in s]

for network in network_file_path:
    
    # Kirchhoff matrix building
    enm = Enm('network')
    enm.read_network(network, sep = '\t')
    enm.gnm_analysis(normalized = False)
    K = enm.gnm.getKirchhoff()

    # Pass the Kirchhoff matrix to hit/commute time 
    hc = IT_HitCommute(K)

    # Report hitting time
    H = hc.buildHitTimes(K)
    hitting_time = H.mean(axis = 1)

    # Report commute time
    C = hc.buildCommuteTimes()
    commute_time = C.mean(axis = 1)

    # PRS Calculation
    enm.get_sensor_effector(use_threshold = True)
    enm.cluster_matrix(enm.prs_mat)

    # Merge results
    prs = pd.DataFrame(enm.df)
    prs = prs.reset_index(drop = True)
    hitting_time = pd.DataFrame(hitting_time)
    hitting_time.columns = ['hitting_time']
    commute_time = pd.DataFrame(commute_time)
    commute_time.columns = ['commute_time']
    data1 = prs.join(hitting_time, on = None, how = 'left', lsuffix = '', rsuffix = '', sort = False)
    data2 = data1.join(commute_time, on = None, how = 'left', lsuffix = '', rsuffix = '', sort = False)
    data2.rename(columns={'orf_name':'aas'},inplace = True)
    data2.drop(['deg', 'trans', 'eigenvec_centr', 'closeness_centr', 'smallest_eigenvec'], axis = 1, inplace = True)

    # Output
    out_path = network.split(os.sep)[-1]
    seq_name = out_path.split("_", 1)[0]
    cur_file_parent_path = os.path.abspath(network + "\\..")
    dyna_out_path = os.path.join(cur_file_parent_path, f"{seq_name}_node_dynamics_info.csv")
    data2.to_csv(dyna_out_path, index = False)
