'''
    Author: Minhong Zhu

'''

import subprocess
import argparse
import os
import pandas as pd
import joblib

from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import GradientBoostingClassifier as GBDT
from sklearn.ensemble import RandomForestClassifier as RF

from funcs import compute_dy, compute_mod, compute_net, load_inputs_fasta
from model.predict import predict


def main():
    # argparse模块：命令行选项、参数和子命令解析器。
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta', help="input can be folder of fasta files or a single fasta file")
    parser.add_argument('--data_dir', default='./data', help='data directory, default:./data')
    parser.add_argument('--out_dir', default='./data', help='output directory of intermediate results')
    parser.add_argument('--method', default='ensemble', help='choose model to predict')
    parser.add_argument('--net', action='store_true', help='False if the network has already been calculated')
    parser.add_argument('--dy', action='store_true', help='False if the dynamic features have already been calculated')
    parser.add_argument('--mod', action='store_true', help='False if the modularity features have already been calculated')
    args = parser.parse_args()
    args.pccn_dir = r'{0}/pccn'.format(args.data_dir)

    # load input fasta file(s)
    load_inputs_fasta(args)

    # compute conservation network
    compute_net(args)

    # compute dynamic features & modularity features 
    dy_feats = compute_dy(args)
    mod_feats = compute_mod(args)
    all_feats = pd.merge(mod_feats, dy_feats, on='Var')
    all_feats.to_csv(r'{0}/all_feats.csv'.format(args.out_dir), index=False)

    predict(feat=all_feats, method=args.method, out_dir=args.out_dir)

    # build

main()

"""
    run: python main.py --input_fasta ./input --conpute_net to skip network construction

    make sure that the output of network construction is put in ./data/pccn/procon

"""



    