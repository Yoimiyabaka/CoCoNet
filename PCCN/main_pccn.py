# -*- coding = utf-8 -*- 
# @Author : ltl
# @Time : 2023/5/15 11:20
# @File : main.py
# @Software : PyCharm
from pathlib import Path
import sys
import cal_procon as cl
from output_procon_analysis import AnalysisMutationGroup,ProConNetwork
import filter_fanmily_protein as ffp
import analysis_procon as ap
import os
from ProCon.myProCon import ProbabilityCalculator, MutualInformation, TripletFinder
import importlib
import get_aligned_clustalw2 as gac
import logging
import argparse

file_root=Path(__file__).resolve().parents[0]
sys.path.append(str(file_root))

log=logging.getLogger('cov2')

parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', help='pccn working directory, default:./data')
args = parser.parse_args()

folder_path=r'{0}/data_fasta'.format(args.data_dir)     #store fasta
data_dir_matching=r'{0}/data_filter'.format(args.data_dir)      #filtered and anglined file
type_file=r"{0}/procon".format(args.data_dir)     #cal_procon output
data_dir_pcn=r'{0}/procon'.format(args.data_dir)      #pcn output
analysis_res_path=r'{0}/data_analysis'.format(args.data_dir)      #analysis output
file_names=sorted(os.listdir(folder_path))

for file_name in file_names[:]:#'file_names[a:b]'set the range of file you want to generate

    log.info('start to generate the net of {}'.format(file_name))

    # Find protein family and filter it
    ffp.filter_family_protein(folder_path,file_name,data_dir_matching)
    
    # Align the protein family file
    gac.get_aligned_clustalw2(data_dir_matching,file_name)

    align_fasta = r"{0}/{1}_protein-matching-IPR042578.filter.fasta.fasta".format(data_dir_matching,file_name)

    # Calculator Probability and output as type.csv
    pc = ProbabilityCalculator.ProbabilityCalculator(align_fasta)
    cl.cal_procon(pc,folder_path,file_name,data_dir_matching,type_file)

    # Analysis procon
    ap.analysis_procon(folder_path,file_name,type_file,analysis_res_path)

    # Output the net
    mutation_groups = AnalysisMutationGroup(folder_path,analysis_res_path,file_name)
    pcn = ProConNetwork(file_name,type_file,data_dir_pcn,mutation_groups, threshold=100)
    pcn.generate_all_node_top_info()
    pcn.output_for_gephi()

    # Reload the package
    importlib.reload(ProbabilityCalculator)

    log.info('the net of {} has already generated'.format(file_name))




