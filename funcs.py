import pandas as pd
import os
import pdb

from dyna import dyna
from module import module
from pccn import pccn
from energy import energy
from utils.files import move_dir_to_dir, move_file_to_dir, validate_input

def compute_net(args):
    if not args.net:
        net = pccn(args)
        net.calculate()
    
def compute_mod(args):
    # compute modularity node features
    if not args.mod:
        mod = module(args)
        mod.aggregate()
        mod_feats = mod.feats
    else:
        mod_file = r'{0}/module_feature.csv'.format(args.out_dir)
        if not os.path.exists(mod_file):
            raise FileNotFoundError(f"File not found: {mod_file}")
        mod_feats = pd.read_csv(mod_file)
    return mod_feats

def compute_dy(args):
    # compute dynamic node features
    if not args.dy:
        dy = dyna(args)
        dy.calculate()
        dy.aggregate(mutation=True)
        dy_feats = dy.dy_feats
    else:
        # make sure that features are stored in out_dir
        dy_file = r'{0}/dynamic_feature.csv'.format(args.out_dir)  
        if not os.path.exists(dy_file):
            raise FileNotFoundError(f"File not found: {dy_file}")
        dy_feats = pd.read_csv(dy_file)
    return dy_feats


def load_inputs_fasta(args):
    # validate inputs
    inputs = args.input_fasta
    input_type = validate_input(inputs)
    # move inputs to data_fasta directory if needed
    data_fasta_dir = r'{0}/data_fasta'.format(args.pccn_dir)
    if input_type == 'Folder' and inputs != data_fasta_dir:
        move_dir_to_dir(inputs, data_fasta_dir)
    elif input_type == 'Fasta file':
        pre_dir = os.path.dirname(inputs)
        if pre_dir != data_fasta_dir:
            move_file_to_dir(inputs, data_fasta_dir)