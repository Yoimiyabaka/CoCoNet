import subprocess
import pandas as pd
import os


class dyna:
    def __init__(self, args):
        self.out_dir = args.out_dir
        self.input_dir = r'{0}/procon'.format(args.pccn_dir)     # relative path of procon w.r.t root dir of Dyna 
        self.dy_feat_names = ['Var', 'eff', 'hitting_time', 'sens', 'commute_time']
        self.dy_feats = pd.DataFrame(columns=self.dy_feat_names)    # store dynamic features
    
    def calculate(self):
        working_dir = r'./Dyna'
        file_dir = r'.{0}'.format(self.input_dir)   # change relative dir path for subprocess
        subprocess.run(['python', 'main_dyna.py', '--dir', file_dir], cwd=working_dir)

    def aggregate(self, mutation=False):
        input_folder = sorted(os.listdir(self.input_dir))
        for sub_dir in input_folder:
            if sub_dir.endswith("_threshold_100"):
                var_name = list(sub_dir.split('_'))[0]
                nodes_dynamic_file = r'{0}/{1}/{2}_node_dynamics_info.csv'.format(self.input_dir, sub_dir, var_name)
                nodes_file =  r'{0}/{1}/{2}_network_node_info.csv'.format(self.input_dir, sub_dir, var_name)
                self.aggr_dy_feats(nodes_dynamic_file, nodes_file, var_name, mutation)
        
        out_file = r'{0}/dynamic_feature.csv'.format(self.out_dir)
        self.dy_feats.to_csv(out_file, index=False)
    
    def aggr_dy_feats(self, node_dynamic_file, nodes_file, var_name, mutation):
        '''
            aggregate dynamic features
            mutation: if True, will only calculate the average features of mutant res
        '''
        feat_dict = {'Var': var_name.split('.fasta')[0]}
        nodes_dy_info = pd.read_csv(node_dynamic_file)
        nodes_dy_info.rename(columns={'aas':'Id'},inplace = True)
        nodes_info = pd.read_csv(nodes_file)
        if mutation:
            nodes_dy_info = pd.merge(nodes_dy_info, nodes_info[['Id', 'is_mutation']], on='Id')
            nodes_dy_info = nodes_dy_info.iloc[list(nodes_dy_info['is_mutation']),:]
            nodes_dy_info.drop(['is_mutation'], axis=1, inplace=True)
        nodes_score = nodes_dy_info.iloc[:,1:].mean(axis=0).to_dict()
        feat_dict.update(nodes_score)
        pre_feats = pd.DataFrame(feat_dict, index=[0])
        self.dy_feats = pd.concat([self.dy_feats, pre_feats], axis=0)


