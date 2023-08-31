import pandas as pd
import argparse

import os


class energy:
    def __init__(self, args) -> None:
        """
            this class is for computing change of energy given the mutation

        """
        self.input_dir = r'{0}/procon'.format(args.pccn_dir)    # input_dir: the root_dir of input data
        self.energy_dir = args.data_dir
        self.out_dir = args.out_dir
        self.feat_names = ['Var', 'energy']
        self.feats = pd.DataFrame(columns=self.feat_names)

    def aggregate(self):
        input_folder = sorted(os.listdir(self.input_dir))
        for sub_dir in input_folder:
            if sub_dir.endswith("_threshold_100"):
                var_name = list(sub_dir.split('_'))[0]
                nodes_info_file = r'{0}/{1}/{2}_network_node_info.csv'.format(self.input_dir, sub_dir, var_name)
                energy_file = r'{0}/energy.csv'.format(self.energy_dir)
                self.compute_energy(nodes_info_file, energy_file, var_name)

        out_file = r'{0}/energy_feature.csv'.format(self.out_dir)
        self.feats.to_csv(out_file, index=False)

    def compute_energy(self, nodes_file, energy_file, var_name):
        feat_dict = {'Var': var_name.split('.fasta')[0]}
        energy = pd.read_csv(energy_file)
        node = pd.read_csv(nodes_file)

        df = pd.merge(node, energy[['mutation', 'change']], left_on='Id', right_on='mutation', how='left')
        result = df.iloc[list(df['is_mutation']),:]['change'].sum(skipna=True)
        feat_dict['energy'] = result
        pre_feats = pd.DataFrame(feat_dict, index=[0])
        self.feats = pd.concat([self.feats, pre_feats], axis=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pccn_dir', default='./data/pccn')
    parser.add_argument('--data_dir', default='./data')
    parser.add_argument('--out_dir', default='./data')
    args = parser.parse_args()

    en = energy(args)
    en.aggregate()

