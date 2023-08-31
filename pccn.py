import subprocess


class pccn:
    def __init__(self, args):
        self.pccn_dir = r'.{0}'.format(args.pccn_dir)

    def calculate(self):
        # call pccn in subprocess
        woring_directory = './PCCN'
        subprocess.run(['python', 'main_pccn.py', '--data_dir', self.pccn_dir], cwd=woring_directory)

        