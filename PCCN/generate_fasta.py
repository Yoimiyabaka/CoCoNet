# -*- coding = utf-8 -*- 
# @Author : ltl
# @Time : 2023/5/17 13:35
# @File : test.py
# @Software : PyCharm
import pandas as pd

f=pd.read_csv(r'./data/procon/Lineages_Fasta_20230518.csv')
print(type(f))
f.set_index('number', inplace=True)
f.index = pd.to_numeric(f.index)
for row in f.itertuples(index=False):
    file = open(r".\data\data_fasta\{}.fasta".format(row.Lineages), "w")
    file.write(">{}\n".format(row.Lineages))
    file.write(row.Sequence)
    file.close()


'''
folder_path='/PCCN/data/data_fasta'
file_names=os.listdir(folder_path)
for file_name in file_names[:]:
    with open('/PCCN/data/data_fasta/{}'.format(file_name),'r') as f:
        fasta = f.readlines()
        first=fasta[0]
        if first[0]!='>':
            seq=''
            for line in fasta[:]:
                seq=seq+line.strip()
            with open('/PCCN/data/data_fasta/{}'.format(file_name),'w') as f:
                f.write('>QRZ20975.1 spike E484K [Expression vector SARS-CoV-2-S-E484K]\n')
                f.write(seq)
    if file_name[-5:]!='fasta':
        os.rename('/PCCN/data/data_fasta/{}'.format(file_name),'/PCCN/data/data_fasta/{}'.format(file_name[:-3]+'fasta'))
'''