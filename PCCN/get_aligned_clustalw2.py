# -*- coding = utf-8 -*- 
# @Author : ltl
# @Time : 2023/5/21 9:16
# @File : get_aligned_clustalw2.py
# @Software : PyCharm
import subprocess
#from Bio.Align.Applications import MuscleCommandline
#from Bio import AlignIO
#from Bio.Align.Applications import ClustalwCommandline
def get_aligned_clustalw2(data_dir_matching,file_name):
    infile_name=r"{0}\{1}_protein-matching-IPR042578.filter.fasta".format(data_dir_matching,file_name)
    outfile_name=r"{0}\{1}_protein-matching-IPR042578.filter.fasta.fasta".format(data_dir_matching,file_name)
    #cline = ClustalwCommandline(r"D:\PCCN\clustalw\clustalw2",infile=infile_name,outfile=outfile_name,outfmt='fasta')
    # Run ClustalW command using subprocess
    clustalw_cmd = r".\clustalw\clustalw2"
    output_format = "fasta"
    command = [clustalw_cmd, "-INFILE=" + infile_name, "-OUTFILE=" + outfile_name, "-OUTPUT=" + output_format]
    subprocess.call(command)
