## Ways to output
if you have the fasta sequences gather in an excel file as : `Lineages_Fasta_20230518.csv`
run :`generate_fasta.py` to generate individual fasta files
else: make sure your fasta files has located in the `/PCCN/data/data_fasta`,and then run the file:
`main.py`
[tips]:about 410 fasta files has already prepared in `data/data_fasta` , so you could just run the `mian.py`and
do nothing.

## Environment requirements
1. biopython==1.79
2. matplotlib
3. networkx
4. numpy
5. pandas
6. requests
7. scipy
8. seaborn
9. sklearn
10. statannotations

## programs explanation
1. Data Preparation.
    * S protein FASTA seq`data/data_fasta`.
    * Variant of concern `data/summary - 20211201.xlsx`.
2. Find [protein family](http://www.ebi.ac.uk/interpro/result/InterProScan/iprscan5-R20210917-073330-0879-28690888-p2m/)
   .
    * Fasta seqs: `data/protein-matching-IPR042578.fasta`.
    * Detail information `data/protein-matching-IPR042578.json` .
3. Seqs of same protein family are too much, so filter seqs with same name and save`filter_fanmily_protein.py`.
    * Fasta seqs: `data/data_filter/{fasta_file_name}_protein-matching-IPR042578.filter.fasta`.
    * Detail information `data/data_filter/{fasta_file_name}_protein-matching-IPR042578.filter.csv`.
4. `get_aligned_clustalw2.py`:`clustalw2` MSA, and export as `data/data_filter{fasta_file_name}_protein-matching-IPR042578.filter.fasta.aligned`.
5. `cal_procon.py` calculates the conservation score, and result is  `data/procon/{fasta_name}_type{}.txt`.
    * Due to the default output format is hard to analysis, it will be converted as csv
      file`data/procon/{file_name}_type/type{}_parse.csv`.
6. `analysis_procon.py` finds the conservation of the variants, and save as `data/procon/{fasta_name}_analysis.json`.
7. `output_procon_analysis` will construct the network and analysis the network.
output:
1)`cache file`
2)`network_node/edge_info.csv`
3)`node_info.csv`
4)`node_info_top.csv`
8.cache file is in the `{fasta_name}_threshold_100`
9.`allosteric_analysis.py`,`combine_with_other_tools.py`has no use here. However you can use it to anaylse 
your output.
