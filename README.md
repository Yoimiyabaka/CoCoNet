# CoCoNet
## 1.Introduction
The specific mutations that occur during the evolution of highly variable viruses will enhance their immune escape ability, which brings great challenges to the development of vaccines and drugs against them. Coronavirus Disease 2019 (Covid-19) virus is a kind of highly variable virus, which has evolved from the original strain to multi-lineage variant strains. We developed a machine learning model-based prediction software (CoCoNet) for the infectiousness of Covid-19 virus mutation lineages. By building a spike protein (S) of the conservative network topological index and calculate mutation module, CoCoNet predict Covid - 19 virus mutation spectrum the infection rate of accuracy is 0.72. Our approach investigates the lineage mutation of COVID-19 from a systems biology perspective and suggests the need for further consideration based on the binding affinity of S protein to its receptor.
## 2.Usage
### Arguments description
- --input_fasta：the fasta format files or folders inputed
- --data_dir：Intermediate data storage path，default: "./data"
- --out_dir：Result storage path，default: "./data"
- --net：Whether to perform coconservative network construction，(default: True)
- --dy：Whether to calculate module topology characteristics or load from the existing results，(default: True)
- --mod：Whether to compute network dynamics features or load from existing results，(default: True)
- --method: usage model，(ensemble or single)(default: ensemble)

### Working flow

```shell
    # following are command lines
    # predict single variant

    python main.py --input_fasta var.fasta


    # predict a group of variants, stored in “input” directory

    python main.py --input_fasta ./input


    # skip network computing

    python main.py --input_fasta var.fasta --net


    # skip module feature computing

    python main.py --input_fasta var.fasta --mod 


    # skip dynamic feature computing

    python main.py --input_fasta var.fasta --dy


    # skip two of them

    python main.py --input_fasta var.fasta --net --dy
```

### Tips
When can you skip some of these steps?
- net: If the current mutation network data (xxx.fasta_threshold_100) has been saved in the './data/pccn/' directory, it can be skipped.
- dy: If the 'out_dir' directory already has dynamic feature.csv, you can skip it.
- mod: If the 'out_dir' directory saved module feature.csv, can be skipped.
## 3.Example


