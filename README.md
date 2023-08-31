# CoCoNet
## 1.Introduction
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

When can you skip some of these steps?
- net: If the current mutation network data (xxx.fasta_threshold_100) has been saved in the './data/pccn/' directory, it can be skipped.
- dy: If the 'out_dir' directory already has dynamic feature.csv, you can skip it.
- mod: If the 'out_dir' directory saved module feature.csv, can be skipped.
