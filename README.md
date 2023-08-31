# CoCoNet
### 参数说明
- --input_fasta：输入的fasta格式文件或文件夹
- --data_dir：中间数据存放路径，default: "./data"
- --out_dir：结果存放路径，default: "./data"
- --net：是否执行共保守网络构建(default: True)
- --dy：是否计算模块拓扑特征或从已有结果加载(default: True)
- --mod：是否计算网络动力学特征或从已有结果加载(default: True)
- --method: 使用模型(ensemble or single)(default: ensemble)

### 运行软件

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

什么时候可以跳过其中某些步骤?
- net: 若./data/pccn/目录下已保存当前突变网络数据(xxx.fasta_threshold_100)，可跳过
- dy: 若out_dir中已保存dynamic_feature.csv，可跳过
- mod: 若out_dir中已保存module_feature.csv，可跳过
