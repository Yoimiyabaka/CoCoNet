# -*- coding:utf-8 -*-
import os
import json
import logging
import logconfig
import pandas as pd
from Bio import SeqIO
import itertools


logconfig.setup_logging()
log = logging.getLogger("cov2")


def load_procon_res(res_path):
    d1, d2 = [pd.read_csv(i) for i in res_path]
    # type1
    r1 = {}
    for i, row in d1.iterrows():
        k = int(row.position[:-1])
        r1[k] = row.to_dict()
    yield r1

    # type2
    r2 = {}
    for i, row in d2.iterrows():
        k = tuple(sorted([int(row.site1[:-1]), int(row.site2[:-1])]))
        r2[k] = row.to_dict()
    yield r2
    '''
    # type3
    r3 = {}
    for i, row in d3.iterrows():
        k = tuple(sorted([int(row.site1[:-1]), int(row.site2[:-1]), int(row.site3[:-1])]))
        r3[k] = row.to_dict()

    yield r3
    '''

def check_aa(seq, aa):
    position = aa[1:-1]
    if position.isdigit():
        position = int(position)
    else:
        print("error", aa)
        return False
    origin = aa[0].upper()
    if seq[position - 1].upper() == origin:
        return True
    else:
        print("error", aa)
        return False


def analysis_procon(folder_path,fasta_file,type_file,analysis_res_path):
    # load
    variation_file = "./data/总结 - 20220709 - filter.xlsx"
    log.info("load")
    variation = pd.read_excel(variation_file, sheet_name=1, index_col=0)
    variation["Year and month first detected"] = pd.to_datetime(variation["Year and month first detected"])
    variation["Year and month first detected"] = variation["Year and month first detected"].dt.strftime('%B %d, %Y')
    log.debug("columns: %s", variation.columns)
    log.debug(variation)

    # remove no clear: Impact on transmissibility	Impact on immunity	Impact on severity
    evidence_columns = ["Impact on transmissibility", "Impact on immunity", "Impact on severity"]
    # print(np.unique(variation[evidence_columns].values.reshape((1, -1))))
    is_vital = variation[evidence_columns].applymap(lambda x: "increase" in x.lower()).any(axis=1)
    variation = variation[is_vital]
    log.debug(f"count {variation.shape}")

    # check
    fasta = next(SeqIO.parse('{0}/{1}'.format(folder_path,fasta_file), "fasta")).seq
    t_fasta = {}
    for i, row in variation.iterrows():
        aas = row.loc["Spike mutations of interest"].split(" ")
        aas = map(lambda x: x.strip(), aas)
        # filter
        aas = filter(lambda x: check_aa(fasta, x), aas)
        aas = list(aas)
        t_fasta[i] = row.to_dict()
        t_fasta[i]["aas"] = aas
    fasta = t_fasta

    # load procon result
    result_dir = r"{0}/{1}_type".format(type_file,fasta_file)
    result_files = [os.path.join(result_dir, "{0}_type{1}_parse.csv".format(fasta_file,i)) for i in range(1, 3)]
    log.debug("parse to dict")
    pr = load_procon_res(result_files)

    log.info("parse type1")
    pr1 = next(pr)
    log.debug("type1: %s", pr1)
    for i, row in fasta.items():
        aas = row["aas"]
        fasta[i]["type1"] = []
        fasta[i]["t1"] = {}
        for aa in aas:
            pst = int(aa[1:-1])
            res = pr1.get(pst, None)
            if res is None:
                res = {"rank": "", "position": aa[1:-1] + aa[0], "information": "", "rate": ""}
            res["aa"] = aa
            fasta[i].get("type1", ).append(res)
            # fasta[i]["t1"][aa] = res
    log.debug("type1 result: %s", fasta)

    log.info("parse type2")
    pr2 = next(pr)
    log.debug("type2: %s", pr2)
    for i, row in fasta.items():
        aas = row["aas"]
        fasta[i]["type2"] = []
        # fasta[i]["t2"] = {}
        for aa2 in itertools.combinations(aas, 2):
            pst2 = tuple(sorted([int(a[1:-1]) for a in aa2]))
            res = pr2.get(pst2, )
            if res is not None:
                fasta[i].get("type2").append(res)
            # fasta[i]["t2"][tuple(aa2)] = res
    log.debug("type2 result: %s", fasta)

    '''
    log.info("parse type3")
    pr3 = next(pr)
    for i, row in fasta.items():
        aas = row["aas"]
        fasta[i]["type3"] = []
        # fasta[i]["t3"] = {}
        for aa3 in itertools.combinations(aas, 3):
            pst3 = tuple(sorted([int(a[1:-1]) for a in aa3]))
            res = pr3.get(pst3, None)
            if res is not None:
                fasta[i].get("type3").append(res)
            # fasta[i]["t3"][pst3] = res
    log.debug("type3 result: %s", fasta)
    '''

    log.info("save")
    analysis_res_path = r"{0}/{1}_analysis.json".format(analysis_res_path,fasta_file)
    with open(analysis_res_path, "w") as f:
        f.write(json.dumps(fasta, indent="\t"))
