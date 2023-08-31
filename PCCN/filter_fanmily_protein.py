# -*- coding:utf-8 -*-
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import json
import logging
import logconfig


logconfig.setup_logging()
log = logging.getLogger("cov2")


def load_fasta(fasta_file):
    seqs = SeqIO.parse(fasta_file, "fasta")
    seqs = [[i.id.split("|")[0], i.id, i.seq] for i in seqs]
    columns = ["id", "desc", "seq", ]
    seqs = pd.DataFrame(seqs, columns=columns)
    return seqs


def load_json(json_file):
    with open(json_file) as f:
        info = json.load(f)
    tinfo = []
    for i in info:
        metadata = i["metadata"]
        metadata.update(metadata.pop("source_organism"))
        tinfo.append(metadata)
    info = pd.DataFrame(tinfo)
    # info = info.rename({"accession": "id"})
    return info


def filter_family_protein(folder_path,query_fasta_file,data_dir):
    fasta_file = r"./data/protein-matching-IPR042578.fasta"
    json_file = r"./data/protein-matching-IPR042578.json"
    fasta = load_fasta(fasta_file)
    info = load_json(json_file)
    log.debug("fasta columns: %s", fasta.columns)
    log.debug("fasta shape %s ", fasta.shape)
    log.debug("info columns: %s", info.columns)
    log.debug("info shape %s ", info.shape)
    # filter
    log.debug("info no duplicated: %s", info.shape[0] - info.fullName.duplicated().sum())
    data = pd.merge(fasta, info, left_on="id", right_on="accession")
    log.debug("merge:\n%s", data)
    # data = data[data.fullName.duplicated()]
    data = data.drop_duplicates("fullName")
    log.debug("drop duplicates:\n%s", data)

    # output csv and fasta
    data.to_csv(r"{0}/{1}_protein-matching-IPR042578.filter.csv".format(data_dir,query_fasta_file))
    res_fasta = list(SeqIO.parse('{0}/{1}'.format(folder_path,query_fasta_file), "fasta"))
    import pdb
    # pdb.set_trace()
    for i, row in data.iterrows():
        rec = SeqRecord(
            row.seq,
            # Seq(row.seq),
            id=row.id,
            # description="|".join([row.fullName, row.source_database])
            description="| " + " | ".join([row.source_database, row.fullName, ])
        )
        res_fasta.append(rec)
    SeqIO.write(res_fasta, r"{0}/{1}_protein-matching-IPR042578.filter.fasta".format(data_dir,query_fasta_file), "fasta")
