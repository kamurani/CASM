"""Loading data from dbPTM benchmark set https://awi.cuhk.edu.cn/dbPTM/download.php"""

import os
from pathlib import Path


KINASE_FAMILIES = [
    "CDK",
    "MAPK",
    "PKA",
    "PKC",
    "CK2",
    "CAMKL",
    "GSK",
    "AKT",
    "CAMK2",
    "CK1",
    "RSK",
    "GRK",
    "PKG",
    "DYRK",
    "MAPKAPK",
    "DMPK",
    "PKD",
    "PDK1",
    "SGK",
    "RAD53",
    "DAPK",
    "PKN",
    "CAMK1",
    "MLCK",
    "NDR",
]

KINASE_FAMILY_DICT = dict(enumerate(KINASE_FAMILIES))

KINASE_TO_INDEX = {item:idx for idx, item in enumerate(KINASE_FAMILIES)}

#print(f"length: {len(KINASE_FAMILY_DICT)}")

root_dir = "../../DATA/dbPTM/kinases/"

from Bio import SeqIO

# Go through all files 
count = {} 
for kinase in KINASE_FAMILIES:
    count[kinase] = {}
    for polarity in ["pos", "neg"]:

        count[kinase][polarity] = 0
        fn = f"{kinase}_{polarity}.fasta"
        fp = Path(root_dir) / fn 

        # Load fasta file
        
        for record in SeqIO.parse(fp, "fasta"):
            count[kinase][polarity] += 1

        print(kinase, polarity, count[kinase][polarity], sep="\t")
     
        #filename = {}
        #fasta_fns = []
        


        # We want a dict that has each site; and has a vector describing the kinase families (with -1 for unknown i.e. it doesn't show up)
        # We can then choose what we do with the -1s 
        #    site, pos, neg, 

