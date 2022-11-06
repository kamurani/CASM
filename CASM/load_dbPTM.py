"""Loading data from dbPTM benchmark set https://awi.cuhk.edu.cn/dbPTM/download.php"""

import os
from pathlib import Path
from Bio import SeqIO

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

aa1to3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}


KINASE_FAMILY_DICT = dict(enumerate(KINASE_FAMILIES))

KINASE_TO_INDEX = {item:idx for idx, item in enumerate(KINASE_FAMILIES)}

#print(f"length: {len(KINASE_FAMILY_DICT)}")

root_dir = "../../DATA/dbPTM/kinases/"



# Go through all files 
count = {} 
chain = 'A'

# Dict of all sites 
sites = {}
for kinase in KINASE_FAMILIES:
    count[kinase] = {}
    for polarity in ["pos", "neg"]:

        count[kinase][polarity] = 0
        fn = f"{kinase}_{polarity}.fasta"
        fp = Path(root_dir) / fn 

        # Load fasta file
        
        
        for record in SeqIO.parse(fp, "fasta"):


            header  = record.id
            seq     = record.seq
            entry_name, pos = header.rsplit("_", 1) 
            
             
            res: str = seq[len(seq) // 2]   # middle character
            res = aa1to3[res]               # 3 letter code
            
            mod_rsd: str = ':'.join([chain, res, pos])
            pos = int(pos)

            d = {
                "pos": pos, 
                "res": res, 
                "mod_rsd": mod_rsd, # Node ID 
            }

            # Index using this tuple 
            site = (entry_name, pos) 
            if site not in sites:
                sites[site] = {"pos": [], "neg": [], "data": d} 

            sites[site][polarity].append(kinase) 

            count[kinase][polarity] += 1




        #print(kinase, polarity, count[kinase][polarity], sep="\t")


# TEST 
entry_name = "TERT_HUMAN"
site = (entry_name, 550)

mod_rsd = sites[site]['data']['mod_rsd']
print(mod_rsd)
#print(sites.keys())
exit(1)

#print(f"{entry_name} {mod_rsd}", f"POS: {sites[site]['pos']}", f"NEG: {sites[site]['neg']}")
for key, item in sites.items():
    print(item["pos"])
    #print(len(item["pos"]), len(item["neg"]), sep="\t")

#filename = {}
#fasta_fns = []



        # We want a dict that has each site; and has a vector describing the kinase families (with -1 for unknown i.e. it doesn't show up)
        # We can then choose what we do with the -1s 
        #    site, pos, neg, 

