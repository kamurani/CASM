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

verbose = True

KINASE_FAMILY_DICT = dict(enumerate(KINASE_FAMILIES))

NUM_KINASE_FAMILIES = len(KINASE_FAMILY_DICT)

KINASE_FAMILY_TO_INDEX = {item:idx for idx, item in enumerate(KINASE_FAMILIES)}

#print(f"length: {len(KINASE_FAMILY_DICT)}")

root_dir = "../../DATA/dbPTM/kinases/"

# Get dict 
uniprot_entry2acc = {}
entry_name_fp = "../../DATA/dbPTM/uniprot_entry_names.tsv"
with open(entry_name_fp) as f:
    next(f)
    for line in f:
        (entry_name, acc_id) = line.split()
        uniprot_entry2acc[entry_name] = acc_id


"""
Convert entry name to uniprot ID
"""
def get_acc_id(entry_name: str):

    if entry_name in uniprot_entry2acc:
        return uniprot_entry2acc[entry_name]
    return None

"""
Return dictionary of all sites
"""
def get_sites():
    # Go through all files 
    count = {} 
    chain = 'A'

    bad_entry_names = []

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
                
                acc_id = get_acc_id(entry_name)
                if acc_id is None: 
                    bad_entry_names.append(entry_name)
                    #print(f"No Uniprot ID mapping from '{entry_name}'")
                    continue 

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
                site = (acc_id, pos) 

                # index using actual node ID 
                #site = (acc_id, mod_rsd)
                if site not in sites:
                    sites[site] = {"pos": [], "neg": [], "mod_rsd": mod_rsd} # "data": d

                if kinase not in sites[site][polarity]:
                    sites[site][polarity].append(kinase) # Should contain unique values

                count[kinase][polarity] += 1

    return sites



def print_entry_names(sites):

    names = [entry_name for (entry_name, pos) in sites.keys()] 
    names = list(set(names)) # only unique values (redundant)

    for n in names:
        print(n)


#sites = get_sites()
# for i in list(sites.items())[0:10]:
#     print(i)


        #print(kinase, polarity, count[kinase][polarity], sep="\t")

def _test_sites():

    #TEST 
    entry_name = "TERT_HUMAN"
    acc_id = get_acc_id(entry_name)
    site = (acc_id, 550)

    mod_rsd = sites[site]['mod_rsd']
    print(mod_rsd)
    print(sites.keys())


#print(f"{entry_name} {mod_rsd}", f"POS: {sites[site]['pos']}", f"NEG: {sites[site]['neg']}")
# for key, item in sites.items():
#     print(item["pos"])
    #print(len(item["pos"]), len(item["neg"]), sep="\t")

#filename = {}
#fasta_fns = []



        # We want a dict that has each site; and has a vector describing the kinase families (with -1 for unknown i.e. it doesn't show up)
        # We can then choose what we do with the -1s 
        #    site, pos, neg, 

