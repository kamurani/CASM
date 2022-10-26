"""Generate training examples"""



"""

KINASE / SUBSTRATE INPUT


Positive examples:

- pairs from known kinase phosphorylation sites 

Negative examples:

- known psite with known kinase; 
    - find kinase family that is (most) different to the kinase family of the known kinases.
        - use distance metric of similarity i.e. family tree of kinases, MSA distance score
        - note that this similarity metric (global sequence) will be influenced by other domains
        not necessarily implicated in substrate binding for phosphorylation 
            - in future, maybe use structural similarity for kinases to get distance score?
    - pick (randomly?) from this family to get pairs that we think are least likely to interact to give 
    negative pairs. 



CLUSTERING

- all known sites, keep similar structures together.

"""



from typing import Dict
import click as c
import pandas as pd

import csv

def load_kin_sub_list(
    path_to_dataset: str,
) -> pd.DataFrame:
    df = pd.read_csv(
        path_to_dataset, 
        sep='\t', 
        header=0, 
        skiprows=3,     # Header starts on LINE 4
        ) 
    
    return df

"""
TODO: create a class that keep this somewhere for easy lookup. function that returns "unknown" if 
no plddt score available for given protein / site. 
"""
def get_plddt_dict(
    filename: str, 

    filter: str = None, # TODO, possibly get only subset or something in this function etc.
) -> Dict[str, Dict]:

    d = {}
    with open(filename) as f:

        tsv_file = csv.reader(f, delimiter="\t")

        for line in tsv_file:
            
            acc, site, plddt, organism = line[0:4]

            info = dict(
                plddt=float(plddt),
                organism=organism,
            )
            if acc not in d:
                d[acc] = {}

            d[acc][site] = info
               
    return d



@c.command()
@c.argument(
    'KIN_SUB_DATASET', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
@c.argument(
    'PLDDT_DATASET', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
def main(
    kin_sub_dataset,
    plddt_dataset,
):
    print(kin_sub_dataset)
    df: pd.DataFrame = load_kin_sub_list(kin_sub_dataset)

    print(len(df))


    # Filter organisms 

    # human kin , human sub (FOR NOW)
    
    organism = "human"

    df = df[df["KIN_ORGANISM"] == organism]
    df = df[df["SUB_ORGANISM"] == organism]

    

    # this is horrific, shouldn't have used a lambda here
    def get_filter_func(d: dict):

        return lambda row: d[row["SUB_ACC_ID"]][row["SUB_MOD_RSD"]]["plddt"] if (row["SUB_ACC_ID"] in d and row["SUB_MOD_RSD"] in d[row["SUB_ACC_ID"]]) else -1 #"UNKNOWN"
    """  
    def func(row):

        acc = row["SUB_ACC_ID"]
        site = row["SUB_MOD_RSD"]
        
        plddt = d[acc][site]["plddt"]
        return plddt 
    """

    # join plddt score 
    func = get_filter_func(get_plddt_dict(filename=plddt_dataset))

    d = get_plddt_dict(filename=plddt_dataset)
    #df['SUB_MOD_PLDDT'] = df.apply(lambda row: d.get(row["SUB_ACC_ID"], {}).get(row["SUB_MOD_RSD"], {}).get("plddt", "UNKNOWN"), axis=1)
    df['SUB_MOD_PLDDT'] = df.apply(func, axis=1)
    df['SUB_MOD_PLDDT'] = df['SUB_MOD_PLDDT'].astype(float)

    plddt_threshold = 60
    dff = df
    dff = dff.loc[dff["SUB_MOD_PLDDT"] >= plddt_threshold]


    print(df[["KIN_ACC_ID", "KINASE", "SUB_ACC_ID", "SUB_MOD_RSD", "SUB_MOD_PLDDT"]])
    print(dff)
    
    
    
    # join rmsd of atp coord prediction for kinase
    #df['KIN_ATP_LOC_RMSD'] = 
    

    # see all kinases for same site; join kinase family; join kinase super / sub families 

    # manually check if multiple unique kinase families (at same site) are similar or different

    # sanity check this, then make a negative example generator 




    # THEN:

    # graphein data loader for kinases, substrates 

    # start training GAT/GCN classifier with negative / positive example pairs 

    


if __name__ == "__main__":
    main()