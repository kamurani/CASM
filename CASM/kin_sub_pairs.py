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

from CASM.kinase.kinase import AlphaFillLoader
from CASM.load_tree import get_distance_matrix, get_kinase_name_map

from typing import Dict, List, Union
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

"""
Creates dictionary for kinase ATP binding site coordinate lookup
"""
def get_atp_site_dict(
    filename: str,
):
    with open(filename) as f:

        tsv_file = csv.reader(f, delimiter="\t")

        d = {}
        for line in tsv_file:
            
            acc = line[0]
            try:
                rmsd, atom, ligand, asym_id, x, y, z = line[1:8]
                coords = (x, y, z)
            except: 
                # UNKNOWN
                rmsd = -1
                coords = None 

            data = dict(
                rmsd=float(rmsd),
                coords=coords,
            )
            d[acc] = data
               
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
@c.argument(
    'ATP_SITE_DATASET', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
@c.argument(
    'ALPHAFILL_DIR', nargs=1,
    type=c.Path(exists=True, file_okay=False, dir_okay=True),
)
@c.argument(
    'KINASE_TREE', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
@c.argument(
    'KINASE_TABLE', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
def main(
    kin_sub_dataset,
    plddt_dataset,
    atp_site_dataset,
    alphafill_dir,
    kinase_tree,
    kinase_table,
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
    kinase_atp_sites: dict = get_atp_site_dict(filename=atp_site_dataset)


    #df['SUB_MOD_PLDDT'] = df.apply(lambda row: d.get(row["SUB_ACC_ID"], {}).get(row["SUB_MOD_RSD"], {}).get("plddt", "UNKNOWN"), axis=1)
    df['SUB_MOD_PLDDT'] = df.apply(func, axis=1)
    df['SUB_MOD_PLDDT'] = df['SUB_MOD_PLDDT'].astype(float)

    
    
    
    
    # join rmsd of atp coord prediction for kinase
    loader = AlphaFillLoader(
        verbose=False,
        cif_dir=alphafill_dir,
    )
    
    # Realtime way using ``loader``
    #df['KIN_ATP_LOC_RMSD'] = df.apply(lambda row: loader.get_coordinates(row['KIN_ACC_ID']).get('rmsd', -1) if loader.get_coordinates(row['KIN_ACC_ID']) else -1, axis=1)
    
    # Using file instead...
    df['KIN_ATP_LOC_RMSD'] = df.apply(lambda row: kinase_atp_sites.get(row["KIN_ACC_ID"], dict(rmsd=-1)).get('rmsd'), axis=1)
    df['KIN_ATP_LOC_RMSD'] = df['KIN_ATP_LOC_RMSD'].astype(float)


    # Join same sites together so we can see if there's >1 kinase for same site
    #df.sort_values(['SUB_ACC_ID','SUB_MOD_RSD'])

    df = df.groupby(['SUB_ACC_ID', 'SUB_MOD_RSD']).agg({'KIN_ACC_ID': lambda x: x.tolist()})
    """
        df.merge(
            grouped,
            left_on=
        )
    """
    

    # Link kinase family name
    
    d = get_kinase_name_map(
        kinase_table,
        "UniprotID",
    )

    xname = get_kinase_name_map(
        kinase_table,
        "xName",
    )

    df['KIN_GROUP_NAME'] = df.apply(lambda row: [d[k]["Group"] if k in d else "UNKNOWN" for k in list(row['KIN_ACC_ID'])], axis=1)
    df['NUM_KINASES'] = df.apply(lambda row: len(set(row['KIN_GROUP_NAME'])), axis=1)
    df['KIN_XNAME'] = df.apply(lambda row: [d[k]["xName"] if k in d else "UNKNOWN" for k in list(row['KIN_ACC_ID'])], axis=1)


    MATRIX = get_distance_matrix(filepath=kinase_tree, kinase_table=kinase_table)

    """
    Get furthest ``n`` kinases away from current kinase.  

    # TODO: if multiple given, maybe get furthest ``n`` from ALL of them somehow? i.e. average distance to all? idk
    """
    def _get_furthest_kinases_single(
        acc_id: str, 
        n: int,
    ):
        if acc_id not in MATRIX:
            return []

        kin_pairs = MATRIX[acc_id]

        # Sort in descending order of pairwise distance
        lists = sorted(kin_pairs.items(), key=lambda item: item[1], reverse=True) # list of tuples, sorted by dict's value
        x, y = zip(*lists)

        return x[0:n]

    def get_furthest_kinases(
        acc_id: Union[str, List[str]],
        n: int, 

    ):
        output = [] 
        if type(acc_id) == str:
            acc_id = [str]

        
        for i in acc_id: 

            output = list(set(output + _get_furthest_kinases_single(i, n=n)))


    



        


        



    N = 20

    # For now, just use the first kinase in list.   
    # TODO: see what happens if we do all of them, then get set of returned list (unique kinases)
    df['FURTHEST_KINASES'] = df.apply(lambda row: get_furthest_kinases(list(row['KIN_ACC_ID']), n=N), axis=1)
    
    df = df.sort_values(['NUM_KINASES'], ascending=False)

    df.to_csv("df_dump.csv", sep="\t", index=True)



    return


    # see all kinases for same site; join kinase family; join kinase super / sub families 
    

    # manually check if multiple unique kinase families (at same site) are similar or different

    # sanity check this, then make a negative example generator 

    # FILTERING

    mod_rsd_plddt_threshold: float  = 60
    kin_atp_rmsd_threshold: float   = 2.0
    dff = df

    dff = dff.loc[dff["SUB_MOD_PLDDT"] >= mod_rsd_plddt_threshold]      # phosphosite in AF2 model has pLDDT score >= threshold 
    dff = dff.loc[dff["KIN_ATP_LOC_RMSD"] >= -1]                        # kinase ATP location is known (``-1`` is error value)
    dff = dff.loc[dff["KIN_ATP_LOC_RMSD"] <= kin_atp_rmsd_threshold ]   # kinase ATP location prediction is above threshold RMSD (AlphaFill model local backbone)


    # Selected fields for printing and sanity checking
    print(df[["KIN_ACC_ID", "KINASE", "SUB_ACC_ID", "SUB_MOD_RSD", "SUB_MOD_PLDDT", "KIN_ATP_LOC_RMSD"]])
    #print(dff)


    # THEN:

    # graphein data loader for kinases, substrates 

    # start training GAT/GCN classifier with negative / positive example pairs 

    


if __name__ == "__main__":
    main()