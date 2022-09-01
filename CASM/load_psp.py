"""Load data from PhosphositePlus (PSP) dataset."""

"""
Agnostic as to how the structure was determined.  Simply receives a directory of already-existing PDB files. 

"""

from operator import inv
import os
from pathlib import Path
import pdb
import pickle
from typing import Callable, List, Union
from subgraphs import get_motif_subgraph
from utils.residue import aa1to3

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


import click as c
import re

import graphein 

from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1
from utils.residue import aa3to1, aa1to3

from tqdm import tqdm


from definitions import GRAPH_NODE_FEATURES


def mod_rsd2node_id(
    mod_rsd: str, 
    chain_id: str = 'A', # TODO by default, assume that we are using AF2 i.e. always one chain, no multimeric proteins. 
) -> str:

    # {RES}{POS}-{PTM}
    p = re.compile("([a-zA-Z]{1})([0-9]+)-([a-zA-Z])")
    
    match = p.search(mod_rsd)
    res = aa1to3(match.group(1))
    pos = match.group(2)
    ptm = match.group(3)

    try:
        return ':'.join([chain_id, res, pos])
    except:
        raise ValueError(f"{mod_rsd} not in valid MOD_RSD format")

def load_psp_list(
    path_to_psp: str,
) -> pd.DataFrame:
    df = pd.read_csv(
        path_to_psp, 
        sep='\t', 
        header=0, 
        skiprows=3,     # Header starts on LINE 4
        ) 
    
    return df

"""Returns a function that can filter a dataframe of residue IDs"""
"""
Note: assumes that data field input (in returned function) is of form A:RES:000
"""

def get_residue_filter(
    residues: Union[List[str], str],
    invert: bool = False,
) -> Callable:
    """
    :param residues: Either a string containing 1-letter codes. 
    :type residues: Union[List[str], str] 
    :param invert: Return True if the input is NOT in the specified list of residues.  Defaults to ``False``. 
    :type invert: bool


    To allow for all residues (i.e. apply no filtering), an empty list can be supplied with ``invert`` set to ``True``.
    
    """
    if type(residues) == str: # String containing 1-letter codes
        residues = residues.upper()
    else:   # List
        residues = "".join(aa1to3(x.upper()) for x in residues)


    return lambda x: not invert if (aa3to1(x.split(':')[1]) in residues) else invert

    A, B = True, False 
    if invert: A, B = B, A


def filter_psp_list(
    df: pd.DataFrame, 
    residues: Union[List, str],


) -> pd.DataFrame:

    pass


"""
Print list of ACC_IDs for structures that were not found in AF database
"""
def print_dataset_pdb_matches(
    df: pd.DataFrame, 
    pdb_dir, 
    fields: List[str] = ["ACC_ID"],
    delim: str = '\t', 
    invert: bool = True, # if true, shows structures that are missing
):
    df_dict = df.to_dict('records')

    for row in df_dict: # tqdm

        p = Path(pdb_dir) 

        p = p / get_pdb_filename(acc_id=row["ACC_ID"])
        p = p.resolve() # get absolute path

        output = "\t".join([row[f] for f in fields])

        pdb_path = str(p)

        # TODO: more elaborate check, e.g. an actual valid PDB file (not just a file with the correct name)
        
        condition = os.path.isfile(pdb_path) != invert # XOR operator
        if condition:
            print(output)
        

def graph_dump(
    filename: str, 
    out_dir: Path, 
    
):
    pass


def get_graph_filename(
    acc_id: str,
    mod_rsd: str,
    radius: float,
    extension: str = "pickle"  # JSON, .dat, etc. 

) -> str:

    return f"AF-{acc_id}-{mod_rsd}-R{radius}.{extension}" 

"""
Get AlphaFold filename from database download
"""
def get_pdb_filename(
    acc_id: str, 

    file_extension: str = "pdb", # could be .pdb.gz 

    model_version: int = 3, 
    ignore_fragments: bool = True, 

) -> str: 

    if ignore_fragments:
        fragment = 1
        return f"AF-{acc_id}-F{fragment}-model_v{model_version}.{file_extension}"
    else:
        raise NotImplementedError(f"Multiple fragments for AF structure not implemented. ")



@c.command()
@c.argument(
    'PTM_DATASET', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
@c.argument('PDB_DIR', nargs=1,
    type=c.Path(exists=True, file_okay=False, dir_okay=True),
)
@c.argument('OUT_DIR', nargs=1,
    type=c.Path(exists=False, file_okay=False, dir_okay=True),
)
@c.option(
    "--organism",
    help="Filter which organism to include.  Includes all organisms by default.",

    type=c.STRING,
    default="ALL",
    show_default=True,
    
)

@c.option(
    "-S",
    "-s",
    "--ser",
    "--SER",
    "--Ser",
    "s",

    is_flag=True,

)
@c.option(
    "-T",
    "-t",
    "--THR",
    "--Thr",
    "--thr",    
    "t",        # dest
    is_flag=True,

)
@c.option(
    "-Y",
    "-y",
    "--tyr",
    "--Tyr",
    "--TYR",
    "y",
    is_flag=True,

)
@c.option(
    "--ptm",
    # TODO
)
@c.option(
    "-r",
    "--radius",
    help="The threshold radius of the motif",
    type=c.FLOAT,
    default=10.0,
    show_default=True,
)
@c.option(
    "--num",
    "--num-rows",
    "num_rows",
    help="Only use the first N rows of the dataset; or last N rows (negative number)",
    default="ALL",
    show_default=True,
    type=c.STRING, 
)
@c.option(
    "--show-matches",
    is_flag=True, 
)
@c.option(
    "--force",
    help="Override existing files", 
    is_flag=True, 
)
def main(
    ptm_dataset,
    pdb_dir, 
    out_dir,

    # options
    organism,
    s,
    t,
    y,
    
    ptm,
    radius,

    num_rows,
    show_matches, 

    force, 
):

    """TODO: prompt user for confirming output directory, filename scheme, etc. before continuing"""


    """TODO: check for overwrite (i.e. grpah files with same name already present; can be overridden with --force"""

    organism = organism.lower()

    # Num rows
    if num_rows.upper() == "ALL":
        num_rows = 0
    else:
        try:
            num_rows = int(num_rows)
        except:
            raise ValueError(f"Specified number of rows '{num_rows}' not an integer.")


    residues = ""
    if s: residues += 'S' 
    if t: residues += 'T'
    if y: residues += 'Y' 


    

    df: pd.DataFrame = load_psp_list(ptm_dataset)
    
    #df = filter_psp_list(df, )

    # Filter organisms if specified
    if not organism.upper() == "ALL":

        # if no results
        if organism in df.ORGANISM.unique():
            df = df[df["ORGANISM"] == organism]
        else:
            raise ValueError(f"Specified organism '{organism}' not in dataset.")


    # Get first (or last) N rows
    if num_rows > 0 and num_rows    <= len(df): df = df.head(num_rows)      # N:    First N rows
    elif num_rows < 0 and -num_rows <= len(df): df = df.tail(-num_rows)     # -N:   Last N rows

    # convert MOD_RSD to node ID
    df['MOD_RSD_ID'] = df.MOD_RSD.apply(mod_rsd2node_id)

    df[['MOD_RSD_POS', 'MOD_RSD_PTM']] = df.MOD_RSD.apply(lambda x: pd.Series(str(x).split('-')))

    # check that it is expected residue (if not, skip structure; sequence info in UniProt may have been updated)

    # If no residues are selected, include all. 
    # Filter to include selected residue types.
    filt = get_residue_filter(residues, invert=(not residues))
    df = df[df['MOD_RSD_ID'].apply(filt)]


    #print(f"Num unique proteins: {len(df['ACC_ID'].unique())}")

    if show_matches: 
        

        df = df.drop_duplicates(subset='ACC_ID', keep='first')
        print_dataset_pdb_matches(
            df=df,
            pdb_dir=pdb_dir,

        )
        return

    """Generate and save graphs to output directory for each psite"""


    p = Path(out_dir)

    p = p / organism  
    
    


    p.mkdir(parents=True, exist_ok=True)

    fn: str = get_graph_filename('5555', 'S100', '8.0')
    filepath = p / fn 


    # Create graph dictionary 


    

    


    # Check if filepath exists
    if os.path.isfile(filepath):
        print("ALREADY A FILE.")

    else: 
        pass
        with filepath.open("w", encoding="utf-8") as f:
            
            pickle.dump(g_dict, f)




    print(p.is_file)



    df_dict = df.to_dict('records')

    for row in tqdm(df_dict):
        
        acc_id      = row["ACC_ID"]
        mod_rsd_id  = row["MOD_RSD_ID"]
        mod_rsd     = row["MOD_RSD_POS"]
        mod_rsd_ptm = row["MOD_RSD_PTM"]


        p = p / organism / get_graph_filename(acc_id, mod_rsd, radius)
        g_out_path = p
        


        #r = radius * radius

        #print(g_out_path)

        """Create subgraph, dump subgraph"""

        p = Path(pdb_dir) 

        p = p / get_pdb_filename(acc_id=acc_id)
        p = p.resolve() # get absolute path

        pdb_path = str(p)
        if not os.path.isfile(pdb_path):
            print(f"No structure for {acc_id}")
            continue
            

        g = get_motif_subgraph(
            pdb_path=pdb_path,
            radius=radius,
            mod_rsd=mod_rsd_id,
        )

        g_dict = {
            "graph": g, 
            "kinase": "UNKNOWN",
        }


        #print(g, mod_rsd)
        #graph_dump(g, )



    # TODO: create path if doesn't exist (i.e. no such 'human' subdirectory; create as required.)
    
    #graph_dump()




    """outfile =  open(out_path, 'wb')
    pickle.dump(data, outfile)
    outfile.close()"""

    #print(df)
    


if __name__ == "__main__":
    main()

