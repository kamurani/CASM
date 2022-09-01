"""Load data from PhosphositePlus (PSP) dataset."""

"""
Agnostic as to how the structure was determined.  Simply receives a directory of already-existing PDB files. 

"""

from pathlib import Path
from typing import Callable, List, Union
from utils.residue import aa1to3

import pandas as pd
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



@c.command()
@c.argument(
    'PTM_DATASET', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
@c.argument('PDB_DIR', nargs=1,
    type=c.Path(exists=True, file_okay=False, dir_okay=True),
)
@c.argument('OUT_DIR', nargs=1,
    type=c.Path(exists=True, file_okay=False, dir_okay=True),
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
):

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

    print(f"Length: {len(df)}")
    return
    # convert MOD_RSD to node ID
    df['MOD_RSD_ID'] = df.MOD_RSD.apply(mod_rsd2node_id)



    df[['MOD_RSD_POS', 'MOD_RSD_PTM']] = df.MOD_RSD.apply(lambda x: pd.Series(str(x).split('-')))

    # check that it is expected residue (if not, skip structure; sequence info in UniProt may have been updated)

    # If no residues are selected, include all. 
    # Filter to include selected residue types.
    filt = get_residue_filter(residues, invert=(not residues))
    df = df[df['MOD_RSD_ID'].apply(filt)]


    print(f"Num unique proteins: {len(df['ACC_ID'].unique())}")


    """Generate and save graphs to output directory for each psite"""


    p = Path(out_dir)


    df_dict = df.to_dict('records')

    return

    for row in tqdm(df_dict):
        
        acc_id      = row["ACC_ID"]
        mod_rsd     = row["MOD_RSD_POS"]
        mod_rsd_ptm = row["MOD_RSD_PTM"]


        Path()

        p = p / organism / get_graph_filename(acc_id, mod_rsd, radius)
        g_out_path = p
        


        #r = radius * radius

        #print(g_out_path)

        """Create subgraph, dump subgraph"""

        g = get_motif_subgraph 


        graph_dump(g, )



    # TODO: create path if doesn't exist (i.e. no such 'human' subdirectory; create as required.)
    
    #graph_dump()




    """outfile =  open(out_path, 'wb')
    pickle.dump(data, outfile)
    outfile.close()"""

    #print(df)
    


if __name__ == "__main__":
    main()

