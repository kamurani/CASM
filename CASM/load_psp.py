"""Load data from PhosphositePlus (PSP) dataset."""

"""
Agnostic as to how the structure was determined.  Simply receives a directory of already-existing PDB files. 

"""

from email.policy import default
from utils.residue import aa1to3

import pandas as pd
import click as c
import re

import graphein 

from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1
from utils.residue import aa3to1, aa1to3


from definitions import GRAPH_NODE_FEATURES


def mod_rsd2node_id(
    mod_rsd: str, 
    chain_id: str = 'A', # TODO by default, assume that we are using AF2 i.e. always one chain, no multimeric proteins. 
) -> str:

    re.compile()

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

@c.command()
@c.argument('PTM_DATASET', nargs=1)
@c.argument('PDB_DIR', nargs=1)
@c.argument('OUT_DIR', nargs=1)

@c.option(
    "--organism",
    help="Filter which organism to include.  Includes all organisms by default.",

    type=c.STRING,
    default="ALL",
    show_default=True,

    
)
def main(
    ptm_dataset,
    pdb_dir, 
    out_dir,

    # options
    organism,
):
    ptm_dataset
    df: pd.DataFrame = load_psp_list(ptm_dataset)

    
    # Filter organisms if specified
    if not organism == "ALL":

        # if no results
        if organism in df.ORGANISM.unique():
            print("YES")
            return

        df = df[df["ORGANISM"] == organism]

    df = df[df[""]]

    """Generate and save graphs to output directory for each psite"""
    

    # convert MOD_RSD to node ID

    # check that it is expected residue (if not, skip structure; sequence info in UniProt may have been updated)

    


if __name__ == "__main__":
    main()

"""    

@c.command()
@c.argument('sites', nargs=1)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
@c.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show extensive program output."
)
@c.option(
    "--debug",
    "-d",
    is_flag=True,
    help="Show extensive program output for debugging."
)
@c.option(
    "--quiet",
    "-q",
    is_flag=True,
    help="Suppress program output."
)
@c.option(
    "--dry-run",
    "--dryrun",
    "-n",
    "is_dryrun",
    is_flag=True,
    help="Print out what the program would do without loading the graphs.",
)
@c.option(
    "--unique",
    "-u",
    is_flag=True,
    help="Only construct graphs for unique motifs.  Duplicate entries (i.e. with different kinases) are ignored."
)

@c.option(
    # TODO: support multiple formats selected at once; i.e. saves more than one file.
    "-o",
    "--graph-format",
    "--output-format",
    "--format",
    type=c.Choice(['NetworkX', 'StellarGraph', 'nx', 'sg'], case_sensitive=False),
    help="Save graphs as NetworkX or StellarGraph instances with feature preprocessing. ",
    default="NetworkX",
    show_default=True,
)
@c.option(
    "-N",
    "--num-psites",
    help="Only consider the first N motifs in a dataset.  Graph construction will continue until N graphs are made, or the end of the dataset is reached.",
    type=c.INT,
    default=-1, 
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
    "--rsa",
    "--rsa-threshold",
    help="The RSA threshold of the motif",
    type=c.FLOAT,
    default=0.0,
    show_default=True,
)
@c.option(
    "--node-features",
    "--nf",
    is_flag=False,
    default=','.join(GRAPH_NODE_FEATURES), show_default=True,
    metavar="<node_features>",
    type=c.STRING,
    help="Which node features to include in the constructed graphs."
)
@c.option(
    "--edge-features",
    "--ef",
    is_flag=False,
)
@c.option(
    "--config",
    "-c",
    help="Path to config.yml file used to specify how to construct the graphs.",
    # TODO: check if path right here?
    default="config.yml",
    show_default=True,
)
def main(
    # POSITIONAL ARGUMENTS:
    phosphosite,  
      
    structures,
    graphs,

    # FLAGS:
    verbose,
    debug,
    quiet,
    is_dryrun,
    unique,

    # PARAMETERS:
    graph_format,
    num_psites,
    radius,
    rsa,
    node_features,
    edge_features,
    config,
):

    print("test")

    path = "/home/cim/STRUCTURAL_MOTIFS/CASM/datasets/Phosphorylation_site_dataset"

    print(phosphosite)
    df: pd.DataFrame = load_psp_list(phosphosite)

    print(df)

"""