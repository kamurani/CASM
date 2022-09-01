
import networkx as nx

from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	


from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa as rsa_function


def get_motif_subgraph(
    pdb_path, #TODO type
    mod_rsd: str,
    radius: float = 10.0, 
    rsa: float = 0.0,

) -> nx.Graph:

    if type(pdb_path) is not str: pdb_path = str(pdb_path)
    
    edge_fns = [
        add_aromatic_interactions,
        add_hydrophobic_interactions,
        add_aromatic_sulphur_interactions,
        add_cation_pi_interactions,
        add_disulfide_interactions,
        add_hydrogen_bond_interactions,
        add_ionic_interactions,
        add_peptide_bonds
    ]
    config = ProteinGraphConfig(
        edge_construction_functions=edge_fns, 
        graph_metadata_functions=[rsa_function], 
        dssp_config=DSSPConfig(),
    )
    g = construct_graph(config, pdb_path=pdb_path)


    return g