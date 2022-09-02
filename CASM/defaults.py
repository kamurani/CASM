"""Defaults"""
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.subgraphs import extract_subgraph_from_point
from graphein.protein.edges.distance import *	


from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa as rsa_function

def get_default_graph_config():

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
        #pdb_dir=pdb_dir,
    )

    return config