
import networkx as nx

from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.subgraphs import extract_subgraph_from_point
from graphein.protein.edges.distance import *	


from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa as rsa_function

'''
Given a graph ``g`` get a subgraph from radius and known phos site
'''
def get_protein_subgraph_radius(g, site, r=10.0):

    # get centre point   
    try:
        x_y_z = node_coords(g, site)
    except ValueError:
        raise ValueError("Specified phospho site isn't in correct format.")       
    
    # Get subgraph
    s_g = extract_subgraph_from_point(g, centre_point=x_y_z, radius=r)
    
    return s_g

def get_motif_subgraph(
    g: nx.Graph, 
    mod_rsd: str, # ID of form A:RES:123
    radius: float = 10.0, 
    rsa: float = 0.0,

) -> nx.Graph:

    #if type(pdb_path) is not str: pdb_path = str(pdb_path)
    
    try:
        g_site: Dict = g.nodes(data=True)[mod_rsd]
    except:
        raise KeyError(f"Node '{mod_rsd}' not in graph.")
                
    g_site_rsd: str = str(g_site['residue_name'])
    g_site_pos: int = int(g_site['residue_number'])

    # Subgraph (radius)
    s_g = get_protein_subgraph_radius(g=g, site=mod_rsd, r=radius)

    print(s_g.nodes())

    # Subgraph (rsa)


    return g