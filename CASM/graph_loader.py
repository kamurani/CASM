


"TODO: add pipeline feature for subgraph steps"
"Get graphs from a set directory given acc_id"
import os
from pathlib import Path
from typing import Callable, List, Union

import networkx as nx

from CASM.utils.a import get_pdb_filename, is_list_of_str
from CASM.subgraphs import get_motif_subgraph

from graphein.protein.graphs import construct_graph



class GraphLoader:
    def __init__(
        self,
        config, 
        pdb_dir: Union[str, Path],
        extension: str = "pdb", # XXX unused
        ignore: bool = False, # ignore missing directories 

        functions: List[Callable] = None, # TODO: specify subgraph selecting functions to be applied 
    ):

        # TODO: store `pdb_dir` in the config , more elegant
        # TODO: check if exists, and contains expected filename extensions
        if type(pdb_dir) is str: 
            pdb_dir = Path(pdb_dir)

        self.pdb_dir = pdb_dir
        if not ignore and not os.path.isdir(str(pdb_dir)):
            raise ValueError(f"No such directory '{str(pdb_dir)}'") 

        self.config = config
        self.extension = extension

    
    
    "Get PDB path as a string"
    def get_pdb_path(self, acc_id: str):

        p = self.pdb_dir / get_pdb_filename(acc_id, file_extension=self.extension)
        p = p.resolve() # get absolute path
        pdb_path = str(p)
        return pdb_path

    "Get graph object from an acc_id, or list of them"
    def get_graph(self, acc_id: Union[str, List[str]]): 

        if type(acc_id) is str: 
            pdb_path: str = self.get_pdb_path(acc_id) 
            
            if not os.path.isfile(pdb_path):
                return None 

            g = construct_graph(config=self.config, pdb_path=pdb_path)
            return g 
        
        elif is_list_of_str(acc_id):
            
            graphs = {}
            for i in acc_id: 

                g: nx.Graph = self.get_graph(i)
                graphs[i] = g
            return graphs 

        else: 
            raise ValueError(f"acc_id must be string or list of strings.")

    

class MotifLoader(GraphLoader):
    def __init__(self, config, pdb_dir: Union[str, Path], extension: str = "pdb", ignore: bool = False, functions: List[Callable] = None,
        radius: float = 10.0, # radius threshold (Å)
        rsa: float = 0.0, 

    ):
        super().__init__(config, pdb_dir, extension, ignore, functions)

        self.radius = radius 
        self.rsa = rsa 

    # TODO: similar to GraphLoader, handle `List[str]` of acc_ids ?
    def get_motif(self, acc_id: str, mod_rsd: str):

        g = self.get_graph(acc_id) 
        g.name = f"{acc_id}@{mod_rsd}"

        # **Reasons for failure**
        # * `mod_rsd` not in `g` (e.g. RSA threshold too high)
        
        try: 
            s_g = get_motif_subgraph(
                g, 
                mod_rsd=mod_rsd, 
                radius=self.radius,
                rsa=self.rsa,
            )
        
        except:
            s_g = None

            # DEBUG 
            print("Subgraph creation failed.")

        return s_g