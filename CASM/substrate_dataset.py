"""Dataset to load substrate graphs with phosphosite node index, and a label (1-hot encoding of kinase family)"""
# Uses dbPTM 

import logging as log
import os
from pathlib import Path
from typing import Callable, Dict, Generator, List, Optional, Tuple
from urllib.error import HTTPError

import networkx as nx
from tqdm import tqdm

from graphein.ml.conversion import GraphFormatConvertor
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graphs_mp, construct_graph
from graphein.protein.utils import (
    download_alphafold_structure,
    download_pdb,
    download_pdb_multiprocessing,
)
from graphein.utils.utils import import_message

try:
    import torch
    from torch_geometric.data import Data, Dataset, InMemoryDataset
except ImportError:
    import_message(
        "graphein.ml.datasets.torch_geometric_dataset",
        "torch_geometric",
        conda_channel="pyg",
        pip_install=True,
    )

import pandas as pd


"""
Convert positive examples into 1-hot encodings 
"""
def convert_kinase_list_1hot():
    pass 
    # TODO

"""
NOTE: assumes all graphs are already loaded (use ProteinGraphDataset to do this beforehand)
"""
class PhosphositeDataset(Dataset):

    def __init__(
        self, 
        name: Optional[str] = "Phosphosite",
        root: Optional[str] = None, 
        transform: Optional[Callable] = None, 
        pre_transform: Optional[Callable] = None, 
        pre_filter: Optional[Callable] = None,

        label_encoding: Optional[str] = "1-hot", # Will create multiple examples for multilabel cases; as opposed to multiple labels in same example 

    ):
        self.name = name 



        self.label_encoding = label_encoding
        if self.label_encoding in ['1-hot', 'one-hot']: 
            #TODO
            pass 


        super().__init__(root, transform, pre_transform, pre_filter)  
    

    def get(self, idx: int):

        node = self.examples[idx]['mod_rsd']
        site = torch.load(
            os.path.join(self.processed_dir, f"{self.examples[idx]['sub']}.pt")
        )

        node: str = mod_rsd 

        # This may fail; e.g. 'off by 1' error in uniprot structure perhaps 

        centre_node = site.node_id.index(node) 
        site.node_index = torch.tensor([centre_node], dtype=torch.long)

        # Label 
        label = torch.tensor(self.examples[idx]['label']).type(torch.LongTensor)

        
