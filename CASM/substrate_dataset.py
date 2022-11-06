"""Dataset to load substrate graphs with phosphosite node index, and a label (1-hot encoding of kinase family)"""
# Uses dbPTM 

import logging as log
import os
from pathlib import Path
import random
from typing import Callable, Dict, Generator, List, Optional, Tuple, Union
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

import torch.nn.functional as F

import pandas as pd

from CASM.load_dbPTM import get_sites, NUM_KINASE_FAMILIES, KINASE_FAMILIES, KINASE_FAMILY_TO_INDEX

"""
Get 1-hot encoding (tensor) of a kinase family
"""
def get_1hot_kinase(k: str):
    if k not in KINASE_FAMILY_TO_INDEX:
        raise KeyError(f"'{k}' not a kinase family")

    idx = KINASE_FAMILY_TO_INDEX[k]
    return F.one_hot(torch.tensor(idx), num_classes=NUM_KINASE_FAMILIES)

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
        sites_dict: dict = get_sites() 
        self.name = name 
        self.label_encoding = label_encoding
        if self.label_encoding in ['1-hot', 'one-hot']: 
            #TODO
            pass 

        # NOTE: WE DO NOT PROCESS HERE; ASSUMES ALL NECESSARY .pt FILES ARE IN /processed !!!
        self.substrates = list(set([
            acc_id 
            for (acc_id, pos) in sites_dict.keys()
            #if os.path.exists(Path(self.processed_dir) / f"{acc_id}.pt")
        ]))
        super().__init__(root, transform, pre_transform, pre_filter) 

        
        self.substrates = list(set([
            acc_id 
            for (acc_id, pos) in sites_dict.keys()
            if os.path.exists(Path(self.processed_dir) / f"{acc_id}.pt")
        ]))

        

        examples: list = [] 

        # Filter sites according to file
        fp = "./dbPTM_no_include"
        count = 0 
        with open(fp) as f:
            for line in f:
                acc, node = line.split()
                pos = int(node.split(':')[-1])

                if (acc, pos) in sites_dict:
                    #print(f"Removing {(acc, pos)} ...")
                    del sites_dict[(acc, pos)]
                    count += 1
                else:
                    print(f"{(acc, pos)} not in sites_dict to begin with.")

        print(f"Removed {count} sites from examples (node not in graph)")        
        # Filter sites 
        for (acc_id, pos), site in sites_dict.items():
    
            if acc_id not in self.substrates:   # Only include ones that exist already as .pt files
                continue
            
            mod_rsd: str = site['mod_rsd']

            # Check that node is in the graph structure 
            check_node_id_ok = True
            check_node_id_ok = False
            if check_node_id_ok:
                try:
                    fn = f"{acc_id}.pt"
                    site_pt = torch.load(
                        os.path.join(self.processed_dir, fn)
                    )
                    centre_node = site_pt.node_id.index(mod_rsd) # Will return value error if MOD_RSD is not in the graph
                except:
                    print(f"Excluding {acc_id} {mod_rsd}") # (node not in {fn} )")
                    continue # don't include in examples
            

            kinase_families = site['pos']
            for k in kinase_families:
                
                #print(f"family: {k}")

                
                label = get_1hot_kinase(k)
                data = {
                    "acc_id": acc_id, 
                    "mod_rsd": mod_rsd,
                    "label": label,
                    "kin": k,
                }

                examples.append(data)

        # DEBUG: print random sample 
        # sample = random.sample(examples, 10)
        # print(sample)
        # exit(1)
        
        self.examples = dict(enumerate(examples))
        

    def len(self) -> int:
        """Returns length of data set (number of examples)."""
        return len(self.examples)

    def process(self):
        pass

    def get(self, idx: int):
        
        example     = self.examples[idx]
        acc_id      = example['acc_id']
        node: str   = example['mod_rsd']

        site = torch.load(
            os.path.join(self.processed_dir, f"{acc_id}.pt")
        )


        # This may fail; e.g. 'off by 1' error in uniprot structure perhaps 
        centre_node = site.node_id.index(node) 
        site.node_index = torch.tensor([centre_node], dtype=torch.long)

        # normalise to get RSA 
        m = max(site.asa)
        site.asa = [a / m for a in site.asa]
        #site.asa = [site.asa]

        site.b_factor = [b / 100 for b in site.b_factor]

        # Label 
        label = torch.tensor(self.examples[idx]['label']).type(torch.LongTensor)

        # Metadata 
        metadata = {
            "acc_id":acc_id, 
            "mod_rsd":node, 

        }

        return site, label, metadata

    @property
    def processed_file_names(self) -> Union[str, List[str], Tuple]:
        
        return [
            f"{s}.pt" 
            for s in self.substrates
        ]

    @property
    def raw_file_names(self) -> Union[str, List[str], Tuple]:
        return [
            f"{s}.pdb" 
            for s in self.substrates
        ]