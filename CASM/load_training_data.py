"""Graphein dataset download method"""

import pandas as pd
import numpy as np

import torch
import torch.nn as nn
import pytorch_lightning as pl
from tqdm.notebook import tqdm
import networkx as nx
import torch_geometric
from torch_geometric.data import Data
from torch_geometric.utils import from_networkx
from sklearn.preprocessing import LabelBinarizer
from sklearn.metrics import f1_score

import math
import glob

import warnings
warnings.filterwarnings("ignore")

import torch
from graphein.ml import ProteinGraphDataset
from graphein.ml.conversion import GraphFormatConvertor
import graphein.protein as gp

###

from CASM.kin_sub_pairs import get_filtered_dataset

fp = "../datasets/Kinase_Substrate_Dataset"
df = get_filtered_dataset(
    path_to_dataset=fp
)

#print(df)

kinases = df.KIN_ACC_ID.unique()
substrates = df.SUB_ACC_ID.unique()

acc_ids = list(kinases) + list(substrates)


print(f"Length: {len(acc_ids)}")

chain = 'A'
chain_selection_map = {}
for a in acc_ids:
    chain_selection_map[a] = chain 

protein_graph_config = gp.ProteinGraphConfig() # TODO: add RSA etc. 


def filter_func(data: torch_geometric.data.Data):

    # TODO
    # Filter based on whether we include in dataset 
    return True




radius = 10.0
def get_subgraph(g: nx.Graph):


    # TODO: 
    # get subgraph around node somehow 
    # store it in ``g`` beforehand...?
    return g

ds = ProteinGraphDataset(
    root = "../alphafold_structures",
    uniprot_ids=acc_ids,

    #graph_label_map=g_lab_map, # NO LABELS FOR NOW
    #node_label_map=node_lab_map,

    #chain_selection_map=chain_selection_map,
    graphein_config=protein_graph_config,

    graph_transformation_funcs=[get_subgraph],

    graph_format_convertor=GraphFormatConvertor(
            src_format="nx", dst_format="pyg"
        ),

    
    
)