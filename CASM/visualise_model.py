"""Visualise model's predictions"""

from typing import List
import pandas as pd
import torch
from torch.utils.data import random_split
from torch_geometric.data import DataLoader

import math
import os

#from CASM.load_training_data import dataset
from CASM.model import GCNN, AttGNN

import torch_optimizer as optim 
from torch.optim.lr_scheduler import MultiStepLR, ReduceLROnPlateau
from metrics import *

from torch_geometric.data import Batch, Data


import torch.nn as nn
import networkx as nx
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree

def main():

    device = torch.device("cuda:1") if torch.cuda.is_available() else torch.device("cpu")

    path = "../GCN.pth"
    path = "../saved_models/GCN_SAVE.pth"

    model = GCNN(
        num_features_pro=23,
    )   
    model.load_state_dict(torch.load(path))

    model.eval()

    # Should be 1
    # TODO; make list of each and cycle through to generate heatmap

    kin = "KIN_P49841.pt"
    sub = "SUB_O00429_S40.pt"

    kin1 = "P49841"
    kin2 = "Q8IV63"


    

    import glob 
    import re

    processed_dir = "../GRAPHS1/processed/"

    kin_filenames = glob.glob("KIN_*", root_dir=processed_dir)
    kin_ids = []
    for k in kin_filenames:
        m = re.search(r'KIN_(.*?).pt', k)
        kin_ids.append(m.group(1))

    # Sorted by acc id order
    kin_ids = sorted(kin_ids)
    # Mapping index to acc_id 
    kin_acc2idx = {val:key for key, val in enumerate(kin_ids)}


    site_filenames = glob.glob("SUB_*", root_dir=processed_dir)
    site_ids = []
    for s in site_filenames:
        m = re.search(r'SUB_(.*?)_(.*?).pt', s)
        sub, mod_rsd = m.group(1), m.group(2)
        site_ids.append((sub, mod_rsd))
      
    # Sort by protein 
    site_ids = sorted(site_ids, key=lambda tup: tup[0])
    site2idx = {val:key for key, val in enumerate(site_ids)}

    #s = "O00429 S40"
    #sub, mod_rsd = s.split('_')



    rows = []
    for (sub, mod_rsd) in site_ids[0:150]:
        # Generate batch for one site, that has all kinases 
        site: Data = torch.load(os.path.join(processed_dir, f"SUB_{sub}_{mod_rsd}.pt"))
        kinases = []
        for k in kin_ids:
            kinase: Data = torch.load(os.path.join(processed_dir, f"KIN_{k}.pt"))
            kinases.append(kinase)

        sites: List[Data] = [site] * len(kinases)





        #kin1 = torch.load(os.path.join(processed_dir, f"KIN_{kin1}.pt"))
        #kin2 = torch.load(os.path.join(processed_dir, f"KIN_{kin2}.pt"))

        #sub1 = torch.load(os.path.join(processed_dir, f"SUB_{sub}_{mod_rsd}.pt"))
        #sub2 = torch.clone(sub1)

        kinases: Batch = Batch.from_data_list(kinases)
        sites: Batch = Batch.from_data_list(sites)

        kinases = kinases.to(device)
        sites = sites.to(device)

        output = model(kinases, sites)

        output = output.flatten().tolist()

        rows.append(output)

    #print()
    import plotly.express as px
    fig = px.imshow(
        rows
    )

    fig.write_html("plot.html")

    return 



    # Test loader
    train_ratio = 0.8
    size = len(dataset)
    training_data, test_data = random_split(dataset, [math.floor(train_ratio * size), size - math.floor(train_ratio * size)])
    print(size, len(training_data), len(test_data))
    test_dataloader = DataLoader(test_data, batch_size=2, num_workers=0, shuffle=False)

    

    

    #for kin, sub, (kin_name, sub_name, mod_rsd), label in test_dataloader:
        
        # kin = kin.to(device)
        # sub = sub.to(device)
        # output = model(kin, sub)

        # print(f"KIN {kin_name} SITE {sub_name} {mod_rsd} LABEL {label} PREDICTION {output}")
        # break

    

    
    

if __name__ == "__main__":
    main()