"""Generate embeddings from model for unlabelled (or labelled) data"""

import pandas as pd
import torch
from torch.utils.data import random_split
from torch_geometric.data import DataLoader

import math

from CASM.model_2 import GCNN2
from CASM.load_dbPTM import KINASE_FAMILIES, get_sites

import torch_optimizer as optim 
from torch.optim.lr_scheduler import MultiStepLR, ReduceLROnPlateau
from metrics import *


import torch.nn as nn
import networkx as nx
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree

import numpy as np


device = torch.device("cuda:1") if torch.cuda.is_available() else torch.device("cpu")

from CASM.substrate_dataset import PhosphositeDataset

model_path = "../saved_models/GCN_MODEL2_SAVE.pth" 
model_path = "../saved_models/GCN_MODEL2_SAVE.pth.pth" 
#model_path = "../saved_models/"

from CASM.model_2 import GCNN2 

def get_features(name):
    def hook(model, input, output):
        features[name] = output.detach()
    return hook 

model = GCNN2(
    features=["asa", "b_factor"],
    use_residue_encoding="1-hot",
)

site_dict = get_sites()

model.load_state_dict(torch.load(model_path))

model.eval()

model.pro_fc1.register_forward_hook(get_features('feats'))


print(model)

dataset = PhosphositeDataset(
    root="../GRAPHS3_MEILER/"
)
batch_size = 4
loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

PREDS = []
FEATS = []
features = {}

kinases = []


#df = pd.DataFrame(columns=[1])
#df[1] = df[1].astype(object)

for idx, (site, label, metadata) in enumerate(loader):

    site = site.to(device)

    acc_id = metadata['acc_id']

    for row in label:
        idx = int(torch.argmax(row))
        kin = KINASE_FAMILIES[idx]

        kinases.append(kin)


    preds = model(site)

    PREDS.append(preds.detach().cpu().numpy())
    FEATS.append(features['feats'].cpu().numpy())
    

with open("kinase_family_labels_model2.txt", 'w') as f:

    for k in kinases:
        f.write("%s\n" % k)

exit(1)
PREDS = np.concatenate(PREDS)
FEATS = np.concatenate(FEATS)

#np.save('predictions_m2_run1.npy', PREDS)
#np.save('embeddings_m2_fc1_128.npy', FEATS)

print('- preds shape:', PREDS.shape)
print('- feats shape:', FEATS.shape)



feats = FEATS
df[1] = feats
df['kinase_family'] = kinases

df.to_csv("model2_embedding_df.tsv", sep='\t')
