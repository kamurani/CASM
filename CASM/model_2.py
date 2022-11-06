"""MODEL 2: Phosphosite input, 1-hot encoding of kinase family output"""

from typing import Union
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric
from torch_geometric.nn import GCNConv, GATConv, global_max_pool as gmp, global_add_pool as gap,global_mean_pool as gep,global_sort_pool
from torch_geometric.utils import dropout_adj
from torch.optim.lr_scheduler import MultiStepLR

import numpy as np


from load_dbPTM import KINASE_FAMILY_DICT


class GCNN2(nn.Module):
    def __init__(
        self, 
        n_output: int = len(KINASE_FAMILY_DICT), 
        num_features_pro: int = None, # Number of features is 1280 if using ESM per-residue embeddings.
        
        output_dim: int = 128, # Length after the first FC layer 
        dropout: float = 0.2,

        features: list = [],

        embedding_method: str = "node",         # or "gep" i.e. global pooling
        use_residue_encoding: str = "meiler", #"1-hot",    # how to encode the residues (NOTE: assumes these features have already been created in the sample data)
        output_activation: str = "softmax",     # TODO: switch between this and individual sigmoid
    ):

        if num_features_pro is None:
            em_len = 20 if embedding_method == "1-hot" else 7 
            num_features_pro = 20 + len(features)

        self.features = features

        super(GCNN2, self).__init__()

        print('GCNN Loaded')
        print(f"Using {use_residue_encoding} node embeddings ...")

        self.embedding_method = embedding_method.lower()
        self.use_residue_encoding = use_residue_encoding.lower()
        self.output_activation = output_activation.lower()
        # for protein 1
        self.n_output = n_output
        self.pro_conv = GCNConv(num_features_pro, num_features_pro)

        # Layer 2 
        # TODO

        self.pro_fc1 = nn.Linear(num_features_pro, output_dim)

        self.relu = nn.LeakyReLU()
        self.dropout = nn.Dropout(dropout)

        self.sigmoid = nn.Sigmoid()
        self.softmax = nn.Softmax(dim=1)

        self.fc1 = nn.Linear(output_dim, 128)
        self.fc2 = nn.Linear(128, 64)
        self.out = nn.Linear(64, self.n_output)

    def forward(self, phosphosite: Union[torch_geometric.data.Data, torch_geometric.data.Batch]):

        
        node_index  = phosphosite.node_index
        num_nodes   = phosphosite.num_nodes
        psite_batch = phosphosite.batch
        edge_index  = phosphosite.edge_index
        
        # Get node features as required 
        
        if self.use_residue_encoding in ["one-hot", "1-hot"]:
            residue_encoding = phosphosite.amino_acid_one_hot
        elif self.use_residue_encoding in ["meiler", "meiler-embedding"]:
            residue_encoding = phosphosite.meiler 
        elif self.use_residue_encoding in ["esm", "esm-embedding"]:
            pass
            residue_encoding = phosphosite.esm_embedding # TODO: check this 

        else:
            raise NotImplementedError(f"Residue encoding '{self.use_residue_encoding}' not implemented.")

        # normalise 
        #print(phosphosite.asa)
        #asa = F.normalize(torch.tensor(phosphosite.asa[0]).float(), p=2, dim=0)

        #print(asa)
        
        #b_factor = [[b / 100 for b in phosphosite.b_factor[0]]] # Get pLDDT as percentage

        #b_factor = phosphosite.b_factor
        #asa = phosphosite.asa


        #self.features 

        #features = [eval(f"phosphosite.{f}") for f in self.features]
        features = []
        if "asa" in self.features:
            features.append(phosphosite.asa)
        if "b_factor" in self.features:
            features.append(phosphosite.b_factor)
        
        features.append(residue_encoding)
        #features = [residue_encoding, asa, b_factor]   # , pro2_data.asa, pro2_data.phi, pro2_data.psi]
        
        # Concatenate along 1st dim (dim=0) to get all nodes for all batches in one tensor
        for i, f in enumerate(features):
                
            features[i] = torch.as_tensor(np.concatenate([np.array(i, dtype=np.float32) for i in f]))
            f = features[i]
            # Vertically stack if only one feature long
            if f.shape == torch.Size([num_nodes]):
                features[i] = f.view(-1, 1)

        # Concatenate all features together
        x = torch.cat(features, dim=1)


        x = x.float()

        

        # First graph convolution layer 
        x = self.pro_conv(x, edge_index)
        x = self.relu(x) 

        # Pooling
        if self.embedding_method in ["gep", "global"]:
            x = gep(x, psite_batch)
        elif self.embedding_method in ["node", "site"]:
            x = x[node_index]
        else: 
            raise NotImplementedError(f"Embedding method '{self.embedding_method}' not implemented.")
        
    
        x = self.pro_fc1(x)
        x = self.relu(x)
        x = self.dropout(x)

        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.dropout(x)
        out = self.out(x)

        if self.output_activation in ["softmax"]:
            out = self.softmax(out)
        else: 
            raise NotImplementedError(f"Output activation '{self.output_activation}' not implemented.")

	
        return out
        

net = GCNN2(
    num_features_pro=9, # 22, 
)
print(net)

"""# GAT"""

class AttGNN(nn.Module):
    def __init__(self, n_output=1, num_features_pro= 1024, output_dim=128, dropout=0.2, heads = 1 ):
        super(AttGNN, self).__init__()

        print('AttGNN Loaded')

        self.hidden = 8
        self.heads = 1
        
        # for protein 1
        self.pro1_conv1 = GATConv(num_features_pro, self.hidden* 16, heads=self.heads, dropout=0.2)
        self.pro1_fc1 = nn.Linear(128, output_dim)


        # for protein 2
        self.pro2_conv1 = GATConv(num_features_pro, self.hidden*16, heads=self.heads, dropout=0.2)
        self.pro2_fc1 = nn.Linear(128, output_dim)

        self.relu = nn.LeakyReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(dropout)
        # combined layers
        self.fc1 = nn.Linear(2 * output_dim, 256)
        self.fc2 = nn.Linear(256, 64)
        self.out = nn.Linear(64, n_output)
        


    def forward(self, pro1_data, pro2_data):

        # get graph input for protein 1 
        pro1_x, pro1_edge_index, pro1_batch = pro1_data.x, pro1_data.edge_index, pro1_data.batch
        # get graph input for protein 2
        pro2_x, pro2_edge_index, pro2_batch = pro2_data.x, pro2_data.edge_index, pro2_data.batch
         
        
        x = self.pro1_conv1(pro1_x, pro1_edge_index)
        x = self.relu(x)
        
	# global pooling
        x = gep(x, pro1_batch)  
       
        # flatten
        x = self.relu(self.pro1_fc1(x))
        x = self.dropout(x)



        xt = self.pro2_conv1(pro2_x, pro2_edge_index)
        xt = self.relu(self.pro2_fc1(xt))
	
	# global pooling
        xt = gep(xt, pro2_batch)  

        # flatten
        xt = self.relu(xt)
        xt = self.dropout(xt)

	
	# Concatenation
        xc = torch.cat((x, xt), 1)

        # add some dense layers
        xc = self.fc1(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        xc = self.fc2(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        out = self.out(xc)
        out = self.sigmoid(out)
        return out

net_GAT = AttGNN()
print(net_GAT)