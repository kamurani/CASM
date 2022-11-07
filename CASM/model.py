# Building model
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric
from torch_geometric.nn import GCNConv, GATConv, global_max_pool as gmp, global_add_pool as gap,global_mean_pool as gep,global_sort_pool
from torch_geometric.utils import dropout_adj
from torch.optim.lr_scheduler import MultiStepLR

import numpy as np




class GCNN(nn.Module):
    def __init__(self, n_output=1, num_features_pro=1024, output_dim=128, dropout=0.2, get_embedding = False):
        super(GCNN, self).__init__()

        print('GCNN Loaded')

        # for protein 1
        self.n_output = n_output
        self.get_embedding = get_embedding
        self.pro1_conv1 = GCNConv(num_features_pro, num_features_pro)
        self.pro1_fc1 = nn.Linear(num_features_pro, output_dim)

        # for protein 2
        self.pro2_conv1 = GCNConv(num_features_pro, num_features_pro)
        self.pro2_fc1 = nn.Linear(num_features_pro, output_dim)

        self.relu = nn.LeakyReLU()
        self.dropout = nn.Dropout(dropout)
        self.sigmoid = nn.Sigmoid()

        # combined layers
        self.fc1 = nn.Linear(2 * output_dim, 256)
        self.fc2 = nn.Linear(256 ,64)
        self.out = nn.Linear(64, self.n_output)

    def forward(self, pro1_data, pro2_data):

        #get graph input for protein 1 
        #print(pro1_data)
        if type(pro1_data) == torch_geometric.data.data.Data:
                num_nodes = pro1_data.num_nodes
                features = [pro1_data.amino_acid_one_hot, pro1_data.asa, pro1_data.phi, pro1_data.psi]

                for i, f in enumerate(features):
                        features[i] = torch.as_tensor(f)

                p1 = torch.cat(features, dim=1)
                print(p1)
        else:
                
                #print(pro1_data.node_id)
                
                # Get node index for SITE
                node_index = pro2_data.node_index
                
                
                num_nodes = pro1_data.num_nodes
                features = [pro1_data.amino_acid_one_hot, pro1_data.asa, pro1_data.phi, pro1_data.psi]
                
                # Concatenate along 1st dim (dim=0) to get all nodes for all batches in one tensor
                for i, f in enumerate(features):
                        
                        features[i] = torch.as_tensor(np.concatenate([np.array(i, dtype=np.float32) for i in f]))
                        f = features[i]
                        # Vertically stack if only one feature long
                        if f.shape == torch.Size([num_nodes]):
                                features[i] = f.view(-1, 1)


                
        
                # aa = pro1_data.amino_acid_one_hot
                # p1 = torch.as_tensor(aa)
                # for f in features:
                #         print(f)

                p1 = torch.cat(features, dim=1)
                p1 = p1.float()


                
                
                #features = [p1.amino_acid_one_hot, p1.asa, p1.phi, p1.psi]

                """

                # WE NEED TO TURN LIST OF LISTS INTO NP ARRAY OF ARRAYS

                for i, f in enumerate(features):
                        
                        #f = np.vstack(f).astype(np.float)

                        #features[i] = torch.as_tensor(np.array([np.array(xi, dtype=np.float) for xi in f], dtype=np.float))
                        features[i] = torch.as_tensor(f)

                print(type(features))
                print(features[0][0])
                
                #exit(1)
                #p1 = torch.cat(features, dim=1)
                """

                #print(p1.amino_acid_one_hot)
                #p1 = torch.cat((p1.amino_acid_one_hot, p1.asa, p1.coords, p1.phi, p1.psi), dim=1)

                pro1_edge_index, pro1_batch = pro1_data.edge_index, pro1_data.batch



                # get graph input for protein 2
                pro2_edge_index, pro2_batch = pro2_data.edge_index, pro2_data.batch

                # sub_node_ids = pro2_data.node_id
                # print("node index:")
                # print(node_index)
                # print("node ids:")
                # print(sub_node_ids)
                # exit(1)
                
                num_nodes = pro2_data.num_nodes
                features = [pro2_data.amino_acid_one_hot, pro2_data.asa, pro2_data.phi, pro2_data.psi]
                
                # Concatenate along 1st dim (dim=0) to get all nodes for all batches in one tensor
                for i, f in enumerate(features):
                        
                        features[i] = torch.as_tensor(np.concatenate([np.array(i, dtype=np.float32) for i in f]))
                        f = features[i]
                        # Vertically stack if only one feature long
                        if f.shape == torch.Size([num_nodes]):
                                features[i] = f.view(-1, 1)

                p2 = torch.cat(features, dim=1)
                p2 = p2.float()

                #print(f"pro1: {p1.shape}")
                #print(f"{p1.dtype}")
                #print(f"pro2: {p2.shape}")



        x = self.pro1_conv1(p1, pro1_edge_index)
        x = self.relu(x)
        
	# global pooling
        # print("before pool: ")
        # print(x.shape)
        # print(x)
        
        x = gep(x, pro1_batch)
        
        # TODO: same for kinase? nearest node to ATP site?
        #x = x[node_index]    
        
        # print("after pool: ")
        # print(x.shape)
        # print(x)
        # exit(1)

        # flatten
        x = self.relu(self.pro1_fc1(x))
        x = self.dropout(x)



        xt = self.pro2_conv1(p2, pro2_edge_index)
        xt = self.relu(xt)

	# global pooling
        #xt = gep(xt, pro2_batch)  

        xt = xt[node_index]

        

        # flatten
        xt = self.relu(self.pro2_fc1(xt))
        if self.get_embedding: 
                
                return xt

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
        

net = GCNN(
    num_features_pro=7,
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