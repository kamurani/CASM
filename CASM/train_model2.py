"""Training a GCN on dataset"""



"""
MODEL 0 a
- simply use kinase labels 


MODEL 0 b
- kinase / sub structure input: yes/no output
"""

"""

"""


"""
MODEL 1
"""

"""
QUESTIONS: 


    - do train / test split (not manually) with random seed, BUT takes into account classes with limited number of examples (e.g. kinase with 1 example)


    - SEMIsupervised: incorporate clustering somehow (so takes into account separation / dissimilarity of substrate protein structure at psite; AND kinase labels )

    - label with kinase superset i.e. family? instead of specific ID


    - apply 'clustering' to shift around the points that were pre-calculated from supervised learning ? or vice versa?
    - i.e. use one approach to get rough embeddings; then 'nudge' the embeddings using the other method.


    - for 'autoencoder' component with encoder / decoder, the target output doesn't need to be all the node features etc. 
    - it only needs to contain each amino acid (could be 1-hot encoding); as a way to go from a vector to a graph correctly.  we only need the adjacency matrix,
    edge distances, and which AA. 
"""


"""
TODO:
- graph loader to generate graphs on demand for training (avoid loading all in memory at once)

SUPERVISED APPROACH: KINASE / SUBSTRATE PAIRS > GNN > CLASSIFIER 

    i. substrate / kinase input: 0 or 1 output
    ii.  substrate input: 1-hot output for all kinases (model must know structure implicitly somehow)

        IDEA: 
            what if we incorporate the input (substrate psite graph) and the kinase structure knowledge, by training
            it with the kinase structure (the graph itself somehow (autoencoder?), embedding, etc.) AS THE EXPECTED
            OUTPUT?

            e.g. substrate / kinase training examples: 
                input: substrate psite embeddings; output: substrate embedding
                this way the neural network will have a more 'continuous' output of a corresponding kinase active site structure.  

                this can be statistically evaluated to find the "closest" kinase (out of the training data output exmamples) either as a reconstructed
                graph (from the output vector), or the vector itself.  i.e. give log likelihoods for each kinase as usual. 


            Also, incorporate LSTM language model of protein sequence for amino acid features (and use whole protein in graph creation), 
            but "centre" the GAT/GCN on the phosphosite residue such that only local features to that region are used in downstream neural network layers. 
            Picking only 2 or 3 layers for the GCN would mean that the rest of the protein has no impact; but the language model still might serve well for 
            this local region, while initially 'seeing' the whole protein ? (will have to empirically test this with different methods e.g. 1-hot amino acid encoding etc.)

            
            To get the 'labels' as vectors, may have to first apply a (different) GNN to just get a way of turning a kinase into a vector. 

            To get the 'labels' as structure, will have to use an autoencoding approach in reverse? 

            To get the 'labels' as vectors, what if we use the SAME GCN way of turning a SUBSTRATE motif into a a vector (which then goes into the FC 
            classifier labels or whatever), but use THOSE weights to generate the 'label' in real time? i.e. input kinase and substrate; both get turned into respective 
            vectors using weights; kinase then gets used as the 'label' for the downstream part of the network.  or will that completely not work at all?


            What if we use a Seq2Vec approach to get embeddings of the kinase; and this is what we give as the 'labels'. 


a. Generate kinase / substrate pairs (names / descriptors in array) with train / test split 
    - take into account kinases with 1 sub (don't use this in test)
    - randomly generate false pairs that are known to be non-interacting 

b. Get actual graphs for feeding into network from (a.)




Kinase
protein id known
1. get active site location 
2. generate graph (using threshold distance)

Sub
protein id known, 
"""


"""
# TODO
- use cross entropy loss (not MSE)
- balance classes 


# ALTERNATIVE
- train on average structure of kinase family 
- use pos / neg examples from dbPTM 
"""

import pandas as pd
import torch
from torch.utils.data import random_split
from torch_geometric.data import DataLoader

import math

from CASM.model_2 import GCNN2

import torch_optimizer as optim 
from torch.optim.lr_scheduler import MultiStepLR, ReduceLROnPlateau
from metrics import *


import torch.nn as nn
import networkx as nx
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree


device = torch.device("cuda:1") if torch.cuda.is_available() else torch.device("cpu")

from CASM.substrate_dataset import PhosphositeDataset

dataset = PhosphositeDataset(
    root="../GRAPHS3_MEILER/"
)
size = len(dataset)
print(f"Length: {size}")

n_epochs_stop = 6
epochs_no_improve = 0
early_stop = False

min_loss = 100


train_ratio = 0.8

batch_size = 4#4 # 64?
num_workers = 0


# Train / test split
training_data, test_data = random_split(dataset, [math.floor(train_ratio * size), size - math.floor(train_ratio * size)])

# Count examples 
def count_examples(dataset):

    count = {
        0: 0,
        1: 0,
    }
    for kin, sub, metadata, label in dataset:

        label = label.flatten().tolist()
        label = int(label[0])
        count[label] += 1

    return count

# for ti in ["training_data", "test_data"]:
#     t = count_examples(eval(ti))
#     total = t[1] + t[0]
#     print(f"{ti}: {t[1]} pos, {t[0]} neg \t {t[1]/(t[0]+t[1]):.2f}:{t[0]/(t[0]+t[1]):.2f} POS:NEG ratio")

  



print(f"Training: {len(training_data)}")
print(f"Testing: {len(test_data)}")

# Kinase / substrate dataset 
train_dataloader = DataLoader(training_data, batch_size=batch_size, num_workers=num_workers, shuffle=True)
test_dataloader = DataLoader(test_data, batch_size=batch_size, num_workers=num_workers, shuffle=True)

print(f"trainloader: {len(train_dataloader)}")
print(f"test loader: {len(test_dataloader)}")

# TODO: validation data loader

def predict(model, device, loader):
    model.eval()
    predictions = torch.Tensor()
    labels = torch.Tensor()
    with torch.no_grad():
        for site, label in loader:
            site = site.to(device)

            #print(torch.Tensor.size(kin.x), torch.Tensor.size(sub.x))
            output = model(site)
            predictions = torch.cat((predictions, output.cpu()), 0)
            labels = torch.cat((labels, label.view(-1,1).cpu()), 0)

    labels = labels.numpy()
    predictions = predictions.numpy()
    return labels.flatten(), predictions.flatten()

def train(model, device, train_dataloader, optimizer, epoch):

    print(f"Training on {len(train_dataloader)} samples....")
    model.train()
    
    # Loss function 

    #loss_func = nn.MSELoss()
    # TODO: assign weight to class (if unbalanced)

    weight = None 
    #weight = torch.tensor([0.9])
    loss_func = nn.CrossEntropyLoss(
        weight=weight,
    )

    predictions_tr = torch.Tensor()
    scheduler = MultiStepLR(optimizer, milestones=[1,5], gamma=0.5)
    labels_tr = torch.Tensor()

    for count, (site, label) in enumerate(train_dataloader):
    
        #print(f"kinase: {type(kinase)}")
        site = site.to(device)
        site = site.to(device)

        optimizer.zero_grad()
        output = model(site)
        predictions_tr = torch.cat((predictions_tr, output.cpu()))

        labels_tr = torch.cat((labels_tr, label.view(-1,1).cpu()), 0)
        loss = loss_func(output, label.view(-1,1).float().to(device))
        loss.backward()
        optimizer.step()

    scheduler.step()
    labels_tr = labels_tr.detach().numpy()
    predictions_tr = predictions_tr.detach().numpy()
    acc_tr = get_accuracy(labels_tr, predictions_tr , 0.5)
    print(f'Epoch {epoch-1} / 30 [==============================] - train_loss : {loss} - train_accuracy : {acc_tr}')


model = GCNN2(
    num_features_pro=7 + 2,
)
#model = model.float()
model.to(device)
num_epochs = 150
best_accuracy = 0
lr = 0.001 # hyperparam
#lr = 0.0005

optimizer = torch.optim.Adam(model.parameters(), lr=lr)

for epoch in range(num_epochs):
    
    train(model, device, train_dataloader, optimizer, epoch+1)
    G, P = predict(model, device, test_dataloader)

    loss = get_mse(G,P)
    accuracy = get_accuracy(G,P, 0.5)
    print(f'Epoch {epoch}/ {num_epochs} [==============================] - val_loss : {loss} - val_accuracy : {accuracy}')


    if(accuracy > best_accuracy):
        best_accuracy = accuracy
        best_acc_epoch = epoch
        torch.save(model.state_dict(), "../saved_models/GCN_MODEL2_SAVE.pth") #path to save the model
        print("Model")
    if(loss < min_loss):
        epochs_no_improve = 0
        min_loss = loss
        min_loss_epoch = epoch
    elif loss > min_loss:
        epochs_no_improve += 1

    if epoch > 5 and epochs_no_improve == n_epochs_stop:
        print('Early stopping!' )
        early_stop = True
        break

print(f'min_val_loss : {min_loss} for epoch {min_loss_epoch} ............... best_val_accuracy : {best_accuracy} for epoch {best_acc_epoch}')
print("Model saved")





