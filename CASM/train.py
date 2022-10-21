"""Training a GCN on dataset"""







"""
MODEL 1
"""

"""
TODO:
- graph loader to generate graphs on demand for training (avoid loading all in memory at once)

SUPERVISED APPROACH: KINASE / SUBSTRATE PAIRS > GNN > CLASSIFIER 


a. Generate kinase / substrate pairs (names / descriptors in array) with train / test split 

b. Get actual graphs for feeding into network from (a.)




Kinase
protein id known
1. get active site location 
2. generate graph (using threshold distance)

Sub
protein id known, 
"""