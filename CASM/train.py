"""Training a GCN on dataset"""







"""
MODEL 1
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
                this way the neural network will 

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