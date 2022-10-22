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

# s
