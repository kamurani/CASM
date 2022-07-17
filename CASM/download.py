



# for each entry in phos dataset: (of specific organism)  
# use accession id to download structure (check if {acc}.pdb already exists in dir first) 
# maybe do this for 1,000 samples randomly picked (at first); then later we can download them all. 

# save a dict with (id, list of psites).  


# for each id, 
#   for each psite,
#   construct graph around psite 
# save this list of graphs in another directory (intermediate processing step to avoid too much RAM used) 

# create dataloader for ML, using directory of graphs (if memory not enough)
# train stellargraph classification task using laplacian matrix distance, to generate embeddings

# project embeddings with method incl. UMAP, tSNE

