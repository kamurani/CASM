# CASM


**C**lustering **A**bstracted **S**tructural **M**otifs


Currently using Phosphosite structural motifs for *Mus musculus* and *Homo sapiens*. 

[PhosphoSitePlus](https://www.phosphosite.org/staticDownloads)

Predicted AF2 PDB structures for Swiss-Prot database download [here](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v3.tar) (542,380 predicted structures)


## TODO

- retrieve from PSP instead of phosphoELM
- apply cutoff for pLDDT score to remove "sequence motifs" from dataset 
- retrieve embedding from modified residue (MODRES) to get a more comparable representation hopefully 
- add more node features / edge features so that interactions can be seen by the algorithm more

- use atom graphs, as opposed to residue and compare clusterings 


- idea from N Warren about using GAT somewhere?


### Improvements
- improve structure graph retrieval / storage to not exhaust our memory :( 
- 




## TODO (future) 

- take into account what the structure is like if other PTMs are present 
- apply clustering to more PTMs


## Journal


### Hyperparameters 

- clustering residue types separately or together (e.g. are `S`, `T`, `Y` SMs comparable?)

### Sequence motif bias 

Sometimes the motif is essentially a sequence motif if it is in a disorderd region of the sequence.  The contribution on the dataset from this can be eliminated by using cutoff metrics for "disorded region".  

1. Check modified residue (ROI), if it has pLDDT < 70, ignore. 
2. More nuanced: look at surrounding residues too; maybe have cutoff by just removing poor confidence residues  from the structural motif (maybe also taking into account if they are sequence-adjacent; versus a similarly low-confidence residue but it is quite far apart -- i.e. it just happens to be nearby but it should not be consistently nearby)
3. etc. 

Keep in mind however that it is a good "sanity check" of whether our structural method actually sees *anything* meaningful; i.e. it should still be able to pick up on sequence motifs and classify these well; but note that this does not guarantee that more "structural" motifs are adequately interpreted.  It shows that the model "isn't wrong", but doesn't necessarily show that it "is right". 


### Unsupervised methods

- unsupervised graph representation learning (learn embeddings)
- use graph-graph distance as a metric to learn by 
- experiment with using different edges (i.e. will have multiple adjacency matrices) 
- use PyMol (or similar) `fit` or `super` function to get similarity metric between two "bubbles" of given radius and confidence (pLDDT) cutoff 
- then do hierarchical / k-means clustering 
- tSNE / UMAP visualisation 
- then apply labels if known (colour by) 
- perhaps use knowledge of kinase ontologies if family relationships are known and meaningful


### Supervised methods 

- use identity of upstream kinase (if known) to enforce algorithm to "pay attention" to certain features that may determine which kinase is implicated (i.e. forces it to take into account the actual relevance of the structural motif)  


### Semi-supervised methods 

Perhaps something similar to [northstar](https://github.com/northstaratlas/northstar) where existing classes are known, but new clusterings can be found also
