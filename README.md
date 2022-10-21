# CASM


**C**lustering **A**bstracted **S**tructural **M**otifs


Currently using Phosphosite structural motifs for *Mus musculus* and *Homo sapiens*. 

[PhosphoSitePlus](https://www.phosphosite.org/staticDownloads)

Predicted AF2 PDB structures for Swiss-Prot database download [here](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v3.tar) (542,380 predicted structures)


## TODO

- use GraphSAGE for representation learning of all nodes in structural motifs 
- use embeddings from psite nodes to do clustering / classification learning? 



## Data processing

PSP dataset

```
# Unique proteins 
cat ../datasets/Phosphorylation_site_dataset | tail -n +5 | cut -f3 | sort | uniq | wc -l 
```


Setup 

```
export PYTHONPATH="${PYTHONPATH}:/home/cim/STRUCTURAL_MOTIFS/CASM/"
```

## Unspervised 

GraphCL -- contrast learning. 

- can be used at node-level, graph-level, or **both**. 

## SCISM 

- use dataset of known interactions (substrate-kinase pair) for training / testing / validation
- use atom or residue graphs in pairs; generate pair of embeddings; concatenate into a vector; apply FC layers of NN to output; repeat for every possible pair;  or apply a GAT/GCN on the graphs themselves (no vectorisation during learning?) 
- visualise CAM for identifying important residues 


alternative idea:
- don't feed in kinase structures; rather, generate representations of kinase structure (similar to what SeqVec does), and this is internally represented in the NN somehow.  this info is used when a protein graph (substrate) is fed in; the network will output scores for each kinase. 


alternative alternative idea:
- just have 127 separate neural networks (each for a given kinase); you feed in a psite and it spits out a score. 


IDEA: transfer learning
- pretrain on large unlabelled dataset (get representation of proteins); then finetune with smaller labelled dataset.  but use same model obviously; it has 'seen' more info 

## TODO

- retrieve from PSP instead of phosphoELM
- apply cutoff for pLDDT score to remove "sequence motifs" from dataset 
- retrieve embedding from modified residue (`MOD_RSD`) to get a more comparable representation hopefully 
- compare this with `RMSD` of the aligned motifs (with `MOD_RSD` of each enforced to be superimposed)
- 
- add more node features / edge features so that interactions can be seen by the algorithm more

- use atom graphs, as opposed to residue and compare clusterings 
- consider all PTMs, not just phosphorylation


- idea from N Warren about using GAT somewhere?

### Training 

- unsupervised representation learning, apply clustering, no labels used 
- supervised: use labels for classification task; THEN get embeddings (or then just run the model on ALL sites, known and unknown)
- semi-supervised: allow clustering for unknown samples and classification for known samples (but together)

### Improvements
- improve structure graph retrieval / storage to not exhaust our memory :( 
- 

- use environment variables e.g. `$STRUCTURE_DIR` and `$PTM_DATASET_PATH` for running workflows easier without continuously specifying paths to things.



### Generalisable to PDB not just AF2?

- [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html) for cross referencing sequence positions


- Deal with uncurated sequences (e.g. A0A6J1BG53 (A0A6J1BG53_9ROSI) ) from alphafold 


## TODO (future) 

- take into account what the structure is like if other PTMs are present 
- apply clustering to more PTMs


- incorporate kinase structure (like when doing protein docking??)
- look at GNN approaches in the literature that did docking type things; could actually incorporate all the known kinase structures AS WELL AS the phosphosites; and then get a score for each kinase from a GNN model (e.g. create residue- or atom- network graph, then the model has a representation of each kinase as well; and it gives a log-flattened score for each (e.g. .9 for one kinase, 0.01 for the others)


## Workflow 

- have a separate command for 'filter dataset' and 'load graphs', allowing control (and config files) to be decoupled for each of the two processes
- may require a temporary dataset file that is created when 'filter dataset' is run; the path to this temporary file may be another environment variable (with some default value)


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
